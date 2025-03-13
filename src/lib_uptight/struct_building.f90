! This file is part of uptight.
!
! uptight is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! uptight is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with uptight. If not, see <https://www.gnu.org/licenses/>.
!
! The module contains routines for generating more advanced structures than wha
! is possible withing the supercell format in the original data_file. 
! The module allows: 
!    Easy creation of hetrostructures
!    Creation of nanowires and nanocrystals with more than 4 side facets
!    Hydrogen passivation of dangling bonds
!    Generation of gen file 
!    Reading an gen file into the UPTIGHT program
!
MODULE struct_building

  !===========================================================================
  USE constants, only : pi
  USE precision
  USE type_defs
  USE lists
  USE input_output
  USE errors
  USE checks, only : check_ion
  USE vectors, only : norm
  USE neighbours
  USE exceptions
  USE rcm_module
  !===========================================================================
  
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: init_basis, make_basis
  PUBLIC :: init_structure
  PUBLIC :: write_gen, read_gen, subs_dg_ions
  !PUBLIC :: save_struct
  !PUBLIC :: struct_build, new_lists, load_struct_def, inside

  REAL ( dp ), PUBLIC, PARAMETER :: nn_max_dist= 4.0d0

CONTAINS


!=================================================================================
SUBROUTINE write_gen(gen_filename, basis, neig_map, rotation, saveH)

  CHARACTER (LEN = 200) :: gen_filename
  TYPE( ion_basis ) :: basis
  TYPE(nearest_neighbours), DIMENSION(:), POINTER :: neig_map
  LOGICAL :: rotation, saveH

  REAL ( dp ), DIMENSION(3) :: xyzstart, xyzcoord
  INTEGER :: i, j, i_basis, nr_atomtypes, type_index, hydro_count
  CHARACTER (LEN = 3), DIMENSION(100) :: atomtypes
  INTEGER, DIMENSION(100) :: atomtype
  INTEGER :: file_num, natoms, n_n
  TYPE (nearest_neighbours), POINTER :: current_near
  REAL (dp) :: teta, phi
  REAL (dp), DIMENSION(3,3) :: Rota
  

  WRITE (*,*) "SAVING GEN FILE"
  ! ------------------------------------------------------------------------
  ! open file and read header
  ! ------------------------------------------------------------------------
  CALL open_file( gen_filename, file_num, "write", format_flag = .TRUE. )

  !------------------------------------------------------------------------
  ! Defines rotation matrix
  !------------------------------------------------------------------------
  Rota(1,1) = sqrt(3.d0)/3.d0*sqrt(2.d0)/2.d0
  Rota(1,2) = sqrt(3.d0)/3.d0*sqrt(2.d0)/2.d0  
  Rota(1,3) = -sqrt(6.d0)/3.d0
  Rota(2,1) = -sqrt(2.d0)/2.d0
  Rota(2,2) = sqrt(2.d0)/2.d0  
  Rota(2,3) = 0.d0
  Rota(3,1) = sqrt(6.d0)/3.d0*sqrt(2.d0)/2.d0
  Rota(3,2) = sqrt(6.d0)/3.d0*sqrt(2.d0)/2.d0  
  Rota(3,3) = sqrt(3.d0)/3.d0  

  !------------------------------------------------------------------------
  !  SAVE .GEN HEADER
  !------------------------------------------------------------------------ 
  natoms = basis%n_basis
  IF ( saveH ) natoms = natoms + basis%n_dg_bond

  IF (basis%periodic) then
     WRITE (file_num, '( I12, A4)') natoms, "   S"
  ELSE
     WRITE (file_num, '( I12, A4)') natoms, "   C"
  ENDIF     

  DO j = 1, size(basis%atomtypes)
     WRITE (file_num,'(a,2x)',advance='NO') basis%atomtypes(j)
  END DO
  WRITE (file_num,*)

  xyzstart = MATMUL(basis%prim,basis%start)/norm(basis%prim)

  hydro_count = 0

  DO i_basis = 1, basis%n_basis

     current_near => neig_map(i_basis)

     xyzcoord = MATMUL( basis%prim, basis%coord( 1:3, i_basis ) )

     if (rotation)  xyzcoord=MATMUL(Rota,xyzcoord)

     IF( saveH ) THEN
        n_n = SIZE(current_near%ind)
     ELSE
        n_n = 0
        DO j = 1, SIZE(current_near%ind)
           IF(current_near%ind(j).LE.basis%n_basis) n_n = n_n+1
        END DO
     END IF

     WRITE (file_num,'(I5, I6, 3F14.6, I6)',ADVANCE = 'NO') &
         basis%mat(i_basis), basis%type(i_basis) , xyzstart + xyzcoord, n_n-1

     DO j = 1, n_n
        
        IF (current_near%ind(j) .NE. i_basis) THEN
           
           IF ( is_dg_bond(j,basis,current_near) ) THEN
              hydro_count = hydro_count + 1
              WRITE (file_num,'(I6)',ADVANCE = 'NO') basis%n_basis+hydro_count
           ELSE
              WRITE (file_num,'(I6)',ADVANCE = 'NO') current_near%ind(j)
           END IF
           
        END IF
        
     END DO
     
     WRITE (file_num, '()', ADVANCE = 'YES')
     
  END DO

  !------------------------------------------------------------------------
  ! SAVE HYDROGEN
  !------------------------------------------------------------------------   
  IF (saveH) THEN
     hydro_count = 0
     DO i_basis = 1, basis%n_basis
        current_near => neig_map(i_basis)
        DO j = 1, SIZE(current_near%ind)
           
           IF (current_near%ind(j) .NE. i_basis) THEN 
              IF ( is_dg_bond(j,basis,current_near) ) THEN
                 hydro_count = hydro_count + 1
                 xyzcoord = MATMUL( basis%prim, basis%coord( 1:3, current_near%ind(j) ) &
                      + current_near%vec(:,j))
                 
                 if (rotation)  xyzcoord=MATMUL(Rota,xyzcoord)
                 
                 WRITE (file_num,'(I5, I6, 3F14.6, I6)',ADVANCE = 'NO') &
                      basis%mat(i_basis), nr_atomtypes , xyzstart + xyzcoord, 1
                 WRITE (file_num,'(I6)',ADVANCE = 'NO') i_basis
                 WRITE (file_num, '()', ADVANCE = 'YES')
              END IF
           END IF
        END DO
     END DO
  END IF
  !------------------------------------------------------------------------
  ! UNIT CELL
  !------------------------------------------------------------------------ 
  WRITE (file_num, '( 3F14.6 )') xyzstart
  WRITE (file_num, '( 3F14.6 )') basis%prim(:,1)
  WRITE (file_num, '( 3F14.6 )') basis%prim(:,2)
  WRITE (file_num, '( 3F14.6 )') basis%prim(:,3)

  CLOSE(file_num)

  WRITE (*,*) "GEN FILE SAVED"

END SUBROUTINE write_gen
!=================================================================================
subroutine parse_symbols(symbols, numsym, atomsym)

  ! INPUT ------------------------------------------------------------------------
  character(100) :: symbols

  ! OUTPUT -----------------------------------------------------------------------
  integer :: numsym
  character(TYPELEN), DIMENSION(:), POINTER :: atomsym

  ! Locals------------------------------------------------------------------------
  character(100) :: tmp_symbol
  character(TYPELEN) :: tmp
  character(1) :: ch1, ch2
  integer :: length_symbols, curr, letter,i

  ! Trim aout the trailing spaces
  symbols=adjustl(symbols)
  length_symbols=len_trim(symbols)

  numsym=count_fields(symbols)

  if(numsym.gt.100) then
     write(*,*) 'Error in gen file: too many symbols'
     call throw_init_exception(ERR_GENERAL)
  endif
  ALLOCATE( atomsym(numsym))

  read(symbols,*) atomsym(1:numsym)

  ! Adjust the type string and add 1 space to 1-ch species P, N, ....
  do i=1,numsym 
    do letter = 2, TYPELEN
      tmp = atomsym(i)
      ch1 = tmp(letter-1:letter-1)
      ch2 = tmp(letter:letter)
      if ( is_upper(ch1) .and. is_upper(ch2) ) then
        tmp(letter+1:TYPELEN) = tmp(letter:TYPELEN-1)
        tmp(letter:letter) = " "
      end if  
      atomsym(i) = tmp
    end do  
  end do

end subroutine parse_symbols
!=================================================================================

integer function count_fields(string)

  character(*) :: string
  integer :: j
  
  count_fields = 0

  j=1
  DO WHILE(j .LE. LEN(string))

     IF ( string(j:j) .EQ. " " ) THEN

        DO WHILE( string(j:j) .EQ. " " )
           j = j + 1
           IF ( j .GT. LEN(string) ) exit

        END DO

     ELSE

        count_fields = count_fields + 1 
        DO WHILE( string(j:j) .NE. " " )
           j = j + 1
           IF ( j .GT. LEN(string) ) exit

        END DO        

     END IF

  ENDDO
  
end function count_fields


logical function is_upper(ch)
  character(1) :: ch
  is_upper = .false.
  if (iachar(ch).ge.65 .and. iachar(ch).le.90) is_upper = .true.

end function is_upper



!=================================================================================

SUBROUTINE read_gen(file_num, n_atoms, nn_list, mat, atomtype, coord)

  INTEGER :: file_num,n_atoms 
  INTEGER, DIMENSION(:,:), INTENT( INOUT ) :: nn_list 
  INTEGER, DIMENSION(:), POINTER :: atomtype 
  INTEGER, DIMENSION(:), POINTER :: mat
  REAL ( dp ), DIMENSION(:,:), POINTER :: coord
  CHARACTER(500) :: line

  ! Local Variables --------------------------------------------------------------
  !
  INTEGER ::  nr_nn, i_basis, j


  DO i_basis = 1, n_atoms

     READ(file_num,'(a)') line

     if( count_fields(line) .eq. 5) THEN 
        
       READ(line,*)  mat(i_basis), atomtype(i_basis), coord(1:3,i_basis)

     else 

       READ(line,*)  mat(i_basis), atomtype(i_basis), coord(1:3,i_basis), &
                     nr_nn,  nn_list(1:nr_nn,i_basis) 

    end if

  END DO

  
END SUBROUTINE read_gen

! ========================================================================
SUBROUTINE reorder(str)
  TYPE( TStructure) :: str

  INTEGER, DIMENSION(:), ALLOCATABLE :: adj, xadj, perm, iperm
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: nn_list
  INTEGER, DIMENSION(:), ALLOCATABLE :: mat 
  INTEGER, DIMENSION(:), ALLOCATABLE :: atomtype
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: coord 

  INTEGER :: i, k, n_atoms, n_adj, err

  n_atoms = size(str%nn_list,2)
  n_adj = count(str%nn_list.ne.0)

  allocate(xadj(n_atoms+1), stat=err)
  allocate(adj(n_adj), stat=err)
  allocate(perm(n_atoms), stat=err)
  allocate(iperm(n_atoms), stat=err)
  ALLOCATE( coord(3, str%n_atoms ), STAT = err )
  ALLOCATE( nn_list(MAX_NUM_NEIGHBOURS, str%n_atoms), STAT = err )
  ALLOCATE( mat(str%n_atoms), STAT = err )
  ALLOCATE( atomtype(str%n_atoms), STAT = err )
  if (err.ne.0) then 
     call throw_init_exception(ERR_GENERAL)
  endif

  print *, 'filling adj'
  n_adj = 1
  xadj(1) = 1
  do i = 1, n_atoms
     do k = 1, MAX_NUM_NEIGHBOURS
         if (str%nn_list(k,i) .eq. 0) exit
         adj(n_adj) = str%nn_list(k,i)
         n_adj = n_adj + 1
     enddo
     xadj(i+1) = n_adj
  end do
  n_adj = n_adj - 1
   
  perm = 0
  
  print *, 'rcm'
  call genrcm(n_atoms, n_adj, xadj, adj, perm)

  do i = 1, n_atoms
    iperm(perm(i)) = i
  enddo  

  print *, 'reordering structure'
  do i = 1, n_atoms
     coord(:,i) = str%coord(:,perm(i))   
     mat(i) = str%mat(perm(i))
     atomtype(i) = str%atomtype(perm(i))
     do k = 1, MAX_NUM_NEIGHBOURS
        if (str%nn_list(k,perm(i)).eq.0) exit
        nn_list(k,i) = iperm(str%nn_list(k,perm(i)))
     enddo   
  enddo 

  str%coord = coord
  str%nn_list = nn_list   
  str%mat = mat
  str%atomtype = atomtype

  deallocate(xadj,adj,perm,iperm,coord,nn_list,mat,atomtype) 

END SUBROUTINE reorder





! ======================================================================================
!  This subroutine reads material information
!  
! id  name struct  type  nr_p  mixing  comp(1)  fraction  database    Ev     gap1 gap2 
! 1   GaN    wz   binary   1   CRY     GaN       1.000     GaN.etb  -0.71033  0.0  0.0
!
!id name  str type    nr_p mix comp1 comp2 x    1-x                    Ev  Ev    gap1 gap2 
! 1 InGaN wz  ternary  2   VCA  GaN   InN  0.9  0.1 GaN.etb InN.etb  -0.71 -0.23  0.0  0.0! 
!
! 
!
!
! ====================================================================================== 
SUBROUTINE read_materials(file_num, nr_materials, materials)

  INTEGER :: file_num, nr_materials
  TYPE( material_data ), DIMENSION(:), POINTER :: materials


  ! Local Variables --------------------------------------------------------------
  ! 
  CHARACTER(500) :: line
  CHARACTER(10) :: name, type, cry, vca
  INTEGER :: nr_parents
  INTEGER :: j, k, mat, err, nargs, nex
 

  READ(file_num,'(a)') line
  nargs = tokenize(line)
  READ(line,*) name, nr_materials
  
  DO j = 1, nr_materials

     READ(file_num,'(a)') line

     nargs = tokenize(line)

     ! Read line up to vca
     READ(line,*) mat, name, cry, type, nr_parents, vca
     
     ALLOCATE(materials(mat)%nr_parents(1), STAT= err)
     IF (err.ne.0) call throw_init_exception(ERR_ALLOC_ERR)

     ! assign to appropriate field
     materials(mat)%nr_parents(1) = nr_parents
     materials(mat)%name = trim(name)
     materials(mat)%cry = trim(cry)
     materials(mat)%type = trim(type)
     materials(mat)%extra_cpl = .false.
        
     SELECT CASE(trim(vca))
     !+++++ PURE CRYSTAL ++++++++++++++++++++++++++++++++++++++++++++++
     CASE('CRY')
        materials(mat)%alloy_random=.false.
        materials(mat)%alloy       =.false.
        materials(mat)%average_cat =.false.
      
        ALLOCATE(materials(mat)%parent(nr_parents), STAT= err)
        ALLOCATE(materials(mat)%bowing(2), STAT= err)
        IF (err.ne.0) call throw_init_exception(ERR_ALLOC_ERR)

        materials(mat)%parent(1:nr_parents)%type=trim(type)
        materials(mat)%parent(1:nr_parents)%e_v = 0.d0
        materials(mat)%bowing=0.d0
        materials(mat)%gap_corr=0.d0

        if (nargs .eq. 6+3*nr_parents+2) then
           READ(line,*) mat, name, cry, type, nr_parents, vca, &
                materials(mat)%parent(1:nr_parents)%name, &
                materials(mat)%parent(1:nr_parents)%content, & 
                materials(mat)%parent(1:nr_parents)%data_file, &
                materials(mat)%gap_corr(1:2)
        
        elseif (nargs .eq. 6+4*nr_parents+2) then
           READ(line,*) mat, name, cry, type, nr_parents, vca, &
                materials(mat)%parent(1:nr_parents)%name, &
                materials(mat)%parent(1:nr_parents)%content, & 
                materials(mat)%parent(1:nr_parents)%data_file, &
                materials(mat)%parent(1:nr_parents)%e_v, &
                materials(mat)%gap_corr(1:2)
        else                 
              write(*,*) 'ERROR: material line has too many fields', &
                        & nargs,'!=',6+4*nr_parents+2
              call throw_init_exception(ERR_INPUT)              
        end if      
     
     !+++++ ALLOY VCA ++++++++++++++++++++++++++++++++++++++++++++++
     CASE('VCA')

        materials(mat)%alloy       =.true.        
        materials(mat)%alloy_random=.false.
        materials(mat)%average_cat =.true.
        
        if (trim(type).eq.'binary') then
           materials(mat)%extra_cpl = .true.
           nex = 1
        else  
           materials(mat)%extra_cpl = .false.
           nex = 0
        end if
        ALLOCATE(materials(mat)%parent(nr_parents+nex), STAT= err)
        ALLOCATE(materials(mat)%bowing(nr_parents), STAT= err)                
        IF (err.ne.0) call throw_init_exception(ERR_ALLOC_ERR)
        
        materials(mat)%parent(1:nr_parents)%type=trim(type)
        materials(mat)%parent(1:nr_parents)%e_v = 0.d0
        materials(mat)%bowing=0.d0
        materials(mat)%gap_corr=0.d0

        if (nargs .eq. 6+4*nr_parents+2+nex) then
           READ(line,*) mat, name, cry, type, nr_parents, vca, &
                materials(mat)%parent(1:nr_parents)%name, &
                materials(mat)%parent(1:nr_parents)%content, &
                materials(mat)%parent(1:nr_parents+nex)%data_file, &
                materials(mat)%bowing(1:nr_parents), &
                materials(mat)%gap_corr(1:2)

        elseif (nargs .eq. 6+5*nr_parents+2+nex) then
           READ(line,*) mat, name, cry, type, nr_parents, vca, &
                materials(mat)%parent(1:nr_parents)%name, &
                materials(mat)%parent(1:nr_parents)%content, &
                materials(mat)%parent(1:nr_parents+nex)%data_file, &
                materials(mat)%parent(1:nr_parents)%e_v, &
                materials(mat)%bowing(1:nr_parents), &
                materials(mat)%gap_corr(1:2)
        
        else                 
              write(*,*) 'ERROR: material line has too many fields', &
                        & nargs,'>',6+4*nr_parents+2+nex
              call throw_init_exception(ERR_INPUT)              
        end if      
     

     !+++++ ALLOY RANDOM  ++++++++++++++++++++++++++++++++++++++++++++++
     CASE('RND')

        materials(mat)%alloy       =.true.       
        materials(mat)%alloy_random=.true.
        materials(mat)%average_cat = .false.
        
        if (trim(type).eq.'binary') then
           materials(mat)%extra_cpl = .true.
           nex = 1
        else  
           materials(mat)%extra_cpl = .false.
           nex = 0
        end if
        ALLOCATE(materials(mat)%parent(nr_parents+nex), STAT= err)
        ALLOCATE(materials(mat)%bowing(nr_parents), STAT= err)
        IF (err.ne.0) call throw_init_exception(ERR_ALLOC_ERR)

        materials(mat)%parent(1:nr_parents)%type=trim(type)
        materials(mat)%parent(1:nr_parents)%content = 1.d0
        materials(mat)%parent(1:nr_parents)%e_v = 0.d0
        materials(mat)%bowing(1:nr_parents)=0.d0
        materials(mat)%gap_corr=0.d0

        if (nargs .eq. 6+4*nr_parents+2+nex) then
           READ(line,*) mat, name, cry, type, nr_parents, vca, &
                materials(mat)%parent(1:nr_parents)%name, &
                materials(mat)%parent(1:nr_parents)%content, &
                materials(mat)%parent(1:nr_parents+nex)%data_file, &
                materials(mat)%bowing(1:nr_parents), &
                materials(mat)%gap_corr(1:2)
 
        elseif (nargs .eq. 6+5*nr_parents+2+nex) then
           READ(line,*) mat, name, cry, type, nr_parents, vca, &
                materials(mat)%parent(1:nr_parents)%name, &
                materials(mat)%parent(1:nr_parents)%content, &
                materials(mat)%parent(1:nr_parents+nex)%data_file, &
                materials(mat)%parent(1:nr_parents)%e_v, &
                materials(mat)%bowing(1:nr_parents), &
                materials(mat)%gap_corr(1:2)
     
        else 
              write(*,*) 'ERROR: material line has too many fields', &
                        & tokenize(line),'>',6+5*nr_parents+2+nex
              call throw_init_exception(ERR_INPUT)              
        
        end if
              
     CASE('default')

        write(*,*) 'ERROR: unrecognized material type',trim(vca)
        call throw_init_exception(ERR_HAM_UNMAT)

     END SELECT
  
   END DO


END SUBROUTINE read_materials


! ======================================================================================
SUBROUTINE read_interfaces(file_num, nr_interfaces, interfaces)

  INTEGER, intent(in) :: file_num
  INTEGER, intent(out) :: nr_interfaces
  TYPE( material_data ), DIMENSION(:,:), POINTER :: interfaces
  
  ! Local Variables --------------------------------------------------------------
  ! 
  CHARACTER(500) :: line
  CHARACTER(10) :: name, type, cry, vca
  INTEGER :: nr_par1, nr_par2, nr_parents
  INTEGER :: j, k, mat1, mat2, err, nargs, nex
 
  READ(file_num,'(a)') line
  nr_interfaces = 0

  DO
     
     READ(file_num,'(a)', end=110) line

     nargs = tokenize(line)
     
     if (nargs == 0) exit


     READ(line,*) mat1, mat2, vca, nr_par1, nr_par2
       
     nr_parents = nr_par1 + nr_par2

     if (nargs .lt. 3*nr_parents + 4) then 
       write(*,*) 'ERROR: interface line has too many fields', &
                        & nargs,' != ',4+3*nr_parents
       call throw_init_exception(ERR_INPUT)              
     endif

     ALLOCATE(interfaces(mat1,mat2)%nr_parents(2), STAT= err)
     ALLOCATE(interfaces(mat1,mat2)%parent(nr_parents), STAT= err)
     IF (err.ne.0) call throw_init_exception(ERR_ALLOC_ERR)
     ALLOCATE(interfaces(mat1,mat2)%bowing(nr_parents), STAT= err)
     IF (err.ne.0) call throw_init_exception(ERR_ALLOC_ERR)

     interfaces(mat1,mat2)%name = "none"
     interfaces(mat1,mat2)%cry = "none"
     interfaces(mat1,mat2)%type = "none"
     interfaces(mat1,mat2)%extra_cpl = .false.
     interfaces(mat1,mat2)%bowing=0.d0
     interfaces(mat1,mat2)%gap_corr=0.d0

     READ(line,*) mat1, mat2, vca, &
          interfaces(mat1,mat2)%nr_parents, &
          interfaces(mat1,mat2)%parent(1:nr_parents)%name, &
          interfaces(mat1,mat2)%parent(1:nr_parents)%content, &
          interfaces(mat1,mat2)%parent(1:nr_parents)%data_file, &
          interfaces(mat1,mat2)%parent(1:nr_parents)%e_v
     
     if (trim(vca)=="VCA") then
        interfaces(mat1,mat2)%alloy        = .true.       
        interfaces(mat1,mat2)%alloy_random = .true.
        interfaces(mat1,mat2)%average_cat  = .true.
     elseif (trim(vca)=="RND") then
        interfaces(mat1,mat2)%alloy        = .true.       
        interfaces(mat1,mat2)%alloy_random = .true.
        interfaces(mat1,mat2)%average_cat  = .false.
     elseif (trim(vca)=="CRY") then
        interfaces(mat1,mat2)%alloy        = .false.       
        interfaces(mat1,mat2)%alloy_random = .false.
        interfaces(mat1,mat2)%average_cat  = .false. 
     end if

     nr_interfaces = nr_interfaces + 1

  END DO   

110 continue

END SUBROUTINE read_interfaces

! ======================================================================================
SUBROUTINE init_structure(verbose, str, materials, nr_mat, interfaces, nr_int)
  INTEGER :: verbose
  TYPE( TStructure) :: str
  TYPE( material_data ), DIMENSION(:), POINTER :: materials
  TYPE( material_data ), DIMENSION(:,:), POINTER :: interfaces
  INTEGER, INTENT(OUT) :: nr_mat, nr_int

  INTEGER :: err, i, j, k, offset, file_num
  LOGICAL :: out_flag  
  character(100) :: symbols
  
  ! ------------------------------------------------------------------------
  ! open file and read header
  ! ------------------------------------------------------------------------
  out_flag=.true.
  if (verbose .eq. 0) out_flag=.false.

  CALL open_file( str%gen_filename, file_num, "read", &
                            format_flag = .TRUE., output_flag = out_flag )
  
  READ (file_num,*) str%n_atoms, str%C_S
  READ (file_num,'(a)') symbols

  ! ------------------------------------------------------------------------
  ! Allocate work space
  ! ------------------------------------------------------------------------
  ALLOCATE( str%coord(3, str%n_atoms ), STAT = err )
  ALLOCATE( str%nn_list(24, str%n_atoms), STAT = err )
  ALLOCATE( str%mat(str%n_atoms), STAT = err )
  ALLOCATE( str%atomtype(str%n_atoms), STAT = err )
  ALLOCATE( str%indexa(str%n_atoms), STAT = err )
  ALLOCATE( str%inv_indexa(str%n_atoms), STAT = err )  
  IF(err.NE.0) CALL alloc_error( 'struct_bulding', 'make_struct', 'work' )

  ! ------------------------------------------------------------------------
  ! parse the symbol line
  ! ------------------------------------------------------------------------
  call parse_symbols(symbols, str%nr_atomtypes, str%atomtypes)

  !IF (verbose.gt.0) THEN
  !   WRITE(*,*) '(read gen) found',str%nr_atomtypes,'symbol(s):'
  !   WRITE(*,'(100(x,A2," "))') str%atomtypes(1:str%nr_atomtypes)
  !END IF

  ! ------------------------------------------------------------------------
  ! Find type nr for hydrogen if non then hydro_type = 0
  ! ------------------------------------------------------------------------
  str%hydro_type = 0
  DO i = 1,str%nr_atomtypes
     IF (index(str%atomtypes(i),trim(hydrosat)) .GT. 0) str%hydro_type = i
  END DO

  str%nn_list = 0

  ! ------------------------------------------------------------------------
  ! Read coordinates 
  ! ------------------------------------------------------------------------
  CALL read_gen(file_num, str%n_atoms, str%nn_list, &
                          str%mat, str%atomtype, str%coord)

  IF (index(str%C_S,"S").gt.0) THEN                      
     READ (file_num,*) str%start(1:3)
     READ (file_num,*) str%prim(:,1)
     READ (file_num,*) str%prim(:,2)
     READ (file_num,*) str%prim(:,3)
  END IF 
  ! ------------------------------------------------------------------------
  ! Count dangling bonds (hydrogens) 
  ! ------------------------------------------------------------------------
  str%n_danglings = 0
  DO i = 1, str%n_atoms 
     IF (str%atomtype(i) .EQ. str%hydro_type) &
                     str%n_danglings = str%n_danglings + 1  
  END DO

  !-------------------------------------------------------------------------
  ! Parse Material lines, crystals and parent materials
  !-------------------------------------------------------------------------
  nr_mat= maxval(str%mat)

  ALLOCATE(materials(nr_mat), STAT= err)
  IF(err.NE.0) CALL alloc_error('struct_building','make_struct','materials')   

  ALLOCATE(interfaces(nr_mat, nr_mat), STAT= err)
  IF(err.NE.0) CALL alloc_error('struct_building','make_struct','interfaces')

  IF (verbose.gt.0) WRITE (*,*) "(make str) Read-in",nr_mat,"materials information"
  call read_materials(file_num,nr_mat,materials)
  call read_interfaces(file_num, nr_int, interfaces)
  IF (verbose.gt.0) WRITE (*,*) "(make str) Read-in ",nr_int," interface information"
  IF (verbose.gt.0) WRITE (*,*) "(make str) gen file read"
  
  close(file_num)
  
  if (verbose.gt.0) WRITE (*,*) "(read gen) structure read"

  ! REORDER THE ATOMIC STRUCTURE
  ! need to fix pot_atoms 
  !CALL reorder(str)

  ! ------------------------------------------------------------------------
  ! Init indexa and inv_indexa used in refine_neighbour_list  
  ! Move Hydrogens to the end
  ! indexa(i) maps the original position of atom i to its ordered position:
  ! e.g., original position 1, ordered position 235
  ! => indexa(1) = 235
  ! inv_index(i) does the opposite
  ! => inv_indexa(235) = 1
  ! ------------------------------------------------------------------------
  offset = 0
  k = str%n_atoms - str%n_danglings
  DO i = 1, str%n_atoms
     IF (str%atomtype(i) .EQ. str%hydro_type) THEN
        offset = offset + 1
        j = k + offset
        str%indexa(i) = j
        str%inv_indexa(j) = i
        if( j.gt.str%n_atoms) then
           write(*,*) i,j,k,offset
           call throw_init_exception(ERR_NN_LIST)
        end if
     ELSE   
        j = i - offset
        str%indexa(i) = j
        str%inv_indexa(j) = i
        if( j.gt.str%n_atoms) then
           write(*,*) i,j,k,offset 
           call throw_init_exception(ERR_NN_LIST)
        endif
     ENDIF
  END DO

END SUBROUTINE init_structure

! ======================================================================================
SUBROUTINE make_basis(verbose, str, basis)

  INTEGER :: verbose
  TYPE( TStructure) :: str
  TYPE( ion_basis ) :: basis
  ! Local Variables --------------------------------------------------------------
  !
  REAL ( dp ), DIMENSION(3,3) :: inv_prim
  REAL ( dp ), DIMENSION( 3 ) :: matrix_start
  INTEGER :: i, j, k, l, p_bound, offset, file_num
  INTEGER :: i_basis, i_mat, i_ion, err, tmp
  INTEGER :: INFO, IPIV(3), WORK(5)
  INTEGER, PARAMETER :: LWORK=5 

  ! ------------------------------------------------------------------------
  ! Allocate arrays in basis type
  ! ------------------------------------------------------------------------
  call create_basis(basis, str%n_atoms)

  ALLOCATE( basis%atomtypes( size(str%atomtypes) ), STAT = err )
  
  IF ( err .NE. 0 ) call throw_init_exception(ERR_ALLOC_ERR) 

  basis%atomtypes = str%atomtypes
  ! ------------------------------------------------------------------------
  ! Resets  "basis%n_basis" 
  ! to the actual number of atoms, not counting dg_bonds
  ! ------------------------------------------------------------------------
  basis%n_dg_bond = str%n_danglings 
  basis%n_basis = str%n_atoms - basis%n_dg_bond

  !write(*,*)  'hydro type =',str%hydro_type
  !write(*,*)  'n. atoms   =',basis%n_basis 
  !write(*,*)  'n. dg bonds=',basis%n_dg_bond
  !write(*,*)  'n. total   =',basis%n_basis + basis%n_dg_bond

  ! ------------------------------------------------------------------------
  ! Set periodic cell
  ! ------------------------------------------------------------------------
  basis%start(1:3)=0.d0
  basis%prim(:,1)=(/1.d0, 0.d0, 0.d0/)
  basis%prim(:,2)=(/0.d0, 1.d0, 0.d0/)
  basis%prim(:,3)=(/0.d0, 0.d0, 1.d0/)

  IF (INDEX(str%C_S, "S") .GT. 0) THEN 
     basis%start(1:3) = str%start(1:3)
     basis%prim(1:3,1:3) = str%prim(1:3,1:3)
     basis%periodic = .true. 
     basis%periodic_BC = "xyz"
  ELSE
     basis%periodic = .false. 
     basis%periodic_BC = "no"
  END IF

  !----------------------------------------------------------------------
  ! Transform coordinates in the lattice basis (fractional coordinates)
  !----------------------------------------------------------------------
  inv_prim = basis%prim

  CALL DGETRF( 3, 3, inv_prim, 3, IPIV, INFO )
  CALL DGETRI( 3, inv_prim, 3, IPIV, WORK, LWORK, INFO )
  
  basis%start = MATMUL(inv_prim,basis%start)*norm(basis%prim)

  basis%rec_latt = 2*pi*transpose(inv_prim)
    

  !============================================================
  ! Initialize the type "basis"
  ! Move Hydrogens to the end
  !============================================================
  offset = 0
  k = basis%n_basis

  DO i = 1, basis%n_basis + basis%n_dg_bond

     IF (str%atomtype(i) .EQ. str%hydro_type) THEN

        offset = offset + 1
        j = k + offset

        if( j.gt.str%n_atoms) then
           write(*,*) i,j,k,offset 
           call throw_init_exception(ERR_NN_LIST)
        endif

        basis%coord(:,j) = MATMUL(inv_prim,str%coord(:,i))
        basis%type(j) = str%atomtype(i)
        basis%mat(j) = str%mat(i)

     ELSE

        j = i - offset

        if( j.gt.str%n_atoms) then
           write(*,*) i,j,k,offset 
           call throw_init_exception(ERR_NN_LIST)
        endif
        
        basis%coord(:,j) = MATMUL(inv_prim,str%coord(:,i))
        basis%type(j) = str%atomtype(i)
        basis%mat(j) = str%mat(i)

     END IF

  END DO
  
  
END SUBROUTINE make_basis

! ==========================================================================
! Substitute H dangling ions with appropriate ion
! Assumes material is crystal binary, binary alloy or ternary (2-ions basis)
! 
! Must be called after make_neighbours_list and init_mat_ion !!
! ==========================================================================
SUBROUTINE subs_dg_ions( basis, materials, neig_map )

  !=========================================================================
  
  ! input arguments :
  TYPE(ion_basis) :: basis
  TYPE(material_data), DIMENSION(:), POINTER :: materials
  TYPE(nearest_neighbours), DIMENSION(:), POINTER :: neig_map
  
  !_________________________________________________________________________

  TYPE (nearest_neighbours), POINTER :: current_near
  INTEGER :: ia, ib, i_ion, i_n, i_mat 

  DO ia = 1, basis%n_basis
  
     current_near => neig_map(ia) 
   
     DO i_n = 1, SIZE( current_near%ind )
          
        ib = current_near%ind(i_n)

        IF( INDEX( basis%atomtypes(basis%type(ib)),trim(hydrosat) ) .GT. 0 ) THEN
           
           call check_ion(basis,ia,materials,i_mat,i_ion)
      
           SELECT CASE (materials(i_mat)%type)
     
           CASE('simple')

              basis%type(ib)= pos( basis%atomtypes, materials(i_mat)%ion(1)%name )

              basis%n_st(ib)= 1
              
           CASE('binary','ternary')

              IF(i_ion.EQ.1) THEN
                 basis%type(ib)= pos( basis%atomtypes, materials(i_mat)%ion(2)%name )
              ELSE
                 basis%type(ib)= pos( basis%atomtypes, materials(i_mat)%ion(1)%name )
              ENDIF
              basis%n_st(ib)= 1
              
           END SELECT

        END IF
     
     END DO
     
  END DO
  
END SUBROUTINE subs_dg_ions

!=========================================================================
SUBROUTINE init_basis(basis, materials)
  ! input arguments :
  TYPE(ion_basis) :: basis
  TYPE(material_data), DIMENSION(:), POINTER :: materials

  INTEGER :: i_a, i_mat, i_ion, err

  do i_a = 1, basis%n_basis

     call check_ion(basis,i_a,materials,i_mat,i_ion)     
     basis%valence(i_a) = materials( i_mat )%ion( i_ion )%valence
     basis%ion(i_a) = i_ion     
  end do

  do i_a = basis%n_basis + 1, basis%n_basis + basis%n_dg_bond

     basis%valence(i_a) = 1

  end do


END SUBROUTINE init_basis
!=========================================================================
! UNFINISHED SUBROUTINE TO CONVERT ETB TO DFTB SK FILES:
SUBROUTINE write_sk_files(mat)

  TYPE(material_data), DIMENSION(:), POINTER :: mat
  
  INTEGER :: i, k, l, nr_materials, a, b
  TYPE(ion_orbit), DIMENSION(:), POINTER :: ion
  TYPE(pair_coupling), DIMENSION(:), POINTER :: pair
  CHARACTER(100) :: file_name
  CHARACTER(1) :: imat

  nr_materials = SIZE(mat)
  
  do i = 1, nr_materials
  
     write(imat,'(i1.1)') i

     do k=1, sum(mat(i)%nr_parents)

        ion => mat(i)%parent(k)%ion
        pair => mat(i)%parent(k)%pair

        ! loop on ions for each pair
        file_name = TRIM(ion(a)%name)//imat//"_" &
                   //TRIM(ion(b)%name)//imat//".sk"

        do l=1,SIZE(ion)

           !ion(i)%name

        end do 



     end do

  end do

END SUBROUTINE write_sk_files

!=========================================================================

 ! tockenize a string into space-separated fields
 ! return the number of fields
 function tokenize(str) result(nf)
    character(256) :: str
    integer :: nf
 
    integer :: i, pos1, pos2
    !character(15) :: word(20)
 
    nf = 0
    pos1 = 1
    !word = ""

    do i = 1, 256

       pos2 = index(trim(str(pos1:)),' ')

       if (pos2 == 0) then
         nf = nf + 1
         !word(nf) = str(pos1:)
         exit
       endif
       if (pos2 == 1) then
          pos1 = pos1 + 1
          cycle
       endif
       nf = nf + 1
       !word(nf) = trim(str(pos1:pos1+pos2-2))
       pos1 = pos2+pos1
    enddo

  end function tokenize

  integer function pos(vec, str)
    character(TYPELEN), DIMENSION(:) :: vec
    character(TYPELEN) :: str

    integer i

    do i =1, size(vec)
       if (vec(i) == str) then 
         pos=i
         exit
       end if
    end do

  end function pos

END MODULE struct_building
