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
module savemofile

use constants
use globals, only : MST, LST
use precision
use upt_param
use type_defs
use input_output
use checks
use sort
use errors

implicit none
private

public :: writemofile, write_eigenstates, append_eigenstate, write_eigenvalues
public :: read_eigenstates, read_eigenvalues, read_old_eigenstates
public :: write_cube, write_jvxl

private :: cspin

  integer, parameter :: BASE = 35
  integer, parameter :: RANGE = 90
  integer, dimension(0:12), parameter :: pow2 = & 
   (/ 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096 /)
  !*  0  1  2  3  4   5   6   7    8    9    10    11    12      
  integer, dimension(0:23), parameter :: edgeVerteces = &
   (/  0, 1, 1, 2, 2, 3, 3, 0, 4, 5, &
  !* 0     1     2     3     4  
     5, 6, 6, 7, 7, 4, 0, 4, 1, 5, 2, 6, 3, 7 /)
  !* 5     6     7     8     9     10    11 
  integer, dimension(0:11), parameter :: edgeTypeTab = &
       (/ 1, 3, 1, 3, 1, 3, 1, 3, 2, 2, 2, 2 /)
  !*********        1=along X, 2=along Y, 3=along Z

  integer, dimension(0:11), parameter :: edgeVertexPlanes = &
   (/ 2, 1, 2, 2, 2, 1, 2, 2, 2, 1, 1, 2 /)
  !** 0  1  2  3  4  5  6  7  8  9  10 11 **
  ! For High to Low scan: edges 1, 5, 9, and 10 are from plane 1(0)

  integer, dimension(0:11), parameter :: edgeVertexPointers = &
   (/ 0, 1, 3, 0, 4, 5, 7, 4, 0, 1, 2, 3 /) 
  ! For High to Low scan

  integer, dimension(0:255), parameter :: insideMaskTable = &
   (/ 0, 265, 515, 778,1030,1295,1541,1804,2060,2309,2575,2822,3082,3331,3593,3840, 400, 153, 915, 666,1430,1183, &
   1941,1692,2460,2197,2975,2710,3482,3219,3993,3728, 560, 825,  51, 314,1590,1855,1077,1340,2620,2869,2111,2358, &
   3642,3891,3129,3376, 928, 681, 419, 170,1958,1711,1445,1196,2988,2725,2479,2214,4010,3747,3497,3232,1120,1385, &
   1635,1898, 102, 367, 613, 876,3180,3429,3695,3942,2154,2403,2665,2912,1520,1273,2035,1786, 502, 255,1013, 764, &
   3580,3317,4095,3830,2554,2291,3065,2800,1616,1881,1107,1370, 598, 863,  85, 348,3676,3925,3167,3414,2650,2899, &
   2137,2384,1984,1737,1475,1226, 966, 719, 453, 204,4044,3781,3535,3270,3018,2755,2505,2240,2240,2505,2755,3018, &
   3270,3535,3781,4044, 204, 453, 719, 966,1226,1475,1737,1984,2384,2137,2899,2650,3414,3167,3925,3676, 348,  85, &
   863, 598,1370,1107,1881,1616,2800,3065,2291,2554,3830,4095,3317,3580, 764,1013, 255, 502,1786,2035,1273,1520,  &
   2912,2665,2403,2154,3942,3695,3429,3180, 876, 613, 367, 102,1898,1635,1385,1120,3232,3497,3747,4010,2214,2479, &
   2725,2988,1196,1445,1711,1958, 170, 419, 681, 928,3376,3129,3891,3642,2358,2111,2869,2620,1340,1077,1855,1590, &
   314,  51, 825, 560,3728,3993,3219,3482,2710,2975,2197,2460,1692,1941,1183,1430, 666, 915, 153, 400,3840,3593,  &
   3331,3082,2822,2575,2309,2060,1804,1541,1295,1030, 778, 515, 265,   0 /) 

 
contains

subroutine writemofile(upt)
type(OUPT), TARGET :: upt
character(MST) :: file_name
type(ion_basis), POINTER :: basis
integer, dimension(:), pointer :: n_st, n_dg
type(material_data), dimension(:), pointer :: mat_data
integer :: i,j,k,l,i_mat,i_ion, file_num, n_spin, id
real(dp), dimension(3) :: xyzcoord
real(dp) :: sum

basis => upt%basis
n_st => upt%basis%n_st
n_dg => upt%basis%n_dg
mat_data => upt%materials
n_spin = upt%n_spin
file_name = 'output.mo'


call open_file( file_name, file_num, operation= "write", &
              format_flag=.TRUE., replace_flag=.TRUE. )

write( file_num, '(a8)') '[HEADER]'
write( file_num, '(a9)') 'Version 2'
write( file_num, '(a)') 'Type DENSITY' 
write( file_num, '(a)') 'ValenceOnly FALSE' 
write( file_num, '(a)') 'EnergyUnits EV'
write( file_num, '(a)') ''
write( file_num, '(a)') ''
write( file_num, '(a7)') '[ATOMS]'

do i=1, basis%n_basis
  xyzcoord = matmul(basis%prim,basis%coord(:,i))
  write( file_num,*) basis%atomtypes(basis%type(i)), xyzcoord/0.529771d0
end do

write( file_num, '(a)') ''
write( file_num, '(a5)') '[STO]'

do i=1, basis%n_basis

  call check_ion(basis, i, mat_data, i_mat, i_ion )

!  write( file_num,'(i5,2x,a)') i,'0 0 0 1 1.565085 0.998181645' 
  do j=1, n_st(i)
    select case(TRIM(mat_data(i_mat)%ion(i_ion)%state(j)) )
    case('s')
       write( file_num,'(i5,2x,a)') i,'0 0 0 1 1.565085 0.998181645' 
       
    case('px')
       write( file_num,'(i5,2x,a)') i,'1 0 0 1 1.658972 1.21151527' 
       
    case('py')
       write( file_num,'(i5,2x,a)') i,'0 1 0 1 1.658972 1.21151527' 
       
    case('pz')
       write( file_num,'(i5,2x,a)') i,'0 0 1 1 1.658972 1.3115127' 
       
    case('s*')
       write( file_num,'(i5,2x,a)') i,'0 0 0 2 1.891185 1.10637709' 
       
    case('dz2r2')
       write( file_num,'(i5,2x,a)') i,'-2 0 0 0 1.891185 1.10637709'  
       
    case('dxy')
       write( file_num,'(i5,2x,a)') i,'1 1 0 0 1.891185 1.10637709' 
       
    case('dzx')
       write( file_num,'(i5,2x,a)') i,'1 0 1 0 1.891185 1.10637709' 
          
    case('dx2y2')
       write( file_num,'(i5,2x,a)') i,'0 -2 0 0 1.891185 1.10637709' 
       
    case('dyz')
       write( file_num,'(i5,2x,a)') i,'0 1 1 0 1.891185 1.10637709' 

    end select
          
   end do
end do

write( file_num, '(a)') ''

do i=1, upt%num_vb

  write( file_num, '(a3,i2,a1)') '[MO',i,']'
  write( file_num, *) upt%eigen_values(i)
  write( file_num, *) '2'

  k = 0; id = 0
  do j=1, basis%n_basis   ! upt%ham%nrow - upt%basis%n_dg_bond

     sum = 0
     do l= 1, n_st(j) 
        k = k + 1
        id = id + 1
        ! Add spin up and spin down
        sum = abs(upt%eigen_vectors(k,i)) + &
              abs(upt%eigen_vectors(k+(n_spin-1)*n_st(j),i))
 
        write( file_num, '(i5,f10.6)') id, abs(upt%eigen_vectors(k,i))
     end do

     k = k + n_st(j) !+ n_spin * n_dg(j)
     !write( file_num, '(i5,f10.6)') j, sum
                
  end do

end do


do i=upt%num_vb+1, upt%num_vb+upt%num_cb
  
  write( file_num, '(a3,i2,a1)') '[MO',i,']'
  write( file_num, *) upt%eigen_values(i)
  write( file_num, *) '0'

  k = 0; id = 0
  do j=1, basis%n_basis   ! upt%ham%nrow - upt%basis%n_dg_bond

     sum = 0
     do l= 1, n_st(j)
        k = k + 1
        id = id + 1
        ! Add spin up and spin down
        sum = abs(upt%eigen_vectors(k,i)) + &
              abs(upt%eigen_vectors(k+(n_spin-1)*n_st(j),i))

        write( file_num, '(i5,f10.6)') id, sum !abs(upt%eigen_vectors(k,i))
     end do

     k = k + n_st(j) ! + n_spin * n_dg(j)
     !write( file_num, '(i5,f10.6)') j, sum
                
  end do
  

  !do j=1, upt%ham%nrow - upt%basis%n_dg_bond
  !      write( file_num, '(i5,f10.6)') j, &
  !              abs(upt%eigen_vectors(j,i))*abs(upt%eigen_vectors(j,i))
  !end do

end do

close(file_num)

write(*,*) "(mo) check: ",id, k, size(upt%eigen_vectors,1)


end subroutine writemofile
! -----------------------------------------------------------------------------


subroutine write_eigenvalues(upt)

  type(OUPT), TARGET :: upt
  character(LST) :: file_name 

  type(ion_basis), POINTER :: basis
  integer, dimension(:), pointer :: n_st, n_dg  
  type(material_data), dimension(:), pointer :: mat_data
  integer :: file_num, i_k, i, n, i_mat, i_ion, i_spin, n_spin 
  real(dp) :: phi
  character(2) :: id

  basis => upt%basis
  n_st => upt%basis%n_st
  n_dg => upt%basis%n_dg
  mat_data => upt%materials
  n_spin = upt%n_spin

  file_name = trim(upt%out_path)//'eigv.dat'
  i_k = 1

  CALL open_file( file_name, file_num, operation = "write", &
       format_flag = .TRUE., replace_flag = .TRUE. )
  
  write( file_num, '(a4)', advance='NO' ) '# N  '
  do i =1, upt%num_vb
    write(id,'(i2.2)') i
    write( file_num, '(a4,a2)', advance='NO' ) '  Ev',id
  enddo
  do i =1, upt%num_cb
    write(id,'(i2.2)') i
    write( file_num, '(a4,a2)', advance='NO' ) '  Ec',id
  enddo
  write( file_num, *)
  
  write( file_num, '(i3)', advance='NO' ) i_k

  do i = 1, upt%num_vb + upt%num_cb 
 
     WRITE( file_num, '(( x1 f13.8 ))', advance='NO' ) &
          upt%eigen_values(i)      

  enddo

  write(file_num,*) 
  
  CLOSE ( file_num )

end subroutine write_eigenvalues
  
! -----------------------------------------------------------------------------

subroutine write_eigenstates(upt)

  type(OUPT), TARGET :: upt

  character(LST) :: file_name 
  integer :: file_num, i

  !file_name = trim(upt%out_path)//trim(upt%state_file)
  file_name = trim(upt%state_file)
  
  CALL open_file( file_name, file_num, operation = "write", &
       stream_flag = .TRUE., replace_flag = .TRUE. )

  ! +1 for electrons, -1 for holes, as in lanczos
  do i = 1, upt%num_vb
    WRITE(file_num) -1, upt%eigen_values(i), upt%eigen_vectors(:,i)
  enddo
  do i = upt%num_vb + 1, upt%num_vb + upt%num_cb
    WRITE(file_num) +1, upt%eigen_values(i), upt%eigen_vectors(:,i)
  enddo

  CLOSE(file_num)
  
  !write(*,*) '(check) nrow=',upt%ham%nrow

end subroutine write_eigenstates

! -----------------------------------------------------------------------------

subroutine append_eigenstate(file_name, eigvec, eigval, particle)
  character(*), intent(in) :: file_name
  complex(dp), dimension(:), target :: eigvec
  real(dp), intent(in) :: eigval
  integer, intent(in) :: particle

  integer :: file_num, err

  CALL open_file( file_name, file_num, operation = "write", &
       stream_flag = .TRUE., replace_flag = .FALSE., output_flag = .FALSE. )
       
  !ENDFILE(file_num, IOSTAT=err)

  write(file_num) particle, eigval, eigvec
  

  close(file_num)

end subroutine append_eigenstate

! -----------------------------------------------------------------------------

character(1) function cspin(i_spin)
  integer :: i_spin

  if (i_spin .eq. 1) then 
     cspin = 'a' 
  else
     cspin = 'b'
  end if

end function cspin

! -----------------------------------------------------------------------------

subroutine read_eigenstates(upt, file_name, ev, ec)

  type(OUPT) :: upt
  character(len=*), intent(in) :: file_name
  integer, intent(out) :: ev, ec

  integer :: file_num, iostatus, i, n 
  real(dp) :: eigval
  integer :: particle_type
  complex(dp) :: matel
  complex(dp), dimension(:), allocatable :: tmp_vec

  integer num_ev

  ! data format:
  ! particle|eigenvalue|eigenvector
  ! INTEGER|REAL(dp)|COMPEX(dp),dimension(nrow)
  !

  num_ev = upt%num_vb + upt%num_cb

  if (.not.associated(upt%eigen_vectors)) then
     write(*,*) "(upt) allocating eigenvectors:",upt%ham%nrow,"x",num_ev
     allocate(upt%eigen_vectors(upt%ham%nrow,num_ev))
     upt%eigen_vectors = (0.d0,0.d0)
  end if

  if (.not.associated(upt%eigen_values)) then
     write(*,*) "(upt) allocating eigenvalue array of size",num_ev
     allocate(upt%eigen_values(num_ev))
     upt%eigen_values = (0.d0,0.d0)
  end if

  if (.not.associated(upt%particles)) then
     !write(*,*) "(upt) allocating particles type array of size",num_ev
     allocate(upt%particles(num_ev))
     upt%particles = 0
  end if

  CALL open_file(file_name, file_num, operation = "read", &
       stream_flag = .TRUE. )

  allocate(tmp_vec(upt%ham%nrow))

  ev = 0
  ec = 0

  do 

    read(file_num,iostat=iostatus) particle_type, eigval, tmp_vec

    if (iostatus .gt. 0) then
      write(*,*) 'An error occurred reading eigenvectors from file', file_name
      exit
    else if (iostatus .lt. 0) then
      !write(*,*) 'Read ', ev + ec, ' states'
      exit
    else
      ! calculate eigenvalue, do some sanity check, identify state
      if (particle_type .eq. -1) then
        if (ev .lt. upt%num_vb) then
          ev = ev + 1
          upt%eigen_vectors(:, ev) = tmp_vec
          upt%eigen_values(ev) = eigval
          upt%particles(ev) = -1
        end if
      else if (particle_type .eq. +1) then
        if (ec .lt. upt%num_cb) then
          ec = ec + 1
          upt%eigen_vectors(:, upt%num_vb + ec) = tmp_vec
          upt%eigen_values(upt%num_vb + ec) = eigval
          upt%particles(upt%num_vb + ec) = +1
        end if
      end if
    end if

  enddo

  CLOSE ( file_num )
  
  write(*,*) '(upt) read',ev,' valence eigenvectors'
  write(*,*) '(upt) read',ec,' conduction eigenvectors'

  deallocate(tmp_vec)  

end subroutine read_eigenstates
! -----------------------------------------------------------------------------
! -----------------------------------------------------------------------------
subroutine read_eigenvalues(upt,file_name,ev,ec)

  type(OUPT) :: upt
  character(LST) :: file_name 
  integer, intent(out) :: ev, ec

  character(6) :: header
  integer :: file_num, i_k, i, n 
  real(dp), dimension(:), allocatable :: tmp_vec

  integer num_ev

  CALL open_file( file_name, file_num, operation = "read", &
       format_flag = .TRUE., replace_flag = .FALSE. )
  
  ev = 0
  ec = 0

  READ( file_num, '(a4)', advance='no') header
  do 
     READ( file_num, '(a6)', advance='no', EOR=100 ) header
     if( index(header,"Ev").ne.0) ev=ev+1
     if( index(header,"Ec").ne.0) ec=ec+1     
  enddo

100 continue

  write(*,*) "found: ",ev," valence states",ec," conduction states"

  if (upt%num_vb .lt. ev .or. upt%num_cb .lt. ec) then
     if (associated(upt%eigen_values)) deallocate(upt%eigen_values)
     if (associated(upt%eigen_vectors)) deallocate(upt%eigen_vectors)
     if (associated(upt%particles)) deallocate(upt%particles)
     upt%num_vb = ev
     upt%num_cb = ec
  endif
 
  ! Additional eigenpairs need to be computed: set start point
  upt%start_vb = ev + 1
  upt%start_cb = ec + 1

  num_ev = upt%num_vb + upt%num_cb

  if (.not.associated(upt%eigen_values)) then
     write(*,*) "(upt) allocating eigenvalues:",num_ev
     allocate(upt%eigen_values(num_ev))
     upt%eigen_values = 0.d0
  end if

  if (.not.associated(upt%particles)) then
     write(*,*) "(upt) allocating particles:",num_ev
     allocate(upt%particles(num_ev))
     upt%particles = 0
  end if

  allocate( tmp_vec(num_ev+1) )

  READ( file_num, * ) tmp_vec(1:ev+ec+1)

  ! store eigenvalues (valence first then conduction)
  ! space is left for additional eigenvalues to be computed 
  do i = 1, ev
     upt%eigen_values(i) = tmp_vec(i+1)
  end do
  
  upt%particles(1:ev) = -1
  
  do i = 1, ec 
     upt%eigen_values(upt%num_vb+i) = tmp_vec(ev+i+1)
  end do
  
  upt%particles(upt%num_vb+1:upt%num_vb+ec) =  1  

  CLOSE ( file_num )

  deallocate(tmp_vec) 


end subroutine read_eigenvalues
! -----------------------------------------------------------------------------
! -----------------------------------------------------------------------------

subroutine read_old_eigenstates(upt,file_name)

  type(OUPT) :: upt
  character(LST) :: file_name 

  integer :: file_num, i_k, i, n 
  real(dp) :: phi

  CALL open_file( file_name, file_num, operation = "read", &
       format_flag = .TRUE., replace_flag = .FALSE. )

  READ( file_num, *) 
  
  do i = 1, upt%num_vb + upt%num_cb

     READ( file_num, *) upt%eigen_vectors(1:upt%ham%nrow, i)
     READ( file_num, *) 
  enddo

  CLOSE ( file_num )


end subroutine read_old_eigenstates
! -----------------------------------------------------------------------------

subroutine write_cube(upt)

type(OUPT), TARGET :: upt

character(LST) :: file_name
type(ion_basis), POINTER :: basis

real(dp), dimension(:,:,:), allocatable :: rho
real(dp), dimension(:,:), allocatable :: xyzcoord
real(dp), dimension(:), allocatable :: atomcharge

real(dp) :: cutoff
integer, dimension(3) :: Np
real(dp), dimension(3) :: mincoord
integer :: i, j, k, l, n, natoms, mo, file_num
logical :: is_exist
integer :: n_basis, n_spin
character(2) :: mo_c

real(dp) :: step, loca

basis => upt%basis
n_spin = upt%n_spin
n_basis = basis%n_basis
step = upt%grid_step

allocate( xyzcoord(3, basis%n_basis + basis%n_dg_bond ) )
allocate( atomcharge( basis%n_basis + basis%n_dg_bond ) )
atomcharge=0.d0

call getNp(basis,step,xyzcoord,mincoord,Np)

write(*,*) 'number of points:',Np

allocate( rho(Np(1),Np(2),Np(3)) )
write(*,*) 'Memory required:', SIZE(rho)*8,'bytes'

do mo = 1, upt%num_vb + upt%num_cb

   write(*,*) '(save cube) computing MO =',mo
   
   call getatomcharges(mo,basis,n_spin,upt%eigen_vectors,atomcharge)   

   write(*,*) '(save cube) normalization=',sum(atomcharge)
 
   loca = localization(atomcharge, 0.5d0)

   write(*,*) '(save cube) localization index (50%)=',loca

   call compute_rho(basis,xyzcoord,mincoord,step,Np,atomcharge,rho)

   call find_cutoff(rho,0.10d0,cutoff)
   write(*,*) "(save cube) iso (10%)=", cutoff 
   call find_cutoff(rho,0.50d0,cutoff)
   write(*,*) "(save cube) iso (50%)=", cutoff 
   call find_cutoff(rho,0.80d0,cutoff)
   write(*,*) "(save cube) iso (80%)=", cutoff 


   write(*,*) '(save cube) norm rho=',get_norm(rho)

   if (mo.le.upt%num_vb) then
     write(mo_c,'(i2.2)') mo
     file_name = trim(upt%out_path)//'mo_vb_'//mo_c//'.cube'
   else
     write(mo_c,'(i2.2)') mo-upt%num_vb
     file_name = trim(upt%out_path)//'mo_cb_'//mo_c//'.cube'
   endif

   !Check if the mo files alreasy exists or not. If yes, not rewrite.
   !Thus the mo files are still of the first k-point
   inquire(file=file_name, exist=is_exist)

   if (.not.is_exist) then

     call open_file( file_name, file_num, operation= "write", &
       format_flag=.TRUE., replace_flag=.TRUE. )

     write( file_num, *) 'CUBE FILE MO'//mo_c
     write( file_num, *) 'x, y, z'
     write( file_num, '(i7,3f12.5)') basis%n_basis, mincoord/ATU
     write( file_num, '(i5,3f12.5)') Np(1), step/ATU, 0.d0, 0.d0
     write( file_num, '(i5,3f12.5)') Np(2), 0.d0, step/ATU, 0.d0
     write( file_num, '(i5,3f12.5)') Np(3), 0.d0, 0.d0, step/ATU 

     do n = 1, basis%n_basis
       write( file_num, '(i5,4e13.5)')  atomZ(basis%atomtypes(basis%type(n))), &
                                    atomcharge(n), xyzcoord(:,n)/ATU
     end do
     !write( file_num, '(i5,i5)') 1,1

     do i = 1, Np(1)
       do j = 1, Np(2)
         do k = 1, Np(3)

           write( file_num, '(e13.5)', ADVANCE='NO') rho(i,j,k)

           if( mod(k-1,6).eq.5 ) write(file_num,*)

         end do
         write(file_num,*)
       end do
     end do

     close(file_num)
   
   endif

end do


deallocate(atomcharge,xyzcoord,rho)


end subroutine write_cube

! -----------------------------------------------------------------------------
subroutine getNp(basis,step,xyzcoord,mincoord,Np)

  type(ion_basis), POINTER :: basis
  real(dp) :: step
  real(dp) :: xyzcoord(:,:)
  real(dp), dimension(3) :: mincoord
  integer, dimension(3) :: Np 

  real(dp), dimension(3) :: side, maxcoord
  integer :: n_basis,k
  real(dp) :: r_min, tau

  n_basis = basis%n_basis

  mincoord = (/ 1.d6, 1.d6, 1.d6 /)
  maxcoord = -mincoord 

   ! conversion of Hubbard U into tau (16/5 U_H) and into 1/Ang
  tau = 3.2d0 * 0.24d0 / ATU
  ! compute radius of influence 
  r_min = 7.d0/tau   

  ! Copy atomic coordinates into xyzcoord with the same ordering 
  ! of the orbital basis set. And find bounding box.
  k = 0
  do k = 1, n_basis
     !k = k + 1
     xyzcoord(:,k) = matmul(basis%prim,basis%coord(:,k))
     
     if(xyzcoord(1,k).lt.mincoord(1)) mincoord(1) = xyzcoord(1,k)
     if(xyzcoord(2,k).lt.mincoord(2)) mincoord(2) = xyzcoord(2,k)
     if(xyzcoord(3,k).lt.mincoord(3)) mincoord(3) = xyzcoord(3,k)
     if(xyzcoord(1,k).gt.maxcoord(1)) maxcoord(1) = xyzcoord(1,k)
     if(xyzcoord(2,k).gt.maxcoord(2)) maxcoord(2) = xyzcoord(2,k)
     if(xyzcoord(3,k).gt.maxcoord(3)) maxcoord(3) = xyzcoord(3,k)
     
  enddo

  mincoord = mincoord - r_min
  maxcoord = maxcoord + r_min
  
  ! construct cube file
  side = maxcoord - mincoord 
  Np = nint(side/step)
  !side(:) = Np(:) * step


end subroutine getNp

! -----------------------------------------------------------------------------
subroutine getatomcharges(mo,basis,n_spin,eigenvec,atomcharge)

  integer :: mo
  type(ion_basis), POINTER :: basis
  integer :: n_spin
  complex(dp), dimension(:,:), pointer :: eigenvec
  real(dp) :: atomcharge(*)

  integer :: last, current, n_basis, n
  integer, dimension(:), pointer :: n_st, n_dg

  n_st => basis%n_st
  n_dg => basis%n_dg
  n_basis = basis%n_basis

  last = 0
  do n = 1, n_basis
     ! Compute atomcharges
     current = last + n_spin * n_st( n )
     
     atomcharge(n) = SUM( ABS( eigenvec(last+1:current, mo ) )**2, 1 )
     
     last = current !+ n_spin * n_dg( n )
     
  end do


end subroutine getatomcharges

! -----------------------------------------------------------------------------
subroutine compute_rho(basis,xyzcoord,mincoord,step,Np,atomcharge,rho)
  type(ion_basis), POINTER :: basis
  real(dp) :: xyzcoord(:,:)
  real(dp), dimension(3) :: mincoord
  real(dp) :: step
  integer, dimension(3) :: Np
  real(dp) :: atomcharge(*)
  real(dp) :: rho(Np(1),Np(2),Np(3))
  

  integer :: n_basis
  real(dp), dimension(3) :: xmin, xmax
  integer, dimension(3) ::  imin, imax
  integer :: i,j,k,n
  real(dp) :: r_min, tau, dist, dx, dy, dz
 
  real(dp), parameter :: one_ov_eight_pi = 1/8.d0/pi

  ! conversion of Hubbard U into tau (16/5 U_H) and into 1/Ang
  tau = 3.2d0 * 0.24d0 / ATU
  ! compute radius of influence 
  r_min = 7.d0/tau  

  n_basis = basis%n_basis

  rho = 0.d0
  do n = 1, n_basis
     
     do i = 1,3
        imin(i) = nint ( (xyzcoord(i,n) - r_min - mincoord(i))/step ) + 1
        imax(i) = nint ( (xyzcoord(i,n) + r_min - mincoord(i))/step ) + 1
        imin(i) = max( 1, imin(i) )
        imax(i) = min( Np(i), imax(i) )
     enddo
     
     do i = imin(1), imax(1)
        do j = imin(2), imax(2)
           do k = imin(3), imax(3)
              
              dx = abs( step*(i-1) + mincoord(1) - xyzcoord(1,n) )
              dy = abs( step*(j-1) + mincoord(2) - xyzcoord(2,n) )
              dz = abs( step*(k-1) + mincoord(3) - xyzcoord(3,n) )              
              dist = sqrt( dx**2 + dy**2 + dz**2 )
              
              rho(i,j,k) = rho(i,j,k) + &
                    atomcharge(n) * tau**3 * exp(- tau * dist) * one_ov_eight_pi
                            
           end do
        end do
     end do
     
  end do
  
  rho = rho * step**3

end subroutine compute_rho
! -----------------------------------------------------------------------------

subroutine write_jvxl(upt)

type(OUPT), TARGET :: upt

character(LST) :: file_name
type(ion_basis), POINTER :: basis

real(dp), dimension(3) :: mincoord
real(dp), dimension(3) :: side
real(dp), dimension(:,:,:), allocatable :: rho
real(dp), dimension(:,:), allocatable :: xyzcoord
real(dp), dimension(:), allocatable :: atomcharge
integer, dimension(:,:,:), allocatable :: indexPlanes
character(1), dimension(:), allocatable :: planes
character(1), dimension(:), allocatable :: edges
character(1), dimension(12) :: string
integer, dimension(:), pointer :: n_st, n_dg
real(dp) :: cutoff
integer, dimension(3) :: Np
integer :: i, j, k, m, n, natoms, mo, file_num, last, current
integer :: n_basis, n_spin
character(2) :: mo_c

real(dp) :: step
integer, parameter :: NS = 1 
REAL(dp), PARAMETER :: fractions(5) = 0.3d0 !(/0.1d0, 0.3d0, 0.5d0, 0.7d0, 0.9d0/)
integer, parameter :: EB = BASE
integer, parameter :: ER = RANGE
integer, parameter :: CB = BASE
integer, parameter :: CR = RANGE
integer, parameter :: NCol = -1
integer :: NPlane
integer :: NEdge 
integer :: insideMask,pt,ptX, count, Nyz
integer, dimension(3,8) :: cubeVertexOffset
integer, dimension(3) :: offset
real(dp), dimension(8) :: vertexValues

cubeVertexOffset(:,1)= (/ 0,0,0 /) !0 pt
cubeVertexOffset(:,2)= (/ 1,0,0 /) !1 pt + yz
cubeVertexOffset(:,3)= (/ 1,0,1 /) !2 pt + yz + 1
cubeVertexOffset(:,4)= (/ 0,0,1 /) !3 pt + 1
cubeVertexOffset(:,5)= (/ 0,1,0 /) !4 pt + z
cubeVertexOffset(:,6)= (/ 1,1,0 /) !5 pt + yz + z
cubeVertexOffset(:,7)= (/ 1,1,1 /) !6 pt + yz + z + 1
cubeVertexOffset(:,8)= (/ 0,1,1 /) !7 pt + z  + 1


basis => upt%basis
n_spin = upt%n_spin
n_basis = basis%n_basis
step = upt%grid_step

allocate( xyzcoord(3, basis%n_basis + basis%n_dg_bond ) )
allocate( atomcharge( basis%n_basis + basis%n_dg_bond ) )
atomcharge=0.d0

call getNp(basis,step,xyzcoord,mincoord,Np)
Nyz = Np(2) * Np(3)


allocate( rho(Np(1),Np(2),Np(3)) )
allocate( planes(0:Np(1)*Np(2)*Np(3)) )
allocate( edges(Np(1)*Np(2)*Np(3)) )
allocate( indexPlanes(3,Nyz,2) )
write(*,*) 'Memory required:', SIZE(rho)*8+SIZE(planes)+SIZE(edges),'bytes'

planes = ''
edges = ''

do mo = 1, upt%num_vb + upt%num_cb

   !write(*,*) '(save jvxl) computing MO =',mo
 
   call getatomcharges(mo,basis,n_spin,upt%eigen_vectors,atomcharge)   

   call compute_rho(basis,xyzcoord,mincoord,step,Np,atomcharge,rho)   

   if (mo.le.upt%num_vb) then
     write(mo_c,'(i2.2)') mo
     file_name = trim(upt%out_path)//'mo_vb_'//mo_c//'.jvxl'
   else
     write(mo_c,'(i2.2)') mo-upt%num_vb
     file_name = trim(upt%out_path)//'mo_cb_'//mo_c//'.jvxl'
   endif

   call open_file( file_name, file_num, operation= "write", &
        format_flag=.TRUE., replace_flag=.TRUE. )

   write( file_num, '(a)') '#JVXL VERSION 1.4'
   write( file_num, '(a)') 'CUBE FILE MO'//mo_c
   write( file_num, '(a)') ' x, y, z'
   write( file_num, '(i7,3f12.5)') -basis%n_basis, mincoord/ATU
   write( file_num, '(i5,3f12.5)') Np(1), step/ATU, 0.d0, 0.d0
   write( file_num, '(i5,3f12.5)') Np(2), 0.d0, step/ATU, 0.d0
   write( file_num, '(i5,3f12.5)') Np(3), 0.d0, 0.d0, step/ATU 

   do n = 1, basis%n_basis
   write( file_num, '(i5,4e13.5)')  atomZ(basis%atomtypes(basis%type(n))), &
                                    atomcharge(n), xyzcoord(:,n)/ATU
   end do

   write( file_num, '(5(i4),a)') -NS,EB,ER,CB,CR, " Jmol voxel format version 1.4" 

   do n = 1, NS

      planes = ''
      edges = ''
      planes(0) = '1'

      ! Find cutoff for a given density fraction (60%)
      call find_cutoff(rho,fractions(n),cutoff)
      
      write(*,*) 'iso',nint(fractions(n)*100),'% cutoff=',cutoff
      
      ! walk across the voxels and find cutting surfaces
      ! count the number of cutting surfaces
      count = 0
      Nplane = 0
      do i = 1, Np(1)
         do j = 1, Np(2)
            do k = 1, Np(3)
               
               count = count + 1
               ! Those outside the surface are marked with '1'
               if(rho(i,j,k).gt.cutoff) then
                  planes(count) = '1'
               endif
               
               if (planes(count).ne.planes(count-1)) then
                  Nplane = Nplane + 1
               endif
               
            end do
         end do
      end do
      Nplane = Nplane + 1
      
      !write(*,*) 'Nplane=',Nplane
      
      ! walk backward through the edges with the marching cube algorithm
      insideMask = 0
      indexPlanes = -1 
      NEdge = 0
      ptX = (Np(1))*Nyz - Np(3) - 1   ! ptX starting point
      
      do i = Np(1)-1, 1, -1
         ptX = ptX - Nyz 
         pt  = ptX
         
         !write(*,*) 'Plane=',i,ptX
         
         indexPlanes(:,:,1) = indexPlanes(:,:,2)
         indexPlanes(:,:,2) = -1

         do j = Np(2)-1, 1, -1
            pt = pt - 1
            do k = Np(3)-1, 1, -1
               pt = pt - 1
               
               insideMask = 0
               do m = 8, 1, -1 
                  
                  offset(:) = cubeVertexOffset(:,m)
                  vertexValues(m) = rho(i+offset(1),j+offset(2),k+offset(3))
                  
                  if (vertexValues(m).gt.cutoff) then
                     insideMask = insideMask + pow2(m-1)
                  endif
                  
               end do
               
               ! check if all points are out or all in
               if (insideMask.eq.0 .or. insideMask.eq.255) cycle
               
               ! encode the surface position
               call processOneCubical(insideMask,vertexValues,pt,Nyz,Np(3), &
                    cutoff, indexPlanes,edges,NEdge)
               
               !do m = 1, 12
               !   if(string(m).ne.'') then 
               !      NEdge = NEdge + 1
               !      edges(NEdge) = string(m)
               !   end if
               !enddo
               
            end do
         end do
         
         !write(*,*) 'NEdge=',NEdge
         
      end do
      
      !write(*,*) 'write header. NEdge=', NEdge
      
      ! Write jvxl file format:
      write( file_num, '(a)') "# surface description"      
      write( file_num, '(E11.2,2x,i7,2x,i7,2x,i4,a,a)') cutoff, NPlane, NEdge, NCol, &
           ' 0.0 0.0 0.0 0.0 rendering:isosurface ID', & 
           ' "isosurface1" fill noMesh noDots notFrontOnly frontlit'
      
      Nplane=0
      count = 0
      do i = 1, Np(1)*Np(2)*Np(3)
         count = count + 1
         if (planes(i).ne.planes(i-1)) then
            call write_fmt_int(file_num,count)
            count = 0
            Nplane = Nplane + 1
         endif
      enddo
      Nplane = Nplane + 1
      call write_fmt_int(file_num,count+1)
      write(file_num,*)
      
      write(*,*) 'writing jvxl. Nplane=',Nplane

      do i = 1, NEdge
         write(file_num,'(a1)',advance='NO') edges(i)
      enddo

      write(file_num,*)

   enddo

   write(file_num,*)
   write(file_num,'(a)') '#-------end of jvxl file data-------'
   write(file_num,'(a)') '# created by TiberCAD '
   write(file_num,'(a)') '# use jmol script: '
   write(file_num,'(a)') '#$ isosurface color yellow "mo_01.jvxl"'   

   close(file_num)
      
end do


deallocate(atomcharge,xyzcoord,rho,planes,edges,indexPlanes)


end subroutine write_jvxl

! ------------------------------------------------------------------------------

!* The key to the algorithm is that we have a catalog that
!* maps the inside-vertex mask to an edge mask, and then
!* each edge is associated with a specific vertex.
!*
!* Each cube vertex may be associated with from 0 to 3 edges,
!* depending upon where it lies in the overall cube of data.
!*
!* When scanning X from low to high, the "leading vertex" is
!* vertex 1 and edgeVertexPlanes[1]. Edges 0, 1, and 9 are
!* associated with vertex 1, and others are associated similarly.
!*
!* When scanning X from high to low, the "leading vertex" is
!* vertex 0 and edgeVertexPlanes[1]. Edges 0, 3, and 8 are
!* associated with vertex 0, and others are associated similarly.
!*
!* edgePointIndexes[iEdge] tracks the vertex index for this
!* specific cubical so that triangles can be created properly.
!*
!*                                                 streaming data offsets
!*     Y
!*         4 --------4--------- 5                   +z --------4--------- +yz+z
!*         /|                   /|                   /|                   /|
!*        / |                  / |                  / |                  / |
!*       /  |                 /  |                 /  |                 /  |
!*      7   8                5   |                7   8                5   |
!*     /    |               /    9               /    |               /    9
!*    /     |              /     |              /     |              /     |
!*   7 --------6--------- 6      |          +z+1 --------6--------- +yz+z+1|
!*   |      |             |      |             |      |             |      |
!*   |      0 ---------0--|----- 1    X        |      0 ---------0--|----- +yz  X
!*   |     /              |     /              |     /              |     /
!*  11    /               10   /              11    /               10   /
!*   |   3                |   1                |   3                |   1
!*   |  /                 |  /                 |  /                 |  /
!*   | /                  | /                  | /                  | /
!*   3 ---------2-------- 2                   +1 ---------2-------- +yz+1
!*    Z                                           Z (inner)
!*
!* edgeVertexPlanes(2)   (1) (scanning x high to low)
!*
!* type 0: x-edges: 0 2 4 6
!* type 1: y-edges: 8 9 10 11
!* type 2: z-edges: 1 3 5 7
!*
!* Data stream offsets for vertices, relative to point 0, based on reading
!* loops {for x {for y {for z}}} 0-->n-1
!* y and z are numbers of grid points in those directions:
!*
!*            0    1      2      3      4      5      6        7
!* lnOffset=  0   +yz   +yz+1   +1     +z    +yz+z  +yz+z+1  +z+1
!*
!* iEdge=        11  10   9   8   7   6   5   4   3   2   1   0
!* edgeVertxPnt=  0,  1,  3,  0,  4,  5,  7,  4,  0,  1,  2,  3 
!* edgeVertxPln=  2,  1,  2,  2,  2,  1,  2,  2,  2,  1,  1,  2
!*
!* These are just looked up in a table. After the first set of cubes,
!* we are only adding points 0, 3, 4 or 7. This means that initially
!* we need two data slices, but after that only one (slice 1):
!*
!*            base
!*           offset 0    1      2      3      4      5      6     7
!*  slice[0]  +yz        0     +1                   +z    +z+1
!*  slice[1]        0                 +1     +z                 +z+1
!*
!*  slice:          1    0      0      1      1      0      0     1
!*
!*  We can request reading of two slices (2*nY*nZ data points) first, then
!*  from then on, just nY*nZ points. "Reading" is really just being handed a
!*  pointer into an array. Perhaps that array is already filled completely;
!*  perhaps it is being read incrementally.
!*
!*  As it is now, the JVXL data are just read into an [nX][nY][nZ] array anyway,
!*  so we can continue to do that with NON progressive files.
!*-------------------------------------------------------------------------------------------
  
subroutine processOneCubical(insideMask,vertexValues,pt,Nyz,Nz,cutoff, &
     indexPlanes,edges,Nedge)
  
  integer, intent(in) :: insideMask, pt, Nyz, Nz
  real(dp) :: vertexValues(0:7)
  real(dp), intent(in) :: cutoff
  integer :: indexPlanes(3,Nyz,2)
  character(1), dimension(:) :: edges
  integer :: Nedge

  integer :: iEdge, xEdge, edgeMask
  integer :: iPlane, iPt, iType, index, vertexA, vertexB
  real(dp) :: f, valueA, valueB
  integer :: linearOffset(0:7)

  linearOffset = (/ 0, Nyz, Nyz+1, 1, Nz, Nyz+Nz, Nyz+Nz+1, Nz+1 /)

  edgeMask = insideMaskTable(insideMask)
  
  do iEdge = 11, 0, -1

     xEdge = pow2(iEdge)
     
     if (iand(edgeMask,xEdge).eq.0) cycle        
     
     iPlane = edgeVertexPlanes(iEdge)
     iPt = mod(pt + linearOffset(edgeVertexPointers(iEdge)), Nyz) + 1
     iType = edgeTypeTab(iEdge)
     index = indexPlanes(iType,iPt,iPlane)
    
     if(index .ge. 0) cycle

     vertexA = edgeVerteces(2*iEdge) 
     vertexB = edgeVerteces(2*iEdge + 1) 

     valueA = vertexValues(vertexA)
     valueB = vertexValues(vertexB)

     f = (cutoff - valueA) / (valueB - valueA);
     
     NEdge = NEdge + 1
     
     edges(NEdge) = jvxlFractionAsChar(f)

     indexPlanes(iType,iPt,iPlane) = NEdge 

  enddo


end subroutine processOneCubical


! ------------------------------------------------------------------------------

function jvxlFractionAsChar(fraction) result(ch)
  real(dp) :: fraction
  integer :: ich
  character(1) :: ch

  if (fraction > 0.9999_dp) then
     fraction = 0.9999_dp
  else if (isnan(fraction)) then
     fraction = 1.0001_dp
  endif

  ich = floor(fraction * RANGE + BASE)
  if (ich < BASE) ch = char(BASE)
  if (ich == 92) ch = char(33)
  ch = char(ich)
 
end function jvxlFractionAsChar
! -----------------------------------------------------------------------------
subroutine find_cutoff(rho,fraction,cutoff)
  
  real(dp), dimension(:,:,:), allocatable :: rho
  real(dp), intent(in) :: fraction
  real(dp), intent(out) :: cutoff

  !locals
  real(dp) :: cutoff1,cutoff2,partial,norm
  integer :: ind

  norm = get_norm(rho)

  cutoff1 = 1.d-10;  cutoff2 = 5.d-4
  partial = 1e6
  ind = 0
  do while( abs(partial-fraction).gt.1.d-4 .and. ind.lt.30 )
     
     cutoff = (cutoff1 + cutoff2)/2
     partial = get_partial(rho,cutoff) * norm
     
     !write(*,*) 'cutoff=',cutoff,'partial=',partial
     
     if (partial.gt.fraction) then
        cutoff1 = cutoff 
     else
        cutoff2 = cutoff 
     endif
     
     ind = ind + 1
  end do
 
  !//////

end subroutine find_cutoff

! -----------------------------------------------------------------------------
function get_norm(rho) result(norm)

  real(dp), dimension(:,:,:), allocatable :: rho 
  
  real(dp) :: norm
  integer :: i,j,k,np1,np2,np3
  np1 = size(rho,1)
  np2 = size(rho,2)
  np3 = size(rho,3)

  norm = 0.d0
  do i = 1, Np1
     do j = 1, Np2
        do k = 1, Np3
          
           norm = norm + rho(i,j,k)
           
        end do
     end do
  end do

end function get_norm

! -----------------------------------------------------------------------------
function get_partial(rho, cutoff) result(partial)

  real(dp), dimension(:,:,:), allocatable :: rho 
  real(dp), intent(in) :: cutoff

  real(dp) :: partial
  integer :: i,j,k,np1,np2,np3

  np1 = size(rho,1)
  np2 = size(rho,2)
  np3 = size(rho,3)

  partial = 0.d0
  do i = 1, Np1
     do j = 1, Np2
        do k = 1, Np3
          
           if(rho(i,j,k).gt.cutoff) then
              partial = partial + rho(i,j,k)
           endif
           
        end do
     end do
  end do
  
end function get_partial
! -----------------------------------------------------------------------------
subroutine write_fmt_int(file_num,int)
  integer :: file_num,int

  if (int.lt.10) then
     write(file_num,'(1x,i1)',advance='NO') int
  elseif (int.lt.100) then
     write(file_num,'(1x,i2)',advance='NO') int  
  elseif (int.lt.1000) then
     write(file_num,'(1x,i3)',advance='NO') int 
  elseif (int.lt.10000) then
     write(file_num,'(1x,i4)',advance='NO') int 
  elseif (int.lt.100000) then
     write(file_num,'(1x,i5)',advance='NO') int 
  elseif (int.lt.1000000) then
     write(file_num,'(1x,i6)',advance='NO') int 
  endif

end subroutine write_fmt_int
! -----------------------------------------------------------------------------
function atomZ(name)
  
  character(2) :: name
  integer :: atomZ
  
  select case (trim(name))
     
   case('Hx')
      atomZ=1
   case('H')
      atomZ=1
   case('He')
      atomZ=2
   
   case('Li')
      atomZ=3
   case('Be')
      atomZ=4 
   case('B')
      atomZ=5 
   case('C')
      atomZ=6 
   case('N')
      atomZ=7
   case('O')
      atomZ=8
   case('F')
      atomZ=9
   case('Ne')
      atomZ=10

   case('Na')
      atomZ=11
   case('Mg')
      atomZ=12 
   case('Al')
      atomZ=13 
   case('Si')
      atomZ=14
   case('P')
      atomZ=15
   case('S')
      atomZ=16
   case('Cl')
      atomZ=17
   case('Ar')
      atomZ=18

   case('K')
      atomZ=19
   case('Ca')
      atomZ=20
   case('Sc')
      atomZ=21
   case('Ti')
      atomZ=22
   case('V')
      atomZ=23
   case('Cr')
      atomZ=24 
   case('Mn')
      atomZ=25
   case('Fe')
      atomZ=26
   case('Co')
      atomZ=27
   case('Ni')
      atomZ=28
   case('Cu')
      atomZ=29
   case('Zn')
      atomZ=30
   case('Ga')
      atomZ=31
   case('Ge')
      atomZ=32
   case('As')
      atomZ=33
   case('Se')
      atomZ=34
   case('Br')
      atomZ=35
   case('Kr')
      atomZ=36


   case('Rb')
      atomZ=37
   case('Sr')
      atomZ=38
   case('Y')
      atomZ=39
   case('Zr')
      atomZ=40
   case('Nb')
      atomZ=41
   case('Mo')
      atomZ=42
   case('Tc')
      atomZ=43
   case('Ru')
      atomZ=44
   case('Rh')
      atomZ=45
   case('Pd')
      atomZ=46
   case('Ag')
      atomZ=47
   case('Cd')
      atomZ=48
   case('In')
      atomZ=49
   case('Sn')
      atomZ=50
   case('Sb')
      atomZ=51
   case('Te')
      atomZ=52
   case('I')
      atomZ=53
   case('Xe')
      atomZ=54

   case('Cs')
      atomZ=55
   case('Ba')
      atomZ=56
   case('Lu')
      atomZ=71
   case('Hf')
      atomZ=72
   case('Ta')
      atomZ=73
   case('W')
      atomZ=74
   case('Re')
      atomZ=75
   case('Os')
      atomZ=76
   case('Ir')
      atomZ=77
   case('Pt')
      atomZ=78
   case('Au')
      atomZ=79
   case('Hg')
      atomZ=80
   case('Tl')
      atomZ=81
   case('Pb')
      atomZ=82
   case('Bi')
      atomZ=83
   case('Po')
      atomZ=84
   case('At')
      atomZ=85
   case('Rn')
      atomZ=86

   case('La')
      atomZ=57
   case('Ce')
      atomZ=58
   case('Nd')
      atomZ=60
   case('Sm')
      atomZ=62
   case('Eu')
      atomZ=63
   case('Gd')
      atomZ=64
   case('Tb')
      atomZ=65
   case('Dy')
      atomZ=66
   case('Ho')
      atomZ=67
   case('Er')
      atomZ=68
   case('Tm')
      atomZ=69
   case('Yb')
      atomZ=70

   end select

 end function atomZ

 ! Localization computes the fraction of atoms on which a given cutoff
 ! of the total probability is located (say 80%)
 function localization(atomcharge,cutoff) result(loca)
   real(dp), dimension(:), target :: atomcharge
   real(dp), intent(in) :: cutoff
   real(dp) :: loca

   integer :: nat, i
   real(dp), dimension(:), pointer :: charges
   integer, dimension(:), pointer :: indx

   nat=size(atomcharge)

   allocate( indx(nat) )

   charges=>atomcharge(:)

   ! Sort in ascending order
   call osort_index(nat, charges, indx)

   loca = 0.d0

   do i = nat, 1, -1
     loca = loca + atomcharge(indx(i))  
     if (loca .ge. cutoff) exit 
   end do  
   
   loca = (nat-i+1)*1.d0/nat

 end function localization


end module savemofile
