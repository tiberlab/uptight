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
!=============================================================================
!
!            Module "vectors" - (c) Jerome Gleize - 2002
!
!=============================================================================
!
! => norm
! => sort_vec
! => sort_vec_by_norm
! => vec_count
! => vec_values
! => cross_product
! => triple_product
!
!=============================================================================

MODULE vectors

  !===========================================================================

  USE precision 

  !===========================================================================

  IMPLICIT NONE
  PRIVATE

  PUBLIC norm, sort_vec, sort_by_norm, vec_count, vec_values
  PUBLIC cross_product, triple_product

  !===========================================================================
  !
  ! NOTE : some routines are generic to handle different argument types
  !        with the same routine.
  !
  !===========================================================================

  INTERFACE norm

     MODULE PROCEDURE norm_int, norm_real, norm_dp, norm_dp_mat 

  END INTERFACE

  !__________________________________________________________________________

  INTERFACE sort_vec
  
     MODULE PROCEDURE sort_vec_int, sort_vec_real, sort_vec_dp
     
  END INTERFACE

  !__________________________________________________________________________

  INTERFACE sort_by_norm
     
     MODULE PROCEDURE sort_by_norm_int, sort_by_norm_real, sort_by_norm_dp
     
  END INTERFACE

  !___________________________________________________________________________

  INTERFACE vec_count
     
     MODULE PROCEDURE vec_count_int, vec_count_real, vec_count_dp
     
  END INTERFACE
  
  !___________________________________________________________________________
  
  INTERFACE vec_values
     
     MODULE PROCEDURE vec_val_int, vec_val_real, vec_val_dp
     
  END INTERFACE
  
  !___________________________________________________________________________
  
  INTERFACE cross_product
     
     MODULE PROCEDURE cross_product_int, cross_product_real, cross_product_dp
     
  END INTERFACE
  
  !___________________________________________________________________________
  
  INTERFACE triple_product
     
     MODULE PROCEDURE triple_product_int, triple_product_real, &
          triple_product_dp
     
  END INTERFACE
 
  !===========================================================================

CONTAINS


  !===========================================================================
  !
  ! Functions "norm_*" : calculate the norm of a vector "vec"
  !                      in cartesian orthogonal axes.
  !
  !===========================================================================
  
  FUNCTION norm_int( vec )

    !=========================================================================

    ! input arguments :

    INTEGER, DIMENSION( : ), INTENT( IN ) :: vec

    !_________________________________________________________________________

    ! output result :
   
    REAL(dp) :: norm_int
    
    !=========================================================================
    
    norm_int = SQRT( DBLE( SUM( vec )**2 ) )

    !=========================================================================
    
  END FUNCTION norm_int
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  FUNCTION norm_real( vec )

    !=========================================================================

    ! input arguments :

    REAL(sp), DIMENSION( : ), INTENT( IN ) :: vec

    !_________________________________________________________________________

    ! output result :
   
    REAL(sp) :: norm_real

    !=========================================================================
    
    norm_real = SQRT( SUM( vec**2 ) )

    !=========================================================================
    
  END FUNCTION norm_real
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  FUNCTION norm_dp( vec )
    
    !=========================================================================
    
    ! input arguments :
    
    REAL(dp), DIMENSION( : ), INTENT( IN ) :: vec
    
    !_________________________________________________________________________
    
    ! output result :
    
    REAL(dp) :: norm_dp
    
    !=========================================================================
    
    norm_dp = SQRT( SUM( vec**2 ) )

    !=========================================================================
    
  END FUNCTION norm_dp


 FUNCTION norm_dp_mat( vec )
    
    !=========================================================================
    
    ! input arguments :
    
    REAL(dp), DIMENSION( :,: ), INTENT( IN ) :: vec
    
    !_________________________________________________________________________
    
    ! output result :
    
    REAL(dp) , DIMENSION( SIZE(vec(:,1)) ) :: norm_dp_mat
    
    !=========================================================================
    
    norm_dp_mat = SQRT( SUM( vec**2,1 ) )

    !=========================================================================
    
  END FUNCTION norm_dp_mat


  !===========================================================================
  !
  ! Functions "sort_vec_*" : sort a vector by increasing values.
  !
  !===========================================================================
  !
  ! The function returns the sorted vector.
  ! An auxiliary copy "vec_aux" is used during sorting, to avoid erasing
  ! the original vector "vec".
  !
  !===========================================================================
  
  FUNCTION sort_vec_int( vec )
    
    !=========================================================================
    
    ! input arguments :
    
    INTEGER, DIMENSION( : ), INTENT( IN ) :: vec
    
    !_________________________________________________________________________

    ! result output :

    INTEGER, DIMENSION( SIZE( vec ) ) :: sort_vec_int

    !_________________________________________________________________________
    
    ! local variables :
    
    INTEGER, DIMENSION( SIZE( vec ) ) :: vec_aux
    
    LOGICAL, DIMENSION( SIZE( vec ) ) :: vec_test

    INTEGER :: mini, maxi, i_sort, n_test
    
    !========================================================================
    
    ! initialization :
    
    vec_aux = vec                 ! auxiliary vector
    
    maxi = MAXVAL( vec_aux )      ! maximum value of "vec"
    
    i_sort = 0                    ! sorted values counter
    
    !=========================================================================
    
    ! loop on sorted values

    DO WHILE ( i_sort .LT. SIZE( vec ) )
       
       !______________________________________________________________________
       
       ! find the current minimum value of "vec_aux"
       
       mini = MINVAL( vec_aux )

       ! label the corresponding elements of "vec_aux"
       
       vec_test = vec_aux .EQ. mini
       n_test = COUNT( vec_test )
       
       ! store these elements in the sorted vector
       
       sort_vec_int( i_sort + 1 : i_sort + n_test ) = PACK( vec_aux, vec_test )
       
       ! update the sorted values counter 
       
       i_sort = i_sort + n_test
       
       ! nullify the current minimum elements of "vec_aux"

       WHERE( vec_test ) vec_aux = maxi + 1
       
       !______________________________________________________________________
       
    END DO
    
    !=========================================================================
    
  END FUNCTION sort_vec_int

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  FUNCTION sort_vec_real( vec )
    
    !=========================================================================
    
    ! input arguments :
    
    REAL(sp), DIMENSION( : ), INTENT( IN ) :: vec
    
    !_________________________________________________________________________

    ! result output :

    REAL(sp), DIMENSION( SIZE( vec ) ) :: sort_vec_real

    !_________________________________________________________________________
    
    ! local variables :

    REAL(sp), DIMENSION( SIZE( vec ) ) :: vec_aux
    
    LOGICAL, DIMENSION( SIZE( vec ) ) :: vec_test
    
    REAL(sp) :: mini, maxi
    
    INTEGER :: i_sort, n_test

    !========================================================================
     
    ! initialization :
    
    vec_aux = vec                 ! auxiliary vector
    
    maxi = MAXVAL( vec_aux )      ! maximum value of "vec"
    
    i_sort = 0                    ! sorted values counter 
        
    !=========================================================================

    ! loop on sorted values
    
    DO WHILE ( i_sort .LT. SIZE( vec ) )
       
       !______________________________________________________________________
        
       ! find the current minimum value of "vec_aux"
       
       mini = MINVAL( vec_aux )

       ! label the corresponding elements of "vec_aux"
       
       vec_test = equiv( DBLE( vec_aux ), DBLE( mini ), prec, fuzzy )
       n_test = COUNT( vec_test )
              
       ! store these elements in the sorted vector
       
       sort_vec_real( i_sort + 1 : i_sort + n_test ) = &
            PACK( vec_aux, vec_test )
       
       ! update the sorted values counter 
       
       i_sort = i_sort + n_test
       
       ! nullify the current minimum elements of "vec_aux"
       
       WHERE( vec_test ) vec_aux = maxi + 1.0

       !______________________________________________________________________
       
    END DO
    
    !=========================================================================
    
  END FUNCTION sort_vec_real

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  FUNCTION sort_vec_dp( vec )
    
    !=========================================================================
    
    ! input arguments :
    
    REAL(dp), DIMENSION( : ), INTENT( IN ) :: vec
    
    !_________________________________________________________________________

    ! result output :

    REAL(dp), DIMENSION( SIZE( vec ) ) :: sort_vec_dp

    !_________________________________________________________________________
    
    ! local variables :

    REAL(dp), DIMENSION( SIZE( vec ) ) :: vec_aux

    LOGICAL, DIMENSION( SIZE( vec ) ) :: vec_test

    REAL(dp) :: mini, maxi

    INTEGER :: i_sort, n_test

    !========================================================================
         
    ! initialization :
    
    vec_aux = vec                 ! auxiliary vector
    
    maxi = MAXVAL( vec_aux )      ! maximum value of "vec"
    
    i_sort = 0                    ! sorted values counter 
        
    !=========================================================================
    
    ! loop on sorted values
    
    DO WHILE ( i_sort .LT. SIZE( vec ) )
       
       !______________________________________________________________________
        
       ! find the current minimum value of "vec_aux"
       
       mini = MINVAL( vec_aux )

       ! label the corresponding elements of "vec_aux"
          
       vec_test = equiv( vec_aux, mini, prec, fuzzy )
       n_test = COUNT( vec_test )
       
       ! store these elements in the sorted vector
       
       sort_vec_dp( i_sort + 1 : i_sort + n_test ) = PACK( vec_aux, vec_test )
       
       ! update the sorted values counter 
       
       i_sort = i_sort + n_test
       
       ! nullify the current minimum elements of "vec_aux"
              
       WHERE( vec_test ) vec_aux = maxi + 1.0d0

       !______________________________________________________________________
       
    END DO
    
    !=========================================================================
    
  END FUNCTION sort_vec_dp



  !===========================================================================
  !
  ! Functions "sort_by_norm_*" : sort an array of vectors by increasing norm.
  !
  !===========================================================================
  !
  ! The function returns the sorted array of vectors.
  ! An auxiliary copy "vec_aux" is used during sorting, to avoid erasing
  ! the original array of vectors "vec".
  !
  ! => INPUT :
  !
  !  "vec" is an ( n1 x n2 ) array of n2 vectors of n1 components each.
  !
  !===========================================================================

  FUNCTION sort_by_norm_int( vec )
    
    !=========================================================================
    
    ! input arguments :
    
    INTEGER, DIMENSION( :, : ), INTENT( IN ) :: vec

    !_________________________________________________________________________
    
    ! output result :
    
    INTEGER, DIMENSION( SIZE( vec,1 ), SIZE( vec, 2 ) ) :: sort_by_norm_int

    !_________________________________________________________________________
    
    ! local variables :
    
    INTEGER, DIMENSION( SIZE( vec,1 ), SIZE( vec, 2 ) ) :: vec_aux
    
    INTEGER, DIMENSION( SIZE( vec, 2 ) ) :: vec_norm
    
    LOGICAL, DIMENSION( SIZE( vec, 2 ) ) :: vec_norm_test
    
    INTEGER :: mini, maxi, i_sort, i_coord, n_test
    
    !========================================================================

    ! initialization :

    vec_aux = vec                  ! auxiliary array
    
    vec_norm = vec( 1, : )**2 + vec( 2, : )**2 + vec( 3, : )**2  ! norm array
    
    maxi = MAXVAL( vec_norm )      ! maximum value of norm
    
    i_sort = 0                     ! sorted elements counter
    
    !=========================================================================
    
    ! loop on sorted elements

    DO WHILE ( i_sort .LT. SIZE( vec_norm ) )
       
       !______________________________________________________________________
       
       ! find the current minimum norm
       
       mini = MINVAL( vec_norm )
       
       ! label the corresponding vectors in the array "vec"
       
       vec_norm_test = ( vec_norm .EQ. mini )
       n_test = COUNT( vec_norm_test )

       !______________________________________________________________________
           
       ! store these elements in the sorted array

       DO i_coord = 1, SIZE( vec, 1 )
          
          sort_by_norm_int( i_coord, i_sort + 1 : i_sort + n_test ) = &
               PACK( vec_aux( i_coord, : ), vec_norm_test )
                    
       END DO
      
       !______________________________________________________________________
      
       ! update the sorted elements counter 

       i_sort = i_sort + n_test

       ! nullify the current elements of minimum norm
       
       WHERE( vec_norm_test ) vec_norm = maxi + 1
       
       !______________________________________________________________________
       
    END DO

    !=========================================================================

  END FUNCTION sort_by_norm_int

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  FUNCTION sort_by_norm_real( vec )
    
    !=========================================================================

    ! input - output arguments :

    REAL(sp), DIMENSION( :, : ), INTENT( INOUT ) :: vec

    !_________________________________________________________________________

    ! result output :

    REAL(sp), DIMENSION( SIZE( vec,1 ), SIZE( vec, 2 ) ) :: sort_by_norm_real

    !_________________________________________________________________________
    
    ! local variables :
    
    REAL(sp), DIMENSION( SIZE( vec,1 ), SIZE( vec, 2 ) ) :: vec_aux
    
    REAL(sp), DIMENSION( SIZE( vec, 2 ) ) :: vec_norm

    LOGICAL, DIMENSION( SIZE( vec, 2 ) ) :: vec_norm_test
    
    REAL(sp) :: mini, maxi
    
    INTEGER :: i_sort, i_coord, n_test
    
    !========================================================================
    
    ! initialization :
    
    vec_aux = vec                  ! auxiliary array
    
    vec_norm = vec( 1, : )**2 + vec( 2, : )**2 + vec( 3, : )**2  ! norm array
    
    maxi = MAXVAL( vec_norm )      ! maximum value of norm
    
    i_sort = 0                     ! sorted elements counter
    
    !=========================================================================
    
    ! loop on sorted elements

    DO WHILE ( i_sort .LT. SIZE( vec_norm ) )
       
       !______________________________________________________________________
       
       ! find the current minimum norm
       
       mini = MINVAL( vec_norm )
       
       ! label the corresponding vectors in the array "vec"
       
       vec_norm_test = equiv( DBLE( vec_norm ), DBLE( mini ), prec, fuzzy )
       n_test = COUNT( vec_norm_test )

       !______________________________________________________________________
           
       ! store these elements in the sorted array
       
       DO i_coord = 1, SIZE( vec, 1 )
          
          sort_by_norm_real( i_coord, i_sort + 1 : i_sort + n_test ) = &
               PACK( vec_aux( i_coord, : ), vec_norm_test )
                    
       END DO
       
       !______________________________________________________________________
       
       ! update the sorted elements counter 
       
       i_sort = i_sort + n_test
       
       ! nullify the current elements of minimum norm
       
       WHERE( vec_norm_test ) vec_norm = maxi + 1.0
       
       !______________________________________________________________________
       
    END DO
       
    !=========================================================================

  END FUNCTION sort_by_norm_real

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  FUNCTION sort_by_norm_dp( vec )
    
    !=========================================================================

    ! input - output arguments :

    REAL(dp), DIMENSION( :, : ), INTENT( INOUT ) :: vec

    !_________________________________________________________________________

    REAL(dp), &
         DIMENSION( SIZE( vec,1 ), SIZE( vec, 2 ) ) :: sort_by_norm_dp

    !_________________________________________________________________________
    
    ! local variables :
    
    REAL(dp), DIMENSION( SIZE( vec,1 ), SIZE( vec, 2 ) ) :: vec_aux
    
    REAL(dp), DIMENSION( SIZE( vec, 2 ) ) :: vec_norm

    LOGICAL, DIMENSION( SIZE( vec, 2 ) ) :: vec_norm_test

    REAL(dp) :: mini, maxi

    INTEGER :: i_sort, i_coord, n_test
    
    !========================================================================
   
    ! initialization :
    
    vec_aux = vec                  ! auxiliary array
    
    vec_norm = vec( 1, : )**2 + vec( 2, : )**2 + vec( 3, : )**2  ! norm array
    
    maxi = MAXVAL( vec_norm )      ! maximum value of norm
    
    i_sort = 0                     ! sorted elements counter
    
    !=========================================================================

    ! loop on sorted elements

    DO WHILE ( i_sort .LT. SIZE( vec_norm ) )
       
       !______________________________________________________________________
       
       ! find the current minimum norm
       
       mini = MINVAL( vec_norm )
       
       ! label the corresponding vectors in the array "vec"
       
       vec_norm_test = equiv( vec_norm, mini, prec, fuzzy )
       n_test = COUNT( vec_norm_test )

       !______________________________________________________________________
           
       ! store these elements in the sorted array
       
       DO i_coord = 1, SIZE( vec, 1 )
          
          sort_by_norm_dp( i_coord, i_sort + 1 : i_sort + n_test ) = &
               PACK( vec_aux( i_coord, : ), vec_norm_test )
                    
       END DO
       
       !______________________________________________________________________
       
       ! update the sorted elements counter 
       
       i_sort = i_sort + n_test
       
       ! nullify the current elements of minimum norm
       
       WHERE( vec_norm_test ) vec_norm = maxi + 1.0d0
       
       !______________________________________________________________________
       
    END DO

    !=========================================================================

  END FUNCTION sort_by_norm_dp



  !===========================================================================
  !
  ! Functions "vec_count_*" :
  !
  ! count the number of distinct elements in a vector
  !
  !===========================================================================

  FUNCTION vec_count_int( vec )

    !=========================================================================

    ! input arguments :

    INTEGER, DIMENSION( : ), INTENT( IN ) :: vec

    !_________________________________________________________________________

    ! output result :

    INTEGER :: vec_count_int

    !_________________________________________________________________________

    ! local variables :

    INTEGER, DIMENSION( SIZE( vec ) ) :: vec_aux
    
    INTEGER :: mini, maxi

    !=========================================================================

    ! initialization :

    vec_count_int = 0

    vec_aux = vec                ! auxiliary vector

    maxi = MAXVAL( vec )         ! maximum value of "vec"
    
    mini = MINVAL( vec_aux )     ! minimum value of "vec_aux"

    !=========================================================================

    ! loop on vector elements : stop when maximum value is reached

    DO WHILE ( mini .LE. maxi )

       !______________________________________________________________________

       ! increment element counter

       vec_count_int = vec_count_int + 1

       ! nullify the corresponding elements of "vec_aux"

       WHERE ( vec_aux .EQ. mini ) vec_aux = maxi + 1

       ! update the minimum value
       
       mini = MINVAL( vec_aux )

       !______________________________________________________________________

    END DO

    !=========================================================================

  END FUNCTION vec_count_int

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  FUNCTION vec_count_real( vec )

    !=========================================================================

    ! input arguments :

    REAL(sp), DIMENSION(:), INTENT( IN ) :: vec

    !_________________________________________________________________________

    ! output result :

    INTEGER :: vec_count_real

    !_________________________________________________________________________

    ! local variables :

    REAL(sp), DIMENSION( SIZE( vec ) ) :: vec_aux

    REAL(sp) :: mini, maxi

    !=========================================================================
   
    ! initialization :

    vec_count_real = 0

    vec_aux = vec                ! auxiliary vector

    maxi = MAXVAL( vec )         ! maximum value of "vec"
    
    mini = MINVAL( vec_aux )     ! minimum value of "vec_aux"

    !=========================================================================
 
    ! loop on vector elements : stop when maximum value is reached

    DO WHILE ( mini .LE. maxi )

       !______________________________________________________________________
       
       ! increment element counter

       vec_count_real = vec_count_real + 1

       ! nullify the corresponding elements of "vec_aux"

       WHERE ( equiv( DBLE( vec_aux ), DBLE( mini ), prec, fuzzy ) ) &
            vec_aux = maxi + 1.0

       ! update the minimum value
       
       mini = MINVAL( vec_aux )
       
       !______________________________________________________________________
       
    END DO
    
    !=========================================================================
   
  END FUNCTION vec_count_real

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  FUNCTION vec_count_dp( vec )

    !=========================================================================

    ! input arguments :

    REAL(dp), DIMENSION( : ), INTENT( IN ) :: vec

    !_________________________________________________________________________

    ! output result :

    INTEGER :: vec_count_dp

    !_________________________________________________________________________

    ! local variables :

    REAL(dp), DIMENSION( SIZE( vec ) ) :: vec_aux
    
    REAL(dp) :: mini, maxi

    !=========================================================================

    ! initialization :

    vec_count_dp = 0
    
    vec_aux = vec                ! auxiliary vector
    
    maxi = MAXVAL( vec )         ! maximum value of "vec"
    
    mini = MINVAL( vec_aux )     ! minimum value of "vec_aux"
    
    !=========================================================================
 
    ! loop on vector elements : stop when maximum value is reached
    
    DO WHILE ( mini .LE. maxi )

       !______________________________________________________________________
       
       ! increment element counter

       vec_count_dp = vec_count_dp + 1

       ! nullify the corresponding elements of "vec_aux"

       WHERE ( equiv( vec_aux, mini, prec, fuzzy ) ) vec_aux = maxi + 1.0d0

       ! update the minimum value
       
       mini = MINVAL( vec_aux )
       
       !______________________________________________________________________
       
    END DO
        
    !=========================================================================

  END FUNCTION vec_count_dp

  

  !===========================================================================
  !
  ! Functions "vec_values_* : gets all distinct values of a vector "vec"
  !
  !===========================================================================

  FUNCTION vec_val_int( vec, n_value )

    !=========================================================================

    ! input arguments :

    INTEGER, DIMENSION( : ), INTENT( IN ) :: vec
    INTEGER,                 INTENT( IN ) :: n_value

    !_________________________________________________________________________
    
    ! output result :
    
    INTEGER, DIMENSION( n_value ) :: vec_val_int

    !_________________________________________________________________________

    ! local variables :

    INTEGER, DIMENSION( SIZE( vec ) ) :: vec_aux
   
    INTEGER :: mini, maxi

    INTEGER :: i_value

    !=========================================================================

    ! initialization
    
    vec_aux = vec               ! auxiliary array
    
    maxi = MAXVAL( vec_aux )    ! maximum value of "vec"
    mini = MINVAL( vec_aux )    ! minimum value of "vec_aux"

    !=========================================================================
    
    ! loop on distinct values

    DO i_value = 1, n_value

       !______________________________________________________________________
       
       ! store the current minimum
       
       vec_val_int( i_value ) = mini
       
       ! nullify the corresponding elements of "vec_aux"
       
       WHERE ( vec_aux .EQ. mini ) vec_aux = maxi + 1
       
       ! update the minimum value of "vec_aux"
       
       mini = MINVAL( vec_aux )

       !______________________________________________________________________

    END DO

    !=========================================================================

  END FUNCTION vec_val_int

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  FUNCTION vec_val_real( vec, n_value )

    !=========================================================================

    ! input arguments :
    
    REAL(sp), DIMENSION( : ), INTENT( IN ) :: vec
    
    INTEGER,              INTENT( IN ) :: n_value

    !_________________________________________________________________________

    ! output result :

    REAL(sp), DIMENSION( n_value ) :: vec_val_real

    !_________________________________________________________________________

    ! local variables :

    REAL(sp), DIMENSION( SIZE( vec ) ) :: vec_aux

    REAL(sp) :: mini, maxi

    INTEGER :: i_value

    !=========================================================================
 
    ! initialization
    
    vec_aux = vec               ! auxiliary array
    
    maxi = MAXVAL( vec_aux )    ! maximum value of "vec"
    mini = MINVAL( vec_aux )    ! minimum value of "vec_aux"
    
    !=========================================================================
    
    ! loop on distinct values
    
    DO i_value = 1, n_value

       !______________________________________________________________________
         
       ! store the current minimum

       vec_val_real( i_value ) = mini

       ! nullify the corresponding elements of "vec_aux"

       WHERE ( equiv( DBLE( vec_aux ), DBLE( mini ), prec, fuzzy ) ) &
            vec_aux = maxi + 1.0
            
       ! update the minimum value of "vec_aux"

       mini = MINVAL( vec_aux )
       
       !______________________________________________________________________

    END DO

    !=========================================================================
    
  END FUNCTION vec_val_real
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  FUNCTION vec_val_dp( vec, n_value )
    
    !=========================================================================

    ! input arguments :

    REAL(dp), DIMENSION( : ), INTENT( IN ) :: vec

    INTEGER,                          INTENT( IN ) :: n_value

    !_________________________________________________________________________

    ! output result :
  
    REAL(dp), DIMENSION( n_value ) :: vec_val_dp

    !_________________________________________________________________________

    ! local variables :

    REAL(dp), DIMENSION( SIZE( vec ) ) :: vec_aux
    
    REAL(dp) :: mini, maxi

    INTEGER :: i_value

    !=========================================================================
 
    ! initialization
    
    vec_aux = vec               ! auxiliary array
    
    maxi = MAXVAL( vec_aux )    ! maximum value of "vec"
    mini = MINVAL( vec_aux )    ! minimum value of "vec_aux"
    
    !=========================================================================
 
    ! loop on distinct values

    DO i_value = 1, n_value
       
       !______________________________________________________________________

       ! store the current minimum
       
       vec_val_dp( i_value ) = mini

       ! nullify the corresponding elements of "vec_aux"
       
       WHERE ( equiv( vec_aux, mini, prec, fuzzy ) ) vec_aux = maxi + 1.0d0

       ! update the minimum value of "vec_aux"

       mini = MINVAL( vec_aux )

       !______________________________________________________________________

    END DO

    !=========================================================================
    
  END FUNCTION vec_val_dp



  !===========================================================================
  !
  ! Functions "cross_product_*" : calculate the cross product of vectors
  !                               vec_a and vec_b, in an orthogonal frame.
  !
  !===========================================================================
 
  FUNCTION cross_product_int( vec_a, vec_b )

    !_________________________________________________________________________
    
    ! input arguments :

    INTEGER, DIMENSION( 3 ), INTENT( IN ) :: vec_a, vec_b

    !_________________________________________________________________________

    ! output result :
   
    INTEGER, DIMENSION( 3 ) :: cross_product_int

    !_________________________________________________________________________

    cross_product_int( 1 ) = vec_a( 2 ) * vec_b( 3 ) - vec_a( 3 ) * vec_b( 2 )
    cross_product_int( 2 ) = vec_a( 3 ) * vec_b( 1 ) - vec_a( 1 ) * vec_b( 3 )
    cross_product_int( 3 ) = vec_a( 1 ) * vec_b( 2 ) - vec_a( 2 ) * vec_b( 1 )
    
    !_________________________________________________________________________

  END FUNCTION cross_product_int

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  FUNCTION cross_product_real( vec_a, vec_b )

    !_________________________________________________________________________
    
    ! input arguments :

    REAL(sp), DIMENSION( 3 ), INTENT( IN ) :: vec_a, vec_b

    !_________________________________________________________________________

    ! output result :
   
    REAL(sp), DIMENSION( 3 ) :: cross_product_real

    !_________________________________________________________________________

    cross_product_real( 1 ) = vec_a( 2 ) * vec_b( 3 ) - vec_a( 3 ) * vec_b( 2 )
    cross_product_real( 2 ) = vec_a( 3 ) * vec_b( 1 ) - vec_a( 1 ) * vec_b( 3 )
    cross_product_real( 3 ) = vec_a( 1 ) * vec_b( 2 ) - vec_a( 2 ) * vec_b( 1 )
    
    !_________________________________________________________________________

  END FUNCTION cross_product_real

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  FUNCTION cross_product_dp( vec_a, vec_b )

    !_________________________________________________________________________
    
    ! input arguments :

    REAL(dp), DIMENSION( 3 ), INTENT( IN ) :: vec_a, vec_b

    !_________________________________________________________________________

    ! output result :
   
    REAL(dp), DIMENSION( 3 ) :: cross_product_dp

    !_________________________________________________________________________

    cross_product_dp( 1 ) = vec_a( 2 ) * vec_b( 3 ) - vec_a( 3 ) * vec_b( 2 )
    cross_product_dp( 2 ) = vec_a( 3 ) * vec_b( 1 ) - vec_a( 1 ) * vec_b( 3 )
    cross_product_dp( 3 ) = vec_a( 1 ) * vec_b( 2 ) - vec_a( 2 ) * vec_b( 1 )
    
    !_________________________________________________________________________

  END FUNCTION cross_product_dp



  !===========================================================================
  !
  ! Functions "triple_product_*" : 
  !
  ! calculate the triple product of vectors vec_a, vec_b, and vec_c,
  ! in an orthogonal frame.
  !
  !===========================================================================
  
  FUNCTION triple_product_int( vec_a, vec_b, vec_c )
    
    !_________________________________________________________________________
    
    ! input arguments :
    
    INTEGER, DIMENSION( 3 ), INTENT( IN ) :: vec_a, vec_b, vec_c

    !_________________________________________________________________________
    
    ! output result :
   
    INTEGER :: triple_product_int

    !_________________________________________________________________________

    triple_product_int = DOT_PRODUCT( vec_a, cross_product( vec_b, vec_c ) )
    
    !_________________________________________________________________________

  END FUNCTION triple_product_int

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  FUNCTION triple_product_real( vec_a, vec_b, vec_c )
    
    !_________________________________________________________________________
    
    ! input arguments :
    
    REAL(sp), DIMENSION( 3 ), INTENT( IN ) :: vec_a, vec_b, vec_c
    
    !_________________________________________________________________________
    
    ! output result :
   
    REAL(sp) :: triple_product_real
    
    !_________________________________________________________________________

    triple_product_real = DOT_PRODUCT( vec_a, cross_product( vec_b, vec_c ) )
    
    !_________________________________________________________________________

  END FUNCTION triple_product_real
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  FUNCTION triple_product_dp( vec_a, vec_b, vec_c )

    !_________________________________________________________________________
    
    ! input arguments :

    REAL(dp), DIMENSION( 3 ), INTENT( IN ) :: vec_a, vec_b, vec_c

    !_________________________________________________________________________

    ! output result :
   
    REAL(dp) :: triple_product_dp

    !_________________________________________________________________________

    triple_product_dp = DOT_PRODUCT( vec_a, cross_product( vec_b, vec_c ) )
    
    !_________________________________________________________________________
    
  END FUNCTION triple_product_dp


  !===========================================================================

END MODULE vectors

!=============================================================================
