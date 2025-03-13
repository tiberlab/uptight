module precision
  
  implicit none
  private

!
! contains :
!
!_____________________________________________________________________________
!
! => general precision definitions
!
! => the subroutine "macheps"
!
! => the generic subroutine "equiv" for test within a given precision
!
! => EPS  :  a small number
!
! => emach : machine accuracy 
!
! => prec :  precision for real data comparison
!
!=============================================================================
  
  integer, parameter :: sp = selected_real_kind(6,30)
  integer, parameter :: dp = selected_real_kind(14,100)
  integer, parameter :: qp = selected_real_kind(28,300)
  
  ! EPS = Small number

  real(dp), parameter :: EPS=1.d-10
  real(dp), parameter :: PREC = 1.0d-04
  real(dp) :: emach 
  
  public :: qp, dp, sp, EPS, PREC, emach, set_machine_acc, equiv


  INTERFACE equiv
     
     MODULE PROCEDURE equiv_val, equiv_vec

  END INTERFACE
  


contains

  subroutine set_machine_acc

    emach = dlamch()

  end subroutine set_machine_acc

  !===========================================================================
  !
  ! Subroutine "macheps" : compute machine accuracy "emach".
  !
  !===========================================================================
  real(dp) function dlamch()

    character(1) :: cmach
    real(dp) :: racc

    racc = 1.d0
    do while((1.d0 + racc) .gt. 1.d0)      
       racc = racc*0.5d0
    enddo
    
    racc = 2.d0*racc
    dlamch=racc

    return

  end function dlamch
  
  !===========================================================================
  !
  ! Function "equiv" : checks if numerical variables a and b are equal.
  !
  ! The comparison is made according to conditions set by the flag "fuzzy" :
  !
  ! => if fuzzy = false :
  !
  !     a equiv b   if   | a - b | < machine accuracy emach
  !
  ! => if fuzzy = true :
  !
  !     a equiv b   if   | a - b | / min( |a|, |b| ) < deviation / 100
  !
  !=========================================================================
  !
  ! INPUT :
  !
  ! => a, b - values to compare.
  ! => deviation - real : tolerance expressed in percent.
  !
  ! OUTPUT :
  !
  ! => equiv - logical.
  !
  !===========================================================================

  FUNCTION equiv_val( a, b, deviation, fuzzy )

    ! input arguments :

    REAL( dp ), INTENT( IN ) :: a, b, deviation
    LOGICAL :: fuzzy

    ! output result :

    LOGICAL :: equiv_val

    ! local variables :

    REAL( dp ) :: difference
    
    equiv_val = .FALSE.

    !_________________________________________________________________________
    
    IF ( ( ABS( a ) .LE. emach ) .AND. ( ABS( b ) .LE. emach ) ) THEN
       
       ! a and b are zero
       
       equiv_val = .TRUE.
       
       !----------------------------------------------------------------------
       
    ELSE IF ( fuzzy ) THEN
       
       ! relative difference is used for comparison.
       
       IF ( ABS( a ) .LE. emach ) THEN
          
          difference = ABS( b )

       ELSEIF ( ABS( b ) .LE. emach ) THEN
          
          difference = ABS( a  )
          
       ELSE

          difference = ABS( a - b ) / MIN( ABS( a ), ABS( b ) )

       END IF

       IF ( difference .LT. deviation / 100.0D0 ) equiv_val = .TRUE.
       
       !----------------------------------------------------------------------
       
    ELSE
       
       ! absolute difference is used for comparison.
       
       IF ( ABS( a - b ) .LE. prec ) equiv_val = .TRUE.
       
       !----------------------------------------------------------------------
       
    END IF


  END FUNCTION equiv_val

  !=========================================================================

  FUNCTION equiv_vec( a, b, deviation, fuzzy )

    ! input arguments :

    REAL( dp ), DIMENSION( : ), INTENT( IN ) :: a
    REAL( dp ),                 INTENT( IN ) :: b, deviation
    LOGICAL :: fuzzy

    ! output result :

    LOGICAL, DIMENSION( SIZE( a ) ) :: equiv_vec

    ! local variables :

    REAL( dp ), DIMENSION( SIZE( a ) ) :: difference
    
    equiv_vec = .FALSE.
       
    ! relative difference is used for comparison.

    IF ( fuzzy ) THEN

       !----------------------------------------------------------------------
       
       WHERE ( ABS( a ) .LE. emach ) difference = ABS( b )

       IF ( ABS( b ) .LE. emach ) THEN
          
          WHERE ( ABS( a ) .GT. emach ) difference = ABS( a )
          
       ELSE
          
          WHERE ( ABS( a ) .GT. emach ) &
               difference = ABS( a - b ) / MIN( ABS( a ), ABS( b ) )
          
       END IF
       
       WHERE ( difference .LT. deviation / 100.0D0 ) equiv_vec = .TRUE.
       
       !----------------------------------------------------------------------
       
    ELSE

       !----------------------------------------------------------------------
      
       ! absolute difference is used for comparison.
       
       WHERE ( ABS( a - b ) .LE. prec ) equiv_vec = .TRUE.
       
       !----------------------------------------------------------------------
       
    END IF

  END FUNCTION equiv_vec
  
  !=========================================================================
  
end module precision
