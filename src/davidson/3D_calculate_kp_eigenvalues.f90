!=============================================================================
!
! Module mod_calculate_kp_eigenvalues3D
!
!   Subroutine calculate_kp_eigenvalues
!
!=============================================================================

MODULE mod_calculate_kp_eigenvalues3D
  
  !===========================================================================
  
  IMPLICIT NONE
  
  !===========================================================================
  
CONTAINS
  
  !===========================================================================
  
  SUBROUTINE calculate_kp_eigenvalues( dim_qr, num_qr, num_states, &
       E_separate, chargeC, kind_kpC, eigenvalueV, spinorM )
  
    !=========================================================================
    !
    ! Calculates k.p eigenstates
    !
    ! dim_qr ... dimension of quantum region
    ! num_qr ... number of quantum region
    !
    ! chargeC  = 'el' or 'hl'
    ! kind_kpC = '8x8' or '6x6'
    ! 
    ! kind_kpC = '8x8' and chargeC = 'el' : The num_states from E_separate
    !                                       upwards are calculated.
    !
    ! kind_kpC = '8x8' and chargeC = 'hl' : The num_states from E_separate
    !                                       downwards are calculated.
    !
    ! kind_kpC = '6x6' and chargeC = 'hl' : The num_states maximum states
    !                                       are calculated.
    !
    !_________________________________________________________________________
    !
    ! Uses method specified in control_numeric: SchroedingerkpEvSolv3D:
    !
    ! - 'it_jam'     uses diagonalization routine of Jacek:
    !
    !       init      = init_itjam
    !       itere_max = itere_stagesV(stage)
    !       eps_itere = epsilon_it_jam
    !       itsub_max = itsub_stagesV(stage)
    !       eps_itsub = epsilon_it_jam
    !       iwrite    = iwtite_it_jam
    !
    !       stage and epsilon_it_jam must be defined in the upper level:
    !       --> level of precision
    !
    !       (thus far only for 6x6 - holes)
    !
    !
    ! - 'arpack'     Arnoldi:
    !
    !       precision = SchroedingerkpResidual3D
    !       
    !       SchroedingerkpIterations3D --> Number of iterations set
    !       to a default value, does not change during simulation.
    !                           
    !       SchroedingerkpResidual3D --> Precision for arnoldi method
    !       is defined according to precision of outer loop.
    !
    !-------------------------------------------------------------------------
    !
    ! Note:
    !
    ! For extremous eigenvalues --> Arnoldi with option 'LR' or 'SR'
    !
    ! For inner eigenvalues --> Arnoldi with option 'SM'
    !
    !-------------------------------------------------------------------------
    !       
    ! NOTE: 
    !
    ! All eigenvalues around E_separate are calculated due to :
    !
    !                    H --> H-E_separate[eV].
    !
    ! This shift is carried out here explicitly.
    !             
    ! For electrons, however, one only needs eigenvalues above this
    ! edge (holes: below) 
    ! 
    ! --> Number of eigenvalues is multiplied with a factor 
    !     mult_num_ev_arnoldi ( module numeric_control3D )
    !     and only the relevant eigenvalues are taken .
    !                 
    !     mult_num_ev_arnoldi( module numeric_control3D ) here is 
    !     just a default value, it is specified in each k.p region.
    !                  
    ! If the calculated eigenvalues are not sufficient, mult_num_ev_arnoldi
    ! ( in the kp-region ) is increased by add_num_ev_arnoldi
    ! ( module numeric_control3D), otherwise decreased by dec_num_ev_arnoldi
    ! ( module numeric_control3D ).
    !
    ! If the calculated eigenvalues are not sufficient, the calculation is
    ! started anew and the specified number of eigenvalues is calculated ...
    !
    !-------------------------------------------------------------------------
    !
    ! note:
    !
    ! arnoldi seems to have problems with degeneracies:
    ! -->  disturb_arnoldi is added to diagonal elements of the upper block
    !      in order to destroy the degeneracies
    !                
    !=========================================================================

    USE numeric_control3D,           ONLY : stage, itere_stagesV, &
         itsub_stagesV, epsilon_it_jam, init_itjam, iwtite_it_jam, &
         add_num_ev_arnoldi, dec_num_ev_arnoldi
   
    USE sparse_matrix,               ONLY : ija,sa
    USE normalization_factors_sg3D,  ONLY : diag_kpV
    USE gitter3D,                    ONLY : norm_sgV
    USE mod_eigensolver_hermitian,   ONLY : eigensolver_hermitian
    USE eigensolver_iter_routines3D, ONLY : matrix_vector_eigencg, &
         precond_eigencg, inner_product_eigencg
    
    USE arnoldi_complex_routines3D,  ONLY : matvec_arnoldi_complex, &
         lin_equation_solver_arnoldi_c
    
    USE sparse_matrix_rout,          ONLY : deallocate_sparse_matrix
    USE mod_eigensolver_herm,        ONLY : arnoldi_complex
    USE gitter3D,                    ONLY : dim_qr3DV, num_qr_regions3D
    USE quantum_solutions3D,         ONLY : kp_el_3DM, kp_hl_3DM
    USE limits,                      ONLY : min_precision_kp, &
         min_iter_itere_kp, min_iter_itsub_kp
    
    USE mod_eigensolver_davidson,    ONLY : eigensolver_davidson
    USE mod_precond_ilu_davidson,    ONLY : transform_sparse_to_ilu3D, &
         E_sigma, setup_david_precond_jac, dealloc_david_precond_jac
    
    USE control_numeric,             ONLY : SchroedingerkpEvSolv3D, &
         SchroedingerkpIterations3D, SchroedingerkpResidual3D

    !=========================================================================

    ! Worked previously only for WNT (Windows) 
    ! Now it also should work on Unix.
    
    USE mod_test_time, ONLY : CallRoutineFirstTime
    USE test_time    , ONLY : TestTime
    
    ! Worked previously only for WNT (Windows) 
    ! Now it also should work on Unix.
    
    !=========================================================================

    IMPLICIT NONE

    !=========================================================================

    REAL    ( KIND(0.0D0) ), DIMENSION( : ),       POINTER :: eigenvalueV
    COMPLEX ( KIND(0.0D0) ), DIMENSION( :, :, : ), POINTER :: spinorM
    
    INTEGER,                   INTENT( IN ) :: dim_qr, num_states, num_qr
    REAL      ( KIND(0.0D0) ), INTENT( IN ) :: E_separate
    CHARACTER ( LEN = 3 ),     INTENT( IN ) :: kind_kpC
    CHARACTER ( LEN = 2 ),     INTENT( IN ) :: chargeC
    
    !_________________________________________________________________________

    COMPLEX ( KIND(0.0D0) ), DIMENSION( :, : ), ALLOCATABLE :: hambandM
    COMPLEX ( KIND(0.0D0) ), DIMENSION( :, : ), ALLOCATABLE :: eigenfunctionsM
    COMPLEX ( KIND(0.0D0) ), DIMENSION( : ),    ALLOCATABLE :: dcwork

    REAL    ( KIND(0.0D0) ), DIMENSION( : ),    ALLOCATABLE :: dpwork, en
    REAL    ( KIND(0.0D0) ), DIMENSION( : ),    ALLOCATABLE :: eigenvalV

    COMPLEX ( KIND(0.0D0) ), DIMENSION( :, : ), POINTER :: &
         eigenfunctions_pointM

    REAL    ( KIND(0.0D0) ), DIMENSION( : ),    POINTER :: &
         eigenval_pointV, norm_kpV
    
    COMPLEX   ( KIND(0.0D0) ) :: zw, sigma

    REAL      ( KIND(0.0D0) ) :: norm, eps_itere, eps_itsub
    REAL      ( KIND(0.0D0) ) :: tolerance, resid_lin_eq, percentage
    REAL      ( KIND(0.0D0) ) :: mult_num_ev_arnoldi, shift, E_separate2

    CHARACTER ( LEN = 2 )     :: which_eigvalC

    INTEGER   ( 4 )           :: lwork, info, maxdim, istatus
    INTEGER   ( 4 )           :: n, i, j, k, count, n1, l

    INTEGER :: max_eigval, max_iterations, iproblem
    INTEGER :: num_ev, init, iwrite
    INTEGER :: itere_max, itsub_max, info_itsub, info_itere
    INTEGER :: iqr, n_converged
  
    LOGICAL :: eigenvectorL

    !=========================================================================

    IF ( ( kind_kpC .NE. '8x8' ) .AND. ( kind_kpC .NE. '6x6' ) ) STOP &
         'Error calculate_kp_eigenvalues: kind_kpC not allowed.'

    IF ( dim_qr3DV( num_qr ) .NE. dim_qr ) STOP &
         'Error calculate_kp_eigenvalues: Inconsistency with quantum region.'

    IF ( ( num_qr .LT. 1 ) .OR. ( num_qr .GT. num_qr_regions3D ) ) STOP &
         'Error calculate_kp_eigenvalues: num_qr out of range.'

    PRINT *, 'kp-in', kind_kpC

    !=========================================================================
    
    SELECT CASE ( kind_kpC )
       
       CASE ( '8x8' )
       k = 8
       
       CASE ( '6x6' )
       k = 6
       
    END SELECT

    !=========================================================================

    SELECT CASE( SchroedingerkpEvSolv3D )
       
       CASE ( 'it_jam', 'arpack', 'davids' )
       IF ( dim_qr * k .NE. ija(1) - 2 ) STOP &
            'Error calculate_kp_eigenvalues: Grid and matrix inconsistent.'
       
    CASE DEFAULT
       STOP 'Error calculate_kp_eigenvalues: &
            &SchroedingerkpEvSolv3D in module control_numeric ill-defined.'
       
    END SELECT
    
    PRINT *, SchroedingerkpEvSolv3D

    !=========================================================================
    !
    ! Setup normalization factor norm_kpV ( same dimension as eigenvectors ),
    ! it can be retrieved from norm_sgV as done below.
    !
    ! Note: In 1D this normalzation factor is already defined in 
    !       subroutine setup_kp_mat1D.
    !
    ! This normalization factor is used to normalize the eigenstates.
    !
    !=========================================================================

    ALLOCATE( norm_kpV( k * dim_qr ), STAT = istatus )
    IF ( istatus .NE. 0 ) STOP &
         'Error calculate_kp_eigenvalues: Allocation failed.'

    DO iqr = 1, dim_qr
       norm_kpV( (iqr-1) * k + 1 : iqr * k ) = norm_sgV( num_qr, iqr )
    END DO

    !=========================================================================

    CallRoutineFirstTime = .TRUE.
    CALL testtime( CallRoutineFirstTime, .FALSE. )

    SELECT CASE( SchroedingerkpEvSolv3D )

       !**********************************************************************
       
       ! Eigenvalue solver of J. Majewski
       
       CASE ( 'it_jam' )

       !**********************************************************************

       IF ( ( kind_kpC .NE. '6x6' ) .AND. ( chargeC .NE. 'hl' ) ) STOP &
            'Error calculate_kp_eigenvalues: it_jam only for 6x6 k.p holes.'

       num_ev = num_states
       n      = dim_qr * k
       
       ALLOCATE( eigenvalV(num_ev), eigenfunctionsM(n,num_ev), STAT = istatus )
       IF ( istatus .NE. 0 ) STOP &
            'Error calculate_kp_eigenvalues: Allocation failed.'
       
       !-------- initial values -----
       
       eigenvalV = 0.0D0
       
       DO i = 1, num_states
          
          DO l = 1, dim_qr
             
             eigenfunctionsM( ( l - 1 ) * 6 + 1 : l * 6, i ) = &
                  spinorM( i, l, 3 : 8 )
             
          END DO
          
       END DO
              
       init      = init_itjam
       
       ! itere_max = MAX( itere_stagesV(stage), min_iter_itere_kp )
       ! can be specified via input file
       
       itere_max = SchroedingerkpIterations3D
       
       ! eps_itere = MIN( epsilon_it_jam, min_precision_kp )
       ! can be specified via input file
       
       eps_itere = SchroedingerkpResidual3D
       
       ! itsub_max = MAX( itsub_stagesV(stage), min_iter_itsub_kp )
       ! can be specified via input file
       
       itsub_max = SchroedingerkpIterations3D
       
       ! eps_itsub = MIN( epsilon_it_jam, min_precision_kp )
       ! can be specified via input file

       eps_itere = SchroedingerkpResidual3D
       
       iwrite    = iwtite_it_jam

       ! itere_max = 30
       ! eps_itere = 1.0d-10
       ! itsub_max = 200
       ! eps_itsub = 1.0d-10
       
       !______________________________________________________________________
       
       SELECT CASE( kind_kpC )
       
          !-------------------------------------------------------------------

          CASE ( '8x8' )
          
          DO i = 1, num_ev
             DO l = 1, dim_qr
                eigenfunctionsM( ( l - 1 ) * 8 + 1 : l * 8, i ) = &
                     spinorM( i, l, 1 : 8 )
             END DO
          END DO
          
          !-------------------------------------------------------------------

          CASE ( '6x6' )
          
          DO i = 1, num_ev
             DO l = 1, dim_qr
                eigenfunctionsM( ( l - 1 ) * 6 + 1 : l * 6, i ) = &
                     spinorM( i, l, 3 : 8 ) 
             END DO
          END DO

          !-------------------------------------------------------------------
       
       CASE DEFAULT
       
          STOP 'Error calculate_kp_eigenvalues: For holes only 6x6 or 8x8 k.p.'
       
          !-------------------------------------------------------------------

       END SELECT

       !______________________________________________________________________

       WRITE (*,*) 'iterative_kp in , ndim= ', n
       WRITE (*,*) 'itere_max, itsub_max, eps_itere, eps_itsub=', &
            itere_max, itsub_max, eps_itere, eps_itsub

       CALL eigensolver_hermitian( n, num_ev, eigenvalV, eigenfunctionsM, &
            init, itere_max, itsub_max, eps_itere, eps_itsub, &
            info_itere, info_itsub, iwrite, matrix_vector_eigencg, &
            precond_eigencg, inner_product_eigencg )
       
       !------------------------------------------------------
       ! Note: This iterative solver searches for the minimum,
       !       for holes however one wants the maximum
       !       So the transformation --> (-1) is done
       !       in every matrix vector multiplication and
       !       the eigenvalues are backtransformed afterwards:
       !       eigenvalV = -eigenvalV
       !
       !       Now, they are in the right order.
       !------------------------------------------------------
       
       eigenvalV = -eigenvalV

       WRITE (*,*) 'Eigenvalues: ', eigenvalV

       !______________________________________________________________________

       SELECT CASE( chargeC )
       
          !-------------------------------------------------------------------

          CASE ( 'el' )

          STOP 'Error calculate_kp_eigenvalues: &
               &Iterative for electrons not yet implemented.'

          !-------------------------------------------------------------------

          ! Holes and backtransform and normalize all.

          CASE ( 'hl' )

          !-------------------------------------------------------------------
          
          count       = 0
          eigenvalueV = 0.0D0
          spinorM     = CMPLX( 0.0D0 )

          DO  i = 1, num_ev
             
             !................................................................
             
             eigenvalueV( i ) = eigenvalV( i )

             DO n1 = 1, dim_qr * 6
                eigenfunctionsM( n1, i ) = &
                     diag_kpV( n1 ) * eigenfunctionsM( n1, i )
             END DO

             norm = SUM( CONJG( eigenfunctionsM( :, i ) ) * &
                  eigenfunctionsM( :, i ) * norm_kpV(:) )
             
             DO n1 = 1, dim_qr * 6
                eigenfunctionsM( n1, i ) = &
                     eigenfunctionsM( n1, i ) / SQRT( norm )
             END DO

             !................................................................

             SELECT CASE( kind_kpC )

                CASE ( '8x8' )

                DO l = 1, dim_qr
                   spinorM( i, l, 1 : 8 ) = &
                        eigenfunctionsM( ( l - 1 ) * 8 + 1 : k * 8, i )
                END DO
                
                CASE ( '6x6' )
                
                DO l = 1, dim_qr
                   spinorM( i, l, 3 : 8 ) = &
                        eigenfunctionsM( ( l - 1 ) * 6 + 1 : l * 6, i )
                 END DO
                 
             CASE DEFAULT
                
                STOP 'Error calculate_kp_eigenvalues: &
                     &For holes only 6x6 or 8x8 k.p!'
                
             END SELECT

             !................................................................

          END DO

          !-------------------------------------------------------------------

       CASE DEFAULT

          STOP 'Error calculate_kp_eigenvalues: chargeC ill-defined.'
       
          !-------------------------------------------------------------------

       END SELECT

       !______________________________________________________________________

       DEALLOCATE( eigenvalV, eigenfunctionsM )

       CALL deallocate_sparse_matrix

       !**********************************************************************
       
       ! Arnoldi
       
       CASE ( 'arpack' )

       !**********************************************************************

       SELECT CASE ( chargeC )

          !___________________________________________________________________

          CASE ( 'hl' )

          !___________________________________________________________________

          SELECT CASE ( kind_kpC )
          
             !----------------------------------------------------------------

             CASE ( '6x6' )

             !----------------------------------------------------------------

             n      = dim_qr * k
             num_ev = num_states
             
             ALLOCATE( eigenval_pointV(num_ev), &
                  eigenfunctions_pointM(n,num_ev), STAT = istatus )
             IF ( istatus .NE. 0 ) STOP &
                  'Error calculate_kp_eigenvalues: Allocation failed.'
             
             !-------- initial values -----
             
             eigenval_pointV = 0.0D0
             
             DO i = 1, num_states
                DO l = 1, dim_qr
                   eigenfunctions_pointM( ( l - 1 ) * 6 + 1 : l * 6, i ) = &
                        spinorM( i, l, 3 : 8 )
                END DO
             END DO
             
             eigenfunctions_pointM( 1:n, num_states+1:num_ev ) = CMPLX( 0.0D0 )
             
             !----------------------------------------------------------------

             eigenvectorL   = .TRUE.
             which_eigvalC  = 'LR'
             max_eigval     = num_states
             tolerance      = MIN( SchroedingerkpResidual3D, min_precision_kp )
             max_iterations =     SchroedingerkpIterations3D
             sigma          = 0.0D0       ! not used
             resid_lin_eq   = 1.0D-10     ! not used
             iproblem       = 1

             !----------------------------------------------------------------

             WRITE (*,*) 'Arnoldi-in', num_states
             WRITE (*,*) 'Tolerance, Max no. iterations, E_separate: ', &
                  tolerance, max_iterations, E_separate

             CALL arnoldi_complex( n, eigenvectorL, which_eigvalC, &
                  max_eigval, tolerance, max_iterations, sigma, &
                  resid_lin_eq, iproblem, eigenval_pointV, &
                  eigenfunctions_pointM, lin_equation_solver_arnoldi_c, &
                  matvec_arnoldi_complex )

             WRITE (*,*) 'Arnoldi-out'
             WRITE (*,*) 'Eigenvalues: ', eigenval_pointV

             !----------------------------------------------------------------

             count       = 0
             eigenvalueV = 0.0D0
             spinorM     = CMPLX( 0.0D0 )

             DO  i = 1, num_ev

                !.............................................................

                eigenvalueV( i ) = eigenval_pointV( i )

                DO n1 = 1, dim_qr * 6
                   eigenfunctions_pointM( n1, i ) = &
                        diag_kpV( n1 ) * eigenfunctions_pointM( n1, i )
                END DO

                norm = SUM( CONJG( eigenfunctions_pointM( :, i ) ) * &
                     eigenfunctions_pointM( :, i ) * norm_kpV(:) )

                IF ( norm .LE. 0.0D0) STOP &
                     'Error calculate_kp_eigenvalues2D: Eigenfunction zero.'
                  
                DO n1 = 1, dim_qr * 6
                   eigenfunctions_pointM( n1, i ) = &
                        eigenfunctions_pointM( n1, i ) / SQRT( norm )
                END DO

                !.............................................................

                SELECT CASE( kind_kpC )

                   CASE ( '8x8' )

                   DO l = 1, dim_qr
                      spinorM( i, l, 1 : 8 ) = &
                           eigenfunctions_pointM( ( l-1) * 8 + 1 : l * 8, i )
                   END DO

                   CASE ( '6x6' )
                   
                   DO l = 1, dim_qr
                      spinorM( i, l, 3 : 8 ) = &
                           eigenfunctions_pointM( (l-1) * 6 + 1 : l * 6, i )
                   END DO

                CASE DEFAULT
                   
                   STOP 'Error calculate_kp_eigenvalues: &
                        &For holes only 6x6 or 8x8 k.p!'
                  
                END SELECT

                !.............................................................

             END DO
             
             !----------------------------------------------------------------

             DEALLOCATE( eigenval_pointV, eigenfunctions_pointM )
             CALL deallocate_sparse_matrix

             !----------------------------------------------------------------

             CASE ( '8x8' )

             !----------------------------------------------------------------
             !
             ! H --> H - E_separate
             !
             ! Note:       H = diag_kpV * sparse_mat * diag_kpV
             !
             ! diagonal elements of sparse matrix 
             !                        --> - E_separate / ( diag_kpV**2 )
             !
             !----------------------------------------------------------------

             n           = dim_qr * k
             sa( 1 : n ) = sa( 1 : n ) - E_separate / ( diag_kpV(1:n)**2 )

             !----------------------------------------------------------------
             ! Note:
             !
             ! num_states eigenvalues above E_separate are requested.
             ! So one calculates num_ev = num_states*mult_num_ev_arnoldi.
             ! Then they are sorted and the lowest num_states eigenvalues
             ! above E_separate are taken.
             ! If they are not sufficient, 
             ! mult_num_ev_arnoldi -> mult_num_ev_arnoldi + add_num_ev_arnoldi
             ! and the calculation is started anew.
             !
             ! If more are calculated than necessary:
             ! mult_num_ev_arnoldi = mult_num_ev_arnoldi - dec_num_ev_arnoldi
             !       
             !----------------------------------------------------------------

10           CONTINUE

             mult_num_ev_arnoldi = kp_hl_3DM( num_qr )%mult_num_ev_arnoldi
             num_ev              = num_states * mult_num_ev_arnoldi

             ALLOCATE( eigenval_pointV(num_ev), &
                  eigenfunctions_pointM(n,num_ev), STAT = istatus )
             IF ( istatus .NE. 0 ) STOP &
                  'Error calculate_kp_eigenvalues: Allocation failed.'

             !-------- initial values -----

             eigenval_pointV = 0.0D0
             
             DO i = 1, num_states
                DO l = 1, dim_qr
                   eigenfunctions_pointM( ( l - 1 ) * 8 + 1 : l * 8, i ) = &
                        spinorM( i, l, 1 : 8 )
                END DO
             END DO
             
             eigenfunctions_pointM( 1:n, num_states+1:num_ev ) = CMPLX(0.0D0)

             !----------------------------------------------------------------

             eigenvectorL   = .TRUE.
             which_eigvalC  = 'SM'
             max_eigval     = num_ev
             tolerance      = MIN( SchroedingerkpResidual3D, min_precision_kp )
             max_iterations = SchroedingerkpIterations3D
             sigma          = 0.0D0      ! not used
             resid_lin_eq   = 1.0D-10    ! not used
             iproblem       = 1

             WRITE (*,*) 'Arnoldi-in', num_states
             WRITE (*,*) 'Tolerance, Max no. iterations,E_separate: ', &
                  tolerance, max_iterations, E_separate

             CALL arnoldi_complex( n, eigenvectorL, which_eigvalC, &
                  max_eigval, tolerance, max_iterations, sigma, &
                  resid_lin_eq, iproblem, eigenval_pointV, &
                  eigenfunctions_pointM, lin_equation_solver_arnoldi_c, &
                  matvec_arnoldi_complex )

             WRITE (*,*) 'Arnoldi-out'

             !-------------  sort eigenvalues --------------------------------

             !  arnoldi_complex provides eigenvalues in descending order.
             
             !----------------------------------------------------------------

             DO i = 2, num_ev
                
                IF ( eigenval_pointV(i) .GT. eigenval_pointV(i-1) ) STOP &
                     'Error calculate_kp_eigenvalues3D:&
                     & Eigenvalues provided in wrong order by arnoldi_complex.'
             
             END DO

             !----------------------------------------------------------------

             eigenval_pointV = eigenval_pointV + E_separate

             WRITE (*,*) 'Eigenvalues', eigenval_pointV

             !----------------------------------------------------------------

             count       = 0
             eigenvalueV = 0.0D0
             spinorM     = CMPLX(0.0D0)

             DO  i = 1, num_ev

                IF ( ( eigenval_pointV(i) .LT. E_separate ) .AND. &
                     ( count .LT. num_states ) ) THEN

                   !..........................................................
                   
                   count = count + 1

                   eigenvalueV( count ) = eigenval_pointV( i )

                   DO n1 = 1, dim_qr * 8
                      eigenfunctions_pointM( n1, i ) = &
                           diag_kpV( n1 ) * eigenfunctions_pointM( n1, i )
                   END DO
                   
                   norm = SUM( CONJG( eigenfunctions_pointM( :, i ) ) * &
                        eigenfunctions_pointM( :, i ) * norm_kpV(:) )
                   
                   IF ( norm .LE. 0.0D0) STOP &
                        'Error calculate_kp_eigenvalues3D: Eigenfunction zero.'

                   DO n1 = 1, dim_qr * 8
                      eigenfunctions_pointM( n1, i ) = &
                           eigenfunctions_pointM( n1, i ) / SQRT( norm )
                   END DO

                   !..........................................................

                   SELECT CASE ( kind_kpC )
                      
                      CASE ( '8x8' )

                      DO l = 1, dim_qr
                         spinorM( count, l, 1 : 8 ) = &
                              eigenfunctions_pointM( (l-1) * 8 + 1 : l * 8, i )
                      END DO
                      
                      CASE ( '6x6' )

                      DO l = 1, dim_qr
                         spinorM( count, l, 3 : 8 ) = &
                              eigenfunctions_pointM( (l-1) * 6 + 1 : l * 6, i )
                      END DO
                      
                   CASE DEFAULT
                      
                      STOP 'Error calculate_kp_eigenvalues: &
                           &For holes only 6x6 or 8x8 k.p!'
                      
                   END SELECT

                   !..........................................................

                ELSE IF ( ( eigenval_pointV(i) .LT. E_separate ) .AND. &
                     ( count .GE. num_states ) ) THEN

                   !..........................................................
                   
                   count = count + 1

                   !..........................................................
                   
                END IF

             END DO

             DEALLOCATE( eigenval_pointV, eigenfunctions_pointM )

             !---------- check if enough eigenvalues have been calculated -----

             !IF (count<num_states) THEN
             
             ! WRITE(*,*) 'calculate_kp_eigenvalues3D: &
             ! &Number of calculated eigenvalues &
             ! & not sufficient.'
             
             ! WRITE(*,*) ' calculate anew'
             
             !mult_num_ev_arnoldi = mult_num_ev_arnoldi + add_num_ev_arnoldi
             !GOTO 10
             
             !ELSE IF (count>num_states) THEN
             
             percentage = DBLE( count ) / num_ev
             mult_num_ev_arnoldi= MAX( ( 1.0D0 / percentage ), &
                  mult_num_ev_arnoldi - dec_num_ev_arnoldi + 0.1D0 )
             
             kp_hl_3DM( num_qr )%mult_num_ev_arnoldi = mult_num_ev_arnoldi

             WRITE (*,*) 'Number of eigenvalues to be calculated reduced.'
             WRITE (*,*) 'New number of eigenvalues: ', mult_num_ev_arnoldi
             
             !END IF

             !-------- end of eigenvalue calculation -------------------------

             CALL deallocate_sparse_matrix

             !----------------------------------------------------------------

          END SELECT

          !___________________________________________________________________
          
          CASE ( 'el' )

          !___________________________________________________________________

          SELECT CASE( kind_kpC )

             !----------------------------------------------------------------

             CASE ( '6x6' )

             !----------------------------------------------------------------

             STOP 'Error mod_calculate_kp_eigenvalues3D: &
                  &Electrons within 6x6 k.p not allowed.'

             CASE ( '8x8' )

             !----------------------------------------------------------------
             !
             ! H --> H-E_separate
             !
             ! Note:       H = diag_kpV * sparse_mat * diag_kpV
             !
             ! diagonal elements of sparse matrix --> -E_separate/(diag_kpV**2)
             !
             !----------------------------------------------------------------

             n       = dim_qr * k
             
             DO n1 = 1, n
                sa(n1) = sa(n1) - E_separate / ( diag_kpV(n1)**2 )
             END DO

             !----------------------------------------------------------------
             !
             ! Note:
             !
             ! num_states eigenvalues above E_separate are requested.
             ! So one calculates num_ev = num_states*mult_num_ev_arnoldi.
             ! Then they are sorted and the lowest num_states eigenvalues
             ! above E_separate are taken.
             ! If they are not sufficient, 
             !
             ! mult_num_ev_arnoldi-> mult_num_ev_arnoldi + add_num_ev_arnoldi
             !
             ! and the calculation is started anew.
             !
             ! If more are calculated than necessary:
             ! 
             !  mult_num_ev_arnoldi = mult_num_ev_arnoldi - dec_num_ev_arnoldi
             !       
             !----------------------------------------------------------------

20           CONTINUE

             mult_num_ev_arnoldi = kp_el_3DM( num_qr )%mult_num_ev_arnoldi
             num_ev              = num_states * mult_num_ev_arnoldi

             ALLOCATE( eigenval_pointV(num_ev), &
                  eigenfunctions_pointM(n,num_ev), STAT = istatus )
             IF ( istatus .NE. 0 ) STOP &
                  'Error calculate_kp_eigenvalues: Allocation failed.'

             !-------- initial values -----

             eigenval_pointV = 0.0D0
             
             DO i = 1, num_states
                DO l = 1, dim_qr
                   eigenfunctions_pointM( (l-1) * 8 + 1 : l * 8, i ) = &
                        spinorM( i, l, 1 : 8 )
                END DO
             END DO
             
             eigenfunctions_pointM( 1:n, num_states+1:num_ev ) = CMPLX(0.0D0)

             !-----------------------------

             eigenvectorL   = .TRUE.
             which_eigvalC  = 'SM'
             max_eigval     = num_ev
             tolerance      = MIN( SchroedingerkpResidual3D, min_precision_kp )
             max_iterations = SchroedingerkpIterations3D
             sigma          = 0.0D0       ! not used
             resid_lin_eq   = 1.0D-10     ! not used
             iproblem       = 1

             WRITE (*,*) 'Arnoldi-in',num_states,num_ev
             WRITE (*,*) 'Tolerance, Max no. iterations, E_separate: ', &
                  tolerance, max_iterations, E_separate

             CALL arnoldi_complex( n, eigenvectorL, which_eigvalC, &
                  max_eigval, tolerance, max_iterations, sigma, &
                  resid_lin_eq, iproblem, eigenval_pointV, &
                  eigenfunctions_pointM, lin_equation_solver_arnoldi_c, &
                  matvec_arnoldi_complex )

             !-------------  sort eigenvalues --------------------------------

             !  arnoldi_complex provides eigenvalues in descending order.
             
             !----------------------------------------------------------------

             DO i = 2, num_ev
                
                IF ( eigenval_pointV(i) .GT. eigenval_pointV(i-1) ) STOP &
                     'Error calculate_kp_eigenvalues3D: &
                     & Eigenvalues provided in wrong order by arnoldi_complex.'
                
             END DO

             !----------------------------------------------------------------

             eigenval_pointV = eigenval_pointV + E_separate

             WRITE (*,*) 'Eigenvalues: ', eigenval_pointV

             !----------------------------------------------------------------

             count       = 0
             eigenvalueV = 0.0D0
             spinorM     = CMPLX(0.0D0)

             DO  i = num_ev, 1, -1

                IF ( ( eigenval_pointV(i) .GT. E_separate ) .AND. &
                     ( count .LT. num_states ) ) THEN
                   
                   !..........................................................

                   count = count + 1
                   eigenvalueV( count ) = eigenval_pointV( i )
                   
                   DO n1 = 1, dim_qr * 8
                      eigenfunctions_pointM( n1, i ) = &
                           diag_kpV( n1 ) * eigenfunctions_pointM( n1, i )
                   END DO

                   norm = SUM( CONJG( eigenfunctions_pointM( :, i ) ) * &
                        eigenfunctions_pointM( :, i ) * norm_kpV(:) )
                   
                   IF ( norm .LE. 0.0D0 ) STOP &
                        'Error calculate_kp_eigenvalues3D: Eigenfunction zero.'

                   DO n1 = 1, dim_qr * 8
                      eigenfunctions_pointM( n1, i ) = &
                           eigenfunctions_pointM( n1, i ) / SQRT( norm )
                   END DO
                   
                   !..........................................................

                   SELECT CASE ( kind_kpC )
                      
                      CASE ( '8x8' )
                      
                      DO l = 1, dim_qr
                         spinorM( count, l, 1 : 8 ) = &
                              eigenfunctions_pointM( (l-1) * 8 + 1 : l * 8, i )
                      END DO
                      
                      CASE ( '6x6' )
                      
                      DO l = 1, dim_qr
                         spinorM( count, l, 3 : 8 ) = &
                              eigenfunctions_pointM( (l-1) * 6 + 1 : l * 6, i )
                      END DO
                      
                   CASE DEFAULT
                      
                      STOP 'Error calculate_kp_eigenvalues: &
                           &For holes only 6x6 or 8x8 k.p!'
                      
                   END SELECT

                   !..........................................................

                ELSE IF ( ( eigenval_pointV(i) .GT. E_separate ) .AND. &
                     ( count .GE. num_states ) ) THEN
                   
                   !..........................................................

                   count = count + 1

                   !..........................................................

                END IF

             END DO

             DEALLOCATE( eigenval_pointV, eigenfunctions_pointM )
             
             !---------- Check if enough eigenvalues have been calculated. ---
             
             !IF (count<num_states) THEN
             
             !  WRITE(*,*) 'calculate_kp_eigenvalues3D: &
             ! &Number of calculated eigenvalues &
             ! & not sufficient.'
             
             !  WRITE(*,*) 'Calculate anew.'
             
             !  mult_num_ev_arnoldi =  mult_num_ev_arnoldi + add_num_ev_arnoldi
             !  GOTO 20
             
             !ELSE IF (count>num_states)THEN
             
             percentage = DBLE( count ) / num_ev
             mult_num_ev_arnoldi = MAX( ( 1.0D0 / percentage ), &
                  mult_num_ev_arnoldi - dec_num_ev_arnoldi + 0.1D0 )
             kp_el_3DM( num_qr )%mult_num_ev_arnoldi = mult_num_ev_arnoldi

             WRITE (*,*) 'Number of eigenvalues to be calculated reduced.'
             WRITE (*,*) 'New number of eigenvalues ', mult_num_ev_arnoldi

             !END IF

             !-------- end of eigenvalue calculation -------------------------

             CALL deallocate_sparse_matrix

             !----------------------------------------------------------------

          END SELECT

          !___________________________________________________________________

       END SELECT

       !**********************************************************************
      
       ! Jacobi-Davidson
       
       CASE ( 'davids' )
       
       !**********************************************************************

       SELECT CASE ( chargeC )
                
          !___________________________________________________________________

          CASE ( 'hl' )

          !___________________________________________________________________
                
          SELECT CASE ( kind_kpC )
          
             !----------------------------------------------------------------
             
             CASE ( '6x6' )

             !----------------------------------------------------------------

             n      = dim_qr * k
             num_ev = num_states

             ALLOCATE( eigenval_pointV(num_ev), &
                  eigenfunctions_pointM(n,num_ev), STAT = istatus )
             IF ( istatus .NE. 0 ) STOP &
                  'Error calculate_kp_eigenvalues: Allocation failed.'
             
             eigenval_pointV = 0.0D0
                   
             !-------- initial values -----

             DO i = 1, num_states
                DO l = 1, dim_qr
                   eigenfunctions_pointM( ( l - 1 ) * 6 + 1 : l * 6, i ) = &
                        spinorM( i, l, 3 : 8 ) 
                END DO
             END DO

             eigenfunctions_pointM( 1:n, num_states+1:num_ev ) = CMPLX(0.0Dd0)
                   
             !----------------------------------------------------------------

             tolerance      = MIN( SchroedingerkpResidual3D, min_precision_kp )
             max_iterations = SchroedingerkpIterations3D
             sigma          = CMPLX( E_separate )
             E_sigma        = sigma

             CALL transform_sparse_to_ilu3D( 'kp' )
             CALL setup_david_precond_jac( sigma )
              
             WRITE (*,*) 'Jacobi-Davidson: No. of eigenvalues: ', num_states
             WRITE (*,*) 'Tolerance, Max. no. of iterations, E_separate: ', &
                  tolerance, max_iterations, E_separate
             
             CALL eigensolver_davidson( n, num_ev, eigenval_pointV, &
                  eigenfunctions_pointM, sigma, n, max_iterations, tolerance )

             WRITE (*,*) 'Jacobi-Davidson eigenvalue solver finished:'
             
             CALL dealloc_david_precond_jac
             
             WRITE(*,*) 'Eigenvalues: ',eigenval_pointV
             
             count       = 0
             eigenvalueV = 0.0D0 
             spinorM     = CMPLX(0.0d0)
             
             DO  i = 1, num_ev

                !.............................................................

                eigenvalueV( i ) = eigenval_pointV( i )
                
                DO n1 = 1, dim_qr * 6
                   eigenfunctions_pointM( n1, i ) = &
                        diag_kpV( n1 ) * eigenfunctions_pointM( n1, i )
                END DO
                
                norm = SUM( CONJG( eigenfunctions_pointM( :, i ) ) * &
                     eigenfunctions_pointM( :, i ) * norm_kpV(:) )
                
                IF ( norm .LE. 0.0D0 ) STOP &
                     'Error calculate_kp_eigenvalues: Eigenfunction zero.'
                
                DO n1 = 1, dim_qr * 6
                   eigenfunctions_pointM( n1, i ) = &
                        eigenfunctions_pointM( n1, i ) / SQRT( norm )
                END DO

                !.............................................................

                SELECT CASE ( kind_kpC )

                   CASE ('6x6')

                   DO l = 1, dim_qr
                      spinorM( i, l, 3 : 8 ) = &
                           eigenfunctions_pointM( (l-1) * 6 + 1 : l * 6, i ) 
                   END DO
                        
                CASE DEFAULT
                   
                   STOP 'Error calculate_kp_eigenvalues: &
                        &For holes only 6x6 k.p or 8x8 k.p!'
                        
                END SELECT

                !.............................................................

             END DO
             
             DEALLOCATE(eigenval_pointV,eigenfunctions_pointM)
             CALL deallocate_sparse_matrix

             !----------------------------------------------------------------

             CASE ( '8x8' )

             !----------------------------------------------------------------
             !
             ! Note:
             !
             ! num_states eigenvalues above E_separate are requested.
             ! So one calculates num_ev = num_states*mult_num_ev_arnoldi.
             ! Then they are sorted and the lowest num_states eigenvalues
             ! above E_separate are taken.
             ! If they are not sufficient, 
             !      
             ! mult_num_ev_arnoldi-> mult_num_ev_arnoldi + add_num_ev_arnoldi
             !
             ! and the calculation is started anew.
             !
             ! If more are calculated than necessary:
             !     
             !  mult_num_ev_arnoldi = mult_num_ev_arnoldi - dec_num_ev_arnoldi
             !       
             !----------------------------------------------------------------
             
15           CONTINUE

             n                   = dim_qr * k
             mult_num_ev_arnoldi = kp_hl_3DM( num_qr )%mult_num_ev_arnoldi
             num_ev              = num_states * mult_num_ev_arnoldi

             ALLOCATE( eigenval_pointV( num_ev ), &
                  eigenfunctions_pointM( n, num_ev ), STAT = istatus )
             IF ( istatus .NE. 0 ) STOP &
                  'Error calculate_kp_eigenvalues: Allocation failed.'
                     
             !-------- initial values -----
             
             DO i = 1, num_states
                DO l = 1, dim_qr
                   eigenfunctions_pointM( ( l - 1 ) * 8 + 1 : l * 8, i ) = &
                        spinorM( i, l, 1 : 8 ) 
                END DO
             END DO

             eigenfunctions_pointM( 1:n, num_states+1:num_ev ) = CMPLX(0.0D0)

             !----------------------------------------------------------------
             
             tolerance      = MIN( SchroedingerkpResidual3D, min_precision_kp )
             max_iterations = SchroedingerkpIterations3D

             CALL transform_sparse_to_ilu3D( 'kp' )

             !-------------------------------------------------------
             ! Note: 
             !
             ! Jacobi-Davidson seems to work better if E_separate is shifted
             ! to -2d0  for holes. This is done here.
             !
             ! transform matrix --> E_separate = -2.0d0
             ! matrix = matrix + shift
             ! E_separate_neu = E_separate + shift =  2.0d0
             !
             !-------------------------------------------------------
             
             WRITE (*,*) 'E_separate', E_separate

             shift       = -E_separate - 2.0D0
             E_separate2 =  E_separate + shift
             sa(1:n)     =  sa(1:n)    + shift

             sigma      = CMPLX( E_separate2 )
             E_sigma    = sigma

             CALL setup_david_precond_jac( sigma )

             WRITE (*,*) 'Jacobi-Davidson-in: ', num_states
             WRITE (*,*) 'Tolerance, Max no. of iterations, E_separate2: ', &
                  tolerance, max_iterations, E_separate2
             
             CALL eigensolver_davidson( n, num_ev, eigenval_pointV, &
                  eigenfunctions_pointM, sigma, n_converged, &
                  max_iterations, tolerance )

             !***  correct for shift ***

             eigenval_pointV = eigenval_pointV - shift
             
             WRITE (*,*) 'Jacobi-Davidson-out'

             CALL dealloc_david_precond_jac

             
             !-------------  sort eigenvalues --------------------------------

             !  Jacobi-Davidson provides eigenvalues in descending order.

             DO i = 2, num_ev
                       
                IF ( eigenval_pointV(i) .GT. eigenval_pointV(i-1) ) STOP &
                     'Error calculate_kp_eigenvalues3D: &
                     & Eigenvalues provided in wrong order by Jacobi-Davidson.'
                     
             END DO

             !----------------------------------------------------------------

             WRITE (*,*) 'Eigenvalues: ', eigenval_pointV
             
             !----------------------------------------------------------------

             count       = 0
             eigenvalueV = 0.0D0
             spinorM     = CMPLX( 0.0D0 )

             DO  i = 1, num_ev

                IF ( ( eigenval_pointV(i) .LT. E_separate ) .AND. &
                     ( count .LT. num_states ) ) THEN

                   !...........................................................

                   count = count + 1
                   eigenvalueV( count ) = eigenval_pointV( i )
                   
                   DO n1 = 1, dim_qr * 8
                      eigenfunctions_pointM( n1, i ) = &
                           diag_kpV( n1 ) * eigenfunctions_pointM( n1, i )
                   END DO
                   
                   norm = SUM( CONJG( eigenfunctions_pointM( :, i ) ) * &
                        eigenfunctions_pointM( :, i ) * norm_kpV(:) )
                   
                   IF ( norm .LE. 0.0D0 ) STOP &
                        'Error calculate_kp_eigenvalues3D: Eigenfunction zero.'
                   
                   DO n1 = 1, dim_qr * 8
                      eigenfunctions_pointM( n1, i ) = &
                           eigenfunctions_pointM( n1, i ) / SQRT( norm )
                   END DO

                   !...........................................................
                   
                   SELECT CASE( kind_kpC )
                      
                      CASE ( '8x8' )
                      
                      DO l = 1, dim_qr
                         spinorM( count, l, 1 : 8 ) = &
                              eigenfunctions_pointM( (l-1) * 8 + 1 : l * 8, i )
                      END DO
                      
                   CASE DEFAULT
                      
                      STOP 'Error calculate_kp_eigenvalues: &
                           &For holes only 6x6 or 8x8 k.p!'
                      
                   END SELECT

                   !...........................................................

                ELSE IF ( ( eigenval_pointV(i) .LT. E_separate ) .AND. &
                     ( count .GE. num_states ) ) THEN
                   
                   !...........................................................

                   count = count + 1

                   !...........................................................

                END IF

             END DO

             !----------------------------------------------------------------

             DEALLOCATE( eigenval_pointV, eigenfunctions_pointM )

             !---------- Check if enough eigenvalues have been calculated. ---
             
             IF ( count .LT. num_states ) THEN

                WRITE (*,*) 'calculate_kp_eigenvalues3D: &
                     &Number of calculated eigenvalues not sufficient.'
                
                WRITE(*,*) 'Calculate anew.'
                
                mult_num_ev_arnoldi = mult_num_ev_arnoldi + add_num_ev_arnoldi
                
                GOTO 15
                
             ELSE IF ( count .GT. num_states ) THEN
                
                percentage = DBLE( count ) / num_ev
                mult_num_ev_arnoldi = MAX( ( 1.0D0 / percentage ), &
                     mult_num_ev_arnoldi - dec_num_ev_arnoldi + 0.1D0 )
                kp_hl_3DM( num_qr )%mult_num_ev_arnoldi = mult_num_ev_arnoldi
                
                WRITE (*,*) 'Number of eigenvalues to be calculated reduced.'
                WRITE (*,*) 'New number of eigenvalues: ', mult_num_ev_arnoldi
                
             END IF
             
             !-------- end of eigenvalue calculation -------------------------

             CALL deallocate_sparse_matrix

             !----------------------------------------------------------------

          END SELECT

          !___________________________________________________________________
          
          CASE ( 'el' )

          !___________________________________________________________________

          SELECT CASE ( kind_kpC )

             !----------------------------------------------------------------

             CASE ( '6x6' )

             STOP 'Error mod_calculate_kp_eigenvalues3D: &
                  &Electrons within 6x6 k.p not allowed.'
             
             !----------------------------------------------------------------

             CASE ( '8x8' )
             
             !----------------------------------------------------------------
             !
             ! Note:
             !
             ! num_states eigenvalues above E_separate are requested,
             ! so one calculates num_ev=num_states*mult_num_ev_arnoldi.
             ! Then they are sorted and the lowest num_states eigenvalues
             ! above E_separate are taken.
             ! 
             ! If they are not sufficient, 
             ! 
             ! mult_num_ev_arnoldi -> mult_num_ev_arnoldi+add_num_ev_arnoldi
             !
             ! and the calculation is started anew.
             !
             ! If more are calculated than necessary:
             !    
             !   mult_num_ev_arnoldi=mult_num_ev_arnoldi-dec_num_ev_arnoldi
             !       
             !----------------------------------------------------------------

25           CONTINUE

             n                   = dim_qr * k
             mult_num_ev_arnoldi = kp_el_3DM( num_qr )%mult_num_ev_arnoldi
             num_ev              = num_states * mult_num_ev_arnoldi
             
             ALLOCATE( eigenval_pointV( num_ev ), &
                  eigenfunctions_pointM( n, num_ev ), STAT = istatus )
             IF ( istatus .NE. 0 ) STOP &
                  'Error calculate_kp_eigenvalues: Allocation failed.'
             
             eigenval_pointV = 0.0D0
             
             !-------- initial values -----
             
             DO i = 1, num_states
                DO l = 1, dim_qr
                   eigenfunctions_pointM( ( l - 1 ) * 8 + 1 : l * 8, i ) = &
                        spinorM( i, l, 1 : 8 ) 
                END DO
             END DO

             eigenfunctions_pointM( 1:n, num_states+1:num_ev ) = CMPLX(0.0D0)
             
             !----------------------------------------------------------------
             
             tolerance      = MIN( SchroedingerkpResidual3D, min_precision_kp )
             max_iterations = SchroedingerkpIterations3D
             
             CALL transform_sparse_to_ilu3D( 'kp' )

             !-------------------------------------------------------
             ! transform matrix --> E_separate = 2.0d0
             ! matrix = matrix + shift
             ! E_separate_neu = E_separate + shift=2.0d0
             !-------------------------------------------------------

             WRITE (*,*) 'E_separate', E_separate
   
             shift       = -E_separate + 2.0D0
             E_separate2 =  E_separate + shift
             sa(1:n)     =  sa(1:n)    + shift
             
             sigma       = CMPLX( E_separate2 )
             E_sigma     = sigma

             CALL setup_david_precond_jac( sigma )
             
             WRITE (*,*) 'Jacobi-Davidson-in: ', num_states
             WRITE (*,*) 'Tolerance, Max no. of iterations, E_separate2: ', &
                  tolerance, max_iterations, E_separate2
             
             CALL eigensolver_davidson( n, num_ev, eigenval_pointV, &
                  eigenfunctions_pointM, sigma, n_converged, &
                  max_iterations, tolerance )
             
             ! correct eigenvalues by shift ****

             eigenval_pointV = eigenval_pointV - shift
             
             WRITE (*,*) 'Jacobi-Davidson-out'
             
             CALL dealloc_david_precond_jac
             
             !-------------  sort eigenvalues --------------------------------
             !
             !  Jacobi-Davidson provides eigenvalues in descending order.
             !
             !----------------------------------------------------------------
             
             DO i = 2, num_ev
                IF ( eigenval_pointV(i) .GT. eigenval_pointV(i-1) ) STOP &
                     'Error calculate_kp_eigenvalues3D:&
                     & Eigenvalues provided in wrong order by Jacobi-Davidson.'
             END DO
             
             !----------------------------------------------------------------
             
             WRITE (*,*) 'here_el'
             WRITE (*,*) 'Eigenvalues: ', eigenval_pointV
             
             !----------------------------------------------------------------
             
             count       = 0
             eigenvalueV = 0.0D0
             spinorM     = CMPLX(0.0D0)
             
             DO  i=num_ev,1,-1
                
                IF ( ( eigenval_pointV(i) .GT. E_separate ) .AND. &
                     ( count .LT. num_states ) ) THEN
                   
                   !...........................................................
                   
                   count = count + 1
                   eigenvalueV( count ) = eigenval_pointV( i )
                   
                   DO n1 = 1, dim_qr * 8
                      eigenfunctions_pointM( n1, i ) = &
                           diag_kpV( n1 ) * eigenfunctions_pointM( n1, i )
                   END DO
                   
                   norm = SUM( CONJG( eigenfunctions_pointM( :, i ) ) * &
                        eigenfunctions_pointM( :, i ) * norm_kpV(:) )
                   
                   IF ( norm .LE. 0.0D0 ) STOP &
                        'Error calculate_kp_eigenvalues2D: Eigenfunction zero.'
                     
                   DO n1 = 1, dim_qr * 8
                      eigenfunctions_pointM( n1, i ) = &
                           eigenfunctions_pointM( n1, i ) / SQRT( norm )
                   END DO
                     
                   !...........................................................

                   SELECT CASE ( kind_kpC )
                      
                      CASE ( '8x8' )
                        
                      DO l = 1, dim_qr
                         spinorM( count, l, 1 : 8 ) = &
                              eigenfunctions_pointM( (l-1) * 8 + 1 : l * 8, i )
                      END DO
                      
                   CASE DEFAULT
                    
                      STOP 'Error calculate_kp_eigenvalues: &
                           &For holes only 6x6 k.p or 8x8 k.p!'
                     
                   END SELECT
                   
                   !...........................................................

                ELSE IF ( ( eigenval_pointV(i) .GT. E_separate ) .AND. &
                     ( count .GE. num_states ) ) THEN

                   !...........................................................
     
                   count = count + 1
                     
                END IF
                  
             END DO
               
             !----------------------------------------------------------------

             DEALLOCATE(eigenval_pointV,eigenfunctions_pointM)
               
             !---------- Check if enough eigenvalues have been calculated. ---

             IF ( count .LT. num_states ) THEN

                WRITE (*,*) 'calculate_kp_eigenvalues3D: &
                     &Number of calculated eigenvalues not sufficient.'
                  
                WRITE (*,*) 'Calculate anew.'
                  
                mult_num_ev_arnoldi =  mult_num_ev_arnoldi + add_num_ev_arnoldi
                GOTO 25
                
             ELSE IF ( count .GT. num_states ) THEN
                
                percentage = DBLE( count ) / num_ev
                mult_num_ev_arnoldi = MAX( ( 1.0D0 / percentage ), &
                     mult_num_ev_arnoldi - dec_num_ev_arnoldi + 0.1D0 )
                kp_el_3DM( num_qr )%mult_num_ev_arnoldi = mult_num_ev_arnoldi
                  
                WRITE (*,*) 'Number of eigenvalues to be calculated reduced.'
                WRITE (*,*) 'New number of eigenvalues: ', mult_num_ev_arnoldi
                  
             END IF

             !-------- end of eigenvalue calculation -------------------------

             CALL deallocate_sparse_matrix
               
          END SELECT
          
          !-------------------------------------------------------------------
          
       END SELECT
       
       !______________________________________________________________________

    END SELECT
    
    DEALLOCATE( norm_kpV )

    CallRoutineFirstTime = .FALSE.

    CALL testtime( CallRoutineFirstTime, .FALSE. )
            
    !=========================================================================
    
  END SUBROUTINE calculate_kp_eigenvalues
  
  !===========================================================================

END MODULE mod_calculate_kp_eigenvalues3D

!=============================================================================
  
