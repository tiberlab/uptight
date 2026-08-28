! kmesh_gen.f90
!
! Module sinh k-mesh cho CPA:
!   1) Monkhorst-Pack full 3D BZ mesh (dung cho BZ sum trong Soven self-consistency,
!      cong thuc 3.3 trong CPA_implementation_plan.md)
!   2) K-path tuyen tinh qua cac diem doi xung cao, hard-code L-Gamma-X cho
!      zinc-blende (dung de ve spectral function A(k,E))
!
! Quy uoc: k-vector duoc bieu dien trong toa do Cartesian (don vi 2*pi/a_lattice),
! tinh tu reciprocal lattice vectors rec_latt(3,3) (hang = 1 rec. vector).
!
MODULE kmesh_gen

  USE precision, ONLY : dp

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: generate_mp_mesh
  PUBLIC :: generate_kpath_LGX

CONTAINS

  !===========================================================================
  ! Subroutine generate_mp_mesh
  !
  ! Sinh Monkhorst-Pack mesh Nx x Ny x Nz trong BZ, tra ve toa do Cartesian.
  !
  ! INPUT:
  !   rec_latt(3,3) : reciprocal lattice vectors, hang i la vector b_i
  !                   (don vi 2*pi/a, tuong thich voi basis%rec_latt)
  !   nk(3)         : so diem chia theo tung chieu (Nx, Ny, Nz)
  !
  ! OUTPUT:
  !   kpts(3, n_k)  : danh sach toa do Cartesian cua cac k-diem
  !                   (ALLOCATABLE, duoc allocate trong subroutine nay)
  !   n_k           : tong so k-diem = Nx*Ny*Nz
  !
  ! Cong thuc Monkhorst-Pack (khong shift):
  !   k_frac(i) = ( (2*m_i - N_i - 1) / (2*N_i) ),  m_i = 1,...,N_i
  !   k_cart = k_frac(1)*b1 + k_frac(2)*b2 + k_frac(3)*b3
  !===========================================================================

  SUBROUTINE generate_mp_mesh( rec_latt, nk, kpts, n_k )

    REAL( dp ), DIMENSION(3,3), INTENT( IN )  :: rec_latt
    INTEGER,    DIMENSION(3),   INTENT( IN )  :: nk

    REAL( dp ), DIMENSION(:,:), ALLOCATABLE, INTENT( OUT ) :: kpts
    INTEGER,                                 INTENT( OUT ) :: n_k

    !=========================================================================
    ! Local variables
    !=========================================================================

    INTEGER :: ix, iy, iz, i_k
    REAL( dp ) :: fx, fy, fz

    !=========================================================================

    n_k = nk(1) * nk(2) * nk(3)

    ALLOCATE( kpts(3, n_k) )

    i_k = 0

    DO ix = 1, nk(1)
       fx = REAL( 2*ix - nk(1) - 1, dp ) / REAL( 2*nk(1), dp )

       DO iy = 1, nk(2)
          fy = REAL( 2*iy - nk(2) - 1, dp ) / REAL( 2*nk(2), dp )

          DO iz = 1, nk(3)
             fz = REAL( 2*iz - nk(3) - 1, dp ) / REAL( 2*nk(3), dp )

             i_k = i_k + 1

             kpts(:, i_k) = fx * rec_latt(1,:) + fy * rec_latt(2,:) &
                                                + fz * rec_latt(3,:)

          END DO
       END DO
    END DO

  END SUBROUTINE generate_mp_mesh


  !===========================================================================
  ! Subroutine generate_kpath_LGX
  !
  ! Sinh k-path tuyen tinh L -> Gamma -> X cho zinc-blende, hard-coded.
  !
  ! Toa do fractional chuan (don vi 2*pi/a) trong he Cartesian cua zinc-blende
  ! (rec. lattice cua FCC):
  !   L     = (0.5, 0.5, 0.5)
  !   Gamma = (0.0, 0.0, 0.0)
  !   X     = (1.0, 0.0, 0.0)
  !
  ! Cac toa do nay duoc nhan voi rec_latt truyen vao (tuong tu generate_mp_mesh)
  ! -> k_cart = f(1)*b1 + f(2)*b2 + f(3)*b3
  !
  ! INPUT:
  !   rec_latt(3,3)  : reciprocal lattice vectors
  !   n_per_segment  : so diem chia tren MOI doan (L->Gamma va Gamma->X),
  !                    khong tinh diem dau doan (tru doan dau tien)
  !                    -> tong so diem = 2*n_per_segment + 1
  !
  ! OUTPUT:
  !   kpts(3, n_k)   : toa do Cartesian cac k-diem doc theo path
  !   k_dist(n_k)    : khoang cach luy ke doc theo path (dung de ve truc x)
  !   n_k            : tong so k-diem
  !   label_pos(3)   : chi so (index) trong kpts ung voi L, Gamma, X
  !                    (dung de danh dau tren truc hoanh khi ve hinh)
  !===========================================================================

  SUBROUTINE generate_kpath_LGX( rec_latt, n_per_segment, &
                                  kpts, k_dist, n_k, label_pos )

    REAL( dp ), DIMENSION(3,3), INTENT( IN )  :: rec_latt
    INTEGER,                    INTENT( IN )  :: n_per_segment

    REAL( dp ), DIMENSION(:,:), ALLOCATABLE, INTENT( OUT ) :: kpts
    REAL( dp ), DIMENSION(:),   ALLOCATABLE, INTENT( OUT ) :: k_dist
    INTEGER,                                 INTENT( OUT ) :: n_k
    INTEGER,    DIMENSION(3),                INTENT( OUT ) :: label_pos

    !=========================================================================
    ! Local variables
    !=========================================================================

    REAL( dp ), DIMENSION(3) :: L_frac, G_frac, X_frac
    REAL( dp ), DIMENSION(3) :: L_cart, G_cart, X_cart
    REAL( dp ), DIMENSION(3) :: kvec_prev, kvec_curr

    INTEGER :: i, i_k
    REAL( dp ) :: t

    !=========================================================================
    ! Diem doi xung cao (fractional, don vi 2*pi/a)
    !=========================================================================

    L_frac = (/ 0.5_dp, 0.5_dp, 0.5_dp /)
    G_frac = (/ 0.0_dp, 0.0_dp, 0.0_dp /)
    X_frac = (/ 1.0_dp, 0.0_dp, 0.0_dp /)

    L_cart = L_frac(1)*rec_latt(1,:) + L_frac(2)*rec_latt(2,:) &
                                      + L_frac(3)*rec_latt(3,:)
    G_cart = G_frac(1)*rec_latt(1,:) + G_frac(2)*rec_latt(2,:) &
                                      + G_frac(3)*rec_latt(3,:)
    X_cart = X_frac(1)*rec_latt(1,:) + X_frac(2)*rec_latt(2,:) &
                                      + X_frac(3)*rec_latt(3,:)

    !=========================================================================
    ! Tong so diem: L -> Gamma (n_per_segment doan, n_per_segment+1 diem)
    !               Gamma -> X (n_per_segment doan, khong lap lai Gamma)
    !=========================================================================

    n_k = 2 * n_per_segment + 1

    ALLOCATE( kpts(3, n_k) )
    ALLOCATE( k_dist(n_k) )

    i_k = 0
    k_dist(1) = 0.0_dp

    !-------------------------------------------------------------------------
    ! Doan L -> Gamma
    !-------------------------------------------------------------------------

    DO i = 0, n_per_segment

       t = REAL(i, dp) / REAL(n_per_segment, dp)

       i_k = i_k + 1
       kpts(:, i_k) = (1.0_dp - t) * L_cart + t * G_cart

       IF ( i_k .GT. 1 ) THEN
          kvec_prev = kpts(:, i_k - 1)
          kvec_curr = kpts(:, i_k)
          k_dist(i_k) = k_dist(i_k - 1) + &
               SQRT( SUM( (kvec_curr - kvec_prev)**2 ) )
       END IF

       IF ( i .EQ. 0 )            label_pos(1) = i_k   ! L
       IF ( i .EQ. n_per_segment) label_pos(2) = i_k   ! Gamma

    END DO

    !-------------------------------------------------------------------------
    ! Doan Gamma -> X (bo qua i=0 vi da co Gamma o tren)
    !-------------------------------------------------------------------------

    DO i = 1, n_per_segment

       t = REAL(i, dp) / REAL(n_per_segment, dp)

       i_k = i_k + 1
       kpts(:, i_k) = (1.0_dp - t) * G_cart + t * X_cart

       kvec_prev = kpts(:, i_k - 1)
       kvec_curr = kpts(:, i_k)
       k_dist(i_k) = k_dist(i_k - 1) + &
            SQRT( SUM( (kvec_curr - kvec_prev)**2 ) )

       IF ( i .EQ. n_per_segment ) label_pos(3) = i_k   ! X

    END DO

  END SUBROUTINE generate_kpath_LGX

END MODULE kmesh_gen
