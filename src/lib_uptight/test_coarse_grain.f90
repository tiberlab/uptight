program test_coarse_grain
  ! Regression for the last stage of the reduction pipeline.  It does not
  ! require METIS: blocks and their retained local bases are supplied directly.
  use precision, only : dp
  use upt_param, only : OUPT
  use coarse_grain, only : cg_lift, cg_clear
  implicit none

  type(OUPT) :: upt
  complex(dp) :: reduced(2,1), physical(3,1)
  real(dp), parameter :: tol = 100.0_dp * epsilon(1.0_dp)
  real(dp) :: s

  s = sqrt(0.5_dp)
  allocate(upt%cg_blocks(2))
  upt%cg_blocks(1)%nrow = 2; upt%cg_blocks(1)%nret = 1
  allocate(upt%cg_blocks(1)%rows(2), upt%cg_blocks(1)%q(2,1))
  upt%cg_blocks(1)%rows = (/ 1, 3 /)
  upt%cg_blocks(1)%q(:,1) = cmplx((/ s, s /), 0.0_dp, dp)
  upt%cg_blocks(2)%nrow = 1; upt%cg_blocks(2)%nret = 1
  allocate(upt%cg_blocks(2)%rows(1), upt%cg_blocks(2)%q(1,1))
  upt%cg_blocks(2)%rows = 2
  upt%cg_blocks(2)%q(1,1) = (1.0_dp, 0.0_dp)

  reduced(:,1) = cmplx((/ 2.0_dp, 3.0_dp /), 0.0_dp, dp)
  call cg_lift(upt, reduced, physical)
  if (maxval(abs(physical(:,1) - cmplx((/ 2.0_dp*s, 3.0_dp, 2.0_dp*s /), 0.0_dp, dp))) > tol) stop 1
  call cg_clear(upt)
end program test_coarse_grain
