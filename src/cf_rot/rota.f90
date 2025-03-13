module crystal_field_d

implicit none

integer, parameter :: dp=8
real(dp), parameter :: pi = 3.1415926535897932d0
real(dp), parameter :: eps = 1.d-15

private

public :: d_rotation, H_d_rotation, tesseral_trans


contains

  !======================================================================
  !
  ! 1. Computes  the Euler angles for the Rotation that takes  
  !
  ! the vector r_orig -> r_dest
  !
  ! 2. Computes  the representation of the Rotation Matrix
  !
  !    on the spinor subspace for l=2 ( d- orbitals ) 
  !
  !======================================================================
  subroutine d_rotation(r_orig,r_dest,R)

    implicit none

    real(dp), intent(in) :: r_orig(3), r_dest(3)
    complex(dp), intent(out) :: R(5,5)

    ! Internal variables

    complex(dp) :: j=(0.d0,1.d0)
    real(dp) :: c(3), u(3), r12(3)
    real(dp) :: psi, chi, phi, teta
    real(dp) :: alpha,gamma,beta
    real(dp) :: bcos,bsin,bcos2,bsin2
    complex(dp) :: eia,eig,eima,eimg
    real(dp) :: bcos3,bsin3,bcos4,bsin4
    
    ! Find cosine directors of the origin vector
    !
    r12(:) = r_orig(:)
    c(:)=r12(:)/sqrt(r12(1)*r12(1)+ r12(2)*r12(2)+r12(3)*r12(3))

    psi=acos(c(3))

    if(abs(psi).gt.EPS) then

       if(c(2)>=0) then
          chi= acos(c(1)/sin(psi))
       else
          chi= pi + acos(-c(1)/sin(psi))
       endif

    else

       chi= 0.d0

    endif
    ! 
    ! Find cosine directors of the destination vector
    !
    r12(:) = r_dest(:)
    u(:)=r12(:)/sqrt(r12(1)*r12(1)+ r12(2)*r12(2)+r12(3)*r12(3))
    !
    teta=acos(u(3))
    !
    if(abs(teta).gt.EPS) then

       if(u(2)>=0) then
          phi= acos(u(1)/sin(teta))
       else
          phi= pi + acos(-u(1)/sin(teta))
       endif

    else

       phi= 0.d0

    endif

    ! define Euler angles of rotations of r_orig -> r_dest
    ! R(a,b,g) = R_3(a) R_2(b) R_3(g)
    !
    ! Obtained by combining R_3(-phi)*R_2(-teta)*R_2(psi)*R_3(chi)
    !
    alpha=-phi
    beta=psi-teta
    gamma=chi
    
    !  Define Some practical utility variables
    !
    eia=exp(j*alpha)
    eig=exp(j*gamma)
    eima=exp(-j*alpha)
    eimg=exp(-j*gamma)

    bcos=cos(beta/2.d0)
    bsin=sin(beta/2.d0)
    bcos2=bcos*bcos
    bsin2=bsin*bsin
    bcos3=bcos2*bcos
    bsin3=bsin2*bsin
    bcos4=bcos2*bcos2
    bsin4=bsin2*bsin2
    ! 
    ! Define the l=2 representation of the SO(3)
    ! Assuming the ordering for the angular momentum basis-set:
    ! |2>  |1>  |0>  |-1>  |-2>
    !
    R(1,1)=eia*eia*eig*eig*bcos4
    R(2,1)=2.d0*eia*eia*eig*bcos3*bsin
    R(3,1)=sqrt(6.d0)*eia*eia*bcos2*bsin2
    R(4,1)=2.d0*eia*eia*eimg*bcos*bsin3
    R(5,1)=eia*eia*eimg*eimg*bsin4

    R(1,2)=-2.d0*eia*eig*eig*bcos3*bsin
    R(2,2)=eia*eig*(bcos4-3.d0*bcos2*bsin2)
    R(3,2)=sqrt(6.d0)*eia*(bcos3*bsin-bcos*bsin3)
    R(4,2)=eia*eimg*(3.d0*bcos2*bsin2-bsin4)
    R(5,2)=2.d0*eia*eimg*eimg*bcos*bsin3

    R(1,3)=sqrt(6.d0)*eig*eig*bcos2*bsin2
    R(2,3)=sqrt(6.d0)*eig*(bcos*bsin3-bcos3*bsin)
    R(3,3)=(bcos4-4.d0*bcos2*bsin2+bsin4)
    R(4,3)=sqrt(6.d0)*eimg*(bcos3*bsin-bcos*bsin3)
    R(5,3)=sqrt(6.d0)*eimg*eimg*bcos2*bsin2

    R(1,4)=-2.d0*eima*eig*eig*bcos*bsin3
    R(2,4)=eima*eig*(3.d0*bcos2*bsin2-bsin4)
    R(3,4)=sqrt(6.d0)*eima*(bcos*bsin3-bcos3*bsin)
    R(4,4)=eima*eimg*(bcos4-3.d0*bcos2*bsin2)
    R(5,4)=2.d0*eima*eimg*eimg*bcos3*bsin

    R(1,5)=eima*eima*eig*eig*bsin4
    R(2,5)=-2.d0*eima*eima*eig*bcos*bsin3
    R(3,5)=sqrt(6.d0)*eima*eima*bcos2*bsin2
    R(4,5)=-2.d0*eima*eima*eimg*bcos3*bsin
    R(5,5)=eima*eima*eimg*eimg*bcos4
    
    !R=conjg(transpose(R))

    return

  end subroutine d_rotation
  ! --------------------------------------------------------
  !======================================================================
  !
  ! 1. Computes  the Unitary transformation that combines 
  !
  !    the |l=2,m> basis set into the tesseral basis set
  ! 
  !    The original angular momentum basis set:
  ! 
  !      |2>       |1>        |0>        |-1>        |-2>     
  !   
  !    is rotated into:
  !
  !    1 |dxy>   = i/sqrt(2)(|-2> - |2>); 
  !    2 |dyz>   = i/sqrt(2)(|-1>+|1>); 
  !    3 |dz2>   = |0>;  
  !    4 |dxz>   = 1/sqrt(2)(|-1> - |1>);
  !    5 |dx2y2> = 1/sqrt(2)(|-2> + |2>)
  !
  ! 2. Transform a matrix represented on the momentum basis
  !    Into a matrix represented on the tesseral basis
  !
  !          + 
  !    Rd = U  Rd  U
  !======================================================================
  subroutine tesseral_trans(Rd)

    complex(dp), intent(inout) :: Rd(5,5)

    complex(dp) :: Rx(5,5)
    complex(dp) :: U(5,5)
    complex(dp) :: j=(0.d0,1.d0)   

    U(:,:)= (0.d0,0.d0)

    U(1,1)= -j/sqrt(2.d0)
    U(2,1)= (0.d0,0.d0)
    U(3,1)= (0.d0,0.d0)
    U(4,1)= (0.d0,0.d0)
    U(5,1)= j/sqrt(2.d0)

    U(1,2)= (0.d0,0.d0)
    U(2,2)= j/sqrt(2.d0)
    U(3,2)= (0.d0,0.d0)
    U(4,2)= j/sqrt(2.d0)
    U(5,2)= (0.d0,0.d0)

    U(1,3)= (0.d0,0.d0)
    U(2,3)= (0.d0,0.d0)
    U(3,3)= (1.d0,0.d0)
    U(4,3)= (0.d0,0.d0) 
    U(5,3)= (0.d0,0.d0)

    U(1,4)= (0.d0,0.d0)
    U(2,4)= -1.d0/sqrt(2.d0)
    U(3,4)= (0.d0,0.d0)
    U(4,4)= 1.d0/sqrt(2.d0)
    U(5,4)= (0.d0,0.d0)
    
    U(1,5)= 1.d0/sqrt(2.d0)
    U(2,5)= (0.d0,0.d0)
    U(3,5)= (0.d0,0.d0) 
    U(4,5)= (0.d0,0.d0)
    U(5,5)= 1.d0/sqrt(2.d0)

    Rx=matmul(conjg(transpose(U)),Rd)
    Rd=(0.d0,0.d0)
    Rd=matmul(Rx,U)
  
  end subroutine tesseral_trans
  ! --------------------------------------------------------
  !======================================================================
  !
  ! 1. Performs the Hamiltonian rotation 
  !
  !    Hd == Diagonal representation of Hd
  ! 
  !    Rd == Rotation matrix
  !   
  !    HR == Rotated Hamiltonian
  !
  !           +
  !    HR = Rd  Hd  Rd
  !======================================================================
  subroutine H_d_rotation(Hd,Rd,HR)
    
    complex(dp), intent(in) :: Hd(5,5), Rd(5,5)
    complex(dp), intent(out) :: HR(5,5)

    complex(dp) :: Td(5,5)

    !Td=matmul(conjg(transpose(Rd)),Hd)
    !HR=matmul(Td,Rd)  
    Td=matmul(Hd,conjg(transpose(Rd)))
    HR=matmul(Rd,Td)

  end subroutine H_d_rotation
  ! -------------------------------------------------------


end module crystal_field_d
