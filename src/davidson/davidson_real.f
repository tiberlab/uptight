      SUBROUTINE Jdqz_real ( alpha, beta, eivec, wanted,
     $     n, shift, eps, kmax, jmax, jmin, method, m, l,
     $     maxnmv, maxstep, lock, order, testspace, work, lwork )

c===========================================================================
c
c     Programmer: Diederik R. Fokkema
c     Modified         : M. van Gijzen
c     Modified 05-24-96: M. Kooper: ra and rb, the Schur matrices of A and B, 
c              added, as well as the vectors sa and sb, which contain the
c              innerproducts of ra with Z and rb with Z. This is added to be
c              enable to compute the eigenvectors in EIVEC
c     Modification 08-27-96: Different computation of eigenvectors, MvG
c
c===========================================================================

      IMPLICIT NONE

c===========================================================================
c
c     .. PARAMETERS ..
c
c===========================================================================

      INTEGER gmres, cgstab
      PARAMETER ( gmres = 1, cgstab = 2 )

      INTEGER kmax, jmax, jmin, method, m, l, maxnmv, maxstep, order
      INTEGER testspace, n, lwork

      DOUBLE PRECISION eps, lock

      DOUBLE COMPLEX shift
      DOUBLE COMPLEX work(n,*), eivec(n,*), alpha(*), beta(*)

      LOGICAL wanted

c===========================================================================
c
c     .. LOCAL PARAMETERS ..
c
c===========================================================================

      LOGICAL loop, try, found, ok

      INTEGER ldvs, ldzwork, iseed(4)

c      PARAMETER ( ldvs = 301, ldzwork = 4 * ldvs )
c
c --------------------------------------------------------------
c    SACCONI changed !!!!!!

      PARAMETER ( ldvs = 1501 , ldzwork = 4 * ldvs )

c --------------------------------------------------------------

c     was changed by S. Birner to make it work for more than 100 eigenvalues
c     ldvs must be larger than jmax = MAX( 3 * kmax, 20 )
c     PARAMETER (ldvs= 50, ldzwork=4*ldvs)

      DOUBLE COMPLEX ma(ldvs,ldvs), mb(ldvs,ldvs),
     $     zma(ldvs,ldvs), zmb(ldvs,ldvs),
     $     vsl(ldvs,ldvs), vsr(ldvs,ldvs),
     $     ra(ldvs,ldvs), rb(ldvs, ldvs),
     $     zwork(4*ldvs), aconv(ldvs), bconv(ldvs)

      INTEGER ldqkz
      PARAMETER ( ldqkz = ldvs )

      INTEGER ipivqkz(ldqkz)

      DOUBLE COMPLEX mqkz(ldqkz,ldqkz), invqkz(ldqkz,ldqkz), f(ldqkz)

      INTEGER i, j, k, info, mxmv, step
      INTEGER d, u, v, w, aux, av, bv, q, z, kz, itmp, tp
      INTEGER solvestep

      DOUBLE PRECISION rnrm, rwork( 3 * ldvs ), deps
      DOUBLE PRECISION dtmp

      DOUBLE COMPLEX zalpha, zbeta, targeta, targetb, evcond
      DOUBLE COMPLEX shifta, shiftb

      DOUBLE COMPLEX zero, one
      PARAMETER ( zero = ( 0.0d0, 0.0d0 ), one = ( 1.0d0, 0.0d0 ) )

c===========================================================================
c
c     .. FUNCTIONS ..
c
c===========================================================================

      DOUBLE PRECISION dznrm2

c===========================================================================
c
c     .. DATA ..
c
c===========================================================================

      DATA iseed / 3, 3, 1966, 29 /

c===========================================================================
c
c     ... Executable Statements
c
c===========================================================================

c...     Are there errors in the PARAMETERS ?

      IF ( kmax .GT. jmax )
     $     CALL error_real( 'jdqz: kmax greater than jmax' )

c     changed by S. Birner to allow more than 50 eigenvalues
c     IF ( jmax .GT. 50 )

      IF ( jmax .GT. MAX( 3 * kmax, 20 ) ) 
     $     CALL error_real( 'jdqz: jmax greater than MAX(3*kmax,20)' )

      IF ( method .NE. 1 .AND. method .NE. 2 )
     $     CALL error_real( 'jdqz: illegal choice for solution method' )

      IF ( order .LT. -2 .OR. order .GT. 2 )
     $     CALL error_real
     $     ( 'illegal value for order, must be between -2 and 2' )

c===========================================================================

c...  d = rhs, these pointers refer to the columns of the workspace
   
      d = 1

c...  Workspace for jdqzmv_real
      
      tp = d + 1

c...  u = pointer to Krylov space GMRES( m ) or Bi-CSTAB( l )
      
      u = tp + 1

c...  v = pointer to search space JDQZ with max dimension jmax
      
      IF ( method .EQ. gmres ) THEN
         v = u + m + 1
      ELSE IF ( method .EQ. cgstab ) THEN
	 v = u + 2*l + 6
      END IF

c...  w = pointer to test subspace JDQZ with max dimension jmax
      
      w = v + jmax

c...  av = pointer to subspace AV with max dimension jmax
     
      av = w + jmax

c...  bv = pointer to subspace BV with max dimension jmax
    
      bv = av + jmax

c...  aux =
      
      aux = bv + jmax

c...  q = pointer to search Schur basis in JDQZ with max dimension kmax
     
      q = aux + jmax

c...  z = pointer to test Schur basis in JDQZ with max dimension kmax
     
      z = q + kmax

c...  kz = pointer to matrix K^{-1}Z_k
     
      kz = z + kmax

      IF ( kz + kmax - 1 .GT. lwork )
     $     CALL error_real ( 'qz: memory fault' )

c===========================================================================
c
c     --- initialize loop
c
c===========================================================================

      ok = .TRUE.

      evcond = CMPLX( SQRT( ABS( shift )**2 + ABS( one )**2 ) )
      shifta = shift / evcond
      shiftb = one / evcond

      targeta = shifta
      targetb = shiftb

      zalpha = shifta
      zbeta = shiftb

      step = 0
      deps = DBLE( one )
      mxmv = 0

      solvestep = 0

      j = 0
      k = 0

c===========================================================================
c
c     --- loop
c
c===========================================================================

 100  CONTINUE

      loop = ( k.LT.kmax .AND. step.LT.maxstep )

      IF ( loop ) THEN

	 step = step + 1
	 solvestep = solvestep + 1

	 IF ( j .EQ. 0 ) THEN
	    
            CALL zlarnv( 2, iseed, n, work( 1, v+j ) )
	    CALL zlarnv( 2, iseed, n, work( 1, w+j ) )

	    DO i = 1, n

	       dtmp = DBLE( work( i, v+j ) )
	       work( i, v+j ) = CMPLX( dtmp, 0d0 )
	       dtmp = DBLE( work( i, w+j ) )
	       work( i, w+j ) = CMPLX( dtmp, 0d0 )

	    END DO

	 ELSE

	    mxmv = maxnmv
	    deps = 2.0d0**( -solvestep )

	    IF ( j .LT. jmin ) THEN

	       mxmv = 1
	       CALL zgmres_real( n, work( 1, v+j ), work(1,d),
     $              m, deps, mxmv, zalpha, zbeta, k+1,
     $              work(1,kz), work(1,q), invqkz, ldqkz,
     $              ipivqkz, f, work(1,u), work(1,tp) )
	    
            ELSEIF ( method .EQ. gmres ) THEN
	    
               mxmv = m
	       CALL zgmres_real( n, work( 1, v+j ), work(1,d),
     $              m, deps, mxmv, zalpha, zbeta, k+1,
     $              work(1,kz), work(1,q), invqkz, ldqkz,
     $              ipivqkz, f, work(1,u), work(1,tp) )
               
	    ELSEIF ( method.EQ. cgstab ) THEN

	       CALL zcgstabl_real( n, work( 1, v+j ), work(1,d),
     $              l, deps, mxmv, zalpha, zbeta, k+1,
     $              work(1,kz), work(1,q), invqkz,
     $              ldqkz, ipivqkz, f, work(1,u), 2*l + 6 )

	    END IF
            
	 END IF
         
	 j = j + 1
         
	 CALL zmgs_real( n, j - 1, work(1,v), work( 1, v+j-1 ), 1 )
	 CALL zmgs_real( n, k, work(1,q), work( 1, v+j-1 ), 1 )
         
	 IF ( testspace .EQ. 1 ) THEN

 	    CALL jdqzmv_real( n, work( 1, v+j-1 ), work( 1, w+j-1 ),
     $           work(1,tp), -conjg( shiftb ), conjg( shifta ) )
            
	 ELSEIF ( testspace .EQ. 2 ) THEN

 	    CALL jdqzmv_real( n, work( 1, v+j-1 ), work( 1, w+j-1 ),
     $           work(1,tp), -conjg( zbeta ), conjg( zalpha ) )

	 ELSEIF ( testspace .EQ. 3 ) THEN

 	    CALL jdqzmv_real( n, work( 1, v+j-1 ),
     $           work( 1, w+j-1 ), work(1,tp), shifta, shiftb )

	 END IF

         CALL zmgs_real( n, j-1, work(1,w), work( 1, w+j-1 ), 1 )
         CALL zmgs_real( n, k, work(1,z), work( 1, w+j-1 ), 1 )

         CALL amul_real( n, work( 1, v+j-1 ), work( 1, av+j-1 ) )
         CALL bmul_real( n, work( 1, v+j-1 ), work( 1, bv+j-1 ) )

         CALL makemm_real( n, j, work(1,w), work(1,av), ma, zma, ldvs )
         CALL makemm_real( n, j, work(1,w), work(1,bv), mb, zmb, ldvs )

         CALL zgegs( 'v', 'v', j, zma, ldvs, zmb, ldvs, alpha, beta,
     $        vsl, ldvs, vsr, ldvs, zwork, ldzwork, rwork, info )

         try = .TRUE.

 200     CONTINUE

         IF ( try ) THEN

c
c           --- Sort the Petrov pairs ---
c

            CALL qzsort_real( targeta, targetb, j, zma, zmb, vsl, vsr,
     $           ldvs, alpha, beta, order )

            zalpha = zma( 1, 1 )
            zbeta = zmb( 1, 1 )

            evcond = CMPLX( SQRT( ABS(zalpha)**2 + ABS(zbeta)**2 ) )

c
c            --- compute new q ---
c

            CALL zgemv( 'n', n, j, one, work(1,v), n, vsr(1,1),
     $           1, zero, work( 1, q+k ), 1 )
            CALL zmgs_real( n, k, work(1,q), work( 1, q+k ), 1 )

c
c           --- compute new z ---
c

            CALL zgemv( 'n', n, j, one, work(1,w), n, vsl(1,1),
     $           1, zero, work( 1, z+k ), 1 )
            CALL zmgs_real( n, k, work(1,z), work( 1, z+k ), 1 )

c
c     --- Make new qkz ---
c

            CALL zcopy( n, work( 1, z+k ), 1, work( 1, kz+k ), 1 )
            CALL precon_real( n, work( 1, kz+k ) )
            CALL mkqkz_real( n, k+1, work(1,q), work(1,kz),
     $           mqkz, invqkz, ldqkz, ipivqkz )

c
c     --- compute new (right) residual= beta Aq - alpha Bq and
c     orthogonalize this vector on Z.
c

            CALL jdqzmv_real( n, work( 1, q+k ), work(1,d),
     $           work(1,tp), zalpha, zbeta )
            CALL zmgs_real( n, k, work(1, z), work(1,d), 0 )

            rnrm = dznrm2( n, work(1,d), 1 ) / DBLE( evcond )
            
            IF ( rnrm .LT. lock .AND. ok ) THEN
               targeta = zalpha
               targetb = zbeta
               ok = .FALSE.
            END IF

            found = ( rnrm .LT. eps .AND.
     $           ( j .GT. 1 .OR. k .EQ. kmax - 1 ) )
            try = found

            IF ( found ) THEN

c     --- increase the number of found evs by 1 ---
             
               k = k + 1

c     --- store the eigenvalue
               
               aconv( k ) = zalpha
               bconv( k ) = zbeta

               solvestep = 0

               IF ( k .EQ. kmax ) GOTO 100

               CALL zgemm( 'n', 'n', n, j-1, j, one, work(1,v), n,
     $              vsr(1,2), ldvs, zero, work(1,aux), n )

               itmp = v
               v = aux
               aux = itmp

               CALL zgemm( 'n', 'n', n, j-1, j, one, work(1,av), n,
     $              vsr(1,2), ldvs, zero, work(1,aux), n )
              
               itmp = av
               av = aux
               aux = itmp
               
               CALL zgemm( 'n', 'n', n, j-1, j, one, work(1,bv), n,
     $              vsr(1,2), ldvs, zero, work(1,aux), n ) 
               
               itmp = bv
               bv = aux
               aux = itmp
               
               CALL zgemm( 'n', 'n', n, j-1, j, one, work(1,w), n,
     $              vsl(1,2), ldvs, zero, work(1,aux), n )

               itmp = w
               w = aux
               aux = itmp
               j = j-1
               
               CALL zlacpy( 'a', j, j, zma(2,2), ldvs, ma, ldvs )
               CALL zlacpy( 'a', j, j, ma, ldvs, zma, ldvs )
               CALL zlacpy( 'a', j, j, zmb(2,2), ldvs, mb, ldvs )
               CALL zlacpy( 'a', j, j, mb, ldvs, zmb, ldvs )
               CALL zlaset( 'a', j, j, zero, one, vsr, ldvs )
               CALL zlaset( 'a', j, j, zero, one, vsl, ldvs )
               
               targeta = shifta
               targetb = shiftb
               ok = .TRUE.
               mxmv = 0
               deps = DBLE( one )

            ELSE IF ( j .EQ. jmax ) THEN

               CALL zgemm( 'n', 'n', n, jmin, j, one, work(1,v), n,
     $              vsr, ldvs, zero, work(1,aux), n )
               
               itmp = v
               v = aux
               aux = itmp
               
               CALL zgemm( 'n', 'n', n, jmin, j, one, work(1,av), n,
     $              vsr, ldvs, zero, work(1,aux), n )

               itmp = av
               av = aux
               aux = itmp

               CALL zgemm( 'n', 'n', n, jmin, j, one, work(1,bv), n,
     $              vsr, ldvs, zero, work(1,aux), n )
             
               itmp = bv
               bv = aux
               aux = itmp
               
               CALL zgemm( 'n', 'n', n, jmin, j, one, work(1,w), n,
     $              vsl, ldvs, zero, work(1,aux), n )

               itmp = w
               w = aux
               aux = itmp
               j = jmin

               CALL zlacpy( 'a', j, j, zma, ldvs, ma, ldvs )
               CALL zlacpy( 'a', j, j, zmb, ldvs, mb, ldvs )
               CALL zlaset( 'a', j, j, zero, one, vsr, ldvs )
               CALL zlaset( 'a', j, j, zero, one, vsl, ldvs )

            END IF

            GOTO 200

         END IF

         GOTO 100

      END IF

c     
c...  Did enough eigenpairs converge?

      kmax = k

      IF ( wanted ) THEN

c...  Compute the Schur matrices if the eigenvectors are
c...  wanted, work(1,tp) is used for temporary storage
c...  Compute RA:

         CALL zlaset( 'l', k, k, zero, zero, ra, ldvs )    

         DO i = 1, k
           
            CALL amul_real( n, work(1,q+i-1), work(1,tp) )
            CALL zgemv( 'c', n, i, one, work(1,z), n, work(1,tp), 1, 
     $                  zero, ra(1,i), 1 )
         
         END DO

c...  Compute RB:

         CALL zlaset ('l', k, k, zero, zero, rb, ldvs)    
         
         DO i = 1, k
           
            CALL bmul_real( n, work( 1, q+i-1 ), work(1,tp) )
            CALL zgemv( 'c', n, i, one, work(1,z), n, work(1,tp), 1, 
     $                  zero, rb(1,i), 1 )
         
         END DO

c     --- The eigenvectors RA and RB  belonging to the found eigenvalues
c     are computed. The Schur vectors in VR and VS are replaced by the
c     eigenvectors of RA and RB

         CALL zgegv( 'N','V', k, ra, ldvs, rb, ldvs, alpha, beta, vsl,
     $        ldvs, vsr, ldvs, zwork, ldzwork, rwork, info )

c     --- Compute the eigenvectors belonging to the found eigenvalues
c     of A and put them in EIVEC

         CALL zgemm( 'n', 'n', n, k, k, one, work(1,q), n,
     $              vsr, ldvs, zero, eivec, n )

      ELSE

c     
c...  Store the Schurvectors in eivec:

         CALL zcopy( k*n, work(1,q), 1, eivec, 1 )
         CALL zcopy( k, aconv, 1, alpha, 1 )
         CALL zcopy( k, bconv, 1, beta, 1 )

      END IF

c===========================================================================

      END

c===========================================================================


      SUBROUTINE jdqzmv_real ( n, x, y, work, alpha, beta )

c===========================================================================
c     
c     Coded by Diederik Fokkema
c     
c===========================================================================

      IMPLICIT NONE

c===========================================================================
c
c     .. PARAMETERS ..
c
c===========================================================================

      INTEGER n
      DOUBLE COMPLEX alpha, beta, x(*), y(*), work(*) 

c===========================================================================
c
c     .. LOCAL ..
c
c===========================================================================

      INTEGER i

c===========================================================================
c
c     .. Executable statements ..
c
c===========================================================================

      CALL amul_real( n, x, work )
      CALL bmul_real( n, x, y )

      DO i = 1, n

         y( i ) = beta * work( i ) - alpha * y( i )

      END DO
      
c===========================================================================

      END


c===========================================================================


      SUBROUTINE error_real ( m )

c===========================================================================
c
c     Coded by Diederik R. Fokkema
c
c     $Id: error_real.f,v 1.4 1995/07/26 09:26:26 fokkema Exp $
c
c===========================================================================

      IMPLICIT NONE

c===========================================================================
c
c     .. PARAMETERS ..
c
c===========================================================================

      CHARACTER m*(*)

c===========================================================================
c
ctex@ \begin{manpage}{ERROR} 
ctex@ \subtitle{Name}
ctex@    ERROR --- Type an error message and stop
ctex@
ctex@ \subtitle{Declaration}
ctex@    %declaration
ctex@ \subtitle{PARAMETERs}
ctex@    \variable{m}
ctex@       character string. On entry m must contain the message string.
ctex@
ctex@ \subtitle{Description}
ctex@    This SUBROUTINE types an error message and stops.
ctex@
ctex@ \end{manpage}
ctex@ \begin{verbatim}
ctex@    % actual code
ctex@ \end{verbatim}
c
c===========================================================================

c===========================================================================
c
c     .. LOCAL ..
c
c     None.
c
c===========================================================================

c===========================================================================
c
c     .. CALLED SUBROUTINES
c
c     None.
c
c===========================================================================

c===========================================================================
c
c     .. Executable statements
c
c===========================================================================

      PRINT 10, m
 10   FORMAT ( /, 1x, 'Error: ', a, / )

c
c     --- Stop
c

      STOP

c===========================================================================
     
      END


c===========================================================================


      SUBROUTINE makemm_real( n, k, w, v, m, zm, ldm )

c===========================================================================
c
c     Coded by Diederik Fokkema
c
c     $Id$
c
c     Time-stamp: <95/08/03 23:33:20 caveman>
c
c===========================================================================

      IMPLICIT NONE

c===========================================================================
c
c     .. PARAMETERS ..
c
c===========================================================================

      INTEGER n, k, ldm
      DOUBLE COMPLEX w(n,*), v(n,*), m(ldm,*), zm(ldm,*)

c===========================================================================
c
c     .. LOCAL ..
c
c===========================================================================

      INTEGER i, j
      DOUBLE COMPLEX zdotc

c===========================================================================
c
c     .. Executable statements ..
c
c===========================================================================

      DO i = 1, k

         DO j = 1, k

            IF ( i.EQ.k .OR. j.EQ.k )
     $           m( i, j ) = zdotc( n, w(1,i), 1, v(1,j), 1 )
           
            zm( i, j ) = m( i, j )

         END DO

      END DO

c===========================================================================

      END


c===========================================================================


      SUBROUTINE mkqkz_real( n, k, q, kq, qkz, invqkz, ldqkz, ipiv )

c===========================================================================
c
c     Coded by Diederik Fokkema
c
c     $Id: mkqkz.f,v 1.1 1995/08/05 09:08:21 caveman Exp $
c
c     Time-stamp: <95/08/03 23:52:52 caveman>
c
c===========================================================================

      IMPLICIT NONE

c===========================================================================
c
c     .. PARAMETERS ..
c
c===========================================================================

      INTEGER n, k, ldqkz, ipiv(*)
      DOUBLE COMPLEX q(n,*), kq(n,*), qkz(ldqkz,*), invqkz(ldqkz,*)

c===========================================================================
c
c     .. LOCAL ..
c
c===========================================================================

      INTEGER i, j, info
      DOUBLE COMPLEX zdotc

c===========================================================================
c
c     .. Executable statements ..
c
c===========================================================================

      DO i = 1, k

         DO j = 1, k

            IF ( i.EQ.k .OR. j.EQ.k ) 
     $           qkz( i, j ) = zdotc( n, q(1,i), 1, kq(1,j), 1 )

            invqkz( i, j ) = qkz( i, j )

         END DO

      END DO

      CALL zgetrf( k, k, invqkz, ldqkz, ipiv, info )

c===========================================================================

      END


c===========================================================================


      SUBROUTINE myexc_real( n, s, t, z, q, ldz, ifst, ilst )

c===========================================================================
c
c     Coded by Diederik Fokkema
c
c     $Id$
c
c     Time-stamp: <95/10/31 23:51:12 caveman>
c
c===========================================================================

      IMPLICIT NONE

c===========================================================================
c
c     .. PARAMETERS ..
c
c===========================================================================

      INTEGER n, ldz, ifst, ilst
      DOUBLE COMPLEX s(ldz,*), t(ldz,*), z(ldz,*), q(ldz,*)

c===========================================================================
c
c     .. LOCAL ..
c
c===========================================================================

      LOGICAL tlts
      INTEGER k, m1, m2, m3
      DOUBLE COMPLEX f, f1, f2, c1, c2, r, sn
      DOUBLE PRECISION cs

c===========================================================================
c
c     .. Executable statements ..
c
c===========================================================================

      IF ( n.EQ.1 .OR. ifst.EQ.ilst ) RETURN
      
      IF ( ifst .LT. ilst ) THEN
      
         m1 = 0
         m2 = -1
         m3 = 1
      
      ELSE
      
         m1 = -1
         m2 = 0
         m3 = -1
      
      END IF
      
      DO k = ifst + m1, ilst + m2, m3

         f = MAX( ABS( t( k+1, k+1 ) ), ABS( s( k+1, k+1 ) ) )
         f1 = t( k+1, k+1 ) / f
         f2 = s( k+1, k+1 ) / f
         tlts = .TRUE.

         IF ( ABS( f1 ) .GT. ABS( f2 ) ) tlts = .FALSE.

         c1 = f1 * s( k, k ) - f2 * t( k, k )
         c2 = f1 * s( k, k+1 ) - f2 * t( k, k+1 )

         CALL zlartg( c2, c1, cs, sn, r)
         CALL zrot( k+1, s( 1, k+1 ), 1, s(1,k), 1, cs, sn )
         CALL zrot( k+1, t( 1, k+1 ), 1, t(1,k), 1, cs, sn )
         CALL zrot( n, q( 1, k+1 ), 1, q(1,k), 1, cs, sn )

         IF ( tlts ) THEN
        
            c1 = s( k, k )
            c2 = s( k+1, k ) 

         ELSE

            c1 = t( k, k )
            c2 = t( k+1, k ) 

         END IF

         CALL zlartg( c1, c2, cs, sn, r)
         CALL zrot( n-k+1, s(k,k), ldz, s( k+1, k ), ldz, cs, sn)
         CALL zrot( n-k+1, t(k,k), ldz, t( k+1, k ), ldz, cs, sn)
         CALL zrot( n, z(1,k), 1, z( 1, k+1 ), 1, cs, conjg(sn) )
      
      END DO

c===========================================================================

      END


c===========================================================================


      SUBROUTINE psolve_real( n, x, nq, q, kz, invqkz, ldqkz, ipiv, f )

c===========================================================================
c
c     Coded by Diederik Fokkema
c
c===========================================================================

      IMPLICIT NONE

c===========================================================================
c
c     .. PARAMETERS ..
c
c===========================================================================

      INTEGER n, nq, ldqkz, ipiv(*)
      DOUBLE COMPLEX x(*), q(n,*), kz(n,*),
     $     invqkz(ldqkz,*), f(*)

c===========================================================================
c
c     .. LOCAL .. 
c
c===========================================================================

      INTEGER info
      DOUBLE COMPLEX zero, one
      PARAMETER (zero=(0.0d0,0.0d0), one=(1.0d0,0.0d0))

c===========================================================================
c
c     .. Executable Statements ..
c
c===========================================================================

      CALL precon_real( n, x )
      CALL zgemv( 'c', n, nq, one, q, n, x, 1, zero, f, 1 )
      CALL zgetrs( 'n', nq, 1, invqkz, ldqkz, ipiv, f, ldqkz, info )
      CALL zgemv( 'n', n, nq, -one, kz, n, f, 1, one, x, 1 )

c===========================================================================

      END


c===========================================================================


      SUBROUTINE qzsort_real( ta, tb, k, s, t, z, q, ldz, alpha, beta,
     $     order )

c===========================================================================
c
c     Coded by Diederik Fokkema
c
c     $Id$
c
c     Time-stamp: <95/08/03 23:34:03 caveman>
c
c===========================================================================

      IMPLICIT NONE

c===========================================================================
c
c     .. PARAMETERS ..
c
c===========================================================================

      INTEGER k, ldz, order
      DOUBLE COMPLEX ta, tb, s(ldz,*), t(ldz,*), z(ldz,*),
     $     q(ldz,*), alpha(*), beta(*)

c===========================================================================
c
c     .. LOCAL ..
c
c===========================================================================

      INTEGER i, j, select_real

c===========================================================================
c
c     .. Executable statements ..
c
c===========================================================================

      DO i = 1, k

         DO j = 1, k

            alpha( j ) = s( j, j )
            beta( j )  = t( j, j )

         END DO

         j = select_real( k-i+1, ta, tb, alpha(i), beta(i), order ) +i-1
         CALL myexc_real( k, s, t, z, q, ldz, j, i )

      END DO

c===========================================================================

      END


c===========================================================================


      INTEGER FUNCTION select_real( n, sa, sb, a, b, order )

c===========================================================================
c
c     Coded by Diederik Fokkema
c     Modified Martin van Gijzen, test on division by zero
c
c===========================================================================

      IMPLICIT NONE

c===========================================================================
c
c     .. PARAMETERS ..
c
c===========================================================================

      INTEGER n, order
      DOUBLE COMPLEX sa, sb, a(*), b(*)

c===========================================================================
c
c     .. LOCAL ..
c
c===========================================================================

      INTEGER i, j
      DOUBLE PRECISION dtmp, optval

c===========================================================================
c
c     .. Executable statements ..
c
c===========================================================================

      j = 1

      IF ( order .LE. 0 ) THEN
         optval =  1.0d99
      ELSE
         optval = -1.0d99
      END IF

      IF ( order .EQ. 0 ) THEN

c     
c...  Nearest to target
c     

         DO i = 1, n

            IF ( b(i) .NE. 0.d0 ) THEN
              
               dtmp = ABS( a(i) / b(i) - sa / sb )

               IF ( dtmp .LT. optval ) THEN
                  j = i
                  optval = dtmp
               END IF

            END IF

         END DO

      ELSEIF ( order .EQ. -1 ) THEN

c
c...  Smallest real part
c

         DO i = 1, n

            IF ( b(i) .EQ. 0.d0 ) THEN

               IF ( DBLE( a(i) ) .LT. 0.d0 ) THEN
                  j = i
                  optval = -1.0d99
               END IF

            ELSE

               dtmp = DBLE( a(i) / b(i) )
               IF ( dtmp .LT.optval ) THEN
                  j = i
                  optval = dtmp
               END IF

            END IF

         END DO

      ELSEIF ( order .EQ. 1 ) THEN

c
c...  Largest real part
c

         DO i = 1, n

            IF ( b(i) .EQ. 0.d0 ) THEN
               
               IF ( DBLE( a(i) ) .GT. 0.d0 ) THEN
                  j = i
                  optval = 1.d99
               END IF

            ELSE

               dtmp = DBLE( a(i) / b(i) )
               IF ( dtmp .GT. optval ) THEN
                  j = i
                  optval = dtmp
               END IF

            END IF

         END DO

      ELSEIF ( order .EQ. -2 ) THEN

c
c...  Smallest imaginari part
c

         DO i = 1, n

            IF ( b(i) .EQ. 0.d0 ) THEN
              
               IF ( AIMAG( a(i) ) .LT. 0.d0 ) THEN
                  j = i
                  optval = -1.0d99
               END IF

            ELSE

               dtmp = AIMAG( a(i) / b(i) )
               IF ( dtmp .LT. optval ) THEN
                  j = i
                  optval = dtmp
               END IF

            END IF

         END DO

      ELSEIF ( order .EQ. 2 ) THEN

c
c...  Largest imaginari part
c

         DO i = 1, n

            IF ( b(i) .EQ. 0.d0 ) THEN

               IF ( AIMAG( a(i) ) .GT. 0.d0 ) THEN
                  j = i
                  optval = 1.0d99
               END IF

            ELSE

               dtmp = AIMAG( a(i) / b(i) )
               IF ( dtmp .GT. optval ) THEN
                  j = i
                  optval = dtmp
               END IF

            END IF

         END DO

      ELSE

         CALL error_real ('unknown order in select_real')

      END IF

      select_real = j

c===========================================================================
      
      END


c===========================================================================


      SUBROUTINE zcgstabl_real (n, x, r, l, eps, mxmv,
     $     zalpha, zbeta, nk, kz, qq, invqkz, ldqkz, jpiv, f,
     $     work, lwork)

c===========================================================================
c
c     Programmer: Diederik R. Fokkema
c
c===========================================================================

      IMPLICIT NONE

c===========================================================================
c
c     .. PARAMETERS ..
c
c===========================================================================
      
      INTEGER l, n, mxmv, nk, lwork, ldqkz, jpiv(*)
      DOUBLE PRECISION  eps
      DOUBLE COMPLEX x(*), r(*), kz(n,*), qq(n,*), work(n,*),
     $     zalpha, zbeta, invqkz(ldqkz,*), f(*)

c===========================================================================
c
c     .. LOCAL ..
c
c===========================================================================

      INTEGER mxl
      PARAMETER (mxl = 32)

      INTEGER i, j, k, m, ipiv(mxl), nmv, info

      INTEGER u, q, w, rr, xp, bp
      LOGICAL rcomp, xadd

      DOUBLE PRECISION maxnorm, delta, bpnorm, tr0norm, trlnorm
      DOUBLE COMPLEX varrho, hatgamma

      DOUBLE PRECISION rnrm, rnrm0, eps1, dznrm2
      DOUBLE COMPLEX alpha, beta, omega, gamma, rho0, rho1,
     $     yr(mxl), ys(mxl), z(mxl,mxl), ztmp

      DOUBLE COMPLEX zero, one
      PARAMETER (zero = (0.0d0,0.0d0), one = (1.0d0,0.0d0))

      DOUBLE COMPLEX zdotc

c===========================================================================
c
c     .. Executable statements ..
c
c===========================================================================

      IF ( mxmv .EQ. 0 ) RETURN

      IF ( l .GT. mxl )
     $     CALL error_real ('l exceeds mxl (zcgstabl_real)')

      u  = 1
      q  = u + ( l+1 )
      rr = q + ( l+1 )
      xp = rr + 1
      bp = xp + 1
      w = bp + 1

      IF ( w .GT. lwork )
     $     CALL error_real ('workspace too small (zcgstabl_real)')

c
c     --- set x to zero and compute first residue
c

      CALL zzeros_real( n, x )
      CALL zscal( n, -one, r, 1 )
      CALL zcopy( n, r, 1, work(1,rr), 1 )
      CALL psolve_real( n, r, nk, qq, kz, invqkz, ldqkz, jpiv, f )

c
c     --- initialize loop
c
      nmv = 0

      rnrm0 = dznrm2( n, r, 1 )
      rnrm = rnrm0
      eps1 = eps * rnrm0

      CALL zcopy( n, x, 1, work(1,xp), 1 )
      CALL zcopy( n, r, 1, work(1,bp), 1 )
      maxnorm = 0d0
      bpnorm = rnrm
      rcomp = .FALSE.
      xadd = .FALSE.
      delta = 1.0d-2

      m = 0
      rho0 = one
      alpha = one
      beta = zero
      omega = one

      CALL zzeros_real( n*(l+1), work(1,u) )
      CALL zzeros_real( n*(l+1), work(1,q) )
      CALL zcopy( n, r, 1, work(1,q), 1 )

c
c     --- loop
c

 1000 CONTINUE
      m = m + l

c
c     --- BiCG part
c

      rho0 = -omega * rho0

      DO k = 1, l

         rho1 = zdotc( n, work(1,rr), 1, work( 1, q+k-1 ), 1 )
         beta = alpha * ( rho1 / rho0 )
         rho0 = rho1
         beta = beta

         DO j = 0, k-1

            CALL zxpay_real( n, work( 1, q+j ), 1, (-beta),
     $           work( 1, u+j ), 1 )

         ENDDO

         CALL jdqzmv_real( n, work( 1, u+k-1 ), work( 1, u+k ),
     $        work(1,w), zalpha, zbeta )

         CALL psolve_real( n, work( 1, u+k ), nk, qq, kz,
     $        invqkz, ldqkz, jpiv, f )

         gamma = zdotc( n, work(1,rr), 1, work( 1, u+k ), 1 )
         alpha = rho0 / gamma

         DO j = 0, k-1

            CALL zaxpy (n, (-alpha), work( 1, u+j+1 ),
     $           1, work( 1, q+j ), 1 )

         ENDDO

         CALL jdqzmv_real( n, work( 1, q+k-1 ), work( 1, q+k ),
     $        work(1,w), zalpha, zbeta )

         CALL psolve_real( n, work( 1, q+k ), nk, qq, kz,
     $        invqkz, ldqkz, jpiv, f )

         CALL zaxpy( n, alpha, work(1,u), 1, x, 1 )

         rnrm = dznrm2( n, work(1,q), 1 )
         maxnorm = MAX( maxnorm, rnrm )
         nmv = nmv+2

      ENDDO

c
c     --- MR part + Maintaining the convergence
c

      DO i = 1, l-1

         DO j = 1, i

            ztmp = zdotc( n, work(1,q+i), 1, work( 1, q+j ), 1 )
            z(i,j) = ztmp
            z(j,i) = CONJG( ztmp )

         ENDDO

         yr(i) = zdotc( n, work( 1, q+i ), 1, work(1,q), 1 )
         ys(i) = zdotc( n, work( 1, q+i ), 1, work(1,q+l), 1 )

      ENDDO

      CALL zgetrf( l-1, l-1, z, mxl, ipiv, info )
      CALL zgetrs( 'n', l-1, 1, z, mxl, ipiv, yr, mxl, info )
      CALL zgetrs( 'n', l-1, 1, z, mxl, ipiv, ys, mxl, info )
      CALL zcopy( n, work(1,q), 1, r, 1 )
      CALL zgemv( 'n', n, l-1, (-one), work( 1, q+1 ), n, yr, 1, one,
     $     r, 1 )
      CALL zgemv( 'n', n, l-1, (-one), work( 1, q+1 ), n, ys, 1, one,
     $     work( 1, q+l ), 1 )

      tr0norm = dznrm2( n, r, 1 )
      trlnorm = dznrm2( n, work( 1, q+l ), 1 )
      varrho = zdotc( n, work( 1, q+l ), 1, r, 1 ) / ( tr0norm*trlnorm )
      hatgamma = varrho / ABS( varrho ) * MAX( ABS( varrho ), 7.0d-1 )
      hatgamma = ( tr0norm / trlnorm ) * hatgamma
      yr(l) = zero
      ys(l) = -one

      CALL zaxpy( l, (-hatgamma), ys, 1, yr, 1 )

      omega = yr(l)

      CALL zgemv( 'n', n, l, one, work(1,q), n, yr, 1, one, x, 1 )
      CALL zgemv( 'n', n, l, (-one), work( 1, u+1 ), n, yr, 1, one,
     $     work(1,u), 1 )
      CALL zaxpy( n, (-hatgamma), work(1,q+l), 1, r, 1 )
      CALL zcopy( n, r, 1, work(1,q), 1 )

c
c     --- reliable update
c

      rnrm = dznrm2( n, work(1,q), 1 )
      maxnorm = MAX( maxnorm, rnrm )
      xadd = ( rnrm .LT. delta*rnrm0 .AND. rnrm0 .LT. maxnorm )
      rcomp = ( ( rnrm .LT.delta*maxnorm .AND.
     $     rnrm0 .LT. maxnorm ) .OR. xadd )

      IF ( rcomp ) THEN
         
         CALL jdqzmv_real( n, x, work(1,q), work(1,w),
     $        zalpha, zbeta )
         CALL psolve_real( n, work(1,q), nk, qq, kz,
     $        invqkz, ldqkz, jpiv, f )
         CALL zxpay_real( n, work(1,bp), 1, -one, work(1,q), 1 )
         maxnorm = rnrm
         
         IF ( xadd ) THEN
            CALL zaxpy( n, one, x, 1, work(1,xp), 1 )
            CALL zzeros_real (n, x)
            CALL zcopy( n, work(1,q), 1, work(1,bp), 1 )
            bpnorm = rnrm
         END IF

      END IF
         
      IF ( nmv .LT. mxmv .AND. rnrm .GT. eps1 ) GOTO 1000

      CALL zaxpy( n, one, work(1,xp), 1, x, 1 )

c
c     --- RETURN
c

      mxmv = nmv
      eps = rnrm/rnrm0

      RETURN

c===========================================================================

      END


c===========================================================================


      SUBROUTINE zgmres_real( n, x, r, mxm, eps, mxmv, 
     $     alpha, beta, k, kz, q, invqkz, ldqkz, ipiv, f, v, tp )

c===========================================================================
c
c     Programmer: Diederik R. Fokkema
c
c===========================================================================

      IMPLICIT NONE

c===========================================================================
c
c     .. PARAMETERS ..
c
c===========================================================================

      INTEGER mxm, mxmv, n, k, ldqkz, ipiv(*)
      DOUBLE PRECISION eps
      DOUBLE COMPLEX x(*), r(*), kz(n,*), q(n,*), v(n,*),
     $     alpha, beta, invqkz(ldqkz,*), f(*), tp(*)

c===========================================================================
c
ctex@ \begin{manpage}{ZGMRES} 
ctex@
ctex@ \subtitle{ZGMRES} 
ctex@    ZGMRES -- Generalized Minimal Residual
ctex@    iterative method for solving linear systems $\beta A-\alpha B = -r$.
ctex@    This SUBROUTINE in specilized for use in JDQZ.
ctex@ 
ctex@ \subtitle{Declaration}
ctex@ \function{SUBROUTINE zgmres (n, x, r, mxm, eps, mxmv, a, ka, b, kb,
ctex@   alpha, beta, k, kz, mqkz, zmqkz, ldvs, q,
ctex@   lu, klu, dlu, v)}
ctex@
ctex@ \subtitle{PARAMETERs}
ctex@    \variable{INTEGER n} 
ctex@       On entry, n specifies the dimension of the matrix A.
ctex@       
ctex@    \variable{x} 
ctex@       DOUBLE COMPLEX array of size n. 
ctex@       On exit, x is overwritten by the approximate solution.
ctex@
ctex@    \variable{r}
ctex@       DOUBLE COMPLEX array of size n. Before entry, the array r 
ctex@       must contain the righthandside of the linear problem Ax=r. 
ctex@       Changed on exit.
ctex@ 
ctex@    \variable{INTEGER mxm} 
ctex@       On entry, mxm specifies the degree of the Minimal Residual
ctex@       polynomial.
ctex@
ctex@    \variable{{DOUBLE PRECISION} eps}
ctex@       On entry, eps specifies the stop tolerance. On exit, eps contains
ctex@       the relative norm of the last residual.
ctex@
ctex@    \variable{INTEGER mxmv}
ctex@       On Entry, mxmv specifies the maximum number of matrix 
ctex@       multiplications. On exit, mxmv contains the number of matrix
ctex@       multiplications performed by the method.
ctex@
ctex@    \variable{{DOUBLE COMPLEX} zalpha}
ctex@       On entry, zalpha specifies $\alpha$. Unchanged on exit.
ctex@ 
ctex@    \variable{{DOUBLE COMPLEX} zbeta}
ctex@       On entry, zbeta specifies $\beta$. Unchanged on exit.
ctex@ 
ctex@    \variable{INTEGER k}
ctex@       On entry, k specifies the number of columns of the arrays
ctex@       kz and q.
ctex@
ctex@    \variable{z}
ctex@       DOUBLE COMPLEX array z, of size (n,k). On entry the array z
ctex@       must contain the preconditioned matrix Z.
ctex@
ctex@    \variable{mqkz}
ctex@       DOUBLE COMPLEX array mqkz, of size (ldvs,k). On entry the array 
ctex@       mqkz must contain the matrix Q'*KZ.
ctex@
ctex@    \variable{zmqkz}
ctex@       DOUBLE COMPLEX array zmqkz, of size (ldvs,k). Workspace. Used to
ctex@       copy mqkz.
ctex@
ctex@    \variable{q}
ctex@       DOUBLE COMPLEX array q, of size (n,k). On entry the array q
ctex@       must contain the preconditioned matrix Q.
ctex@
ctex@    \variable{v}
ctex@       DOUBLE COMPLEX array of size (n,mxm+1). Workspace.
ctex@
ctex@ \subtitle{Description}
ctex@    ***
ctex@
ctex@ \subtitle{See Also}
ctex@    ***
ctex@
ctex@ \subtitle{References}
ctex@    ***
ctex@ 
ctex@ \subtitle{Bugs}
ctex@    ***
ctex@
ctex@ \subtitle{Author}
ctex@     Diederik R.\ Fokkema
ctex@
ctex@ \end{manpage}
ctex@ \begin{verbatim}
ctex@    % actual code
ctex@ \end{verbatim}
ctex@
c
c===========================================================================

c===========================================================================
c
c     .. LOCAL ..
c
c===========================================================================
    
      LOGICAL restrt, loop
      INTEGER maxm
      PARAMETER (maxm = 100)

      INTEGER i, m, m1, nmv
      DOUBLE PRECISION rnrm0, rnrm, eps1, c(maxm)
      DOUBLE COMPLEX hh(maxm,maxm-1), rs(maxm), s(maxm), y(maxm), rcs
      DOUBLE COMPLEX zero, one
      PARAMETER (zero = (0.0d0,0.0d0), one = (1.0d0,0.0d0))

      DOUBLE PRECISION  dznrm2

      DOUBLE COMPLEX ztmp, zdotc

c===========================================================================
c
c     .. Executable Statements ..
c
c===========================================================================

      IF ( mxm .GT. maxm-1 )
     $     CALL error_real ('mxm larger than maxm (zgmres_real)')

c
c     --- Initialize first residue
c

      CALL zzeros_real( n, x )
      CALL zscal( n, -one, r, 1 )
      CALL psolve_real( n, r, k, q, kz, invqkz, ldqkz, ipiv, f )

c
c     --- initialize loop
c

      rnrm0  = dznrm2( n, r, 1 )
      rnrm = rnrm0
      eps1  = eps * rnrm

      nmv = 0

      CALL zcopy( n, r, 1, v(1,1), 1 )

c         
c     --- outer loop
c

 1000 restrt = ( nmv .LT. mxmv .AND. rnrm .GT. eps1 )

      IF ( restrt ) THEN

         ztmp = one / rnrm
         CALL zscal( n, ztmp, v(1,1), 1 )
         rs(1) = rnrm

c
c     --- inner loop
c

         m = 0

 2000    loop = ( nmv .LT. mxmv .AND. m .LT. mxm .AND. rnrm .GT. eps1 )
        
         IF ( loop ) THEN
         
            m  = m + 1
            m1 = m + 1
            CALL jdqzmv_real( n, v(1,m), v(1,m1), tp, alpha, beta )
            CALL psolve_real( n, v(1,m1), k, q, kz, invqkz,
     $           ldqkz, ipiv, f )
            nmv = nmv + 1 
 
            DO i = 1, m

               ztmp = zdotc( n, v(1,i), 1, v(1,m1), 1 )
               hh(i,m) = ztmp
               CALL zaxpy( n, (-ztmp), v(1,i), 1, v(1,m1), 1 )

            ENDDO

            ztmp = dznrm2( n, v(1,m1), 1 )
            hh( m1, m ) = ztmp

            CALL zscal( n, (one/ztmp), v(1,m1), 1 )

            DO i = 1, m-1
               CALL zrot( 1, hh(i,m), 1, hh( i+1, m ), 1, c(i), s(i) )
            ENDDO

            CALL zlartg( hh(m,m), hh(m1,m), c(m), s(m), rcs )
           
            hh(m,m) = rcs
            hh(m1,m) = zero
            rs(m1) = zero
            
            CALL zrot(  1, rs(m), 1, rs(m1), 1, c(m), s(m) )
            rnrm = ABS( rs(m1) )
            GOTO 2000

         END IF

c
c     --- compute approximate solution x
c

         CALL zcopy( m, rs, 1, y, 1 )
         CALL ztrsv( 'u', 'n', 'n', m, hh, maxm, y, 1 )
         CALL zgemv( 'n', n, m, one, v, n, y, 1, one, x, 1 )

c
c     --- compute residual for restart
c

         CALL jdqzmv_real( n, x, v(1,2), tp, alpha, beta )
         CALL psolve_real( n, v(1,2), k, q, kz, invqkz,
     $        ldqkz, ipiv, f )
         CALL zcopy( n, r, 1, v(1,1), 1 )
         CALL zaxpy( n, -one, v(1,2), 1, v(1,1), 1 )
         rnrm = dznrm2( n, v(1,1), 1 )

         GOTO 1000

      END IF

c
c     --- RETURN
c

      eps = rnrm / rnrm0
      mxmv = nmv

      RETURN

c===========================================================================

      END


c===========================================================================


      SUBROUTINE zmgs_real (n, k, v, w, job )

c===========================================================================
c
c     Coded by Diederik Fokkema
c     Modified 05-23-96: M. Kooper: job =1 corrected, array YWORK added
c
c===========================================================================

      IMPLICIT NONE

c===========================================================================
c
c     .. PARAMETERS ..
c
c===========================================================================
            
      INTEGER n, k, job
      DOUBLE COMPLEX v(n,*), w(*)

c===========================================================================
c
c     .. LOCAL ..
c
c===========================================================================
      
      INTEGER i
      DOUBLE PRECISION s0, s1, dznrm2
      DOUBLE COMPLEX znrm

c===========================================================================
c
c     .. Executable statements ..
c
c===========================================================================

      s1 = dznrm2( n, w, 1 )

      DO i = 1, k
	 s0 = s1
	 CALL zortho_real( n, v(1,i), w, s0, s1, znrm )
      END DO

      IF ( job .EQ. 0 ) THEN
	 RETURN
      ELSE
         znrm  = 1.0d0 / s1
         CALL zscal( n, znrm, w, 1 )
      END IF

      RETURN

c===========================================================================

      END


c===========================================================================


      SUBROUTINE zortho_real( n, v, w, s0, s1, znrm )

c===========================================================================
c
c     Coded by Diederik Fokkema
c
c===========================================================================

      IMPLICIT NONE

c===========================================================================
c
c     .. PARAMETERS ..
c
c===========================================================================

      INTEGER n
      DOUBLE PRECISION s0, s1
      DOUBLE COMPLEX v(*), w(*), znrm

c===========================================================================
c
c     .. LOCAL ..
c
c===========================================================================

      DOUBLE PRECISION kappa, dznrm2
      DOUBLE COMPLEX ztmp, zdotc
      PARAMETER (kappa=1d2)

c===========================================================================
c
c     .. Executable statements ..
c
c===========================================================================

      znrm = zdotc( n, v, 1, w, 1 )
      CALL zaxpy( n, (-znrm), v, 1, w, 1 )
      s1 = dznrm2( n, w, 1 )

      IF ( s1 .GT. s0 / kappa ) THEN
	 GOTO 100
      ELSE
	 s0 = s1
	 ztmp = zdotc( n, v, 1, w, 1 )
         znrm = znrm + ztmp
	 CALL zaxpy( n, (-ztmp), v, 1, w, 1 )
	 s1 = dznrm2( n, w, 1 )
	 IF ( s1 .GT. s0 / kappa ) THEN
	    GOTO 100
	 ELSE
	    CALL error_real ('zero vector in zmgs_real')
	 END IF

      END IF

 100  CONTINUE

      RETURN

c===========================================================================-

      END


c===========================================================================


      SUBROUTINE zones_real( n, x )

c===========================================================================
c
c     Coded by Diederik Fokkema
c
c     $Id$
c
c     Time-stamp: <95/07/30 21:45:46 caveman>
c
c===========================================================================

c===========================================================================
c
c     .. PARAMETERS ..
c
c===========================================================================

      IMPLICIT NONE
      INTEGER n
      DOUBLE COMPLEX x(*)

c===========================================================================
c
c     .. LOCAL ..
c
c===========================================================================

      INTEGER i

c===========================================================================
c
c     .. Executable statements ..
c
c===========================================================================

      DO i = 1, n
         x(i) = ( 1.0d0, 0.0d0 )
      END DO

c===========================================================================c

      END
      

c===========================================================================


      SUBROUTINE zxpay_real( n, dx, incx, da, dy, incy )

c===========================================================================
c
c     modified by:  D.R. Fokkema
c     01/06/94
c
c     a vector plus constant times a vector.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
c===========================================================================

      IMPLICIT NONE

c===========================================================================
c
c     .. PARAMETERS ..
c
c===========================================================================

      DOUBLE COMPLEX dx(*),dy(*),da
      INTEGER i,incx,incy,ix,iy,m,mp1,n

c===========================================================================
      
      IF ( n .LE. 0 ) RETURN
      IF ( incx .EQ. 1 .AND. incy .EQ. 1 ) GOTO 20

c     
c     code for unequal increments or equal increments
c     not equal to 1
c
      
      ix = 1
      iy = 1

      IF ( incx .LT. 0 ) ix = ( -n + 1 ) * incx + 1
      IF ( incy .LT. 0 ) iy = ( -n + 1 ) * incy + 1
      
      DO 10 i = 1, n
        dy(iy) = da * dy(iy) + dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 CONTINUE

      RETURN

c     
c     code for both increments equal to 1
c     
c     
c     clean-up loop
c     

   20 m = mod( n, 4 )

      IF ( m .EQ. 0 ) GOTO 40

      DO 30 i = 1, m
        dy(i) = da * dy(i) + dx(i)
   30 CONTINUE

      IF ( n .LT. 4 ) RETURN

   40 mp1 = m + 1

      DO 50 i = mp1, n, 4
        dy(i) = da * dy(i) + dx(i)
        dy( i + 1 ) = da * dy( i + 1 ) + dx( i + 1 )
        dy( i + 2 ) = da * dy( i + 2 ) + dx( i + 2 )
        dy( i + 3 ) = da * dy( i + 3 ) + dx( i + 3 )
   50 CONTINUE

      RETURN

c===========================================================================
      
      END


c===========================================================================



      SUBROUTINE zzeros_real( n, x )

c===========================================================================
c
c     Coded by Diederik Fokkema
c
c     $Id$
c
c     Time-stamp: <95/07/30 21:48:00 caveman>
c
c===========================================================================

      IMPLICIT NONE

c===========================================================================
c
c     .. PARAMETERS ..
c
c===========================================================================

      INTEGER n
      DOUBLE COMPLEX x(*)

c===========================================================================
c
c     .. LOCAL ..
c
c===========================================================================

      INTEGER i

c===========================================================================
c
c     .. Executable statements ..
c
c===========================================================================

      DO i = 1, n
         x( i ) = ( 0.0d0, 0.0d0 ) 
      END DO

c===========================================================================

      END

