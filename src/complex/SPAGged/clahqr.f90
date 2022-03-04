!*==clahqr.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b CLAHQR computes the eigenvalues and Schur factorization of an upper Hessenberg matrix, using the double-shift/single-shift QR algorithm.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CLAHQR + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clahqr.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clahqr.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clahqr.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CLAHQR( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILOZ,
!                          IHIZ, Z, LDZ, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, N
!       LOGICAL            WANTT, WANTZ
!       ..
!       .. Array Arguments ..
!       COMPLEX            H( LDH, * ), W( * ), Z( LDZ, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    CLAHQR is an auxiliary routine called by CHSEQR to update the
!>    eigenvalues and Schur decomposition already computed by CHSEQR, by
!>    dealing with the Hessenberg submatrix in rows and columns ILO to
!>    IHI.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] WANTT
!> \verbatim
!>          WANTT is LOGICAL
!>          = .TRUE. : the full Schur form T is required;
!>          = .FALSE.: only eigenvalues are required.
!> \endverbatim
!>
!> \param[in] WANTZ
!> \verbatim
!>          WANTZ is LOGICAL
!>          = .TRUE. : the matrix of Schur vectors Z is required;
!>          = .FALSE.: Schur vectors are not required.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix H.  N >= 0.
!> \endverbatim
!>
!> \param[in] ILO
!> \verbatim
!>          ILO is INTEGER
!> \endverbatim
!>
!> \param[in] IHI
!> \verbatim
!>          IHI is INTEGER
!>          It is assumed that H is already upper triangular in rows and
!>          columns IHI+1:N, and that H(ILO,ILO-1) = 0 (unless ILO = 1).
!>          CLAHQR works primarily with the Hessenberg submatrix in rows
!>          and columns ILO to IHI, but applies transformations to all of
!>          H if WANTT is .TRUE..
!>          1 <= ILO <= max(1,IHI); IHI <= N.
!> \endverbatim
!>
!> \param[in,out] H
!> \verbatim
!>          H is COMPLEX array, dimension (LDH,N)
!>          On entry, the upper Hessenberg matrix H.
!>          On exit, if INFO is zero and if WANTT is .TRUE., then H
!>          is upper triangular in rows and columns ILO:IHI.  If INFO
!>          is zero and if WANTT is .FALSE., then the contents of H
!>          are unspecified on exit.  The output state of H in case
!>          INF is positive is below under the description of INFO.
!> \endverbatim
!>
!> \param[in] LDH
!> \verbatim
!>          LDH is INTEGER
!>          The leading dimension of the array H. LDH >= max(1,N).
!> \endverbatim
!>
!> \param[out] W
!> \verbatim
!>          W is COMPLEX array, dimension (N)
!>          The computed eigenvalues ILO to IHI are stored in the
!>          corresponding elements of W. If WANTT is .TRUE., the
!>          eigenvalues are stored in the same order as on the diagonal
!>          of the Schur form returned in H, with W(i) = H(i,i).
!> \endverbatim
!>
!> \param[in] ILOZ
!> \verbatim
!>          ILOZ is INTEGER
!> \endverbatim
!>
!> \param[in] IHIZ
!> \verbatim
!>          IHIZ is INTEGER
!>          Specify the rows of Z to which transformations must be
!>          applied if WANTZ is .TRUE..
!>          1 <= ILOZ <= ILO; IHI <= IHIZ <= N.
!> \endverbatim
!>
!> \param[in,out] Z
!> \verbatim
!>          Z is COMPLEX array, dimension (LDZ,N)
!>          If WANTZ is .TRUE., on entry Z must contain the current
!>          matrix Z of transformations accumulated by CHSEQR, and on
!>          exit Z has been updated; transformations are applied only to
!>          the submatrix Z(ILOZ:IHIZ,ILO:IHI).
!>          If WANTZ is .FALSE., Z is not referenced.
!> \endverbatim
!>
!> \param[in] LDZ
!> \verbatim
!>          LDZ is INTEGER
!>          The leading dimension of the array Z. LDZ >= max(1,N).
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>           = 0:  successful exit
!>           > 0:  if INFO = i, CLAHQR failed to compute all the
!>                  eigenvalues ILO to IHI in a total of 30 iterations
!>                  per eigenvalue; elements i+1:ihi of W contain
!>                  those eigenvalues which have been successfully
!>                  computed.
!>
!>                  If INFO > 0 and WANTT is .FALSE., then on exit,
!>                  the remaining unconverged eigenvalues are the
!>                  eigenvalues of the upper Hessenberg matrix
!>                  rows and columns ILO through INFO of the final,
!>                  output value of H.
!>
!>                  If INFO > 0 and WANTT is .TRUE., then on exit
!>          (*)       (initial value of H)*U  = U*(final value of H)
!>                  where U is an orthogonal matrix.    The final
!>                  value of H is upper Hessenberg and triangular in
!>                  rows and columns INFO+1 through IHI.
!>
!>                  If INFO > 0 and WANTZ is .TRUE., then on exit
!>                      (final value of Z)  = (initial value of Z)*U
!>                  where U is the orthogonal matrix in (*)
!>                  (regardless of the value of WANTT.)
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \date December 2016
!
!> \ingroup complexOTHERauxiliary
!
!> \par Contributors:
!  ==================
!>
!> \verbatim
!>
!>     02-96 Based on modifications by
!>     David Day, Sandia National Laboratory, USA
!>
!>     12-04 Further modifications by
!>     Ralph Byers, University of Kansas, USA
!>     This is a modified version of CLAHQR from LAPACK version 3.0.
!>     It is (1) more robust against overflow and underflow and
!>     (2) adopts the more conservative Ahues & Tisseur stopping
!>     criterion (LAWN 122, 1997).
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE CLAHQR(Wantt,Wantz,N,Ilo,Ihi,H,Ldh,W,Iloz,Ihiz,Z,Ldz,  &
     &                  Info)
      USE S_CCOPY
      USE S_CLADIV
      USE S_CLARFG
      USE S_CSCAL
      USE S_SLABAD
      USE S_SLAMCH
      IMPLICIT NONE
!*--CLAHQR205
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E0,0.0E0) ,                   &
     &                         ONE = (1.0E0,0.0E0)
      REAL , PARAMETER  ::  RZERO = 0.0E0 , RONE = 1.0E0 ,              &
     &                      HALF = 0.5E0 , DAT1 = 3.0E0/4.0E0
      INTEGER , PARAMETER  ::  KEXSH = 10
!
! Dummy argument declarations rewritten by SPAG
!
      LOGICAL , INTENT(IN) :: Wantt
      LOGICAL , INTENT(IN) :: Wantz
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Ilo
      INTEGER , INTENT(IN) :: Ihi
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldh,*) :: H
      INTEGER :: Ldh
      COMPLEX , INTENT(OUT) , DIMENSION(*) :: W
      INTEGER , INTENT(IN) :: Iloz
      INTEGER , INTENT(IN) :: Ihiz
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldz,*) :: Z
      INTEGER , INTENT(IN) :: Ldz
      INTEGER , INTENT(OUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL :: aa , ab , ba , bb , h10 , h21 , rtemp , s , safmax ,      &
     &        safmin , smlnum , sx , t2 , tst , ulp
      REAL :: CABS1
      COMPLEX :: cdum , h11 , h11s , h22 , sc , sum , t , t1 , temp ,   &
     &           u , v2 , x , y
      INTEGER :: i , i1 , i2 , itmax , its , j , jhi , jlo , k , kdefl ,&
     &           l , m , nh , nz
      COMPLEX , DIMENSION(2) :: v
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
!  =========================================================
!
!     .. Parameters ..
!     ..
!     .. Local Scalars ..
!     ..
!     .. Local Arrays ..
!     ..
!     .. External Functions ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Statement Functions ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Statement Function definitions ..
      CABS1(cdum) = ABS(REAL(cdum)) + ABS(AIMAG(cdum))
!     ..
!     .. Executable Statements ..
!
      Info = 0
!
!     Quick return if possible
!
      IF ( N==0 ) RETURN
      IF ( Ilo==Ihi ) THEN
         W(Ilo) = H(Ilo,Ilo)
         RETURN
      ENDIF
!
!     ==== clear out the trash ====
      DO j = Ilo , Ihi - 3
         H(j+2,j) = ZERO
         H(j+3,j) = ZERO
      ENDDO
      IF ( Ilo<=Ihi-2 ) H(Ihi,Ihi-2) = ZERO
!     ==== ensure that subdiagonal entries are real ====
      IF ( Wantt ) THEN
         jlo = 1
         jhi = N
      ELSE
         jlo = Ilo
         jhi = Ihi
      ENDIF
      DO i = Ilo + 1 , Ihi
         IF ( AIMAG(H(i,i-1))/=RZERO ) THEN
!           ==== The following redundant normalization
!           .    avoids problems with both gradual and
!           .    sudden underflow in ABS(H(I,I-1)) ====
            sc = H(i,i-1)/CABS1(H(i,i-1))
            sc = CONJG(sc)/ABS(sc)
            H(i,i-1) = ABS(H(i,i-1))
            CALL CSCAL(jhi-i+1,sc,H(i,i),Ldh)
            CALL CSCAL(MIN(jhi,i+1)-jlo+1,CONJG(sc),H(jlo,i),1)
            IF ( Wantz ) CALL CSCAL(Ihiz-Iloz+1,CONJG(sc),Z(Iloz,i),1)
         ENDIF
      ENDDO
!
      nh = Ihi - Ilo + 1
      nz = Ihiz - Iloz + 1
!
!     Set machine-dependent constants for the stopping criterion.
!
      safmin = SLAMCH('SAFE MINIMUM')
      safmax = RONE/safmin
      CALL SLABAD(safmin,safmax)
      ulp = SLAMCH('PRECISION')
      smlnum = safmin*(REAL(nh)/ulp)
!
!     I1 and I2 are the indices of the first row and last column of H
!     to which transformations must be applied. If eigenvalues only are
!     being computed, I1 and I2 are set inside the main loop.
!
      IF ( Wantt ) THEN
         i1 = 1
         i2 = N
      ENDIF
!
!     ITMAX is the total number of QR iterations allowed.
!
      itmax = 30*MAX(10,nh)
!
!     KDEFL counts the number of iterations since a deflation
!
      kdefl = 0
!
!     The main loop begins here. I is the loop index and decreases from
!     IHI to ILO in steps of 1. Each iteration of the loop works
!     with the active submatrix in rows and columns L to I.
!     Eigenvalues I+1 to IHI have already converged. Either L = ILO, or
!     H(L,L-1) is negligible so that the matrix splits.
!
      i = Ihi
 100  IF ( i<Ilo ) GOTO 99999
!
!     Perform QR iterations on rows and columns ILO to I until a
!     submatrix of order 1 splits off at the bottom because a
!     subdiagonal element has become negligible.
!
      l = Ilo
      DO its = 0 , itmax
!
!        Look for a single small subdiagonal element.
!
         DO k = i , l + 1 , -1
            IF ( CABS1(H(k,k-1))<=smlnum ) EXIT
            tst = CABS1(H(k-1,k-1)) + CABS1(H(k,k))
            IF ( tst==ZERO ) THEN
               IF ( k-2>=Ilo ) tst = tst + ABS(REAL(H(k-1,k-2)))
               IF ( k+1<=Ihi ) tst = tst + ABS(REAL(H(k+1,k)))
            ENDIF
!           ==== The following is a conservative small subdiagonal
!           .    deflation criterion due to Ahues & Tisseur (LAWN 122,
!           .    1997). It has better mathematical foundation and
!           .    improves accuracy in some examples.  ====
            IF ( ABS(REAL(H(k,k-1)))<=ulp*tst ) THEN
               ab = MAX(CABS1(H(k,k-1)),CABS1(H(k-1,k)))
               ba = MIN(CABS1(H(k,k-1)),CABS1(H(k-1,k)))
               aa = MAX(CABS1(H(k,k)),CABS1(H(k-1,k-1)-H(k,k)))
               bb = MIN(CABS1(H(k,k)),CABS1(H(k-1,k-1)-H(k,k)))
               s = aa + ab
               IF ( ba*(ab/s)<=MAX(smlnum,ulp*(bb*(aa/s))) ) EXIT
            ENDIF
         ENDDO
         l = k
!
!           H(L,L-1) is negligible
!
         IF ( l>Ilo ) H(l,l-1) = ZERO
!
!        Exit from loop if a submatrix of order 1 has split off.
!
         IF ( l>=i ) GOTO 200
         kdefl = kdefl + 1
!
!        Now the active submatrix is in rows and columns L to I. If
!        eigenvalues only are being computed, only the active submatrix
!        need be transformed.
!
         IF ( .NOT.Wantt ) THEN
            i1 = l
            i2 = i
         ENDIF
!
         IF ( MOD(kdefl,2*KEXSH)==0 ) THEN
!
!           Exceptional shift.
!
            s = DAT1*ABS(REAL(H(i,i-1)))
            t = s + H(i,i)
         ELSEIF ( MOD(kdefl,KEXSH)==0 ) THEN
!
!           Exceptional shift.
!
            s = DAT1*ABS(REAL(H(l+1,l)))
            t = s + H(l,l)
         ELSE
!
!           Wilkinson's shift.
!
            t = H(i,i)
            u = SQRT(H(i-1,i))*SQRT(H(i,i-1))
            s = CABS1(u)
            IF ( s/=RZERO ) THEN
               x = HALF*(H(i-1,i-1)-t)
               sx = CABS1(x)
               s = MAX(s,CABS1(x))
               y = s*SQRT((x/s)**2+(u/s)**2)
               IF ( sx>RZERO ) THEN
                  IF ( REAL(x/sx)*REAL(y)+AIMAG(x/sx)*AIMAG(y)<RZERO )  &
     &                 y = -y
               ENDIF
               t = t - u*CLADIV(u,(x+y))
            ENDIF
         ENDIF
!
!        Look for two consecutive small subdiagonal elements.
!
         DO m = i - 1 , l + 1 , -1
!
!           Determine the effect of starting the single-shift QR
!           iteration at row M, and see if this would make H(M,M-1)
!           negligible.
!
            h11 = H(m,m)
            h22 = H(m+1,m+1)
            h11s = h11 - t
            h21 = REAL(H(m+1,m))
            s = CABS1(h11s) + ABS(h21)
            h11s = h11s/s
            h21 = h21/s
            v(1) = h11s
            v(2) = h21
            h10 = REAL(H(m,m-1))
            IF ( ABS(h10)*ABS(h21)                                      &
     &           <=ulp*(CABS1(h11s)*(CABS1(h11)+CABS1(h22))) ) GOTO 150
         ENDDO
         h11 = H(l,l)
         h22 = H(l+1,l+1)
         h11s = h11 - t
         h21 = REAL(H(l+1,l))
         s = CABS1(h11s) + ABS(h21)
         h11s = h11s/s
         h21 = h21/s
         v(1) = h11s
         v(2) = h21
!
!        Single-shift QR step
!
 150     DO k = m , i - 1
!
!           The first iteration of this loop determines a reflection G
!           from the vector V and applies it from left and right to H,
!           thus creating a nonzero bulge below the subdiagonal.
!
!           Each subsequent iteration determines a reflection G to
!           restore the Hessenberg form in the (K-1)th column, and thus
!           chases the bulge one step toward the bottom of the active
!           submatrix.
!
!           V(2) is always real before the call to CLARFG, and hence
!           after the call T2 ( = T1*V(2) ) is also real.
!
            IF ( k>m ) CALL CCOPY(2,H(k,k-1),1,v,1)
            CALL CLARFG(2,v(1),v(2),1,t1)
            IF ( k>m ) THEN
               H(k,k-1) = v(1)
               H(k+1,k-1) = ZERO
            ENDIF
            v2 = v(2)
            t2 = REAL(t1*v2)
!
!           Apply G from the left to transform the rows of the matrix
!           in columns K to I2.
!
            DO j = k , i2
               sum = CONJG(t1)*H(k,j) + t2*H(k+1,j)
               H(k,j) = H(k,j) - sum
               H(k+1,j) = H(k+1,j) - sum*v2
            ENDDO
!
!           Apply G from the right to transform the columns of the
!           matrix in rows I1 to min(K+2,I).
!
            DO j = i1 , MIN(k+2,i)
               sum = t1*H(j,k) + t2*H(j,k+1)
               H(j,k) = H(j,k) - sum
               H(j,k+1) = H(j,k+1) - sum*CONJG(v2)
            ENDDO
!
            IF ( Wantz ) THEN
!
!              Accumulate transformations in the matrix Z
!
               DO j = Iloz , Ihiz
                  sum = t1*Z(j,k) + t2*Z(j,k+1)
                  Z(j,k) = Z(j,k) - sum
                  Z(j,k+1) = Z(j,k+1) - sum*CONJG(v2)
               ENDDO
            ENDIF
!
            IF ( k==m .AND. m>l ) THEN
!
!              If the QR step was started at row M > L because two
!              consecutive small subdiagonals were found, then extra
!              scaling must be performed to ensure that H(M,M-1) remains
!              real.
!
               temp = ONE - t1
               temp = temp/ABS(temp)
               H(m+1,m) = H(m+1,m)*CONJG(temp)
               IF ( m+2<=i ) H(m+2,m+1) = H(m+2,m+1)*temp
               DO j = m , i
                  IF ( j/=m+1 ) THEN
                     IF ( i2>j ) CALL CSCAL(i2-j,temp,H(j,j+1),Ldh)
                     CALL CSCAL(j-i1,CONJG(temp),H(i1,j),1)
                     IF ( Wantz ) CALL CSCAL(nz,CONJG(temp),Z(Iloz,j),1)
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
!
!        Ensure that H(I,I-1) is real.
!
         temp = H(i,i-1)
         IF ( AIMAG(temp)/=RZERO ) THEN
            rtemp = ABS(temp)
            H(i,i-1) = rtemp
            temp = temp/rtemp
            IF ( i2>i ) CALL CSCAL(i2-i,CONJG(temp),H(i,i+1),Ldh)
            CALL CSCAL(i-i1,temp,H(i1,i),1)
            IF ( Wantz ) CALL CSCAL(nz,temp,Z(Iloz,i),1)
         ENDIF
!
      ENDDO
!
!     Failure to converge in remaining number of iterations
!
      Info = i
      RETURN
!
!
!     H(I,I-1) is negligible: one eigenvalue has converged.
!
 200  W(i) = H(i,i)
!     reset deflation counter
      kdefl = 0
!
!     return to start of the main loop with new value of I.
!
      i = l - 1
      GOTO 100
!
!
!     End of CLAHQR
!
99999 END SUBROUTINE CLAHQR
