!*==dlahqr.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b DLAHQR computes the eigenvalues and Schur factorization of an upper Hessenberg matrix, using the double-shift/single-shift QR algorithm.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLAHQR + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlahqr.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlahqr.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlahqr.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLAHQR( WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI,
!                          ILOZ, IHIZ, Z, LDZ, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, N
!       LOGICAL            WANTT, WANTZ
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   H( LDH, * ), WI( * ), WR( * ), Z( LDZ, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    DLAHQR is an auxiliary routine called by DHSEQR to update the
!>    eigenvalues and Schur decomposition already computed by DHSEQR, by
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
!>          It is assumed that H is already upper quasi-triangular in
!>          rows and columns IHI+1:N, and that H(ILO,ILO-1) = 0 (unless
!>          ILO = 1). DLAHQR works primarily with the Hessenberg
!>          submatrix in rows and columns ILO to IHI, but applies
!>          transformations to all of H if WANTT is .TRUE..
!>          1 <= ILO <= max(1,IHI); IHI <= N.
!> \endverbatim
!>
!> \param[in,out] H
!> \verbatim
!>          H is DOUBLE PRECISION array, dimension (LDH,N)
!>          On entry, the upper Hessenberg matrix H.
!>          On exit, if INFO is zero and if WANTT is .TRUE., H is upper
!>          quasi-triangular in rows and columns ILO:IHI, with any
!>          2-by-2 diagonal blocks in standard form. If INFO is zero
!>          and WANTT is .FALSE., the contents of H are unspecified on
!>          exit.  The output state of H if INFO is nonzero is given
!>          below under the description of INFO.
!> \endverbatim
!>
!> \param[in] LDH
!> \verbatim
!>          LDH is INTEGER
!>          The leading dimension of the array H. LDH >= max(1,N).
!> \endverbatim
!>
!> \param[out] WR
!> \verbatim
!>          WR is DOUBLE PRECISION array, dimension (N)
!> \endverbatim
!>
!> \param[out] WI
!> \verbatim
!>          WI is DOUBLE PRECISION array, dimension (N)
!>          The real and imaginary parts, respectively, of the computed
!>          eigenvalues ILO to IHI are stored in the corresponding
!>          elements of WR and WI. If two eigenvalues are computed as a
!>          complex conjugate pair, they are stored in consecutive
!>          elements of WR and WI, say the i-th and (i+1)th, with
!>          WI(i) > 0 and WI(i+1) < 0. If WANTT is .TRUE., the
!>          eigenvalues are stored in the same order as on the diagonal
!>          of the Schur form returned in H, with WR(i) = H(i,i), and, if
!>          H(i:i+1,i:i+1) is a 2-by-2 diagonal block,
!>          WI(i) = sqrt(H(i+1,i)*H(i,i+1)) and WI(i+1) = -WI(i).
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
!>          Z is DOUBLE PRECISION array, dimension (LDZ,N)
!>          If WANTZ is .TRUE., on entry Z must contain the current
!>          matrix Z of transformations accumulated by DHSEQR, and on
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
!>           > 0:  If INFO = i, DLAHQR failed to compute all the
!>                  eigenvalues ILO to IHI in a total of 30 iterations
!>                  per eigenvalue; elements i+1:ihi of WR and WI
!>                  contain those eigenvalues which have been
!>                  successfully computed.
!>
!>                  If INFO > 0 and WANTT is .FALSE., then on exit,
!>                  the remaining unconverged eigenvalues are the
!>                  eigenvalues of the upper Hessenberg matrix rows
!>                  and columns ILO through INFO of the final, output
!>                  value of H.
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
!> \ingroup doubleOTHERauxiliary
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>     02-96 Based on modifications by
!>     David Day, Sandia National Laboratory, USA
!>
!>     12-04 Further modifications by
!>     Ralph Byers, University of Kansas, USA
!>     This is a modified version of DLAHQR from LAPACK version 3.0.
!>     It is (1) more robust against overflow and underflow and
!>     (2) adopts the more conservative Ahues & Tisseur stopping
!>     criterion (LAWN 122, 1997).
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE DLAHQR(Wantt,Wantz,N,Ilo,Ihi,H,Ldh,Wr,Wi,Iloz,Ihiz,Z,  &
     &                  Ldz,Info)
      IMPLICIT NONE
!*--DLAHQR211
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Ihi , Ihiz , Ilo , Iloz , Info , Ldh , Ldz , N
      LOGICAL Wantt , Wantz
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION H(Ldh,*) , Wi(*) , Wr(*) , Z(Ldz,*)
!     ..
!
!  =========================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE , TWO
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0)
      DOUBLE PRECISION DAT1 , DAT2
      PARAMETER (DAT1=3.0D0/4.0D0,DAT2=-0.4375D0)
      INTEGER KEXSH
      PARAMETER (KEXSH=10)
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION aa , ab , ba , bb , cs , det , h11 , h12 , h21 , &
     &                 h21s , h22 , rt1i , rt1r , rt2i , rt2r , rtdisc ,&
     &                 s , safmax , safmin , smlnum , sn , sum , t1 ,   &
     &                 t2 , t3 , tr , tst , ulp , v2 , v3
      INTEGER i , i1 , i2 , its , itmax , j , k , l , m , nh , nr , nz ,&
     &        kdefl
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION v(3)
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMCH
      EXTERNAL DLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL DCOPY , DLABAD , DLANV2 , DLARFG , DROT
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , DBLE , MAX , MIN , SQRT
!     ..
!     .. Executable Statements ..
!
      Info = 0
!
!     Quick return if possible
!
      IF ( N==0 ) RETURN
      IF ( Ilo==Ihi ) THEN
         Wr(Ilo) = H(Ilo,Ilo)
         Wi(Ilo) = ZERO
         RETURN
      ENDIF
!
!     ==== clear out the trash ====
      DO j = Ilo , Ihi - 3
         H(j+2,j) = ZERO
         H(j+3,j) = ZERO
      ENDDO
      IF ( Ilo<=Ihi-2 ) H(Ihi,Ihi-2) = ZERO
!
      nh = Ihi - Ilo + 1
      nz = Ihiz - Iloz + 1
!
!     Set machine-dependent constants for the stopping criterion.
!
      safmin = DLAMCH('SAFE MINIMUM')
      safmax = ONE/safmin
      CALL DLABAD(safmin,safmax)
      ulp = DLAMCH('PRECISION')
      smlnum = safmin*(DBLE(nh)/ulp)
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
!     IHI to ILO in steps of 1 or 2. Each iteration of the loop works
!     with the active submatrix in rows and columns L to I.
!     Eigenvalues I+1 to IHI have already converged. Either L = ILO or
!     H(L,L-1) is negligible so that the matrix splits.
!
      i = Ihi
 100  l = Ilo
      IF ( i<Ilo ) GOTO 99999
!
!     Perform QR iterations on rows and columns ILO to I until a
!     submatrix of order 1 or 2 splits off at the bottom because a
!     subdiagonal element has become negligible.
!
      DO its = 0 , itmax
!
!        Look for a single small subdiagonal element.
!
         DO k = i , l + 1 , -1
            IF ( ABS(H(k,k-1))<=smlnum ) EXIT
            tst = ABS(H(k-1,k-1)) + ABS(H(k,k))
            IF ( tst==ZERO ) THEN
               IF ( k-2>=Ilo ) tst = tst + ABS(H(k-1,k-2))
               IF ( k+1<=Ihi ) tst = tst + ABS(H(k+1,k))
            ENDIF
!           ==== The following is a conservative small subdiagonal
!           .    deflation  criterion due to Ahues & Tisseur (LAWN 122,
!           .    1997). It has better mathematical foundation and
!           .    improves accuracy in some cases.  ====
            IF ( ABS(H(k,k-1))<=ulp*tst ) THEN
               ab = MAX(ABS(H(k,k-1)),ABS(H(k-1,k)))
               ba = MIN(ABS(H(k,k-1)),ABS(H(k-1,k)))
               aa = MAX(ABS(H(k,k)),ABS(H(k-1,k-1)-H(k,k)))
               bb = MIN(ABS(H(k,k)),ABS(H(k-1,k-1)-H(k,k)))
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
!        Exit from loop if a submatrix of order 1 or 2 has split off.
!
         IF ( l>=i-1 ) GOTO 200
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
            s = ABS(H(i,i-1)) + ABS(H(i-1,i-2))
            h11 = DAT1*s + H(i,i)
            h12 = DAT2*s
            h21 = s
            h22 = h11
         ELSEIF ( MOD(kdefl,KEXSH)==0 ) THEN
!
!           Exceptional shift.
!
            s = ABS(H(l+1,l)) + ABS(H(l+2,l+1))
            h11 = DAT1*s + H(l,l)
            h12 = DAT2*s
            h21 = s
            h22 = h11
         ELSE
!
!           Prepare to use Francis' double shift
!           (i.e. 2nd degree generalized Rayleigh quotient)
!
            h11 = H(i-1,i-1)
            h21 = H(i,i-1)
            h12 = H(i-1,i)
            h22 = H(i,i)
         ENDIF
         s = ABS(h11) + ABS(h12) + ABS(h21) + ABS(h22)
         IF ( s==ZERO ) THEN
            rt1r = ZERO
            rt1i = ZERO
            rt2r = ZERO
            rt2i = ZERO
         ELSE
            h11 = h11/s
            h21 = h21/s
            h12 = h12/s
            h22 = h22/s
            tr = (h11+h22)/TWO
            det = (h11-tr)*(h22-tr) - h12*h21
            rtdisc = SQRT(ABS(det))
            IF ( det>=ZERO ) THEN
!
!              ==== complex conjugate shifts ====
!
               rt1r = tr*s
               rt2r = rt1r
               rt1i = rtdisc*s
               rt2i = -rt1i
            ELSE
!
!              ==== real shifts (use only one of them)  ====
!
               rt1r = tr + rtdisc
               rt2r = tr - rtdisc
               IF ( ABS(rt1r-h22)<=ABS(rt2r-h22) ) THEN
                  rt1r = rt1r*s
                  rt2r = rt1r
               ELSE
                  rt2r = rt2r*s
                  rt1r = rt2r
               ENDIF
               rt1i = ZERO
               rt2i = ZERO
            ENDIF
         ENDIF
!
!        Look for two consecutive small subdiagonal elements.
!
         DO m = i - 2 , l , -1
!           Determine the effect of starting the double-shift QR
!           iteration at row M, and see if this would make H(M,M-1)
!           negligible.  (The following uses scaling to avoid
!           overflows and most underflows.)
!
            h21s = H(m+1,m)
            s = ABS(H(m,m)-rt2r) + ABS(rt2i) + ABS(h21s)
            h21s = H(m+1,m)/s
            v(1) = h21s*H(m,m+1) + (H(m,m)-rt1r)*((H(m,m)-rt2r)/s)      &
     &             - rt1i*(rt2i/s)
            v(2) = h21s*(H(m,m)+H(m+1,m+1)-rt1r-rt2r)
            v(3) = h21s*H(m+2,m+1)
            s = ABS(v(1)) + ABS(v(2)) + ABS(v(3))
            v(1) = v(1)/s
            v(2) = v(2)/s
            v(3) = v(3)/s
            IF ( m==l ) EXIT
            IF ( ABS(H(m,m-1))*(ABS(v(2))+ABS(v(3)))<=ulp*ABS(v(1))     &
     &           *(ABS(H(m-1,m-1))+ABS(H(m,m))+ABS(H(m+1,m+1))) ) EXIT
         ENDDO
!
!        Double-shift QR step
!
         DO k = m , i - 1
!
!           The first iteration of this loop determines a reflection G
!           from the vector V and applies it from left and right to H,
!           thus creating a nonzero bulge below the subdiagonal.
!
!           Each subsequent iteration determines a reflection G to
!           restore the Hessenberg form in the (K-1)th column, and thus
!           chases the bulge one step toward the bottom of the active
!           submatrix. NR is the order of G.
!
            nr = MIN(3,i-k+1)
            IF ( k>m ) CALL DCOPY(nr,H(k,k-1),1,v,1)
            CALL DLARFG(nr,v(1),v(2),1,t1)
            IF ( k>m ) THEN
               H(k,k-1) = v(1)
               H(k+1,k-1) = ZERO
               IF ( k<i-1 ) H(k+2,k-1) = ZERO
            ELSEIF ( m>l ) THEN
!               ==== Use the following instead of
!               .    H( K, K-1 ) = -H( K, K-1 ) to
!               .    avoid a bug when v(2) and v(3)
!               .    underflow. ====
               H(k,k-1) = H(k,k-1)*(ONE-t1)
            ENDIF
            v2 = v(2)
            t2 = t1*v2
            IF ( nr==3 ) THEN
               v3 = v(3)
               t3 = t1*v3
!
!              Apply G from the left to transform the rows of the matrix
!              in columns K to I2.
!
               DO j = k , i2
                  sum = H(k,j) + v2*H(k+1,j) + v3*H(k+2,j)
                  H(k,j) = H(k,j) - sum*t1
                  H(k+1,j) = H(k+1,j) - sum*t2
                  H(k+2,j) = H(k+2,j) - sum*t3
               ENDDO
!
!              Apply G from the right to transform the columns of the
!              matrix in rows I1 to min(K+3,I).
!
               DO j = i1 , MIN(k+3,i)
                  sum = H(j,k) + v2*H(j,k+1) + v3*H(j,k+2)
                  H(j,k) = H(j,k) - sum*t1
                  H(j,k+1) = H(j,k+1) - sum*t2
                  H(j,k+2) = H(j,k+2) - sum*t3
               ENDDO
!
               IF ( Wantz ) THEN
!
!                 Accumulate transformations in the matrix Z
!
                  DO j = Iloz , Ihiz
                     sum = Z(j,k) + v2*Z(j,k+1) + v3*Z(j,k+2)
                     Z(j,k) = Z(j,k) - sum*t1
                     Z(j,k+1) = Z(j,k+1) - sum*t2
                     Z(j,k+2) = Z(j,k+2) - sum*t3
                  ENDDO
               ENDIF
            ELSEIF ( nr==2 ) THEN
!
!              Apply G from the left to transform the rows of the matrix
!              in columns K to I2.
!
               DO j = k , i2
                  sum = H(k,j) + v2*H(k+1,j)
                  H(k,j) = H(k,j) - sum*t1
                  H(k+1,j) = H(k+1,j) - sum*t2
               ENDDO
!
!              Apply G from the right to transform the columns of the
!              matrix in rows I1 to min(K+3,I).
!
               DO j = i1 , i
                  sum = H(j,k) + v2*H(j,k+1)
                  H(j,k) = H(j,k) - sum*t1
                  H(j,k+1) = H(j,k+1) - sum*t2
               ENDDO
!
               IF ( Wantz ) THEN
!
!                 Accumulate transformations in the matrix Z
!
                  DO j = Iloz , Ihiz
                     sum = Z(j,k) + v2*Z(j,k+1)
                     Z(j,k) = Z(j,k) - sum*t1
                     Z(j,k+1) = Z(j,k+1) - sum*t2
                  ENDDO
               ENDIF
            ENDIF
         ENDDO
!
      ENDDO
!
!     Failure to converge in remaining number of iterations
!
      Info = i
      RETURN
!
!
 200  IF ( l==i ) THEN
!
!        H(I,I-1) is negligible: one eigenvalue has converged.
!
         Wr(i) = H(i,i)
         Wi(i) = ZERO
      ELSEIF ( l==i-1 ) THEN
!
!        H(I-1,I-2) is negligible: a pair of eigenvalues have converged.
!
!        Transform the 2-by-2 submatrix to standard Schur form,
!        and compute and store the eigenvalues.
!
         CALL DLANV2(H(i-1,i-1),H(i-1,i),H(i,i-1),H(i,i),Wr(i-1),Wi(i-1)&
     &               ,Wr(i),Wi(i),cs,sn)
!
         IF ( Wantt ) THEN
!
!           Apply the transformation to the rest of H.
!
            IF ( i2>i ) CALL DROT(i2-i,H(i-1,i+1),Ldh,H(i,i+1),Ldh,cs,  &
     &                            sn)
            CALL DROT(i-i1-1,H(i1,i-1),1,H(i1,i),1,cs,sn)
         ENDIF
!
!           Apply the transformation to Z.
!
         IF ( Wantz ) CALL DROT(nz,Z(Iloz,i-1),1,Z(Iloz,i),1,cs,sn)
      ENDIF
!     reset deflation counter
      kdefl = 0
!
!     return to start of the main loop with new value of I.
!
      i = l - 1
      GOTO 100
!
!
!     End of DLAHQR
!
99999 END SUBROUTINE DLAHQR
