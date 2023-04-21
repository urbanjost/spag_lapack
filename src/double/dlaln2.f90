!*==dlaln2.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b DLALN2 solves a 1-by-1 or 2-by-2 linear system of equations of the specified form.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLALN2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaln2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaln2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaln2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLALN2( LTRANS, NA, NW, SMIN, CA, A, LDA, D1, D2, B,
!                          LDB, WR, WI, X, LDX, SCALE, XNORM, INFO )
!
!       .. Scalar Arguments ..
!       LOGICAL            LTRANS
!       INTEGER            INFO, LDA, LDB, LDX, NA, NW
!       DOUBLE PRECISION   CA, D1, D2, SCALE, SMIN, WI, WR, XNORM
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), X( LDX, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLALN2 solves a system of the form  (ca A - w D ) X = s B
!> or (ca A**T - w D) X = s B   with possible scaling ("s") and
!> perturbation of A.  (A**T means A-transpose.)
!>
!> A is an NA x NA real matrix, ca is a real scalar, D is an NA x NA
!> real diagonal matrix, w is a real or complex value, and X and B are
!> NA x 1 matrices -- real if w is real, complex if w is complex.  NA
!> may be 1 or 2.
!>
!> If w is complex, X and B are represented as NA x 2 matrices,
!> the first column of each being the real part and the second
!> being the imaginary part.
!>
!> "s" is a scaling factor (<= 1), computed by DLALN2, which is
!> so chosen that X can be computed without overflow.  X is further
!> scaled if necessary to assure that norm(ca A - w D)*norm(X) is less
!> than overflow.
!>
!> If both singular values of (ca A - w D) are less than SMIN,
!> SMIN*identity will be used instead of (ca A - w D).  If only one
!> singular value is less than SMIN, one element of (ca A - w D) will be
!> perturbed enough to make the smallest singular value roughly SMIN.
!> If both singular values are at least SMIN, (ca A - w D) will not be
!> perturbed.  In any case, the perturbation will be at most some small
!> multiple of max( SMIN, ulp*norm(ca A - w D) ).  The singular values
!> are computed by infinity-norm approximations, and thus will only be
!> correct to a factor of 2 or so.
!>
!> Note: all input quantities are assumed to be smaller than overflow
!> by a reasonable factor.  (See BIGNUM.)
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] LTRANS
!> \verbatim
!>          LTRANS is LOGICAL
!>          =.TRUE.:  A-transpose will be used.
!>          =.FALSE.: A will be used (not transposed.)
!> \endverbatim
!>
!> \param[in] NA
!> \verbatim
!>          NA is INTEGER
!>          The size of the matrix A.  It may (only) be 1 or 2.
!> \endverbatim
!>
!> \param[in] NW
!> \verbatim
!>          NW is INTEGER
!>          1 if "w" is real, 2 if "w" is complex.  It may only be 1
!>          or 2.
!> \endverbatim
!>
!> \param[in] SMIN
!> \verbatim
!>          SMIN is DOUBLE PRECISION
!>          The desired lower bound on the singular values of A.  This
!>          should be a safe distance away from underflow or overflow,
!>          say, between (underflow/machine precision) and  (machine
!>          precision * overflow ).  (See BIGNUM and ULP.)
!> \endverbatim
!>
!> \param[in] CA
!> \verbatim
!>          CA is DOUBLE PRECISION
!>          The coefficient c, which A is multiplied by.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,NA)
!>          The NA x NA matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of A.  It must be at least NA.
!> \endverbatim
!>
!> \param[in] D1
!> \verbatim
!>          D1 is DOUBLE PRECISION
!>          The 1,1 element in the diagonal matrix D.
!> \endverbatim
!>
!> \param[in] D2
!> \verbatim
!>          D2 is DOUBLE PRECISION
!>          The 2,2 element in the diagonal matrix D.  Not used if NA=1.
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is DOUBLE PRECISION array, dimension (LDB,NW)
!>          The NA x NW matrix B (right-hand side).  If NW=2 ("w" is
!>          complex), column 1 contains the real part of B and column 2
!>          contains the imaginary part.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of B.  It must be at least NA.
!> \endverbatim
!>
!> \param[in] WR
!> \verbatim
!>          WR is DOUBLE PRECISION
!>          The real part of the scalar "w".
!> \endverbatim
!>
!> \param[in] WI
!> \verbatim
!>          WI is DOUBLE PRECISION
!>          The imaginary part of the scalar "w".  Not used if NW=1.
!> \endverbatim
!>
!> \param[out] X
!> \verbatim
!>          X is DOUBLE PRECISION array, dimension (LDX,NW)
!>          The NA x NW matrix X (unknowns), as computed by DLALN2.
!>          If NW=2 ("w" is complex), on exit, column 1 will contain
!>          the real part of X and column 2 will contain the imaginary
!>          part.
!> \endverbatim
!>
!> \param[in] LDX
!> \verbatim
!>          LDX is INTEGER
!>          The leading dimension of X.  It must be at least NA.
!> \endverbatim
!>
!> \param[out] SCALE
!> \verbatim
!>          SCALE is DOUBLE PRECISION
!>          The scale factor that B must be multiplied by to insure
!>          that overflow does not occur when computing X.  Thus,
!>          (ca A - w D) X  will be SCALE*B, not B (ignoring
!>          perturbations of A.)  It will be at most 1.
!> \endverbatim
!>
!> \param[out] XNORM
!> \verbatim
!>          XNORM is DOUBLE PRECISION
!>          The infinity-norm of X, when X is regarded as an NA x NW
!>          real matrix.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          An error flag.  It will be set to zero if no error occurs,
!>          a negative number if an argument is in error, or a positive
!>          number if  ca A - w D  had to be perturbed.
!>          The possible values are:
!>          = 0: No error occurred, and (ca A - w D) did not have to be
!>                 perturbed.
!>          = 1: (ca A - w D) had to be perturbed to make its smallest
!>               (or only) singular value greater than SMIN.
!>          NOTE: In the interests of speed, this routine does not
!>                check the inputs for errors.
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
!  =====================================================================
      SUBROUTINE DLALN2(Ltrans,Na,Nw,Smin,Ca,A,Lda,D1,D2,B,Ldb,Wr,Wi,X, &
     &                  Ldx,Scale,Xnorm,Info)
      IMPLICIT NONE
!*--DLALN2222
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      LOGICAL Ltrans
      INTEGER Info , Lda , Ldb , Ldx , Na , Nw
      DOUBLE PRECISION Ca , D1 , D2 , Scale , Smin , Wi , Wr , Xnorm
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION A(Lda,*) , B(Ldb,*) , X(Ldx,*)
!     ..
!
! =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D0,ONE=1.0D0)
      DOUBLE PRECISION TWO
      PARAMETER (TWO=2.0D0)
!     ..
!     .. Local Scalars ..
      INTEGER icmax , j
      DOUBLE PRECISION bbnd , bi1 , bi2 , bignum , bnorm , br1 , br2 ,  &
     &                 ci21 , ci22 , cmax , cnorm , cr21 , cr22 , csi , &
     &                 csr , li21 , lr21 , smini , smlnum , temp ,      &
     &                 u22abs , ui11 , ui11r , ui12 , ui12s , ui22 ,    &
     &                 ur11 , ur11r , ur12 , ur12s , ur22 , xi1 , xi2 , &
     &                 xr1 , xr2
!     ..
!     .. Local Arrays ..
      LOGICAL rswap(4) , zswap(4)
      INTEGER ipivot(4,4)
      DOUBLE PRECISION ci(2,2) , civ(4) , cr(2,2) , crv(4)
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMCH
      EXTERNAL DLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL DLADIV
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX
!     ..
!     .. Equivalences ..
      EQUIVALENCE (ci(1,1),civ(1))
      EQUIVALENCE (cr(1,1),crv(1))
!     ..
!     .. Data statements ..
      DATA zswap/.FALSE. , .FALSE. , .TRUE. , .TRUE./
      DATA rswap/.FALSE. , .TRUE. , .FALSE. , .TRUE./
      DATA ipivot/1 , 2 , 3 , 4 , 2 , 1 , 4 , 3 , 3 , 4 , 1 , 2 , 4 ,   &
     &     3 , 2 , 1/
!     ..
!     .. Executable Statements ..
!
!     Compute BIGNUM
!
      smlnum = TWO*DLAMCH('Safe minimum')
      bignum = ONE/smlnum
      smini = MAX(Smin,smlnum)
!
!     Don't check for input errors
!
      Info = 0
!
!     Standard Initializations
!
      Scale = ONE
!
      IF ( Na/=1 ) THEN
!
!        2x2 System
!
!        Compute the real part of  C = ca A - w D  (or  ca A**T - w D )
!
         cr(1,1) = Ca*A(1,1) - Wr*D1
         cr(2,2) = Ca*A(2,2) - Wr*D2
         IF ( Ltrans ) THEN
            cr(1,2) = Ca*A(2,1)
            cr(2,1) = Ca*A(1,2)
         ELSE
            cr(2,1) = Ca*A(2,1)
            cr(1,2) = Ca*A(1,2)
         ENDIF
!
         IF ( Nw==1 ) THEN
!
!           Real 2x2 system  (w is real)
!
!           Find the largest element in C
!
            cmax = ZERO
            icmax = 0
!
            DO j = 1 , 4
               IF ( ABS(crv(j))>cmax ) THEN
                  cmax = ABS(crv(j))
                  icmax = j
               ENDIF
            ENDDO
!
!           If norm(C) < SMINI, use SMINI*identity.
!
            IF ( cmax<smini ) THEN
               bnorm = MAX(ABS(B(1,1)),ABS(B(2,1)))
               IF ( smini<ONE .AND. bnorm>ONE ) THEN
                  IF ( bnorm>bignum*smini ) Scale = ONE/bnorm
               ENDIF
               temp = Scale/smini
               X(1,1) = temp*B(1,1)
               X(2,1) = temp*B(2,1)
               Xnorm = temp*bnorm
               Info = 1
               RETURN
            ENDIF
!
!           Gaussian elimination with complete pivoting.
!
            ur11 = crv(icmax)
            cr21 = crv(ipivot(2,icmax))
            ur12 = crv(ipivot(3,icmax))
            cr22 = crv(ipivot(4,icmax))
            ur11r = ONE/ur11
            lr21 = ur11r*cr21
            ur22 = cr22 - ur12*lr21
!
!           If smaller pivot < SMINI, use SMINI
!
            IF ( ABS(ur22)<smini ) THEN
               ur22 = smini
               Info = 1
            ENDIF
            IF ( rswap(icmax) ) THEN
               br1 = B(2,1)
               br2 = B(1,1)
            ELSE
               br1 = B(1,1)
               br2 = B(2,1)
            ENDIF
            br2 = br2 - lr21*br1
            bbnd = MAX(ABS(br1*(ur22*ur11r)),ABS(br2))
            IF ( bbnd>ONE .AND. ABS(ur22)<ONE ) THEN
               IF ( bbnd>=bignum*ABS(ur22) ) Scale = ONE/bbnd
            ENDIF
!
            xr2 = (br2*Scale)/ur22
            xr1 = (Scale*br1)*ur11r - xr2*(ur11r*ur12)
            IF ( zswap(icmax) ) THEN
               X(1,1) = xr2
               X(2,1) = xr1
            ELSE
               X(1,1) = xr1
               X(2,1) = xr2
            ENDIF
            Xnorm = MAX(ABS(xr1),ABS(xr2))
!
!           Further scaling if  norm(A) norm(X) > overflow
!
            IF ( Xnorm>ONE .AND. cmax>ONE ) THEN
               IF ( Xnorm>bignum/cmax ) THEN
                  temp = cmax/bignum
                  X(1,1) = temp*X(1,1)
                  X(2,1) = temp*X(2,1)
                  Xnorm = temp*Xnorm
                  Scale = temp*Scale
               ENDIF
            ENDIF
         ELSE
!
!           Complex 2x2 system  (w is complex)
!
!           Find the largest element in C
!
            ci(1,1) = -Wi*D1
            ci(2,1) = ZERO
            ci(1,2) = ZERO
            ci(2,2) = -Wi*D2
            cmax = ZERO
            icmax = 0
!
            DO j = 1 , 4
               IF ( ABS(crv(j))+ABS(civ(j))>cmax ) THEN
                  cmax = ABS(crv(j)) + ABS(civ(j))
                  icmax = j
               ENDIF
            ENDDO
!
!           If norm(C) < SMINI, use SMINI*identity.
!
            IF ( cmax<smini ) THEN
               bnorm = MAX(ABS(B(1,1))+ABS(B(1,2)),ABS(B(2,1))          &
     &                 +ABS(B(2,2)))
               IF ( smini<ONE .AND. bnorm>ONE ) THEN
                  IF ( bnorm>bignum*smini ) Scale = ONE/bnorm
               ENDIF
               temp = Scale/smini
               X(1,1) = temp*B(1,1)
               X(2,1) = temp*B(2,1)
               X(1,2) = temp*B(1,2)
               X(2,2) = temp*B(2,2)
               Xnorm = temp*bnorm
               Info = 1
               RETURN
            ENDIF
!
!           Gaussian elimination with complete pivoting.
!
            ur11 = crv(icmax)
            ui11 = civ(icmax)
            cr21 = crv(ipivot(2,icmax))
            ci21 = civ(ipivot(2,icmax))
            ur12 = crv(ipivot(3,icmax))
            ui12 = civ(ipivot(3,icmax))
            cr22 = crv(ipivot(4,icmax))
            ci22 = civ(ipivot(4,icmax))
            IF ( icmax==1 .OR. icmax==4 ) THEN
!
!              Code when off-diagonals of pivoted C are real
!
               IF ( ABS(ur11)>ABS(ui11) ) THEN
                  temp = ui11/ur11
                  ur11r = ONE/(ur11*(ONE+temp**2))
                  ui11r = -temp*ur11r
               ELSE
                  temp = ur11/ui11
                  ui11r = -ONE/(ui11*(ONE+temp**2))
                  ur11r = -temp*ui11r
               ENDIF
               lr21 = cr21*ur11r
               li21 = cr21*ui11r
               ur12s = ur12*ur11r
               ui12s = ur12*ui11r
               ur22 = cr22 - ur12*lr21
               ui22 = ci22 - ur12*li21
            ELSE
!
!              Code when diagonals of pivoted C are real
!
               ur11r = ONE/ur11
               ui11r = ZERO
               lr21 = cr21*ur11r
               li21 = ci21*ur11r
               ur12s = ur12*ur11r
               ui12s = ui12*ur11r
               ur22 = cr22 - ur12*lr21 + ui12*li21
               ui22 = -ur12*li21 - ui12*lr21
            ENDIF
            u22abs = ABS(ur22) + ABS(ui22)
!
!           If smaller pivot < SMINI, use SMINI
!
            IF ( u22abs<smini ) THEN
               ur22 = smini
               ui22 = ZERO
               Info = 1
            ENDIF
            IF ( rswap(icmax) ) THEN
               br2 = B(1,1)
               br1 = B(2,1)
               bi2 = B(1,2)
               bi1 = B(2,2)
            ELSE
               br1 = B(1,1)
               br2 = B(2,1)
               bi1 = B(1,2)
               bi2 = B(2,2)
            ENDIF
            br2 = br2 - lr21*br1 + li21*bi1
            bi2 = bi2 - li21*br1 - lr21*bi1
            bbnd = MAX((ABS(br1)+ABS(bi1))                              &
     &             *(u22abs*(ABS(ur11r)+ABS(ui11r))),ABS(br2)+ABS(bi2))
            IF ( bbnd>ONE .AND. u22abs<ONE ) THEN
               IF ( bbnd>=bignum*u22abs ) THEN
                  Scale = ONE/bbnd
                  br1 = Scale*br1
                  bi1 = Scale*bi1
                  br2 = Scale*br2
                  bi2 = Scale*bi2
               ENDIF
            ENDIF
!
            CALL DLADIV(br2,bi2,ur22,ui22,xr2,xi2)
            xr1 = ur11r*br1 - ui11r*bi1 - ur12s*xr2 + ui12s*xi2
            xi1 = ui11r*br1 + ur11r*bi1 - ui12s*xr2 - ur12s*xi2
            IF ( zswap(icmax) ) THEN
               X(1,1) = xr2
               X(2,1) = xr1
               X(1,2) = xi2
               X(2,2) = xi1
            ELSE
               X(1,1) = xr1
               X(2,1) = xr2
               X(1,2) = xi1
               X(2,2) = xi2
            ENDIF
            Xnorm = MAX(ABS(xr1)+ABS(xi1),ABS(xr2)+ABS(xi2))
!
!           Further scaling if  norm(A) norm(X) > overflow
!
            IF ( Xnorm>ONE .AND. cmax>ONE ) THEN
               IF ( Xnorm>bignum/cmax ) THEN
                  temp = cmax/bignum
                  X(1,1) = temp*X(1,1)
                  X(2,1) = temp*X(2,1)
                  X(1,2) = temp*X(1,2)
                  X(2,2) = temp*X(2,2)
                  Xnorm = temp*Xnorm
                  Scale = temp*Scale
               ENDIF
            ENDIF
         ENDIF
!
!        1 x 1  (i.e., scalar) system   C X = B
!
      ELSEIF ( Nw==1 ) THEN
!
!           Real 1x1 system.
!
!           C = ca A - w D
!
         csr = Ca*A(1,1) - Wr*D1
         cnorm = ABS(csr)
!
!           If | C | < SMINI, use C = SMINI
!
         IF ( cnorm<smini ) THEN
            csr = smini
            cnorm = smini
            Info = 1
         ENDIF
!
!           Check scaling for  X = B / C
!
         bnorm = ABS(B(1,1))
         IF ( cnorm<ONE .AND. bnorm>ONE ) THEN
            IF ( bnorm>bignum*cnorm ) Scale = ONE/bnorm
         ENDIF
!
!           Compute X
!
         X(1,1) = (B(1,1)*Scale)/csr
         Xnorm = ABS(X(1,1))
      ELSE
!
!           Complex 1x1 system (w is complex)
!
!           C = ca A - w D
!
         csr = Ca*A(1,1) - Wr*D1
         csi = -Wi*D1
         cnorm = ABS(csr) + ABS(csi)
!
!           If | C | < SMINI, use C = SMINI
!
         IF ( cnorm<smini ) THEN
            csr = smini
            csi = ZERO
            cnorm = smini
            Info = 1
         ENDIF
!
!           Check scaling for  X = B / C
!
         bnorm = ABS(B(1,1)) + ABS(B(1,2))
         IF ( cnorm<ONE .AND. bnorm>ONE ) THEN
            IF ( bnorm>bignum*cnorm ) Scale = ONE/bnorm
         ENDIF
!
!           Compute X
!
         CALL DLADIV(Scale*B(1,1),Scale*B(1,2),csr,csi,X(1,1),X(1,2))
         Xnorm = ABS(X(1,1)) + ABS(X(1,2))
!
      ENDIF
!
!
!     End of DLALN2
!
      END SUBROUTINE DLALN2
