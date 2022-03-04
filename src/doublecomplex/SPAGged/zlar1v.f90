!*==zlar1v.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b ZLAR1V computes the (scaled) r-th column of the inverse of the submatrix in rows b1 through bn of the tridiagonal matrix LDLT - λI.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZLAR1V + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlar1v.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlar1v.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlar1v.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZLAR1V( N, B1, BN, LAMBDA, D, L, LD, LLD,
!                  PIVMIN, GAPTOL, Z, WANTNC, NEGCNT, ZTZ, MINGMA,
!                  R, ISUPPZ, NRMINV, RESID, RQCORR, WORK )
!
!       .. Scalar Arguments ..
!       LOGICAL            WANTNC
!       INTEGER   B1, BN, N, NEGCNT, R
!       DOUBLE PRECISION   GAPTOL, LAMBDA, MINGMA, NRMINV, PIVMIN, RESID,
!      $                   RQCORR, ZTZ
!       ..
!       .. Array Arguments ..
!       INTEGER            ISUPPZ( * )
!       DOUBLE PRECISION   D( * ), L( * ), LD( * ), LLD( * ),
!      $                  WORK( * )
!       COMPLEX*16       Z( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZLAR1V computes the (scaled) r-th column of the inverse of
!> the sumbmatrix in rows B1 through BN of the tridiagonal matrix
!> L D L**T - sigma I. When sigma is close to an eigenvalue, the
!> computed vector is an accurate eigenvector. Usually, r corresponds
!> to the index where the eigenvector is largest in magnitude.
!> The following steps accomplish this computation :
!> (a) Stationary qd transform,  L D L**T - sigma I = L(+) D(+) L(+)**T,
!> (b) Progressive qd transform, L D L**T - sigma I = U(-) D(-) U(-)**T,
!> (c) Computation of the diagonal elements of the inverse of
!>     L D L**T - sigma I by combining the above transforms, and choosing
!>     r as the index where the diagonal of the inverse is (one of the)
!>     largest in magnitude.
!> (d) Computation of the (scaled) r-th column of the inverse using the
!>     twisted factorization obtained by combining the top part of the
!>     the stationary and the bottom part of the progressive transform.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>           The order of the matrix L D L**T.
!> \endverbatim
!>
!> \param[in] B1
!> \verbatim
!>          B1 is INTEGER
!>           First index of the submatrix of L D L**T.
!> \endverbatim
!>
!> \param[in] BN
!> \verbatim
!>          BN is INTEGER
!>           Last index of the submatrix of L D L**T.
!> \endverbatim
!>
!> \param[in] LAMBDA
!> \verbatim
!>          LAMBDA is DOUBLE PRECISION
!>           The shift. In order to compute an accurate eigenvector,
!>           LAMBDA should be a good approximation to an eigenvalue
!>           of L D L**T.
!> \endverbatim
!>
!> \param[in] L
!> \verbatim
!>          L is DOUBLE PRECISION array, dimension (N-1)
!>           The (n-1) subdiagonal elements of the unit bidiagonal matrix
!>           L, in elements 1 to N-1.
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is DOUBLE PRECISION array, dimension (N)
!>           The n diagonal elements of the diagonal matrix D.
!> \endverbatim
!>
!> \param[in] LD
!> \verbatim
!>          LD is DOUBLE PRECISION array, dimension (N-1)
!>           The n-1 elements L(i)*D(i).
!> \endverbatim
!>
!> \param[in] LLD
!> \verbatim
!>          LLD is DOUBLE PRECISION array, dimension (N-1)
!>           The n-1 elements L(i)*L(i)*D(i).
!> \endverbatim
!>
!> \param[in] PIVMIN
!> \verbatim
!>          PIVMIN is DOUBLE PRECISION
!>           The minimum pivot in the Sturm sequence.
!> \endverbatim
!>
!> \param[in] GAPTOL
!> \verbatim
!>          GAPTOL is DOUBLE PRECISION
!>           Tolerance that indicates when eigenvector entries are negligible
!>           w.r.t. their contribution to the residual.
!> \endverbatim
!>
!> \param[in,out] Z
!> \verbatim
!>          Z is COMPLEX*16 array, dimension (N)
!>           On input, all entries of Z must be set to 0.
!>           On output, Z contains the (scaled) r-th column of the
!>           inverse. The scaling is such that Z(R) equals 1.
!> \endverbatim
!>
!> \param[in] WANTNC
!> \verbatim
!>          WANTNC is LOGICAL
!>           Specifies whether NEGCNT has to be computed.
!> \endverbatim
!>
!> \param[out] NEGCNT
!> \verbatim
!>          NEGCNT is INTEGER
!>           If WANTNC is .TRUE. then NEGCNT = the number of pivots < pivmin
!>           in the  matrix factorization L D L**T, and NEGCNT = -1 otherwise.
!> \endverbatim
!>
!> \param[out] ZTZ
!> \verbatim
!>          ZTZ is DOUBLE PRECISION
!>           The square of the 2-norm of Z.
!> \endverbatim
!>
!> \param[out] MINGMA
!> \verbatim
!>          MINGMA is DOUBLE PRECISION
!>           The reciprocal of the largest (in magnitude) diagonal
!>           element of the inverse of L D L**T - sigma I.
!> \endverbatim
!>
!> \param[in,out] R
!> \verbatim
!>          R is INTEGER
!>           The twist index for the twisted factorization used to
!>           compute Z.
!>           On input, 0 <= R <= N. If R is input as 0, R is set to
!>           the index where (L D L**T - sigma I)^{-1} is largest
!>           in magnitude. If 1 <= R <= N, R is unchanged.
!>           On output, R contains the twist index used to compute Z.
!>           Ideally, R designates the position of the maximum entry in the
!>           eigenvector.
!> \endverbatim
!>
!> \param[out] ISUPPZ
!> \verbatim
!>          ISUPPZ is INTEGER array, dimension (2)
!>           The support of the vector in Z, i.e., the vector Z is
!>           nonzero only in elements ISUPPZ(1) through ISUPPZ( 2 ).
!> \endverbatim
!>
!> \param[out] NRMINV
!> \verbatim
!>          NRMINV is DOUBLE PRECISION
!>           NRMINV = 1/SQRT( ZTZ )
!> \endverbatim
!>
!> \param[out] RESID
!> \verbatim
!>          RESID is DOUBLE PRECISION
!>           The residual of the FP vector.
!>           RESID = ABS( MINGMA )/SQRT( ZTZ )
!> \endverbatim
!>
!> \param[out] RQCORR
!> \verbatim
!>          RQCORR is DOUBLE PRECISION
!>           The Rayleigh Quotient correction to LAMBDA.
!>           RQCORR = MINGMA*TMP
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (4*N)
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
!> \ingroup complex16OTHERauxiliary
!
!> \par Contributors:
!  ==================
!>
!> Beresford Parlett, University of California, Berkeley, USA \n
!> Jim Demmel, University of California, Berkeley, USA \n
!> Inderjit Dhillon, University of Texas, Austin, USA \n
!> Osni Marques, LBNL/NERSC, USA \n
!> Christof Voemel, University of California, Berkeley, USA
!
!  =====================================================================
      SUBROUTINE ZLAR1V(N,B1,Bn,Lambda,D,L,Ld,Lld,Pivmin,Gaptol,Z,      &
     &                  Wantnc,Negcnt,Ztz,Mingma,R,Isuppz,Nrminv,Resid, &
     &                  Rqcorr,Work)
      USE F77KINDS                        
      USE S_DISNAN
      USE S_DLAMCH
      IMPLICIT NONE
!*--ZLAR1V237
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0
      COMPLEX(CX16KIND) , PARAMETER  ::  CONE = (1.0D0,0.0D0)
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: B1
      INTEGER , INTENT(IN) :: Bn
      REAL(R8KIND) , INTENT(IN) :: Lambda
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: D
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: L
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: Ld
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: Lld
      REAL(R8KIND) , INTENT(IN) :: Pivmin
      REAL(R8KIND) , INTENT(IN) :: Gaptol
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Z
      LOGICAL , INTENT(IN) :: Wantnc
      INTEGER , INTENT(OUT) :: Negcnt
      REAL(R8KIND) , INTENT(INOUT) :: Ztz
      REAL(R8KIND) , INTENT(INOUT) :: Mingma
      INTEGER , INTENT(INOUT) :: R
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Isuppz
      REAL(R8KIND) , INTENT(INOUT) :: Nrminv
      REAL(R8KIND) , INTENT(OUT) :: Resid
      REAL(R8KIND) , INTENT(OUT) :: Rqcorr
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: dminus , dplus , eps , s , tmp
      INTEGER :: i , indlpl , indp , inds , indumn , neg1 , neg2 , r1 , &
     &           r2
      LOGICAL :: sawnan1 , sawnan2
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
 
!     ..
!     .. Local Scalars ..
!     ..
!     .. External Functions ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
      eps = DLAMCH('Precision')
 
 
      IF ( R==0 ) THEN
         r1 = B1
         r2 = Bn
      ELSE
         r1 = R
         r2 = R
      ENDIF
 
!     Storage for LPLUS
      indlpl = 0
!     Storage for UMINUS
      indumn = N
      inds = 2*N + 1
      indp = 3*N + 1
 
      IF ( B1==1 ) THEN
         Work(inds) = ZERO
      ELSE
         Work(inds+B1-1) = Lld(B1-1)
      ENDIF
 
!
!     Compute the stationary transform (using the differential form)
!     until the index R2.
!
      sawnan1 = .FALSE.
      neg1 = 0
      s = Work(inds+B1-1) - Lambda
      DO i = B1 , r1 - 1
         dplus = D(i) + s
         Work(indlpl+i) = Ld(i)/dplus
         IF ( dplus<ZERO ) neg1 = neg1 + 1
         Work(inds+i) = s*Work(indlpl+i)*L(i)
         s = Work(inds+i) - Lambda
      ENDDO
      sawnan1 = DISNAN(s)
      IF ( .NOT.(sawnan1) ) THEN
         DO i = r1 , r2 - 1
            dplus = D(i) + s
            Work(indlpl+i) = Ld(i)/dplus
            Work(inds+i) = s*Work(indlpl+i)*L(i)
            s = Work(inds+i) - Lambda
         ENDDO
         sawnan1 = DISNAN(s)
      ENDIF
!
      IF ( sawnan1 ) THEN
!        Runs a slower version of the above loop if a NaN is detected
         neg1 = 0
         s = Work(inds+B1-1) - Lambda
         DO i = B1 , r1 - 1
            dplus = D(i) + s
            IF ( ABS(dplus)<Pivmin ) dplus = -Pivmin
            Work(indlpl+i) = Ld(i)/dplus
            IF ( dplus<ZERO ) neg1 = neg1 + 1
            Work(inds+i) = s*Work(indlpl+i)*L(i)
            IF ( Work(indlpl+i)==ZERO ) Work(inds+i) = Lld(i)
            s = Work(inds+i) - Lambda
         ENDDO
         DO i = r1 , r2 - 1
            dplus = D(i) + s
            IF ( ABS(dplus)<Pivmin ) dplus = -Pivmin
            Work(indlpl+i) = Ld(i)/dplus
            Work(inds+i) = s*Work(indlpl+i)*L(i)
            IF ( Work(indlpl+i)==ZERO ) Work(inds+i) = Lld(i)
            s = Work(inds+i) - Lambda
         ENDDO
      ENDIF
!
!     Compute the progressive transform (using the differential form)
!     until the index R1
!
      sawnan2 = .FALSE.
      neg2 = 0
      Work(indp+Bn-1) = D(Bn) - Lambda
      DO i = Bn - 1 , r1 , -1
         dminus = Lld(i) + Work(indp+i)
         tmp = D(i)/dminus
         IF ( dminus<ZERO ) neg2 = neg2 + 1
         Work(indumn+i) = L(i)*tmp
         Work(indp+i-1) = Work(indp+i)*tmp - Lambda
      ENDDO
      tmp = Work(indp+r1-1)
      sawnan2 = DISNAN(tmp)
 
      IF ( sawnan2 ) THEN
!        Runs a slower version of the above loop if a NaN is detected
         neg2 = 0
         DO i = Bn - 1 , r1 , -1
            dminus = Lld(i) + Work(indp+i)
            IF ( ABS(dminus)<Pivmin ) dminus = -Pivmin
            tmp = D(i)/dminus
            IF ( dminus<ZERO ) neg2 = neg2 + 1
            Work(indumn+i) = L(i)*tmp
            Work(indp+i-1) = Work(indp+i)*tmp - Lambda
            IF ( tmp==ZERO ) Work(indp+i-1) = D(i) - Lambda
         ENDDO
      ENDIF
!
!     Find the index (from R1 to R2) of the largest (in magnitude)
!     diagonal element of the inverse
!
      Mingma = Work(inds+r1-1) + Work(indp+r1-1)
      IF ( Mingma<ZERO ) neg1 = neg1 + 1
      IF ( Wantnc ) THEN
         Negcnt = neg1 + neg2
      ELSE
         Negcnt = -1
      ENDIF
      IF ( ABS(Mingma)==ZERO ) Mingma = eps*Work(inds+r1-1)
      R = r1
      DO i = r1 , r2 - 1
         tmp = Work(inds+i) + Work(indp+i)
         IF ( tmp==ZERO ) tmp = eps*Work(inds+i)
         IF ( ABS(tmp)<=ABS(Mingma) ) THEN
            Mingma = tmp
            R = i + 1
         ENDIF
      ENDDO
!
!     Compute the FP vector: solve N^T v = e_r
!
      Isuppz(1) = B1
      Isuppz(2) = Bn
      Z(R) = CONE
      Ztz = ONE
!
!     Compute the FP vector upwards from R
!
      IF ( .NOT.sawnan1 .AND. .NOT.sawnan2 ) THEN
         DO i = R - 1 , B1 , -1
            Z(i) = -(Work(indlpl+i)*Z(i+1))
            IF ( (ABS(Z(i))+ABS(Z(i+1)))*ABS(Ld(i))<Gaptol ) THEN
               Z(i) = ZERO
               Isuppz(1) = i + 1
               EXIT
            ENDIF
            Ztz = Ztz + DBLE(Z(i)*Z(i))
         ENDDO
      ELSE
!        Run slower loop if NaN occurred.
         DO i = R - 1 , B1 , -1
            IF ( Z(i+1)==ZERO ) THEN
               Z(i) = -(Ld(i+1)/Ld(i))*Z(i+2)
            ELSE
               Z(i) = -(Work(indlpl+i)*Z(i+1))
            ENDIF
            IF ( (ABS(Z(i))+ABS(Z(i+1)))*ABS(Ld(i))<Gaptol ) THEN
               Z(i) = ZERO
               Isuppz(1) = i + 1
               EXIT
            ENDIF
            Ztz = Ztz + DBLE(Z(i)*Z(i))
         ENDDO
      ENDIF
 
!     Compute the FP vector downwards from R in blocks of size BLKSIZ
      IF ( .NOT.sawnan1 .AND. .NOT.sawnan2 ) THEN
         DO i = R , Bn - 1
            Z(i+1) = -(Work(indumn+i)*Z(i))
            IF ( (ABS(Z(i))+ABS(Z(i+1)))*ABS(Ld(i))<Gaptol ) THEN
               Z(i+1) = ZERO
               Isuppz(2) = i
               EXIT
            ENDIF
            Ztz = Ztz + DBLE(Z(i+1)*Z(i+1))
         ENDDO
      ELSE
!        Run slower loop if NaN occurred.
         DO i = R , Bn - 1
            IF ( Z(i)==ZERO ) THEN
               Z(i+1) = -(Ld(i-1)/Ld(i))*Z(i-1)
            ELSE
               Z(i+1) = -(Work(indumn+i)*Z(i))
            ENDIF
            IF ( (ABS(Z(i))+ABS(Z(i+1)))*ABS(Ld(i))<Gaptol ) THEN
               Z(i+1) = ZERO
               Isuppz(2) = i
               EXIT
            ENDIF
            Ztz = Ztz + DBLE(Z(i+1)*Z(i+1))
         ENDDO
      ENDIF
!
!     Compute quantities for convergence test
!
      tmp = ONE/Ztz
      Nrminv = SQRT(tmp)
      Resid = ABS(Mingma)*Nrminv
      Rqcorr = Mingma*tmp
!
!
!
!     End of ZLAR1V
!
      END SUBROUTINE ZLAR1V
