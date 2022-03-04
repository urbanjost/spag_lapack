!*==sggbak.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b SGGBAK
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SGGBAK + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sggbak.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sggbak.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sggbak.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SGGBAK( JOB, SIDE, N, ILO, IHI, LSCALE, RSCALE, M, V,
!                          LDV, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          JOB, SIDE
!       INTEGER            IHI, ILO, INFO, LDV, M, N
!       ..
!       .. Array Arguments ..
!       REAL               LSCALE( * ), RSCALE( * ), V( LDV, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SGGBAK forms the right or left eigenvectors of a real generalized
!> eigenvalue problem A*x = lambda*B*x, by backward transformation on
!> the computed eigenvectors of the balanced pair of matrices output by
!> SGGBAL.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] JOB
!> \verbatim
!>          JOB is CHARACTER*1
!>          Specifies the type of backward transformation required:
!>          = 'N':  do nothing, return immediately;
!>          = 'P':  do backward transformation for permutation only;
!>          = 'S':  do backward transformation for scaling only;
!>          = 'B':  do backward transformations for both permutation and
!>                  scaling.
!>          JOB must be the same as the argument JOB supplied to SGGBAL.
!> \endverbatim
!>
!> \param[in] SIDE
!> \verbatim
!>          SIDE is CHARACTER*1
!>          = 'R':  V contains right eigenvectors;
!>          = 'L':  V contains left eigenvectors.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of rows of the matrix V.  N >= 0.
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
!>          The integers ILO and IHI determined by SGGBAL.
!>          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.
!> \endverbatim
!>
!> \param[in] LSCALE
!> \verbatim
!>          LSCALE is REAL array, dimension (N)
!>          Details of the permutations and/or scaling factors applied
!>          to the left side of A and B, as returned by SGGBAL.
!> \endverbatim
!>
!> \param[in] RSCALE
!> \verbatim
!>          RSCALE is REAL array, dimension (N)
!>          Details of the permutations and/or scaling factors applied
!>          to the right side of A and B, as returned by SGGBAL.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of columns of the matrix V.  M >= 0.
!> \endverbatim
!>
!> \param[in,out] V
!> \verbatim
!>          V is REAL array, dimension (LDV,M)
!>          On entry, the matrix of right or left eigenvectors to be
!>          transformed, as returned by STGEVC.
!>          On exit, V is overwritten by the transformed eigenvectors.
!> \endverbatim
!>
!> \param[in] LDV
!> \verbatim
!>          LDV is INTEGER
!>          The leading dimension of the matrix V. LDV >= max(1,N).
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit.
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
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
!> \ingroup realGBcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  See R.C. Ward, Balancing the generalized eigenvalue problem,
!>                 SIAM J. Sci. Stat. Comp. 2 (1981), 141-152.
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE SGGBAK(Job,Side,N,Ilo,Ihi,Lscale,Rscale,M,V,Ldv,Info)
      USE S_LSAME
      USE S_SSCAL
      USE S_SSWAP
      USE S_XERBLA
      IMPLICIT NONE
!*--SGGBAK154
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Job
      CHARACTER :: Side
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Ilo
      INTEGER , INTENT(IN) :: Ihi
      REAL , DIMENSION(*) :: Lscale
      REAL , DIMENSION(*) :: Rscale
      INTEGER :: M
      REAL , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: i , k
      LOGICAL :: leftv , rightv
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
!     ..
!     .. External Functions ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters
!
      rightv = LSAME(Side,'R')
      leftv = LSAME(Side,'L')
!
      Info = 0
      IF ( .NOT.LSAME(Job,'N') .AND. .NOT.LSAME(Job,'P') .AND.          &
     &     .NOT.LSAME(Job,'S') .AND. .NOT.LSAME(Job,'B') ) THEN
         Info = -1
      ELSEIF ( .NOT.rightv .AND. .NOT.leftv ) THEN
         Info = -2
      ELSEIF ( N<0 ) THEN
         Info = -3
      ELSEIF ( Ilo<1 ) THEN
         Info = -4
      ELSEIF ( N==0 .AND. Ihi==0 .AND. Ilo/=1 ) THEN
         Info = -4
      ELSEIF ( N>0 .AND. (Ihi<Ilo .OR. Ihi>MAX(1,N)) ) THEN
         Info = -5
      ELSEIF ( N==0 .AND. Ilo==1 .AND. Ihi/=0 ) THEN
         Info = -5
      ELSEIF ( M<0 ) THEN
         Info = -8
      ELSEIF ( Ldv<MAX(1,N) ) THEN
         Info = -10
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('SGGBAK',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 ) RETURN
      IF ( M==0 ) RETURN
      IF ( LSAME(Job,'N') ) RETURN
!
      IF ( Ilo/=Ihi ) THEN
!
!     Backward balance
!
         IF ( LSAME(Job,'S') .OR. LSAME(Job,'B') ) THEN
!
!        Backward transformation on right eigenvectors
!
            IF ( rightv ) THEN
               DO i = Ilo , Ihi
                  CALL SSCAL(M,Rscale(i),V(i,1),Ldv)
               ENDDO
            ENDIF
!
!        Backward transformation on left eigenvectors
!
            IF ( leftv ) THEN
               DO i = Ilo , Ihi
                  CALL SSCAL(M,Lscale(i),V(i,1),Ldv)
               ENDDO
            ENDIF
         ENDIF
      ENDIF
!
!     Backward permutation
!
      IF ( LSAME(Job,'P') .OR. LSAME(Job,'B') ) THEN
!
!        Backward permutation on right eigenvectors
!
         IF ( rightv ) THEN
            IF ( Ilo/=1 ) THEN
!
               DO i = Ilo - 1 , 1 , -1
                  k = Rscale(i)
                  IF ( k/=i ) CALL SSWAP(M,V(i,1),Ldv,V(k,1),Ldv)
               ENDDO
            ENDIF
!
            IF ( Ihi/=N ) THEN
               DO i = Ihi + 1 , N
                  k = Rscale(i)
                  IF ( k/=i ) CALL SSWAP(M,V(i,1),Ldv,V(k,1),Ldv)
               ENDDO
            ENDIF
         ENDIF
!
!        Backward permutation on left eigenvectors
!
         IF ( leftv ) THEN
            IF ( Ilo/=1 ) THEN
               DO i = Ilo - 1 , 1 , -1
                  k = Lscale(i)
                  IF ( k/=i ) CALL SSWAP(M,V(i,1),Ldv,V(k,1),Ldv)
               ENDDO
            ENDIF
!
            IF ( Ihi/=N ) THEN
               DO i = Ihi + 1 , N
                  k = Lscale(i)
                  IF ( k/=i ) CALL SSWAP(M,V(i,1),Ldv,V(k,1),Ldv)
               ENDDO
            ENDIF
         ENDIF
      ENDIF
!
!
!
!     End of SGGBAK
!
      END SUBROUTINE SGGBAK
