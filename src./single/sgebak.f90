!*==sgebak.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b SGEBAK
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SGEBAK + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgebak.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgebak.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgebak.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SGEBAK( JOB, SIDE, N, ILO, IHI, SCALE, M, V, LDV,
!                          INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          JOB, SIDE
!       INTEGER            IHI, ILO, INFO, LDV, M, N
!       ..
!       .. Array Arguments ..
!       REAL               V( LDV, * ), SCALE( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SGEBAK forms the right or left eigenvectors of a real general matrix
!> by backward transformation on the computed eigenvectors of the
!> balanced matrix output by SGEBAL.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] JOB
!> \verbatim
!>          JOB is CHARACTER*1
!>          Specifies the type of backward transformation required:
!>          = 'N': do nothing, return immediately;
!>          = 'P': do backward transformation for permutation only;
!>          = 'S': do backward transformation for scaling only;
!>          = 'B': do backward transformations for both permutation and
!>                 scaling.
!>          JOB must be the same as the argument JOB supplied to SGEBAL.
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
!>          The integers ILO and IHI determined by SGEBAL.
!>          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.
!> \endverbatim
!>
!> \param[in] SCALE
!> \verbatim
!>          SCALE is REAL array, dimension (N)
!>          Details of the permutation and scaling factors, as returned
!>          by SGEBAL.
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
!>          transformed, as returned by SHSEIN or STREVC.
!>          On exit, V is overwritten by the transformed eigenvectors.
!> \endverbatim
!>
!> \param[in] LDV
!> \verbatim
!>          LDV is INTEGER
!>          The leading dimension of the array V. LDV >= max(1,N).
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
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
!> \ingroup realGEcomputational
!
!  =====================================================================
      SUBROUTINE SGEBAK(Job,Side,N,Ilo,Ihi,Scale,M,V,Ldv,Info)
      USE S_LSAME
      USE S_SSCAL
      USE S_SSWAP
      USE S_XERBLA
      IMPLICIT NONE
!*--SGEBAK137
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ONE = 1.0E+0
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Job
      CHARACTER :: Side
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Ilo
      INTEGER , INTENT(IN) :: Ihi
      REAL , INTENT(IN) , DIMENSION(*) :: Scale
      INTEGER :: M
      REAL , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: i , ii , k
      LOGICAL :: leftv , rightv
      REAL :: s
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
!     .. External Subroutines ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
!     Decode and Test the input parameters
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
      ELSEIF ( Ilo<1 .OR. Ilo>MAX(1,N) ) THEN
         Info = -4
      ELSEIF ( Ihi<MIN(Ilo,N) .OR. Ihi>N ) THEN
         Info = -5
      ELSEIF ( M<0 ) THEN
         Info = -7
      ELSEIF ( Ldv<MAX(1,N) ) THEN
         Info = -9
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('SGEBAK',-Info)
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
            IF ( rightv ) THEN
               DO i = Ilo , Ihi
                  s = Scale(i)
                  CALL SSCAL(M,s,V(i,1),Ldv)
               ENDDO
            ENDIF
!
            IF ( leftv ) THEN
               DO i = Ilo , Ihi
                  s = ONE/Scale(i)
                  CALL SSCAL(M,s,V(i,1),Ldv)
               ENDDO
            ENDIF
!
         ENDIF
      ENDIF
!
!     Backward permutation
!
!     For  I = ILO-1 step -1 until 1,
!              IHI+1 step 1 until N do --
!
      IF ( LSAME(Job,'P') .OR. LSAME(Job,'B') ) THEN
         IF ( rightv ) THEN
            DO ii = 1 , N
               i = ii
               IF ( i<Ilo .OR. i>Ihi ) THEN
                  IF ( i<Ilo ) i = Ilo - ii
                  k = Scale(i)
                  IF ( k/=i ) CALL SSWAP(M,V(i,1),Ldv,V(k,1),Ldv)
               ENDIF
            ENDDO
         ENDIF
!
         IF ( leftv ) THEN
            DO ii = 1 , N
               i = ii
               IF ( i<Ilo .OR. i>Ihi ) THEN
                  IF ( i<Ilo ) i = Ilo - ii
                  k = Scale(i)
                  IF ( k/=i ) CALL SSWAP(M,V(i,1),Ldv,V(k,1),Ldv)
               ENDIF
            ENDDO
         ENDIF
      ENDIF
!
!
!     End of SGEBAK
!
      END SUBROUTINE SGEBAK
