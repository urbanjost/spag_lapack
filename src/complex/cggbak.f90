!*==cggbak.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b CGGBAK
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CGGBAK + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cggbak.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cggbak.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cggbak.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CGGBAK( JOB, SIDE, N, ILO, IHI, LSCALE, RSCALE, M, V,
!                          LDV, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          JOB, SIDE
!       INTEGER            IHI, ILO, INFO, LDV, M, N
!       ..
!       .. Array Arguments ..
!       REAL               LSCALE( * ), RSCALE( * )
!       COMPLEX            V( LDV, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CGGBAK forms the right or left eigenvectors of a complex generalized
!> eigenvalue problem A*x = lambda*B*x, by backward transformation on
!> the computed eigenvectors of the balanced pair of matrices output by
!> CGGBAL.
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
!>          JOB must be the same as the argument JOB supplied to CGGBAL.
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
!>          The integers ILO and IHI determined by CGGBAL.
!>          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.
!> \endverbatim
!>
!> \param[in] LSCALE
!> \verbatim
!>          LSCALE is REAL array, dimension (N)
!>          Details of the permutations and/or scaling factors applied
!>          to the left side of A and B, as returned by CGGBAL.
!> \endverbatim
!>
!> \param[in] RSCALE
!> \verbatim
!>          RSCALE is REAL array, dimension (N)
!>          Details of the permutations and/or scaling factors applied
!>          to the right side of A and B, as returned by CGGBAL.
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
!>          V is COMPLEX array, dimension (LDV,M)
!>          On entry, the matrix of right or left eigenvectors to be
!>          transformed, as returned by CTGEVC.
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
!> \ingroup complexGBcomputational
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
      SUBROUTINE CGGBAK(Job,Side,N,Ilo,Ihi,Lscale,Rscale,M,V,Ldv,Info)
      IMPLICIT NONE
!*--CGGBAK151
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Job , Side
      INTEGER Ihi , Ilo , Info , Ldv , M , N
!     ..
!     .. Array Arguments ..
      REAL Lscale(*) , Rscale(*)
      COMPLEX V(Ldv,*)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      LOGICAL leftv , rightv
      INTEGER i , k
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL CSSCAL , CSWAP , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX
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
         CALL XERBLA('CGGBAK',-Info)
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
                  CALL CSSCAL(M,Rscale(i),V(i,1),Ldv)
               ENDDO
            ENDIF
!
!        Backward transformation on left eigenvectors
!
            IF ( leftv ) THEN
               DO i = Ilo , Ihi
                  CALL CSSCAL(M,Lscale(i),V(i,1),Ldv)
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
               DO i = Ilo - 1 , 1 , -1
                  k = Rscale(i)
                  IF ( k/=i ) CALL CSWAP(M,V(i,1),Ldv,V(k,1),Ldv)
               ENDDO
            ENDIF
!
            IF ( Ihi/=N ) THEN
               DO i = Ihi + 1 , N
                  k = Rscale(i)
                  IF ( k/=i ) CALL CSWAP(M,V(i,1),Ldv,V(k,1),Ldv)
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
                  IF ( k/=i ) CALL CSWAP(M,V(i,1),Ldv,V(k,1),Ldv)
               ENDDO
            ENDIF
!
            IF ( Ihi/=N ) THEN
               DO i = Ihi + 1 , N
                  k = Lscale(i)
                  IF ( k/=i ) CALL CSWAP(M,V(i,1),Ldv,V(k,1),Ldv)
               ENDDO
            ENDIF
         ENDIF
      ENDIF
!
!
!
!     End of CGGBAK
!
      END SUBROUTINE CGGBAK
