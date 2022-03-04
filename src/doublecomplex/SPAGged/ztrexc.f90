!*==ztrexc.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b ZTREXC
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZTREXC + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztrexc.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztrexc.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztrexc.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZTREXC( COMPQ, N, T, LDT, Q, LDQ, IFST, ILST, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          COMPQ
!       INTEGER            IFST, ILST, INFO, LDQ, LDT, N
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         Q( LDQ, * ), T( LDT, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZTREXC reorders the Schur factorization of a complex matrix
!> A = Q*T*Q**H, so that the diagonal element of T with row index IFST
!> is moved to row ILST.
!>
!> The Schur form T is reordered by a unitary similarity transformation
!> Z**H*T*Z, and optionally the matrix Q of Schur vectors is updated by
!> postmultplying it with Z.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] COMPQ
!> \verbatim
!>          COMPQ is CHARACTER*1
!>          = 'V':  update the matrix Q of Schur vectors;
!>          = 'N':  do not update Q.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix T. N >= 0.
!>          If N == 0 arguments ILST and IFST may be any value.
!> \endverbatim
!>
!> \param[in,out] T
!> \verbatim
!>          T is COMPLEX*16 array, dimension (LDT,N)
!>          On entry, the upper triangular matrix T.
!>          On exit, the reordered upper triangular matrix.
!> \endverbatim
!>
!> \param[in] LDT
!> \verbatim
!>          LDT is INTEGER
!>          The leading dimension of the array T. LDT >= max(1,N).
!> \endverbatim
!>
!> \param[in,out] Q
!> \verbatim
!>          Q is COMPLEX*16 array, dimension (LDQ,N)
!>          On entry, if COMPQ = 'V', the matrix Q of Schur vectors.
!>          On exit, if COMPQ = 'V', Q has been postmultiplied by the
!>          unitary transformation matrix Z which reorders T.
!>          If COMPQ = 'N', Q is not referenced.
!> \endverbatim
!>
!> \param[in] LDQ
!> \verbatim
!>          LDQ is INTEGER
!>          The leading dimension of the array Q.  LDQ >= 1, and if
!>          COMPQ = 'V', LDQ >= max(1,N).
!> \endverbatim
!>
!> \param[in] IFST
!> \verbatim
!>          IFST is INTEGER
!> \endverbatim
!>
!> \param[in] ILST
!> \verbatim
!>          ILST is INTEGER
!>
!>          Specify the reordering of the diagonal elements of T:
!>          The element with row index IFST is moved to row ILST by a
!>          sequence of transpositions between adjacent elements.
!>          1 <= IFST <= N; 1 <= ILST <= N.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
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
!> \ingroup complex16OTHERcomputational
!
!  =====================================================================
      SUBROUTINE ZTREXC(Compq,N,T,Ldt,Q,Ldq,Ifst,Ilst,Info)
      USE F77KINDS                        
      USE S_LSAME
      USE S_XERBLA
      USE S_ZLARTG
      USE S_ZROT
      IMPLICIT NONE
!*--ZTREXC135
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Compq
      INTEGER :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      COMPLEX(CX16KIND) , DIMENSION(Ldq,*) :: Q
      INTEGER , INTENT(IN) :: Ldq
      INTEGER , INTENT(IN) :: Ifst
      INTEGER , INTENT(IN) :: Ilst
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: cs
      INTEGER :: k , m1 , m2 , m3
      COMPLEX(CX16KIND) :: sn , t11 , t22 , temp
      LOGICAL :: wantq
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
!     Decode and test the input parameters.
!
      Info = 0
      wantq = LSAME(Compq,'V')
      IF ( .NOT.LSAME(Compq,'N') .AND. .NOT.wantq ) THEN
         Info = -1
      ELSEIF ( N<0 ) THEN
         Info = -2
      ELSEIF ( Ldt<MAX(1,N) ) THEN
         Info = -4
      ELSEIF ( Ldq<1 .OR. (wantq .AND. Ldq<MAX(1,N)) ) THEN
         Info = -6
      ELSEIF ( (Ifst<1 .OR. Ifst>N) .AND. (N>0) ) THEN
         Info = -7
      ELSEIF ( (Ilst<1 .OR. Ilst>N) .AND. (N>0) ) THEN
         Info = -8
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('ZTREXC',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N<=1 .OR. Ifst==Ilst ) RETURN
!
      IF ( Ifst<Ilst ) THEN
!
!        Move the IFST-th diagonal element forward down the diagonal.
!
         m1 = 0
         m2 = -1
         m3 = 1
      ELSE
!
!        Move the IFST-th diagonal element backward up the diagonal.
!
         m1 = -1
         m2 = 0
         m3 = -1
      ENDIF
!
      DO k = Ifst + m1 , Ilst + m2 , m3
!
!        Interchange the k-th and (k+1)-th diagonal elements.
!
         t11 = T(k,k)
         t22 = T(k+1,k+1)
!
!        Determine the transformation to perform the interchange.
!
         CALL ZLARTG(T(k,k+1),t22-t11,cs,sn,temp)
!
!        Apply transformation to the matrix T.
!
         IF ( k+2<=N ) CALL ZROT(N-k-1,T(k,k+2),Ldt,T(k+1,k+2),Ldt,cs,  &
     &                           sn)
         CALL ZROT(k-1,T(1,k),1,T(1,k+1),1,cs,DCONJG(sn))
!
         T(k,k) = t22
         T(k+1,k+1) = t11
!
!
!           Accumulate transformation in the matrix Q.
!
         IF ( wantq ) CALL ZROT(N,Q(1,k),1,Q(1,k+1),1,cs,DCONJG(sn))
!
      ENDDO
!
!
!     End of ZTREXC
!
      END SUBROUTINE ZTREXC
