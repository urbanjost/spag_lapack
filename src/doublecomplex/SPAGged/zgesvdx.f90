!*==zgesvdx.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief <b> ZGESVDX computes the singular value decomposition (SVD) for GE matrices</b>
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZGESVDX + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgesvdx.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgesvdx.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgesvdx.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!     SUBROUTINE ZGESVDX( JOBU, JOBVT, RANGE, M, N, A, LDA, VL, VU,
!    $                    IL, IU, NS, S, U, LDU, VT, LDVT, WORK,
!    $                    LWORK, RWORK, IWORK, INFO )
!
!
!     .. Scalar Arguments ..
!      CHARACTER          JOBU, JOBVT, RANGE
!      INTEGER            IL, INFO, IU, LDA, LDU, LDVT, LWORK, M, N, NS
!      DOUBLE PRECISION   VL, VU
!     ..
!     .. Array Arguments ..
!      INTEGER            IWORK( * )
!      DOUBLE PRECISION   S( * ), RWORK( * )
!      COMPLEX*16         A( LDA, * ), U( LDU, * ), VT( LDVT, * ),
!     $                   WORK( * )
!     ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>  ZGESVDX computes the singular value decomposition (SVD) of a complex
!>  M-by-N matrix A, optionally computing the left and/or right singular
!>  vectors. The SVD is written
!>
!>      A = U * SIGMA * transpose(V)
!>
!>  where SIGMA is an M-by-N matrix which is zero except for its
!>  min(m,n) diagonal elements, U is an M-by-M unitary matrix, and
!>  V is an N-by-N unitary matrix.  The diagonal elements of SIGMA
!>  are the singular values of A; they are real and non-negative, and
!>  are returned in descending order.  The first min(m,n) columns of
!>  U and V are the left and right singular vectors of A.
!>
!>  ZGESVDX uses an eigenvalue problem for obtaining the SVD, which
!>  allows for the computation of a subset of singular values and
!>  vectors. See DBDSVDX for details.
!>
!>  Note that the routine returns V**T, not V.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] JOBU
!> \verbatim
!>          JOBU is CHARACTER*1
!>          Specifies options for computing all or part of the matrix U:
!>          = 'V':  the first min(m,n) columns of U (the left singular
!>                  vectors) or as specified by RANGE are returned in
!>                  the array U;
!>          = 'N':  no columns of U (no left singular vectors) are
!>                  computed.
!> \endverbatim
!>
!> \param[in] JOBVT
!> \verbatim
!>          JOBVT is CHARACTER*1
!>           Specifies options for computing all or part of the matrix
!>           V**T:
!>           = 'V':  the first min(m,n) rows of V**T (the right singular
!>                   vectors) or as specified by RANGE are returned in
!>                   the array VT;
!>           = 'N':  no rows of V**T (no right singular vectors) are
!>                   computed.
!> \endverbatim
!>
!> \param[in] RANGE
!> \verbatim
!>          RANGE is CHARACTER*1
!>          = 'A': all singular values will be found.
!>          = 'V': all singular values in the half-open interval (VL,VU]
!>                 will be found.
!>          = 'I': the IL-th through IU-th singular values will be found.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the input matrix A.  M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the input matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          On entry, the M-by-N matrix A.
!>          On exit, the contents of A are destroyed.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,M).
!> \endverbatim
!>
!> \param[in] VL
!> \verbatim
!>          VL is DOUBLE PRECISION
!>          If RANGE='V', the lower bound of the interval to
!>          be searched for singular values. VU > VL.
!>          Not referenced if RANGE = 'A' or 'I'.
!> \endverbatim
!>
!> \param[in] VU
!> \verbatim
!>          VU is DOUBLE PRECISION
!>          If RANGE='V', the upper bound of the interval to
!>          be searched for singular values. VU > VL.
!>          Not referenced if RANGE = 'A' or 'I'.
!> \endverbatim
!>
!> \param[in] IL
!> \verbatim
!>          IL is INTEGER
!>          If RANGE='I', the index of the
!>          smallest singular value to be returned.
!>          1 <= IL <= IU <= min(M,N), if min(M,N) > 0.
!>          Not referenced if RANGE = 'A' or 'V'.
!> \endverbatim
!>
!> \param[in] IU
!> \verbatim
!>          IU is INTEGER
!>          If RANGE='I', the index of the
!>          largest singular value to be returned.
!>          1 <= IL <= IU <= min(M,N), if min(M,N) > 0.
!>          Not referenced if RANGE = 'A' or 'V'.
!> \endverbatim
!>
!> \param[out] NS
!> \verbatim
!>          NS is INTEGER
!>          The total number of singular values found,
!>          0 <= NS <= min(M,N).
!>          If RANGE = 'A', NS = min(M,N); if RANGE = 'I', NS = IU-IL+1.
!> \endverbatim
!>
!> \param[out] S
!> \verbatim
!>          S is DOUBLE PRECISION array, dimension (min(M,N))
!>          The singular values of A, sorted so that S(i) >= S(i+1).
!> \endverbatim
!>
!> \param[out] U
!> \verbatim
!>          U is COMPLEX*16 array, dimension (LDU,UCOL)
!>          If JOBU = 'V', U contains columns of U (the left singular
!>          vectors, stored columnwise) as specified by RANGE; if
!>          JOBU = 'N', U is not referenced.
!>          Note: The user must ensure that UCOL >= NS; if RANGE = 'V',
!>          the exact value of NS is not known in advance and an upper
!>          bound must be used.
!> \endverbatim
!>
!> \param[in] LDU
!> \verbatim
!>          LDU is INTEGER
!>          The leading dimension of the array U.  LDU >= 1; if
!>          JOBU = 'V', LDU >= M.
!> \endverbatim
!>
!> \param[out] VT
!> \verbatim
!>          VT is COMPLEX*16 array, dimension (LDVT,N)
!>          If JOBVT = 'V', VT contains the rows of V**T (the right singular
!>          vectors, stored rowwise) as specified by RANGE; if JOBVT = 'N',
!>          VT is not referenced.
!>          Note: The user must ensure that LDVT >= NS; if RANGE = 'V',
!>          the exact value of NS is not known in advance and an upper
!>          bound must be used.
!> \endverbatim
!>
!> \param[in] LDVT
!> \verbatim
!>          LDVT is INTEGER
!>          The leading dimension of the array VT.  LDVT >= 1; if
!>          JOBVT = 'V', LDVT >= NS (see above).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK;
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK.
!>          LWORK >= MAX(1,MIN(M,N)*(MIN(M,N)+4)) for the paths (see
!>          comments inside the code):
!>             - PATH 1  (M much larger than N)
!>             - PATH 1t (N much larger than M)
!>          LWORK >= MAX(1,MIN(M,N)*2+MAX(M,N)) for the other paths.
!>          For good performance, LWORK should generally be larger.
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the WORK array, returns
!>          this value as the first entry of the WORK array, and no error
!>          message related to LWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array, dimension (MAX(1,LRWORK))
!>          LRWORK >= MIN(M,N)*(MIN(M,N)*2+15*MIN(M,N)).
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (12*MIN(M,N))
!>          If INFO = 0, the first NS elements of IWORK are zero. If INFO > 0,
!>          then IWORK contains the indices of the eigenvectors that failed
!>          to converge in DBDSVDX/DSTEVX.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>     INFO is INTEGER
!>           = 0:  successful exit
!>           < 0:  if INFO = -i, the i-th argument had an illegal value
!>           > 0:  if INFO = i, then i eigenvectors failed to converge
!>                 in DBDSVDX/DSTEVX.
!>                 if INFO = N*2 + 1, an internal error occurred in
!>                 DBDSVDX
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
!> \date June 2016
!
!> \ingroup complex16GEsing
!
!  =====================================================================
      SUBROUTINE ZGESVDX(Jobu,Jobvt,Range,M,N,A,Lda,Vl,Vu,Il,Iu,Ns,S,U, &
     &                   Ldu,Vt,Ldvt,Work,Lwork,Rwork,Iwork,Info)
      USE F77KINDS                        
      USE S_DBDSVDX
      USE S_DLAMCH
      USE S_DLASCL
      USE S_ILAENV
      USE S_LSAME
      USE S_XERBLA
      USE S_ZGEBRD
      USE S_ZGELQF
      USE S_ZGEQRF
      USE S_ZLACPY
      USE S_ZLANGE
      USE S_ZLASCL
      USE S_ZLASET
      USE S_ZUNMBR
      USE S_ZUNMLQ
      USE S_ZUNMQR
      IMPLICIT NONE
!*--ZGESVDX290
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX(CX16KIND) , PARAMETER  ::  CZERO = (0.0D0,0.0D0)
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Jobu
      CHARACTER :: Jobvt
      CHARACTER :: Range
      INTEGER :: M
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) :: Vl
      REAL(R8KIND) :: Vu
      INTEGER , INTENT(IN) :: Il
      INTEGER , INTENT(IN) :: Iu
      INTEGER , INTENT(INOUT) :: Ns
      REAL(R8KIND) , DIMENSION(*) :: S
      COMPLEX(CX16KIND) , DIMENSION(Ldu,*) :: U
      INTEGER :: Ldu
      COMPLEX(CX16KIND) , DIMENSION(Ldvt,*) :: Vt
      INTEGER :: Ldvt
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: abstol , anrm , bignum , eps , smlnum
      LOGICAL :: alls , inds , lquery , vals , wantu , wantvt
      REAL(R8KIND) , DIMENSION(1) :: dum
      INTEGER :: i , id , ie , ierr , ilqf , iltgk , iqrf , iscl ,      &
     &           itau , itaup , itauq , itemp , itempr , itgkz , iutgk ,&
     &           j , k , maxwrk , minmn , minwrk , mnthr
      CHARACTER :: jobz , rngtgk
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
!     .. Local Arrays ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. External Functions ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments.
!
      Ns = 0
      Info = 0
      abstol = 2*DLAMCH('S')
      lquery = (Lwork==-1)
      minmn = MIN(M,N)
 
      wantu = LSAME(Jobu,'V')
      wantvt = LSAME(Jobvt,'V')
      IF ( wantu .OR. wantvt ) THEN
         jobz = 'V'
      ELSE
         jobz = 'N'
      ENDIF
      alls = LSAME(Range,'A')
      vals = LSAME(Range,'V')
      inds = LSAME(Range,'I')
!
      Info = 0
      IF ( .NOT.LSAME(Jobu,'V') .AND. .NOT.LSAME(Jobu,'N') ) THEN
         Info = -1
      ELSEIF ( .NOT.LSAME(Jobvt,'V') .AND. .NOT.LSAME(Jobvt,'N') ) THEN
         Info = -2
      ELSEIF ( .NOT.(alls .OR. vals .OR. inds) ) THEN
         Info = -3
      ELSEIF ( M<0 ) THEN
         Info = -4
      ELSEIF ( N<0 ) THEN
         Info = -5
      ELSEIF ( M>Lda ) THEN
         Info = -7
      ELSEIF ( minmn>0 ) THEN
         IF ( vals ) THEN
            IF ( Vl<ZERO ) THEN
               Info = -8
            ELSEIF ( Vu<=Vl ) THEN
               Info = -9
            ENDIF
         ELSEIF ( inds ) THEN
            IF ( Il<1 .OR. Il>MAX(1,minmn) ) THEN
               Info = -10
            ELSEIF ( Iu<MIN(minmn,Il) .OR. Iu>minmn ) THEN
               Info = -11
            ENDIF
         ENDIF
         IF ( Info==0 ) THEN
            IF ( wantu .AND. Ldu<M ) THEN
               Info = -15
            ELSEIF ( wantvt ) THEN
               IF ( inds ) THEN
                  IF ( Ldvt<Iu-Il+1 ) Info = -17
               ELSEIF ( Ldvt<minmn ) THEN
                  Info = -17
               ENDIF
            ENDIF
         ENDIF
      ENDIF
!
!     Compute workspace
!     (Note: Comments in the code beginning "Workspace:" describe the
!     minimal amount of workspace needed at that point in the code,
!     as well as the preferred amount for good performance.
!     NB refers to the optimal block size for the immediately
!     following subroutine, as returned by ILAENV.)
!
      IF ( Info==0 ) THEN
         minwrk = 1
         maxwrk = 1
         IF ( minmn>0 ) THEN
            IF ( M>=N ) THEN
               mnthr = ILAENV(6,'ZGESVD',Jobu//Jobvt,M,N,0,0)
               IF ( M>=mnthr ) THEN
!
!                 Path 1 (M much larger than N)
!
                  minwrk = N*(N+5)
                  maxwrk = N + N*ILAENV(1,'ZGEQRF',' ',M,N,-1,-1)
                  maxwrk = MAX(maxwrk,                                  &
     &                     N*N+2*N+2*N*ILAENV(1,'ZGEBRD',' ',N,N,-1,-1))
                  IF ( wantu .OR. wantvt )                              &
     &                 maxwrk = MAX(maxwrk,N*N+2*N+N*ILAENV(1,'ZUNMQR', &
     &                 'LN',N,N,N,-1))
               ELSE
!
!                 Path 2 (M at least N, but not much larger)
!
                  minwrk = 3*N + M
                  maxwrk = 2*N + (M+N)*ILAENV(1,'ZGEBRD',' ',M,N,-1,-1)
                  IF ( wantu .OR. wantvt )                              &
     &                 maxwrk = MAX(maxwrk,2*N+N*ILAENV(1,'ZUNMQR','LN',&
     &                 N,N,N,-1))
               ENDIF
            ELSE
               mnthr = ILAENV(6,'ZGESVD',Jobu//Jobvt,M,N,0,0)
               IF ( N>=mnthr ) THEN
!
!                 Path 1t (N much larger than M)
!
                  minwrk = M*(M+5)
                  maxwrk = M + M*ILAENV(1,'ZGELQF',' ',M,N,-1,-1)
                  maxwrk = MAX(maxwrk,                                  &
     &                     M*M+2*M+2*M*ILAENV(1,'ZGEBRD',' ',M,M,-1,-1))
                  IF ( wantu .OR. wantvt )                              &
     &                 maxwrk = MAX(maxwrk,M*M+2*M+M*ILAENV(1,'ZUNMQR', &
     &                 'LN',M,M,M,-1))
               ELSE
!
!                 Path 2t (N greater than M, but not much larger)
!
!
                  minwrk = 3*M + N
                  maxwrk = 2*M + (M+N)*ILAENV(1,'ZGEBRD',' ',M,N,-1,-1)
                  IF ( wantu .OR. wantvt )                              &
     &                 maxwrk = MAX(maxwrk,2*M+M*ILAENV(1,'ZUNMQR','LN',&
     &                 M,M,M,-1))
               ENDIF
            ENDIF
         ENDIF
         maxwrk = MAX(maxwrk,minwrk)
         Work(1) = DCMPLX(DBLE(maxwrk),ZERO)
!
         IF ( Lwork<minwrk .AND. .NOT.lquery ) Info = -19
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('ZGESVDX',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( M==0 .OR. N==0 ) RETURN
!
!     Set singular values indices accord to RANGE='A'.
!
      IF ( alls ) THEN
         rngtgk = 'I'
         iltgk = 1
         iutgk = MIN(M,N)
      ELSEIF ( inds ) THEN
         rngtgk = 'I'
         iltgk = Il
         iutgk = Iu
      ELSE
         rngtgk = 'V'
         iltgk = 0
         iutgk = 0
      ENDIF
!
!     Get machine constants
!
      eps = DLAMCH('P')
      smlnum = SQRT(DLAMCH('S'))/eps
      bignum = ONE/smlnum
!
!     Scale A if max element outside range [SMLNUM,BIGNUM]
!
      anrm = ZLANGE('M',M,N,A,Lda,dum)
      iscl = 0
      IF ( anrm>ZERO .AND. anrm<smlnum ) THEN
         iscl = 1
         CALL ZLASCL('G',0,0,anrm,smlnum,M,N,A,Lda,Info)
      ELSEIF ( anrm>bignum ) THEN
         iscl = 1
         CALL ZLASCL('G',0,0,anrm,bignum,M,N,A,Lda,Info)
      ENDIF
!
      IF ( M>=N ) THEN
!
!        A has at least as many rows as columns. If A has sufficiently
!        more rows than columns, first reduce A using the QR
!        decomposition.
!
         IF ( M>=mnthr ) THEN
!
!           Path 1 (M much larger than N):
!           A = Q * R = Q * ( QB * B * PB**T )
!                     = Q * ( QB * ( UB * S * VB**T ) * PB**T )
!           U = Q * QB * UB; V**T = VB**T * PB**T
!
!           Compute A=Q*R
!           (Workspace: need 2*N, prefer N+N*NB)
!
            itau = 1
            itemp = itau + N
            CALL ZGEQRF(M,N,A,Lda,Work(itau),Work(itemp),Lwork-itemp+1, &
     &                  Info)
!
!           Copy R into WORK and bidiagonalize it:
!           (Workspace: need N*N+3*N, prefer N*N+N+2*N*NB)
!
            iqrf = itemp
            itauq = itemp + N*N
            itaup = itauq + N
            itemp = itaup + N
            id = 1
            ie = id + N
            itgkz = ie + N
            CALL ZLACPY('U',N,N,A,Lda,Work(iqrf),N)
            CALL ZLASET('L',N-1,N-1,CZERO,CZERO,Work(iqrf+1),N)
            CALL ZGEBRD(N,N,Work(iqrf),N,Rwork(id),Rwork(ie),Work(itauq)&
     &                  ,Work(itaup),Work(itemp),Lwork-itemp+1,Info)
            itempr = itgkz + N*(N*2+1)
!
!           Solve eigenvalue problem TGK*Z=Z*S.
!           (Workspace: need 2*N*N+14*N)
!
            CALL DBDSVDX('U',jobz,rngtgk,N,Rwork(id),Rwork(ie),Vl,Vu,   &
     &                   iltgk,iutgk,Ns,S,Rwork(itgkz),N*2,Rwork(itempr)&
     &                   ,Iwork,Info)
!
!           If needed, compute left singular vectors.
!
            IF ( wantu ) THEN
               k = itgkz
               DO i = 1 , Ns
                  DO j = 1 , N
                     U(j,i) = DCMPLX(Rwork(k),ZERO)
                     k = k + 1
                  ENDDO
                  k = k + N
               ENDDO
               CALL ZLASET('A',M-N,Ns,CZERO,CZERO,U(N+1,1),Ldu)
!
!              Call ZUNMBR to compute QB*UB.
!              (Workspace in WORK( ITEMP ): need N, prefer N*NB)
!
               CALL ZUNMBR('Q','L','N',N,Ns,N,Work(iqrf),N,Work(itauq), &
     &                     U,Ldu,Work(itemp),Lwork-itemp+1,Info)
!
!              Call ZUNMQR to compute Q*(QB*UB).
!              (Workspace in WORK( ITEMP ): need N, prefer N*NB)
!
               CALL ZUNMQR('L','N',M,Ns,N,A,Lda,Work(itau),U,Ldu,       &
     &                     Work(itemp),Lwork-itemp+1,Info)
            ENDIF
!
!           If needed, compute right singular vectors.
!
            IF ( wantvt ) THEN
               k = itgkz + N
               DO i = 1 , Ns
                  DO j = 1 , N
                     Vt(i,j) = DCMPLX(Rwork(k),ZERO)
                     k = k + 1
                  ENDDO
                  k = k + N
               ENDDO
!
!              Call ZUNMBR to compute VB**T * PB**T
!              (Workspace in WORK( ITEMP ): need N, prefer N*NB)
!
               CALL ZUNMBR('P','R','C',Ns,N,N,Work(iqrf),N,Work(itaup), &
     &                     Vt,Ldvt,Work(itemp),Lwork-itemp+1,Info)
            ENDIF
         ELSE
!
!           Path 2 (M at least N, but not much larger)
!           Reduce A to bidiagonal form without QR decomposition
!           A = QB * B * PB**T = QB * ( UB * S * VB**T ) * PB**T
!           U = QB * UB; V**T = VB**T * PB**T
!
!           Bidiagonalize A
!           (Workspace: need 2*N+M, prefer 2*N+(M+N)*NB)
!
            itauq = 1
            itaup = itauq + N
            itemp = itaup + N
            id = 1
            ie = id + N
            itgkz = ie + N
            CALL ZGEBRD(M,N,A,Lda,Rwork(id),Rwork(ie),Work(itauq),      &
     &                  Work(itaup),Work(itemp),Lwork-itemp+1,Info)
            itempr = itgkz + N*(N*2+1)
!
!           Solve eigenvalue problem TGK*Z=Z*S.
!           (Workspace: need 2*N*N+14*N)
!
            CALL DBDSVDX('U',jobz,rngtgk,N,Rwork(id),Rwork(ie),Vl,Vu,   &
     &                   iltgk,iutgk,Ns,S,Rwork(itgkz),N*2,Rwork(itempr)&
     &                   ,Iwork,Info)
!
!           If needed, compute left singular vectors.
!
            IF ( wantu ) THEN
               k = itgkz
               DO i = 1 , Ns
                  DO j = 1 , N
                     U(j,i) = DCMPLX(Rwork(k),ZERO)
                     k = k + 1
                  ENDDO
                  k = k + N
               ENDDO
               CALL ZLASET('A',M-N,Ns,CZERO,CZERO,U(N+1,1),Ldu)
!
!              Call ZUNMBR to compute QB*UB.
!              (Workspace in WORK( ITEMP ): need N, prefer N*NB)
!
               CALL ZUNMBR('Q','L','N',M,Ns,N,A,Lda,Work(itauq),U,Ldu,  &
     &                     Work(itemp),Lwork-itemp+1,ierr)
            ENDIF
!
!           If needed, compute right singular vectors.
!
            IF ( wantvt ) THEN
               k = itgkz + N
               DO i = 1 , Ns
                  DO j = 1 , N
                     Vt(i,j) = DCMPLX(Rwork(k),ZERO)
                     k = k + 1
                  ENDDO
                  k = k + N
               ENDDO
!
!              Call ZUNMBR to compute VB**T * PB**T
!              (Workspace in WORK( ITEMP ): need N, prefer N*NB)
!
               CALL ZUNMBR('P','R','C',Ns,N,N,A,Lda,Work(itaup),Vt,Ldvt,&
     &                     Work(itemp),Lwork-itemp+1,ierr)
            ENDIF
         ENDIF
!
!        A has more columns than rows. If A has sufficiently more
!        columns than rows, first reduce A using the LQ decomposition.
!
      ELSEIF ( N>=mnthr ) THEN
!
!           Path 1t (N much larger than M):
!           A = L * Q = ( QB * B * PB**T ) * Q
!                     = ( QB * ( UB * S * VB**T ) * PB**T ) * Q
!           U = QB * UB ; V**T = VB**T * PB**T * Q
!
!           Compute A=L*Q
!           (Workspace: need 2*M, prefer M+M*NB)
!
         itau = 1
         itemp = itau + M
         CALL ZGELQF(M,N,A,Lda,Work(itau),Work(itemp),Lwork-itemp+1,    &
     &               Info)
 
!           Copy L into WORK and bidiagonalize it:
!           (Workspace in WORK( ITEMP ): need M*M+3*M, prefer M*M+M+2*M*NB)
!
         ilqf = itemp
         itauq = ilqf + M*M
         itaup = itauq + M
         itemp = itaup + M
         id = 1
         ie = id + M
         itgkz = ie + M
         CALL ZLACPY('L',M,M,A,Lda,Work(ilqf),M)
         CALL ZLASET('U',M-1,M-1,CZERO,CZERO,Work(ilqf+M),M)
         CALL ZGEBRD(M,M,Work(ilqf),M,Rwork(id),Rwork(ie),Work(itauq),  &
     &               Work(itaup),Work(itemp),Lwork-itemp+1,Info)
         itempr = itgkz + M*(M*2+1)
!
!           Solve eigenvalue problem TGK*Z=Z*S.
!           (Workspace: need 2*M*M+14*M)
!
         CALL DBDSVDX('U',jobz,rngtgk,M,Rwork(id),Rwork(ie),Vl,Vu,iltgk,&
     &                iutgk,Ns,S,Rwork(itgkz),M*2,Rwork(itempr),Iwork,  &
     &                Info)
!
!           If needed, compute left singular vectors.
!
         IF ( wantu ) THEN
            k = itgkz
            DO i = 1 , Ns
               DO j = 1 , M
                  U(j,i) = DCMPLX(Rwork(k),ZERO)
                  k = k + 1
               ENDDO
               k = k + M
            ENDDO
!
!              Call ZUNMBR to compute QB*UB.
!              (Workspace in WORK( ITEMP ): need M, prefer M*NB)
!
            CALL ZUNMBR('Q','L','N',M,Ns,M,Work(ilqf),M,Work(itauq),U,  &
     &                  Ldu,Work(itemp),Lwork-itemp+1,Info)
         ENDIF
!
!           If needed, compute right singular vectors.
!
         IF ( wantvt ) THEN
            k = itgkz + M
            DO i = 1 , Ns
               DO j = 1 , M
                  Vt(i,j) = DCMPLX(Rwork(k),ZERO)
                  k = k + 1
               ENDDO
               k = k + M
            ENDDO
            CALL ZLASET('A',Ns,N-M,CZERO,CZERO,Vt(1,M+1),Ldvt)
!
!              Call ZUNMBR to compute (VB**T)*(PB**T)
!              (Workspace in WORK( ITEMP ): need M, prefer M*NB)
!
            CALL ZUNMBR('P','R','C',Ns,M,M,Work(ilqf),M,Work(itaup),Vt, &
     &                  Ldvt,Work(itemp),Lwork-itemp+1,Info)
!
!              Call ZUNMLQ to compute ((VB**T)*(PB**T))*Q.
!              (Workspace in WORK( ITEMP ): need M, prefer M*NB)
!
            CALL ZUNMLQ('R','N',Ns,N,M,A,Lda,Work(itau),Vt,Ldvt,        &
     &                  Work(itemp),Lwork-itemp+1,Info)
         ENDIF
      ELSE
!
!           Path 2t (N greater than M, but not much larger)
!           Reduce to bidiagonal form without LQ decomposition
!           A = QB * B * PB**T = QB * ( UB * S * VB**T ) * PB**T
!           U = QB * UB; V**T = VB**T * PB**T
!
!           Bidiagonalize A
!           (Workspace: need 2*M+N, prefer 2*M+(M+N)*NB)
!
         itauq = 1
         itaup = itauq + M
         itemp = itaup + M
         id = 1
         ie = id + M
         itgkz = ie + M
         CALL ZGEBRD(M,N,A,Lda,Rwork(id),Rwork(ie),Work(itauq),         &
     &               Work(itaup),Work(itemp),Lwork-itemp+1,Info)
         itempr = itgkz + M*(M*2+1)
!
!           Solve eigenvalue problem TGK*Z=Z*S.
!           (Workspace: need 2*M*M+14*M)
!
         CALL DBDSVDX('L',jobz,rngtgk,M,Rwork(id),Rwork(ie),Vl,Vu,iltgk,&
     &                iutgk,Ns,S,Rwork(itgkz),M*2,Rwork(itempr),Iwork,  &
     &                Info)
!
!           If needed, compute left singular vectors.
!
         IF ( wantu ) THEN
            k = itgkz
            DO i = 1 , Ns
               DO j = 1 , M
                  U(j,i) = DCMPLX(Rwork(k),ZERO)
                  k = k + 1
               ENDDO
               k = k + M
            ENDDO
!
!              Call ZUNMBR to compute QB*UB.
!              (Workspace in WORK( ITEMP ): need M, prefer M*NB)
!
            CALL ZUNMBR('Q','L','N',M,Ns,N,A,Lda,Work(itauq),U,Ldu,     &
     &                  Work(itemp),Lwork-itemp+1,Info)
         ENDIF
!
!           If needed, compute right singular vectors.
!
         IF ( wantvt ) THEN
            k = itgkz + M
            DO i = 1 , Ns
               DO j = 1 , M
                  Vt(i,j) = DCMPLX(Rwork(k),ZERO)
                  k = k + 1
               ENDDO
               k = k + M
            ENDDO
            CALL ZLASET('A',Ns,N-M,CZERO,CZERO,Vt(1,M+1),Ldvt)
!
!              Call ZUNMBR to compute VB**T * PB**T
!              (Workspace in WORK( ITEMP ): need M, prefer M*NB)
!
            CALL ZUNMBR('P','R','C',Ns,N,M,A,Lda,Work(itaup),Vt,Ldvt,   &
     &                  Work(itemp),Lwork-itemp+1,Info)
         ENDIF
      ENDIF
!
!     Undo scaling if necessary
!
      IF ( iscl==1 ) THEN
         IF ( anrm>bignum ) CALL DLASCL('G',0,0,bignum,anrm,minmn,1,S,  &
     &                                  minmn,Info)
         IF ( anrm<smlnum ) CALL DLASCL('G',0,0,smlnum,anrm,minmn,1,S,  &
     &                                  minmn,Info)
      ENDIF
!
!     Return optimal workspace in WORK(1)
!
      Work(1) = DCMPLX(DBLE(maxwrk),ZERO)
!
!
!     End of ZGESVDX
!
      END SUBROUTINE ZGESVDX
