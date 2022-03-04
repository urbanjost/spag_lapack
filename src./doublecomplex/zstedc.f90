!*==zstedc.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b ZSTEDC
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZSTEDC + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zstedc.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zstedc.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zstedc.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZSTEDC( COMPZ, N, D, E, Z, LDZ, WORK, LWORK, RWORK,
!                          LRWORK, IWORK, LIWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          COMPZ
!       INTEGER            INFO, LDZ, LIWORK, LRWORK, LWORK, N
!       ..
!       .. Array Arguments ..
!       INTEGER            IWORK( * )
!       DOUBLE PRECISION   D( * ), E( * ), RWORK( * )
!       COMPLEX*16         WORK( * ), Z( LDZ, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZSTEDC computes all eigenvalues and, optionally, eigenvectors of a
!> symmetric tridiagonal matrix using the divide and conquer method.
!> The eigenvectors of a full or band complex Hermitian matrix can also
!> be found if ZHETRD or ZHPTRD or ZHBTRD has been used to reduce this
!> matrix to tridiagonal form.
!>
!> This code makes very mild assumptions about floating point
!> arithmetic. It will work on machines with a guard digit in
!> add/subtract, or on those binary machines without guard digits
!> which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or Cray-2.
!> It could conceivably fail on hexadecimal or decimal machines
!> without guard digits, but we know of none.  See DLAED3 for details.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] COMPZ
!> \verbatim
!>          COMPZ is CHARACTER*1
!>          = 'N':  Compute eigenvalues only.
!>          = 'I':  Compute eigenvectors of tridiagonal matrix also.
!>          = 'V':  Compute eigenvectors of original Hermitian matrix
!>                  also.  On entry, Z contains the unitary matrix used
!>                  to reduce the original matrix to tridiagonal form.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The dimension of the symmetric tridiagonal matrix.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] D
!> \verbatim
!>          D is DOUBLE PRECISION array, dimension (N)
!>          On entry, the diagonal elements of the tridiagonal matrix.
!>          On exit, if INFO = 0, the eigenvalues in ascending order.
!> \endverbatim
!>
!> \param[in,out] E
!> \verbatim
!>          E is DOUBLE PRECISION array, dimension (N-1)
!>          On entry, the subdiagonal elements of the tridiagonal matrix.
!>          On exit, E has been destroyed.
!> \endverbatim
!>
!> \param[in,out] Z
!> \verbatim
!>          Z is COMPLEX*16 array, dimension (LDZ,N)
!>          On entry, if COMPZ = 'V', then Z contains the unitary
!>          matrix used in the reduction to tridiagonal form.
!>          On exit, if INFO = 0, then if COMPZ = 'V', Z contains the
!>          orthonormal eigenvectors of the original Hermitian matrix,
!>          and if COMPZ = 'I', Z contains the orthonormal eigenvectors
!>          of the symmetric tridiagonal matrix.
!>          If  COMPZ = 'N', then Z is not referenced.
!> \endverbatim
!>
!> \param[in] LDZ
!> \verbatim
!>          LDZ is INTEGER
!>          The leading dimension of the array Z.  LDZ >= 1.
!>          If eigenvectors are desired, then LDZ >= max(1,N).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK.
!>          If COMPZ = 'N' or 'I', or N <= 1, LWORK must be at least 1.
!>          If COMPZ = 'V' and N > 1, LWORK must be at least N*N.
!>          Note that for COMPZ = 'V', then if N is less than or
!>          equal to the minimum divide size, usually 25, then LWORK need
!>          only be 1.
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal sizes of the WORK, RWORK and
!>          IWORK arrays, returns these values as the first entries of
!>          the WORK, RWORK and IWORK arrays, and no error message
!>          related to LWORK or LRWORK or LIWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array, dimension (MAX(1,LRWORK))
!>          On exit, if INFO = 0, RWORK(1) returns the optimal LRWORK.
!> \endverbatim
!>
!> \param[in] LRWORK
!> \verbatim
!>          LRWORK is INTEGER
!>          The dimension of the array RWORK.
!>          If COMPZ = 'N' or N <= 1, LRWORK must be at least 1.
!>          If COMPZ = 'V' and N > 1, LRWORK must be at least
!>                         1 + 3*N + 2*N*lg N + 4*N**2 ,
!>                         where lg( N ) = smallest integer k such
!>                         that 2**k >= N.
!>          If COMPZ = 'I' and N > 1, LRWORK must be at least
!>                         1 + 4*N + 2*N**2 .
!>          Note that for COMPZ = 'I' or 'V', then if N is less than or
!>          equal to the minimum divide size, usually 25, then LRWORK
!>          need only be max(1,2*(N-1)).
!>
!>          If LRWORK = -1, then a workspace query is assumed; the
!>          routine only calculates the optimal sizes of the WORK, RWORK
!>          and IWORK arrays, returns these values as the first entries
!>          of the WORK, RWORK and IWORK arrays, and no error message
!>          related to LWORK or LRWORK or LIWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (MAX(1,LIWORK))
!>          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.
!> \endverbatim
!>
!> \param[in] LIWORK
!> \verbatim
!>          LIWORK is INTEGER
!>          The dimension of the array IWORK.
!>          If COMPZ = 'N' or N <= 1, LIWORK must be at least 1.
!>          If COMPZ = 'V' or N > 1,  LIWORK must be at least
!>                                    6 + 6*N + 5*N*lg N.
!>          If COMPZ = 'I' or N > 1,  LIWORK must be at least
!>                                    3 + 5*N .
!>          Note that for COMPZ = 'I' or 'V', then if N is less than or
!>          equal to the minimum divide size, usually 25, then LIWORK
!>          need only be 1.
!>
!>          If LIWORK = -1, then a workspace query is assumed; the
!>          routine only calculates the optimal sizes of the WORK, RWORK
!>          and IWORK arrays, returns these values as the first entries
!>          of the WORK, RWORK and IWORK arrays, and no error message
!>          related to LWORK or LRWORK or LIWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit.
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
!>          > 0:  The algorithm failed to compute an eigenvalue while
!>                working on the submatrix lying in rows and columns
!>                INFO/(N+1) through mod(INFO,N+1).
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
!> \date June 2017
!
!> \ingroup complex16OTHERcomputational
!
!> \par Contributors:
!  ==================
!>
!> Jeff Rutter, Computer Science Division, University of California
!> at Berkeley, USA
!
!  =====================================================================
      SUBROUTINE ZSTEDC(Compz,N,D,E,Z,Ldz,Work,Lwork,Rwork,Lrwork,Iwork,&
     &                  Liwork,Info)
      USE F77KINDS                        
      USE S_DLAMCH
      USE S_DLANST
      USE S_DLASCL
      USE S_DLASET
      USE S_DSTEDC
      USE S_DSTEQR
      USE S_DSTERF
      USE S_ILAENV
      USE S_LSAME
      USE S_XERBLA
      USE S_ZLACPY
      USE S_ZLACRM
      USE S_ZLAED0
      USE S_ZSTEQR
      USE S_ZSWAP
      IMPLICIT NONE
!*--ZSTEDC232
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0 ,        &
     &                              TWO = 2.0D0
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Compz
      INTEGER :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: D
      REAL(R8KIND) , DIMENSION(*) :: E
      COMPLEX(CX16KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(IN) :: Lrwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER :: Liwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: eps , orgnrm , p , tiny
      INTEGER :: finish , i , icompz , ii , j , k , lgn , liwmin , ll , &
     &           lrwmin , lwmin , m , smlsiz , start
      LOGICAL :: lquery
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
!     Test the input parameters.
!
      Info = 0
      lquery = (Lwork==-1 .OR. Lrwork==-1 .OR. Liwork==-1)
!
      IF ( LSAME(Compz,'N') ) THEN
         icompz = 0
      ELSEIF ( LSAME(Compz,'V') ) THEN
         icompz = 1
      ELSEIF ( LSAME(Compz,'I') ) THEN
         icompz = 2
      ELSE
         icompz = -1
      ENDIF
      IF ( icompz<0 ) THEN
         Info = -1
      ELSEIF ( N<0 ) THEN
         Info = -2
      ELSEIF ( (Ldz<1) .OR. (icompz>0 .AND. Ldz<MAX(1,N)) ) THEN
         Info = -6
      ENDIF
!
      IF ( Info==0 ) THEN
!
!        Compute the workspace requirements
!
         smlsiz = ILAENV(9,'ZSTEDC',' ',0,0,0,0)
         IF ( N<=1 .OR. icompz==0 ) THEN
            lwmin = 1
            liwmin = 1
            lrwmin = 1
         ELSEIF ( N<=smlsiz ) THEN
            lwmin = 1
            liwmin = 1
            lrwmin = 2*(N-1)
         ELSEIF ( icompz==1 ) THEN
            lgn = INT(LOG(DBLE(N))/LOG(TWO))
            IF ( 2**lgn<N ) lgn = lgn + 1
            IF ( 2**lgn<N ) lgn = lgn + 1
            lwmin = N*N
            lrwmin = 1 + 3*N + 2*N*lgn + 4*N**2
            liwmin = 6 + 6*N + 5*N*lgn
         ELSEIF ( icompz==2 ) THEN
            lwmin = 1
            lrwmin = 1 + 4*N + 2*N**2
            liwmin = 3 + 5*N
         ENDIF
         Work(1) = lwmin
         Rwork(1) = lrwmin
         Iwork(1) = liwmin
!
         IF ( Lwork<lwmin .AND. .NOT.lquery ) THEN
            Info = -8
         ELSEIF ( Lrwork<lrwmin .AND. .NOT.lquery ) THEN
            Info = -10
         ELSEIF ( Liwork<liwmin .AND. .NOT.lquery ) THEN
            Info = -12
         ENDIF
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('ZSTEDC',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 ) RETURN
      IF ( N==1 ) THEN
         IF ( icompz/=0 ) Z(1,1) = ONE
         RETURN
      ENDIF
!
!     If the following conditional clause is removed, then the routine
!     will use the Divide and Conquer routine to compute only the
!     eigenvalues, which requires (3N + 3N**2) real workspace and
!     (2 + 5N + 2N lg(N)) integer workspace.
!     Since on many architectures DSTERF is much faster than any other
!     algorithm for finding eigenvalues only, it is used here
!     as the default. If the conditional clause is removed, then
!     information on the size of workspace needs to be changed.
!
!     If COMPZ = 'N', use DSTERF to compute the eigenvalues.
!
      IF ( icompz==0 ) THEN
         CALL DSTERF(N,D,E,Info)
         GOTO 100
      ENDIF
!
!     If N is smaller than the minimum divide size (SMLSIZ+1), then
!     solve the problem with another solver.
!
      IF ( N<=smlsiz ) THEN
!
         CALL ZSTEQR(Compz,N,D,E,Z,Ldz,Rwork,Info)
!
      ELSE
!
!        If COMPZ = 'I', we simply call DSTEDC instead.
!
         IF ( icompz==2 ) THEN
            CALL DLASET('Full',N,N,ZERO,ONE,Rwork,N)
            ll = N*N + 1
            CALL DSTEDC('I',N,D,E,Rwork,N,Rwork(ll),Lrwork-ll+1,Iwork,  &
     &                  Liwork,Info)
            DO j = 1 , N
               DO i = 1 , N
                  Z(i,j) = Rwork((j-1)*N+i)
               ENDDO
            ENDDO
            GOTO 100
         ENDIF
!
!        From now on, only option left to be handled is COMPZ = 'V',
!        i.e. ICOMPZ = 1.
!
!        Scale.
!
         orgnrm = DLANST('M',N,D,E)
         IF ( orgnrm/=ZERO ) THEN
!
            eps = DLAMCH('Epsilon')
!
            start = 1
            DO
!
!        while ( START <= N )
!
               IF ( start<=N ) THEN
!
!           Let FINISH be the position of the next subdiagonal entry
!           such that E( FINISH ) <= TINY or FINISH = N if no such
!           subdiagonal exists.  The matrix identified by the elements
!           between START and FINISH constitutes an independent
!           sub-problem.
!
                  finish = start
                  DO
                     IF ( finish<N ) THEN
                        tiny = eps*SQRT(ABS(D(finish)))                 &
     &                         *SQRT(ABS(D(finish+1)))
                        IF ( ABS(E(finish))>tiny ) THEN
                           finish = finish + 1
                           CYCLE
                        ENDIF
                     ENDIF
!
!           (Sub) Problem determined.  Compute its size and solve it.
!
                     m = finish - start + 1
                     IF ( m>smlsiz ) THEN
!
!              Scale.
!
                        orgnrm = DLANST('M',m,D(start),E(start))
                        CALL DLASCL('G',0,0,orgnrm,ONE,m,1,D(start),m,  &
     &                              Info)
                        CALL DLASCL('G',0,0,orgnrm,ONE,m-1,1,E(start),  &
     &                              m-1,Info)
!
                        CALL ZLAED0(N,m,D(start),E(start),Z(1,start),   &
     &                              Ldz,Work,N,Rwork,Iwork,Info)
                        IF ( Info>0 ) THEN
                           Info = (Info/(m+1)+start-1)*(N+1)            &
     &                            + MOD(Info,(m+1)) + start - 1
                           GOTO 100
                        ENDIF
!
!              Scale back.
!
                        CALL DLASCL('G',0,0,ONE,orgnrm,m,1,D(start),m,  &
     &                              Info)
!
                     ELSE
                        CALL DSTEQR('I',m,D(start),E(start),Rwork,m,    &
     &                              Rwork(m*m+1),Info)
                        CALL ZLACRM(N,m,Z(1,start),Ldz,Rwork,m,Work,N,  &
     &                              Rwork(m*m+1))
                        CALL ZLACPY('A',N,m,Work,N,Z(1,start),Ldz)
                        IF ( Info>0 ) THEN
                           Info = start*(N+1) + finish
                           GOTO 100
                        ENDIF
                     ENDIF
!
                     start = finish + 1
                     GOTO 20
                  ENDDO
               ENDIF
!
!        endwhile
!
!
!        Use Selection Sort to minimize swaps of eigenvectors
!
               DO ii = 2 , N
                  i = ii - 1
                  k = i
                  p = D(i)
                  DO j = ii , N
                     IF ( D(j)<p ) THEN
                        k = j
                        p = D(j)
                     ENDIF
                  ENDDO
                  IF ( k/=i ) THEN
                     D(k) = D(i)
                     D(i) = p
                     CALL ZSWAP(N,Z(1,i),1,Z(1,k),1)
                  ENDIF
               ENDDO
               EXIT
 20         ENDDO
         ENDIF
      ENDIF
!
 100  Work(1) = lwmin
      Rwork(1) = lrwmin
      Iwork(1) = liwmin
!
!
!     End of ZSTEDC
!
      END SUBROUTINE ZSTEDC
