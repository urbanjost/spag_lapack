!*==sstedc.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b SSTEDC
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SSTEDC + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sstedc.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sstedc.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sstedc.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SSTEDC( COMPZ, N, D, E, Z, LDZ, WORK, LWORK, IWORK,
!                          LIWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          COMPZ
!       INTEGER            INFO, LDZ, LIWORK, LWORK, N
!       ..
!       .. Array Arguments ..
!       INTEGER            IWORK( * )
!       REAL               D( * ), E( * ), WORK( * ), Z( LDZ, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SSTEDC computes all eigenvalues and, optionally, eigenvectors of a
!> symmetric tridiagonal matrix using the divide and conquer method.
!> The eigenvectors of a full or band real symmetric matrix can also be
!> found if SSYTRD or SSPTRD or SSBTRD has been used to reduce this
!> matrix to tridiagonal form.
!>
!> This code makes very mild assumptions about floating point
!> arithmetic. It will work on machines with a guard digit in
!> add/subtract, or on those binary machines without guard digits
!> which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or Cray-2.
!> It could conceivably fail on hexadecimal or decimal machines
!> without guard digits, but we know of none.  See SLAED3 for details.
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
!>          = 'V':  Compute eigenvectors of original dense symmetric
!>                  matrix also.  On entry, Z contains the orthogonal
!>                  matrix used to reduce the original matrix to
!>                  tridiagonal form.
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
!>          D is REAL array, dimension (N)
!>          On entry, the diagonal elements of the tridiagonal matrix.
!>          On exit, if INFO = 0, the eigenvalues in ascending order.
!> \endverbatim
!>
!> \param[in,out] E
!> \verbatim
!>          E is REAL array, dimension (N-1)
!>          On entry, the subdiagonal elements of the tridiagonal matrix.
!>          On exit, E has been destroyed.
!> \endverbatim
!>
!> \param[in,out] Z
!> \verbatim
!>          Z is REAL array, dimension (LDZ,N)
!>          On entry, if COMPZ = 'V', then Z contains the orthogonal
!>          matrix used in the reduction to tridiagonal form.
!>          On exit, if INFO = 0, then if COMPZ = 'V', Z contains the
!>          orthonormal eigenvectors of the original symmetric matrix,
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
!>          WORK is REAL array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK.
!>          If COMPZ = 'N' or N <= 1 then LWORK must be at least 1.
!>          If COMPZ = 'V' and N > 1 then LWORK must be at least
!>                         ( 1 + 3*N + 2*N*lg N + 4*N**2 ),
!>                         where lg( N ) = smallest integer k such
!>                         that 2**k >= N.
!>          If COMPZ = 'I' and N > 1 then LWORK must be at least
!>                         ( 1 + 4*N + N**2 ).
!>          Note that for COMPZ = 'I' or 'V', then if N is less than or
!>          equal to the minimum divide size, usually 25, then LWORK need
!>          only be max(1,2*(N-1)).
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the WORK array, returns
!>          this value as the first entry of the WORK array, and no error
!>          message related to LWORK is issued by XERBLA.
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
!>          If COMPZ = 'N' or N <= 1 then LIWORK must be at least 1.
!>          If COMPZ = 'V' and N > 1 then LIWORK must be at least
!>                         ( 6 + 6*N + 5*N*lg N ).
!>          If COMPZ = 'I' and N > 1 then LIWORK must be at least
!>                         ( 3 + 5*N ).
!>          Note that for COMPZ = 'I' or 'V', then if N is less than or
!>          equal to the minimum divide size, usually 25, then LIWORK
!>          need only be 1.
!>
!>          If LIWORK = -1, then a workspace query is assumed; the
!>          routine only calculates the optimal size of the IWORK array,
!>          returns this value as the first entry of the IWORK array, and
!>          no error message related to LIWORK is issued by XERBLA.
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
!> \date December 2016
!
!> \ingroup auxOTHERcomputational
!
!> \par Contributors:
!  ==================
!>
!> Jeff Rutter, Computer Science Division, University of California
!> at Berkeley, USA \n
!>  Modified by Francoise Tisseur, University of Tennessee
!>
!  =====================================================================
      SUBROUTINE SSTEDC(Compz,N,D,E,Z,Ldz,Work,Lwork,Iwork,Liwork,Info)
      USE S_ILAENV
      USE S_LSAME
      USE S_SGEMM
      USE S_SLACPY
      USE S_SLAED0
      USE S_SLAMCH
      USE S_SLANST
      USE S_SLASCL
      USE S_SLASET
      USE S_SLASRT
      USE S_SSTEQR
      USE S_SSTERF
      USE S_SSWAP
      USE S_XERBLA
      IMPLICIT NONE
!*--SSTEDC205
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0 , TWO = 2.0E0
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Compz
      INTEGER :: N
      REAL , INTENT(INOUT) , DIMENSION(*) :: D
      REAL , DIMENSION(*) :: E
      REAL , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(IN) :: Liwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL :: eps , orgnrm , p , tiny
      INTEGER :: finish , i , icompz , ii , j , k , lgn , liwmin ,      &
     &           lwmin , m , smlsiz , start , storez , strtrw
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
      lquery = (Lwork==-1 .OR. Liwork==-1)
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
         smlsiz = ILAENV(9,'SSTEDC',' ',0,0,0,0)
         IF ( N<=1 .OR. icompz==0 ) THEN
            liwmin = 1
            lwmin = 1
         ELSEIF ( N<=smlsiz ) THEN
            liwmin = 1
            lwmin = 2*(N-1)
         ELSE
            lgn = INT(LOG(REAL(N))/LOG(TWO))
            IF ( 2**lgn<N ) lgn = lgn + 1
            IF ( 2**lgn<N ) lgn = lgn + 1
            IF ( icompz==1 ) THEN
               lwmin = 1 + 3*N + 2*N*lgn + 4*N**2
               liwmin = 6 + 6*N + 5*N*lgn
            ELSEIF ( icompz==2 ) THEN
               lwmin = 1 + 4*N + N**2
               liwmin = 3 + 5*N
            ENDIF
         ENDIF
         Work(1) = lwmin
         Iwork(1) = liwmin
!
         IF ( Lwork<lwmin .AND. .NOT.lquery ) THEN
            Info = -8
         ELSEIF ( Liwork<liwmin .AND. .NOT.lquery ) THEN
            Info = -10
         ENDIF
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('SSTEDC',-Info)
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
!     Since on many architectures SSTERF is much faster than any other
!     algorithm for finding eigenvalues only, it is used here
!     as the default. If the conditional clause is removed, then
!     information on the size of workspace needs to be changed.
!
!     If COMPZ = 'N', use SSTERF to compute the eigenvalues.
!
      IF ( icompz==0 ) THEN
         CALL SSTERF(N,D,E,Info)
         GOTO 100
      ENDIF
!
!     If N is smaller than the minimum divide size (SMLSIZ+1), then
!     solve the problem with another solver.
!
      IF ( N<=smlsiz ) THEN
!
         CALL SSTEQR(Compz,N,D,E,Z,Ldz,Work,Info)
!
      ELSE
!
!        If COMPZ = 'V', the Z matrix must be stored elsewhere for later
!        use.
!
         IF ( icompz==1 ) THEN
            storez = 1 + N*N
         ELSE
            storez = 1
         ENDIF
!
         IF ( icompz==2 ) CALL SLASET('Full',N,N,ZERO,ONE,Z,Ldz)
!
!        Scale.
!
         orgnrm = SLANST('M',N,D,E)
         IF ( orgnrm/=ZERO ) THEN
!
            eps = SLAMCH('Epsilon')
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
                     IF ( m==1 ) THEN
                        start = finish + 1
                        GOTO 20
                     ENDIF
                     IF ( m>smlsiz ) THEN
!
!              Scale.
!
                        orgnrm = SLANST('M',m,D(start),E(start))
                        CALL SLASCL('G',0,0,orgnrm,ONE,m,1,D(start),m,  &
     &                              Info)
                        CALL SLASCL('G',0,0,orgnrm,ONE,m-1,1,E(start),  &
     &                              m-1,Info)
!
                        IF ( icompz==1 ) THEN
                           strtrw = 1
                        ELSE
                           strtrw = start
                        ENDIF
                        CALL SLAED0(icompz,N,m,D(start),E(start),       &
     &                              Z(strtrw,start),Ldz,Work(1),N,      &
     &                              Work(storez),Iwork,Info)
                        IF ( Info/=0 ) THEN
                           Info = (Info/(m+1)+start-1)*(N+1)            &
     &                            + MOD(Info,(m+1)) + start - 1
                           GOTO 100
                        ENDIF
!
!              Scale back.
!
                        CALL SLASCL('G',0,0,ONE,orgnrm,m,1,D(start),m,  &
     &                              Info)
!
                     ELSE
                        IF ( icompz==1 ) THEN
!
!                 Since QR won't update a Z matrix which is larger than
!                 the length of D, we must solve the sub-problem in a
!                 workspace and then multiply back into Z.
!
                           CALL SSTEQR('I',m,D(start),E(start),Work,m,  &
     &                                 Work(m*m+1),Info)
                           CALL SLACPY('A',N,m,Z(1,start),Ldz,          &
     &                                 Work(storez),N)
                           CALL SGEMM('N','N',N,m,m,ONE,Work(storez),N, &
     &                                Work,m,ZERO,Z(1,start),Ldz)
                        ELSEIF ( icompz==2 ) THEN
                           CALL SSTEQR('I',m,D(start),E(start),         &
     &                                 Z(start,start),Ldz,Work,Info)
                        ELSE
                           CALL SSTERF(m,D(start),E(start),Info)
                        ENDIF
                        IF ( Info/=0 ) THEN
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
               IF ( icompz==0 ) THEN
!
!          Use Quick Sort
!
                  CALL SLASRT('I',N,D,Info)
!
               ELSE
!
!          Use Selection Sort to minimize swaps of eigenvectors
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
                        CALL SSWAP(N,Z(1,i),1,Z(1,k),1)
                     ENDIF
                  ENDDO
               ENDIF
               EXIT
 20         ENDDO
         ENDIF
      ENDIF
!
 100  Work(1) = lwmin
      Iwork(1) = liwmin
!
!
!     End of SSTEDC
!
      END SUBROUTINE SSTEDC
