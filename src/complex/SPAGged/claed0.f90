!*==claed0.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b CLAED0 used by sstedc. Computes all eigenvalues and corresponding eigenvectors of an unreduced symmetric tridiagonal matrix using the divide and conquer method.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CLAED0 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/claed0.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/claed0.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/claed0.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CLAED0( QSIZ, N, D, E, Q, LDQ, QSTORE, LDQS, RWORK,
!                          IWORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDQ, LDQS, N, QSIZ
!       ..
!       .. Array Arguments ..
!       INTEGER            IWORK( * )
!       REAL               D( * ), E( * ), RWORK( * )
!       COMPLEX            Q( LDQ, * ), QSTORE( LDQS, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> Using the divide and conquer method, CLAED0 computes all eigenvalues
!> of a symmetric tridiagonal matrix which is one diagonal block of
!> those from reducing a dense or band Hermitian matrix and
!> corresponding eigenvectors of the dense or band matrix.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] QSIZ
!> \verbatim
!>          QSIZ is INTEGER
!>         The dimension of the unitary matrix used to reduce
!>         the full matrix to tridiagonal form.  QSIZ >= N if ICOMPQ = 1.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>         The dimension of the symmetric tridiagonal matrix.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] D
!> \verbatim
!>          D is REAL array, dimension (N)
!>         On entry, the diagonal elements of the tridiagonal matrix.
!>         On exit, the eigenvalues in ascending order.
!> \endverbatim
!>
!> \param[in,out] E
!> \verbatim
!>          E is REAL array, dimension (N-1)
!>         On entry, the off-diagonal elements of the tridiagonal matrix.
!>         On exit, E has been destroyed.
!> \endverbatim
!>
!> \param[in,out] Q
!> \verbatim
!>          Q is COMPLEX array, dimension (LDQ,N)
!>         On entry, Q must contain an QSIZ x N matrix whose columns
!>         unitarily orthonormal. It is a part of the unitary matrix
!>         that reduces the full dense Hermitian matrix to a
!>         (reducible) symmetric tridiagonal matrix.
!> \endverbatim
!>
!> \param[in] LDQ
!> \verbatim
!>          LDQ is INTEGER
!>         The leading dimension of the array Q.  LDQ >= max(1,N).
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array,
!>         the dimension of IWORK must be at least
!>                      6 + 6*N + 5*N*lg N
!>                      ( lg( N ) = smallest integer k
!>                                  such that 2^k >= N )
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array,
!>                               dimension (1 + 3*N + 2*N*lg N + 3*N**2)
!>                        ( lg( N ) = smallest integer k
!>                                    such that 2^k >= N )
!> \endverbatim
!>
!> \param[out] QSTORE
!> \verbatim
!>          QSTORE is COMPLEX array, dimension (LDQS, N)
!>         Used to store parts of
!>         the eigenvector matrix when the updating matrix multiplies
!>         take place.
!> \endverbatim
!>
!> \param[in] LDQS
!> \verbatim
!>          LDQS is INTEGER
!>         The leading dimension of the array QSTORE.
!>         LDQS >= max(1,N).
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
!> \ingroup complexOTHERcomputational
!
!  =====================================================================
      SUBROUTINE CLAED0(Qsiz,N,D,E,Q,Ldq,Qstore,Ldqs,Rwork,Iwork,Info)
      USE S_CCOPY
      USE S_CLACRM
      USE S_CLAED7
      USE S_ILAENV
      USE S_SCOPY
      USE S_SSTEQR
      USE S_XERBLA
      IMPLICIT NONE
!*--CLAED0155
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  TWO = 2.E+0
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER :: Qsiz
      INTEGER :: N
      REAL , INTENT(INOUT) , DIMENSION(*) :: D
      REAL , DIMENSION(*) :: E
      COMPLEX , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      COMPLEX , DIMENSION(Ldqs,*) :: Qstore
      INTEGER :: Ldqs
      REAL , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: curlvl , curprb , curr , i , igivcl , igivnm , igivpt ,&
     &           indxq , iperm , iprmpt , iq , iqptr , iwrem , j , k ,  &
     &           lgn , ll , matsiz , msd2 , smlsiz , smm1 , spm1 ,      &
     &           spm2 , submat , subpbs , tlvls
      REAL :: temp
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
!  =====================================================================
!
!  Warning:      N could be as big as QSIZ!
!
!     .. Parameters ..
!     ..
!     .. Local Scalars ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. External Functions ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      Info = 0
!
!     IF( ICOMPQ .LT. 0 .OR. ICOMPQ .GT. 2 ) THEN
!        INFO = -1
!     ELSE IF( ( ICOMPQ .EQ. 1 ) .AND. ( QSIZ .LT. MAX( 0, N ) ) )
!    $        THEN
      IF ( Qsiz<MAX(0,N) ) THEN
         Info = -1
      ELSEIF ( N<0 ) THEN
         Info = -2
      ELSEIF ( Ldq<MAX(1,N) ) THEN
         Info = -6
      ELSEIF ( Ldqs<MAX(1,N) ) THEN
         Info = -8
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('CLAED0',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 ) RETURN
!
      smlsiz = ILAENV(9,'CLAED0',' ',0,0,0,0)
!
!     Determine the size and placement of the submatrices, and save in
!     the leading elements of IWORK.
!
      Iwork(1) = N
      subpbs = 1
      tlvls = 0
      DO
         IF ( Iwork(subpbs)>smlsiz ) THEN
            DO j = subpbs , 1 , -1
               Iwork(2*j) = (Iwork(j)+1)/2
               Iwork(2*j-1) = Iwork(j)/2
            ENDDO
            tlvls = tlvls + 1
            subpbs = 2*subpbs
            CYCLE
         ENDIF
         DO j = 2 , subpbs
            Iwork(j) = Iwork(j) + Iwork(j-1)
         ENDDO
!
!     Divide the matrix into SUBPBS submatrices of size at most SMLSIZ+1
!     using rank-1 modifications (cuts).
!
         spm1 = subpbs - 1
         DO i = 1 , spm1
            submat = Iwork(i) + 1
            smm1 = submat - 1
            D(smm1) = D(smm1) - ABS(E(smm1))
            D(submat) = D(submat) - ABS(E(smm1))
         ENDDO
!
         indxq = 4*N + 3
!
!     Set up workspaces for eigenvalues only/accumulate new vectors
!     routine
!
         temp = LOG(REAL(N))/LOG(TWO)
         lgn = INT(temp)
         IF ( 2**lgn<N ) lgn = lgn + 1
         IF ( 2**lgn<N ) lgn = lgn + 1
         iprmpt = indxq + N + 1
         iperm = iprmpt + N*lgn
         iqptr = iperm + N*lgn
         igivpt = iqptr + N + 2
         igivcl = igivpt + N*lgn
!
         igivnm = 1
         iq = igivnm + 2*N*lgn
         iwrem = iq + N**2 + 1
!     Initialize pointers
         DO i = 0 , subpbs
            Iwork(iprmpt+i) = 1
            Iwork(igivpt+i) = 1
         ENDDO
         Iwork(iqptr) = 1
!
!     Solve each submatrix eigenproblem at the bottom of the divide and
!     conquer tree.
!
         curr = 0
         DO i = 0 , spm1
            IF ( i==0 ) THEN
               submat = 1
               matsiz = Iwork(1)
            ELSE
               submat = Iwork(i) + 1
               matsiz = Iwork(i+1) - Iwork(i)
            ENDIF
            ll = iq - 1 + Iwork(iqptr+curr)
            CALL SSTEQR('I',matsiz,D(submat),E(submat),Rwork(ll),matsiz,&
     &                  Rwork,Info)
            CALL CLACRM(Qsiz,matsiz,Q(1,submat),Ldq,Rwork(ll),matsiz,   &
     &                  Qstore(1,submat),Ldqs,Rwork(iwrem))
            Iwork(iqptr+curr+1) = Iwork(iqptr+curr) + matsiz**2
            curr = curr + 1
            IF ( Info>0 ) THEN
               Info = submat*(N+1) + submat + matsiz - 1
               RETURN
            ENDIF
            k = 1
            DO j = submat , Iwork(i+1)
               Iwork(indxq+j) = k
               k = k + 1
            ENDDO
         ENDDO
!
!     Successively merge eigensystems of adjacent submatrices
!     into eigensystem for the corresponding larger matrix.
!
!     while ( SUBPBS > 1 )
!
         curlvl = 1
         EXIT
      ENDDO
      DO
         IF ( subpbs>1 ) THEN
            spm2 = subpbs - 2
            DO i = 0 , spm2 , 2
               IF ( i==0 ) THEN
                  submat = 1
                  matsiz = Iwork(2)
                  msd2 = Iwork(1)
                  curprb = 0
               ELSE
                  submat = Iwork(i) + 1
                  matsiz = Iwork(i+2) - Iwork(i)
                  msd2 = matsiz/2
                  curprb = curprb + 1
               ENDIF
!
!     Merge lower order eigensystems (of size MSD2 and MATSIZ - MSD2)
!     into an eigensystem of size MATSIZ.  CLAED7 handles the case
!     when the eigenvectors of a full or band Hermitian matrix (which
!     was reduced to tridiagonal form) are desired.
!
!     I am free to use Q as a valuable working space until Loop 150.
!
               CALL CLAED7(matsiz,msd2,Qsiz,tlvls,curlvl,curprb,        &
     &                     D(submat),Qstore(1,submat),Ldqs,             &
     &                     E(submat+msd2-1),Iwork(indxq+submat),        &
     &                     Rwork(iq),Iwork(iqptr),Iwork(iprmpt),        &
     &                     Iwork(iperm),Iwork(igivpt),Iwork(igivcl),    &
     &                     Rwork(igivnm),Q(1,submat),Rwork(iwrem),      &
     &                     Iwork(subpbs+1),Info)
               IF ( Info>0 ) THEN
                  Info = submat*(N+1) + submat + matsiz - 1
                  RETURN
               ENDIF
               Iwork(i/2+1) = Iwork(i+2)
            ENDDO
            subpbs = subpbs/2
            curlvl = curlvl + 1
            CYCLE
         ENDIF
!
!     end while
!
!     Re-merge the eigenvalues/vectors which were deflated at the final
!     merge step.
!
         DO i = 1 , N
            j = Iwork(indxq+i)
            Rwork(i) = D(j)
            CALL CCOPY(Qsiz,Qstore(1,j),1,Q(1,i),1)
         ENDDO
         CALL SCOPY(N,Rwork,1,D,1)
         EXIT
      ENDDO
!
!
!     End of CLAED0
!
      END SUBROUTINE CLAED0
