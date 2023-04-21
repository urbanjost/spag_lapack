!*==slaed0.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b SLAED0 used by sstedc. Computes all eigenvalues and corresponding eigenvectors of an unreduced symmetric tridiagonal matrix using the divide and conquer method.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SLAED0 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slaed0.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slaed0.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slaed0.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLAED0( ICOMPQ, QSIZ, N, D, E, Q, LDQ, QSTORE, LDQS,
!                          WORK, IWORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            ICOMPQ, INFO, LDQ, LDQS, N, QSIZ
!       ..
!       .. Array Arguments ..
!       INTEGER            IWORK( * )
!       REAL               D( * ), E( * ), Q( LDQ, * ), QSTORE( LDQS, * ),
!      $                   WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLAED0 computes all eigenvalues and corresponding eigenvectors of a
!> symmetric tridiagonal matrix using the divide and conquer method.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] ICOMPQ
!> \verbatim
!>          ICOMPQ is INTEGER
!>          = 0:  Compute eigenvalues only.
!>          = 1:  Compute eigenvectors of original dense symmetric matrix
!>                also.  On entry, Q contains the orthogonal matrix used
!>                to reduce the original matrix to tridiagonal form.
!>          = 2:  Compute eigenvalues and eigenvectors of tridiagonal
!>                matrix.
!> \endverbatim
!>
!> \param[in] QSIZ
!> \verbatim
!>          QSIZ is INTEGER
!>         The dimension of the orthogonal matrix used to reduce
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
!>         On entry, the main diagonal of the tridiagonal matrix.
!>         On exit, its eigenvalues.
!> \endverbatim
!>
!> \param[in] E
!> \verbatim
!>          E is REAL array, dimension (N-1)
!>         The off-diagonal elements of the tridiagonal matrix.
!>         On exit, E has been destroyed.
!> \endverbatim
!>
!> \param[in,out] Q
!> \verbatim
!>          Q is REAL array, dimension (LDQ, N)
!>         On entry, Q must contain an N-by-N orthogonal matrix.
!>         If ICOMPQ = 0    Q is not referenced.
!>         If ICOMPQ = 1    On entry, Q is a subset of the columns of the
!>                          orthogonal matrix used to reduce the full
!>                          matrix to tridiagonal form corresponding to
!>                          the subset of the full matrix which is being
!>                          decomposed at this time.
!>         If ICOMPQ = 2    On entry, Q will be the identity matrix.
!>                          On exit, Q contains the eigenvectors of the
!>                          tridiagonal matrix.
!> \endverbatim
!>
!> \param[in] LDQ
!> \verbatim
!>          LDQ is INTEGER
!>         The leading dimension of the array Q.  If eigenvectors are
!>         desired, then  LDQ >= max(1,N).  In any case,  LDQ >= 1.
!> \endverbatim
!>
!> \param[out] QSTORE
!> \verbatim
!>          QSTORE is REAL array, dimension (LDQS, N)
!>         Referenced only when ICOMPQ = 1.  Used to store parts of
!>         the eigenvector matrix when the updating matrix multiplies
!>         take place.
!> \endverbatim
!>
!> \param[in] LDQS
!> \verbatim
!>          LDQS is INTEGER
!>         The leading dimension of the array QSTORE.  If ICOMPQ = 1,
!>         then  LDQS >= max(1,N).  In any case,  LDQS >= 1.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array,
!>         If ICOMPQ = 0 or 1, the dimension of WORK must be at least
!>                     1 + 3*N + 2*N*lg N + 3*N**2
!>                     ( lg( N ) = smallest integer k
!>                                 such that 2^k >= N )
!>         If ICOMPQ = 2, the dimension of WORK must be at least
!>                     4*N + N**2.
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array,
!>         If ICOMPQ = 0 or 1, the dimension of IWORK must be at least
!>                        6 + 6*N + 5*N*lg N.
!>                        ( lg( N ) = smallest integer k
!>                                    such that 2^k >= N )
!>         If ICOMPQ = 2, the dimension of IWORK must be at least
!>                        3 + 5*N.
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
!> at Berkeley, USA
!
!  =====================================================================
      SUBROUTINE SLAED0(Icompq,Qsiz,N,D,E,Q,Ldq,Qstore,Ldqs,Work,Iwork, &
     &                  Info)
      IMPLICIT NONE
!*--SLAED0176
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Icompq , Info , Ldq , Ldqs , N , Qsiz
!     ..
!     .. Array Arguments ..
      INTEGER Iwork(*)
      REAL D(*) , E(*) , Q(Ldq,*) , Qstore(Ldqs,*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE , TWO
      PARAMETER (ZERO=0.E0,ONE=1.E0,TWO=2.E0)
!     ..
!     .. Local Scalars ..
      INTEGER curlvl , curprb , curr , i , igivcl , igivnm , igivpt ,   &
     &        indxq , iperm , iprmpt , iq , iqptr , iwrem , j , k ,     &
     &        lgn , matsiz , msd2 , smlsiz , smm1 , spm1 , spm2 ,       &
     &        submat , subpbs , tlvls
      REAL temp
!     ..
!     .. External Subroutines ..
      EXTERNAL SCOPY , SGEMM , SLACPY , SLAED1 , SLAED7 , SSTEQR ,      &
     &         XERBLA
!     ..
!     .. External Functions ..
      INTEGER ILAENV
      EXTERNAL ILAENV
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , INT , LOG , MAX , REAL
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      Info = 0
!
      IF ( Icompq<0 .OR. Icompq>2 ) THEN
         Info = -1
      ELSEIF ( (Icompq==1) .AND. (Qsiz<MAX(0,N)) ) THEN
         Info = -2
      ELSEIF ( N<0 ) THEN
         Info = -3
      ELSEIF ( Ldq<MAX(1,N) ) THEN
         Info = -7
      ELSEIF ( Ldqs<MAX(1,N) ) THEN
         Info = -9
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('SLAED0',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 ) RETURN
!
      smlsiz = ILAENV(9,'SLAED0',' ',0,0,0,0)
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
         IF ( Icompq/=2 ) THEN
!
!        Set up workspaces for eigenvalues only/accumulate new vectors
!        routine
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
!
!        Initialize pointers
!
            DO i = 0 , subpbs
               Iwork(iprmpt+i) = 1
               Iwork(igivpt+i) = 1
            ENDDO
            Iwork(iqptr) = 1
         ENDIF
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
            IF ( Icompq==2 ) THEN
               CALL SSTEQR('I',matsiz,D(submat),E(submat),              &
     &                     Q(submat,submat),Ldq,Work,Info)
               IF ( Info/=0 ) GOTO 100
            ELSE
               CALL SSTEQR('I',matsiz,D(submat),E(submat),              &
     &                     Work(iq-1+Iwork(iqptr+curr)),matsiz,Work,    &
     &                     Info)
               IF ( Info/=0 ) GOTO 100
               IF ( Icompq==1 ) CALL SGEMM('N','N',Qsiz,matsiz,matsiz,  &
     &              ONE,Q(1,submat),Ldq,Work(iq-1+Iwork(iqptr+curr)),   &
     &              matsiz,ZERO,Qstore(1,submat),Ldqs)
               Iwork(iqptr+curr+1) = Iwork(iqptr+curr) + matsiz**2
               curr = curr + 1
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
!     into an eigensystem of size MATSIZ.
!     SLAED1 is used only for the full eigensystem of a tridiagonal
!     matrix.
!     SLAED7 handles the cases in which eigenvalues only or eigenvalues
!     and eigenvectors of a full symmetric matrix (which was reduced to
!     tridiagonal form) are desired.
!
               IF ( Icompq==2 ) THEN
                  CALL SLAED1(matsiz,D(submat),Q(submat,submat),Ldq,    &
     &                        Iwork(indxq+submat),E(submat+msd2-1),msd2,&
     &                        Work,Iwork(subpbs+1),Info)
               ELSE
                  CALL SLAED7(Icompq,matsiz,Qsiz,tlvls,curlvl,curprb,   &
     &                        D(submat),Qstore(1,submat),Ldqs,          &
     &                        Iwork(indxq+submat),E(submat+msd2-1),msd2,&
     &                        Work(iq),Iwork(iqptr),Iwork(iprmpt),      &
     &                        Iwork(iperm),Iwork(igivpt),Iwork(igivcl), &
     &                        Work(igivnm),Work(iwrem),Iwork(subpbs+1), &
     &                        Info)
               ENDIF
               IF ( Info/=0 ) GOTO 100
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
         IF ( Icompq==1 ) THEN
            DO i = 1 , N
               j = Iwork(indxq+i)
               Work(i) = D(j)
               CALL SCOPY(Qsiz,Qstore(1,j),1,Q(1,i),1)
            ENDDO
            CALL SCOPY(N,Work,1,D,1)
         ELSEIF ( Icompq==2 ) THEN
            DO i = 1 , N
               j = Iwork(indxq+i)
               Work(i) = D(j)
               CALL SCOPY(N,Q(1,j),1,Work(N*i+1),1)
            ENDDO
            CALL SCOPY(N,Work,1,D,1)
            CALL SLACPY('A',N,N,Work(N+1),N,Q,Ldq)
         ELSE
            DO i = 1 , N
               j = Iwork(indxq+i)
               Work(i) = D(j)
            ENDDO
            CALL SCOPY(N,Work,1,D,1)
         ENDIF
         GOTO 99999
      ENDDO
!
 100  Info = submat*(N+1) + submat + matsiz - 1
!
!
!     End of SLAED0
!
99999 END SUBROUTINE SLAED0
