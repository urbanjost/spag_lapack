!*==dggbal.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b DGGBAL
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DGGBAL + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dggbal.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dggbal.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dggbal.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGGBAL( JOB, N, A, LDA, B, LDB, ILO, IHI, LSCALE,
!                          RSCALE, WORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          JOB
!       INTEGER            IHI, ILO, INFO, LDA, LDB, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), LSCALE( * ),
!      $                   RSCALE( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DGGBAL balances a pair of general real matrices (A,B).  This
!> involves, first, permuting A and B by similarity transformations to
!> isolate eigenvalues in the first 1 to ILO$-$1 and last IHI+1 to N
!> elements on the diagonal; and second, applying a diagonal similarity
!> transformation to rows and columns ILO to IHI to make the rows
!> and columns as close in norm as possible. Both steps are optional.
!>
!> Balancing may reduce the 1-norm of the matrices, and improve the
!> accuracy of the computed eigenvalues and/or eigenvectors in the
!> generalized eigenvalue problem A*x = lambda*B*x.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] JOB
!> \verbatim
!>          JOB is CHARACTER*1
!>          Specifies the operations to be performed on A and B:
!>          = 'N':  none:  simply set ILO = 1, IHI = N, LSCALE(I) = 1.0
!>                  and RSCALE(I) = 1.0 for i = 1,...,N.
!>          = 'P':  permute only;
!>          = 'S':  scale only;
!>          = 'B':  both permute and scale.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrices A and B.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          On entry, the input matrix A.
!>          On exit,  A is overwritten by the balanced matrix.
!>          If JOB = 'N', A is not referenced.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A. LDA >= max(1,N).
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is DOUBLE PRECISION array, dimension (LDB,N)
!>          On entry, the input matrix B.
!>          On exit,  B is overwritten by the balanced matrix.
!>          If JOB = 'N', B is not referenced.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B. LDB >= max(1,N).
!> \endverbatim
!>
!> \param[out] ILO
!> \verbatim
!>          ILO is INTEGER
!> \endverbatim
!>
!> \param[out] IHI
!> \verbatim
!>          IHI is INTEGER
!>          ILO and IHI are set to integers such that on exit
!>          A(i,j) = 0 and B(i,j) = 0 if i > j and
!>          j = 1,...,ILO-1 or i = IHI+1,...,N.
!>          If JOB = 'N' or 'S', ILO = 1 and IHI = N.
!> \endverbatim
!>
!> \param[out] LSCALE
!> \verbatim
!>          LSCALE is DOUBLE PRECISION array, dimension (N)
!>          Details of the permutations and scaling factors applied
!>          to the left side of A and B.  If P(j) is the index of the
!>          row interchanged with row j, and D(j)
!>          is the scaling factor applied to row j, then
!>            LSCALE(j) = P(j)    for J = 1,...,ILO-1
!>                      = D(j)    for J = ILO,...,IHI
!>                      = P(j)    for J = IHI+1,...,N.
!>          The order in which the interchanges are made is N to IHI+1,
!>          then 1 to ILO-1.
!> \endverbatim
!>
!> \param[out] RSCALE
!> \verbatim
!>          RSCALE is DOUBLE PRECISION array, dimension (N)
!>          Details of the permutations and scaling factors applied
!>          to the right side of A and B.  If P(j) is the index of the
!>          column interchanged with column j, and D(j)
!>          is the scaling factor applied to column j, then
!>            LSCALE(j) = P(j)    for J = 1,...,ILO-1
!>                      = D(j)    for J = ILO,...,IHI
!>                      = P(j)    for J = IHI+1,...,N.
!>          The order in which the interchanges are made is N to IHI+1,
!>          then 1 to ILO-1.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (lwork)
!>          lwork must be at least max(1,6*N) when JOB = 'S' or 'B', and
!>          at least 1 when JOB = 'N' or 'P'.
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
!> \ingroup doubleGBcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  See R.C. WARD, Balancing the generalized eigenvalue problem,
!>                 SIAM J. Sci. Stat. Comp. 2 (1981), 141-152.
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE DGGBAL(Job,N,A,Lda,B,Ldb,Ilo,Ihi,Lscale,Rscale,Work,   &
     &                  Info)
      IMPLICIT NONE
!*--DGGBAL181
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Job
      INTEGER Ihi , Ilo , Info , Lda , Ldb , N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION A(Lda,*) , B(Ldb,*) , Lscale(*) , Rscale(*) ,    &
     &                 Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , HALF , ONE
      PARAMETER (ZERO=0.0D+0,HALF=0.5D+0,ONE=1.0D+0)
      DOUBLE PRECISION THREE , SCLFAC
      PARAMETER (THREE=3.0D+0,SCLFAC=1.0D+1)
!     ..
!     .. Local Scalars ..
      INTEGER i , icab , iflow , ip1 , ir , irab , it , j , jc , jp1 ,  &
     &        k , kount , l , lcab , lm1 , lrab , lsfmax , lsfmin , m , &
     &        nr , nrp2
      DOUBLE PRECISION alpha , basl , beta , cab , cmax , coef , coef2 ,&
     &                 coef5 , cor , ew , ewc , gamma , pgamma , rab ,  &
     &                 sfmax , sfmin , sum , t , ta , tb , tc
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      INTEGER IDAMAX
      DOUBLE PRECISION DDOT , DLAMCH
      EXTERNAL LSAME , IDAMAX , DDOT , DLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL DAXPY , DSCAL , DSWAP , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , DBLE , INT , LOG10 , MAX , MIN , SIGN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters
!
      Info = 0
      IF ( .NOT.LSAME(Job,'N') .AND. .NOT.LSAME(Job,'P') .AND.          &
     &     .NOT.LSAME(Job,'S') .AND. .NOT.LSAME(Job,'B') ) THEN
         Info = -1
      ELSEIF ( N<0 ) THEN
         Info = -2
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -4
      ELSEIF ( Ldb<MAX(1,N) ) THEN
         Info = -6
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('DGGBAL',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 ) THEN
         Ilo = 1
         Ihi = N
         RETURN
      ENDIF
!
      IF ( N==1 ) THEN
         Ilo = 1
         Ihi = N
         Lscale(1) = ONE
         Rscale(1) = ONE
         RETURN
      ENDIF
!
      IF ( LSAME(Job,'N') ) THEN
         Ilo = 1
         Ihi = N
         DO i = 1 , N
            Lscale(i) = ONE
            Rscale(i) = ONE
         ENDDO
         RETURN
      ENDIF
!
      k = 1
      l = N
!
      IF ( LSAME(Job,'S') ) GOTO 800
!
 100  lm1 = l - 1
      DO i = l , 1 , -1
         DO j = 1 , lm1
            jp1 = j + 1
            IF ( A(i,j)/=ZERO .OR. B(i,j)/=ZERO ) GOTO 150
         ENDDO
         j = l
         GOTO 200
!
 150     DO j = jp1 , l
            IF ( A(i,j)/=ZERO .OR. B(i,j)/=ZERO ) GOTO 300
         ENDDO
         j = jp1 - 1
!
 200     m = l
         iflow = 1
         GOTO 700
 300  ENDDO
!
 400  DO j = k , l
         DO i = k , lm1
            ip1 = i + 1
            IF ( A(i,j)/=ZERO .OR. B(i,j)/=ZERO ) GOTO 450
         ENDDO
         i = l
         GOTO 500
 450     DO i = ip1 , l
            IF ( A(i,j)/=ZERO .OR. B(i,j)/=ZERO ) GOTO 600
         ENDDO
         i = ip1 - 1
 500     m = k
         iflow = 2
         GOTO 700
 600  ENDDO
      GOTO 800
!
!     Permute rows M and I
!
 700  Lscale(m) = i
      IF ( i/=m ) THEN
         CALL DSWAP(N-k+1,A(i,k),Lda,A(m,k),Lda)
         CALL DSWAP(N-k+1,B(i,k),Ldb,B(m,k),Ldb)
      ENDIF
!
!     Permute columns M and J
!
      Rscale(m) = j
      IF ( j/=m ) THEN
         CALL DSWAP(l,A(1,j),1,A(1,m),1)
         CALL DSWAP(l,B(1,j),1,B(1,m),1)
      ENDIF
!
      IF ( iflow==1 ) THEN
!
!     Permute the matrices A and B to isolate the eigenvalues.
!
!     Find row with one nonzero in columns 1 through L
!
         l = lm1
         IF ( l/=1 ) GOTO 100
!
         Rscale(1) = ONE
         Lscale(1) = ONE
      ELSEIF ( iflow==2 ) THEN
!
!     Find column with one nonzero in rows K through N
!
         k = k + 1
         GOTO 400
      ENDIF
!
 800  Ilo = k
      Ihi = l
!
      IF ( LSAME(Job,'P') ) THEN
         DO i = Ilo , Ihi
            Lscale(i) = ONE
            Rscale(i) = ONE
         ENDDO
         RETURN
      ENDIF
!
      IF ( Ilo==Ihi ) RETURN
!
!     Balance the submatrix in rows ILO to IHI.
!
      nr = Ihi - Ilo + 1
      DO i = Ilo , Ihi
         Rscale(i) = ZERO
         Lscale(i) = ZERO
!
         Work(i) = ZERO
         Work(i+N) = ZERO
         Work(i+2*N) = ZERO
         Work(i+3*N) = ZERO
         Work(i+4*N) = ZERO
         Work(i+5*N) = ZERO
      ENDDO
!
!     Compute right side vector in resulting linear equations
!
      basl = LOG10(SCLFAC)
      DO i = Ilo , Ihi
         DO j = Ilo , Ihi
            tb = B(i,j)
            ta = A(i,j)
            IF ( ta/=ZERO ) ta = LOG10(ABS(ta))/basl
            IF ( tb/=ZERO ) tb = LOG10(ABS(tb))/basl
            Work(i+4*N) = Work(i+4*N) - ta - tb
            Work(j+5*N) = Work(j+5*N) - ta - tb
         ENDDO
      ENDDO
!
      coef = ONE/DBLE(2*nr)
      coef2 = coef*coef
      coef5 = HALF*coef2
      nrp2 = nr + 2
      beta = ZERO
      it = 1
      DO
!
!     Start generalized conjugate gradient iteration
!
!
         gamma = DDOT(nr,Work(Ilo+4*N),1,Work(Ilo+4*N),1)               &
     &           + DDOT(nr,Work(Ilo+5*N),1,Work(Ilo+5*N),1)
!
         ew = ZERO
         ewc = ZERO
         DO i = Ilo , Ihi
            ew = ew + Work(i+4*N)
            ewc = ewc + Work(i+5*N)
         ENDDO
!
         gamma = coef*gamma - coef2*(ew**2+ewc**2) - coef5*(ew-ewc)**2
         IF ( gamma==ZERO ) EXIT
         IF ( it/=1 ) beta = gamma/pgamma
         t = coef5*(ewc-THREE*ew)
         tc = coef5*(ew-THREE*ewc)
!
         CALL DSCAL(nr,beta,Work(Ilo),1)
         CALL DSCAL(nr,beta,Work(Ilo+N),1)
!
         CALL DAXPY(nr,coef,Work(Ilo+4*N),1,Work(Ilo+N),1)
         CALL DAXPY(nr,coef,Work(Ilo+5*N),1,Work(Ilo),1)
!
         DO i = Ilo , Ihi
            Work(i) = Work(i) + tc
            Work(i+N) = Work(i+N) + t
         ENDDO
!
!     Apply matrix to vector
!
         DO i = Ilo , Ihi
            kount = 0
            sum = ZERO
            DO j = Ilo , Ihi
               IF ( A(i,j)/=ZERO ) THEN
                  kount = kount + 1
                  sum = sum + Work(j)
               ENDIF
               IF ( B(i,j)/=ZERO ) THEN
                  kount = kount + 1
                  sum = sum + Work(j)
               ENDIF
            ENDDO
            Work(i+2*N) = DBLE(kount)*Work(i+N) + sum
         ENDDO
!
         DO j = Ilo , Ihi
            kount = 0
            sum = ZERO
            DO i = Ilo , Ihi
               IF ( A(i,j)/=ZERO ) THEN
                  kount = kount + 1
                  sum = sum + Work(i+N)
               ENDIF
               IF ( B(i,j)/=ZERO ) THEN
                  kount = kount + 1
                  sum = sum + Work(i+N)
               ENDIF
            ENDDO
            Work(j+3*N) = DBLE(kount)*Work(j) + sum
         ENDDO
!
         sum = DDOT(nr,Work(Ilo+N),1,Work(Ilo+2*N),1)                   &
     &         + DDOT(nr,Work(Ilo),1,Work(Ilo+3*N),1)
         alpha = gamma/sum
!
!     Determine correction to current iteration
!
         cmax = ZERO
         DO i = Ilo , Ihi
            cor = alpha*Work(i+N)
            IF ( ABS(cor)>cmax ) cmax = ABS(cor)
            Lscale(i) = Lscale(i) + cor
            cor = alpha*Work(i)
            IF ( ABS(cor)>cmax ) cmax = ABS(cor)
            Rscale(i) = Rscale(i) + cor
         ENDDO
         IF ( cmax<HALF ) EXIT
!
         CALL DAXPY(nr,-alpha,Work(Ilo+2*N),1,Work(Ilo+4*N),1)
         CALL DAXPY(nr,-alpha,Work(Ilo+3*N),1,Work(Ilo+5*N),1)
!
         pgamma = gamma
         it = it + 1
         IF ( it>nrp2 ) EXIT
      ENDDO
!
!     End generalized conjugate gradient iteration
!
      sfmin = DLAMCH('S')
      sfmax = ONE/sfmin
      lsfmin = INT(LOG10(sfmin)/basl+ONE)
      lsfmax = INT(LOG10(sfmax)/basl)
      DO i = Ilo , Ihi
         irab = IDAMAX(N-Ilo+1,A(i,Ilo),Lda)
         rab = ABS(A(i,irab+Ilo-1))
         irab = IDAMAX(N-Ilo+1,B(i,Ilo),Ldb)
         rab = MAX(rab,ABS(B(i,irab+Ilo-1)))
         lrab = INT(LOG10(rab+sfmin)/basl+ONE)
         ir = INT(Lscale(i)+SIGN(HALF,Lscale(i)))
         ir = MIN(MAX(ir,lsfmin),lsfmax,lsfmax-lrab)
         Lscale(i) = SCLFAC**ir
         icab = IDAMAX(Ihi,A(1,i),1)
         cab = ABS(A(icab,i))
         icab = IDAMAX(Ihi,B(1,i),1)
         cab = MAX(cab,ABS(B(icab,i)))
         lcab = INT(LOG10(cab+sfmin)/basl+ONE)
         jc = INT(Rscale(i)+SIGN(HALF,Rscale(i)))
         jc = MIN(MAX(jc,lsfmin),lsfmax,lsfmax-lcab)
         Rscale(i) = SCLFAC**jc
      ENDDO
!
!     Row scaling of matrices A and B
!
      DO i = Ilo , Ihi
         CALL DSCAL(N-Ilo+1,Lscale(i),A(i,Ilo),Lda)
         CALL DSCAL(N-Ilo+1,Lscale(i),B(i,Ilo),Ldb)
      ENDDO
!
!     Column scaling of matrices A and B
!
      DO j = Ilo , Ihi
         CALL DSCAL(Ihi,Rscale(j),A(1,j),1)
         CALL DSCAL(Ihi,Rscale(j),B(1,j),1)
      ENDDO
!
!
!     End of DGGBAL
!
      END SUBROUTINE DGGBAL
