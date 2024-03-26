!*==sget24.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b sget24
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE sget24( COMP, JTYPE, THRESH, ISEED, NOUNIT, N, A, LDA,
!                          H, HT, WR, WI, WRT, WIT, WRTMP, WITMP, VS,
!                          LDVS, VS1, RCDEIN, RCDVIN, NSLCT, ISLCT,
!                          RESULT, WORK, LWORK, IWORK, BWORK, INFO )
!
!       .. Scalar Arguments ..
!       LOGICAL            COMP
!       INTEGER            INFO, JTYPE, LDA, LDVS, LWORK, N, NOUNIT, NSLCT
!       REAL               RCDEIN, RCDVIN, THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            BWORK( * )
!       INTEGER            ISEED( 4 ), ISLCT( * ), IWORK( * )
!       REAL               A( LDA, * ), H( LDA, * ), HT( LDA, * ),
!      $                   RESULT( 17 ), VS( LDVS, * ), VS1( LDVS, * ),
!      $                   WI( * ), WIT( * ), WITMP( * ), WORK( * ),
!      $                   WR( * ), WRT( * ), WRTMP( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    sget24 checks the nonsymmetric eigenvalue (Schur form) problem
!>    expert driver SGEESX.
!>
!>    If COMP = .FALSE., the first 13 of the following tests will be
!>    be performed on the input matrix A, and also tests 14 and 15
!>    if LWORK is sufficiently large.
!>    If COMP = .TRUE., all 17 test will be performed.
!>
!>    (1)     0 if T is in Schur form, 1/ulp otherwise
!>           (no sorting of eigenvalues)
!>
!>    (2)     | A - VS T VS' | / ( n |A| ulp )
!>
!>      Here VS is the matrix of Schur eigenvectors, and T is in Schur
!>      form  (no sorting of eigenvalues).
!>
!>    (3)     | I - VS VS' | / ( n ulp ) (no sorting of eigenvalues).
!>
!>    (4)     0     if WR+sqrt(-1)*WI are eigenvalues of T
!>            1/ulp otherwise
!>            (no sorting of eigenvalues)
!>
!>    (5)     0     if T(with VS) = T(without VS),
!>            1/ulp otherwise
!>            (no sorting of eigenvalues)
!>
!>    (6)     0     if eigenvalues(with VS) = eigenvalues(without VS),
!>            1/ulp otherwise
!>            (no sorting of eigenvalues)
!>
!>    (7)     0 if T is in Schur form, 1/ulp otherwise
!>            (with sorting of eigenvalues)
!>
!>    (8)     | A - VS T VS' | / ( n |A| ulp )
!>
!>      Here VS is the matrix of Schur eigenvectors, and T is in Schur
!>      form  (with sorting of eigenvalues).
!>
!>    (9)     | I - VS VS' | / ( n ulp ) (with sorting of eigenvalues).
!>
!>    (10)    0     if WR+sqrt(-1)*WI are eigenvalues of T
!>            1/ulp otherwise
!>            If workspace sufficient, also compare WR, WI with and
!>            without reciprocal condition numbers
!>            (with sorting of eigenvalues)
!>
!>    (11)    0     if T(with VS) = T(without VS),
!>            1/ulp otherwise
!>            If workspace sufficient, also compare T with and without
!>            reciprocal condition numbers
!>            (with sorting of eigenvalues)
!>
!>    (12)    0     if eigenvalues(with VS) = eigenvalues(without VS),
!>            1/ulp otherwise
!>            If workspace sufficient, also compare VS with and without
!>            reciprocal condition numbers
!>            (with sorting of eigenvalues)
!>
!>    (13)    if sorting worked and SDIM is the number of
!>            eigenvalues which were SELECTed
!>            If workspace sufficient, also compare SDIM with and
!>            without reciprocal condition numbers
!>
!>    (14)    if RCONDE the same no matter if VS and/or RCONDV computed
!>
!>    (15)    if RCONDV the same no matter if VS and/or RCONDE computed
!>
!>    (16)  |RCONDE - RCDEIN| / cond(RCONDE)
!>
!>       RCONDE is the reciprocal average eigenvalue condition number
!>       computed by SGEESX and RCDEIN (the precomputed true value)
!>       is supplied as input.  cond(RCONDE) is the condition number
!>       of RCONDE, and takes errors in computing RCONDE into account,
!>       so that the resulting quantity should be O(ULP). cond(RCONDE)
!>       is essentially given by norm(A)/RCONDV.
!>
!>    (17)  |RCONDV - RCDVIN| / cond(RCONDV)
!>
!>       RCONDV is the reciprocal right invariant subspace condition
!>       number computed by SGEESX and RCDVIN (the precomputed true
!>       value) is supplied as input. cond(RCONDV) is the condition
!>       number of RCONDV, and takes errors in computing RCONDV into
!>       account, so that the resulting quantity should be O(ULP).
!>       cond(RCONDV) is essentially given by norm(A)/RCONDE.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] COMP
!> \verbatim
!>          COMP is LOGICAL
!>          COMP describes which input tests to perform:
!>            = .FALSE. if the computed condition numbers are not to
!>                      be tested against RCDVIN and RCDEIN
!>            = .TRUE.  if they are to be compared
!> \endverbatim
!>
!> \param[in] JTYPE
!> \verbatim
!>          JTYPE is INTEGER
!>          Type of input matrix. Used to label output if error occurs.
!> \endverbatim
!>
!> \param[in] ISEED
!> \verbatim
!>          ISEED is INTEGER array, dimension (4)
!>          If COMP = .FALSE., the random number generator seed
!>          used to produce matrix.
!>          If COMP = .TRUE., ISEED(1) = the number of the example.
!>          Used to label output if error occurs.
!> \endverbatim
!>
!> \param[in] THRESH
!> \verbatim
!>          THRESH is REAL
!>          A test will count as "failed" if the "error", computed as
!>          described above, exceeds THRESH.  Note that the error
!>          is scaled to be O(1), so THRESH should be a reasonably
!>          small multiple of 1, e.g., 10 or 100.  In particular,
!>          it should not depend on the precision (single vs. double)
!>          or the size of the matrix.  It must be at least zero.
!> \endverbatim
!>
!> \param[in] NOUNIT
!> \verbatim
!>          NOUNIT is INTEGER
!>          The FORTRAN unit number for printing out error messages
!>          (e.g., if a routine returns INFO not equal to 0.)
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The dimension of A. N must be at least 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is REAL array, dimension (LDA, N)
!>          Used to hold the matrix whose eigenvalues are to be
!>          computed.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of A, and H. LDA must be at
!>          least 1 and at least N.
!> \endverbatim
!>
!> \param[out] H
!> \verbatim
!>          H is REAL array, dimension (LDA, N)
!>          Another copy of the test matrix A, modified by SGEESX.
!> \endverbatim
!>
!> \param[out] HT
!> \verbatim
!>          HT is REAL array, dimension (LDA, N)
!>          Yet another copy of the test matrix A, modified by SGEESX.
!> \endverbatim
!>
!> \param[out] WR
!> \verbatim
!>          WR is REAL array, dimension (N)
!> \endverbatim
!>
!> \param[out] WI
!> \verbatim
!>          WI is REAL array, dimension (N)
!>
!>          The real and imaginary parts of the eigenvalues of A.
!>          On exit, WR + WI*i are the eigenvalues of the matrix in A.
!> \endverbatim
!>
!> \param[out] WRT
!> \verbatim
!>          WRT is REAL array, dimension (N)
!> \endverbatim
!>
!> \param[out] WIT
!> \verbatim
!>          WIT is REAL array, dimension (N)
!>
!>          Like WR, WI, these arrays contain the eigenvalues of A,
!>          but those computed when SGEESX only computes a partial
!>          eigendecomposition, i.e. not Schur vectors
!> \endverbatim
!>
!> \param[out] WRTMP
!> \verbatim
!>          WRTMP is REAL array, dimension (N)
!> \endverbatim
!>
!> \param[out] WITMP
!> \verbatim
!>          WITMP is REAL array, dimension (N)
!>
!>          Like WR, WI, these arrays contain the eigenvalues of A,
!>          but sorted by increasing real part.
!> \endverbatim
!>
!> \param[out] VS
!> \verbatim
!>          VS is REAL array, dimension (LDVS, N)
!>          VS holds the computed Schur vectors.
!> \endverbatim
!>
!> \param[in] LDVS
!> \verbatim
!>          LDVS is INTEGER
!>          Leading dimension of VS. Must be at least max(1, N).
!> \endverbatim
!>
!> \param[out] VS1
!> \verbatim
!>          VS1 is REAL array, dimension (LDVS, N)
!>          VS1 holds another copy of the computed Schur vectors.
!> \endverbatim
!>
!> \param[in] RCDEIN
!> \verbatim
!>          RCDEIN is REAL
!>          When COMP = .TRUE. RCDEIN holds the precomputed reciprocal
!>          condition number for the average of selected eigenvalues.
!> \endverbatim
!>
!> \param[in] RCDVIN
!> \verbatim
!>          RCDVIN is REAL
!>          When COMP = .TRUE. RCDVIN holds the precomputed reciprocal
!>          condition number for the selected right invariant subspace.
!> \endverbatim
!>
!> \param[in] NSLCT
!> \verbatim
!>          NSLCT is INTEGER
!>          When COMP = .TRUE. the number of selected eigenvalues
!>          corresponding to the precomputed values RCDEIN and RCDVIN.
!> \endverbatim
!>
!> \param[in] ISLCT
!> \verbatim
!>          ISLCT is INTEGER array, dimension (NSLCT)
!>          When COMP = .TRUE. ISLCT selects the eigenvalues of the
!>          input matrix corresponding to the precomputed values RCDEIN
!>          and RCDVIN. For I=1, ... ,NSLCT, if ISLCT(I) = J, then the
!>          eigenvalue with the J-th largest real part is selected.
!>          Not referenced if COMP = .FALSE.
!> \endverbatim
!>
!> \param[out] RESULT
!> \verbatim
!>          RESULT is REAL array, dimension (17)
!>          The values computed by the 17 tests described above.
!>          The values are currently limited to 1/ulp, to avoid
!>          overflow.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (LWORK)
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The number of entries in WORK to be passed to SGEESX. This
!>          must be at least 3*N, and N+N**2 if tests 14--16 are to
!>          be performed.
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (N*N)
!> \endverbatim
!>
!> \param[out] BWORK
!> \verbatim
!>          BWORK is LOGICAL array, dimension (N)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          If 0,  successful exit.
!>          If <0, input parameter -INFO had an incorrect value.
!>          If >0, SGEESX returned an error code, the absolute
!>                 value of which is returned.
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
!> \ingroup single_eig
!
!  =====================================================================
      SUBROUTINE sget24(Comp,Jtype,Thresh,Iseed,Nounit,N,A,Lda,H,Ht,Wr, &
     &                  Wi,Wrt,Wit,Wrtmp,Witmp,Vs,Ldvs,Vs1,Rcdein,      &
     &                  Rcdvin,Nslct,Islct,Result,Work,Lwork,Iwork,     &
     &                  Bwork,Info)
      IMPLICIT NONE
!*--SGET24347
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      LOGICAL Comp
      INTEGER Info , Jtype , Lda , Ldvs , Lwork , N , Nounit , Nslct
      REAL Rcdein , Rcdvin , Thresh
!     ..
!     .. Array Arguments ..
      LOGICAL Bwork(*)
      INTEGER Iseed(4) , Islct(*) , Iwork(*)
      REAL A(Lda,*) , H(Lda,*) , Ht(Lda,*) , Result(17) , Vs(Ldvs,*) ,  &
     &     Vs1(Ldvs,*) , Wi(*) , Wit(*) , Witmp(*) , Work(*) , Wr(*) ,  &
     &     Wrt(*) , Wrtmp(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE
      PARAMETER (ZERO=0.0E0,ONE=1.0E0)
      REAL EPSIN
      PARAMETER (EPSIN=5.9605E-8)
!     ..
!     .. Local Scalars ..
      CHARACTER sort
      INTEGER i , iinfo , isort , itmp , j , kmin , knteig , liwork ,   &
     &        rsub , sdim , sdim1
      REAL anorm , eps , rcnde1 , rcndv1 , rconde , rcondv , smlnum ,   &
     &     tmp , tol , tolin , ulp , ulpinv , v , vimin , vrmin , wnorm
!     ..
!     .. Local Arrays ..
      INTEGER ipnt(20)
!     ..
!     .. Arrays in Common ..
      LOGICAL SELval(20)
      REAL SELwi(20) , SELwr(20)
!     ..
!     .. Scalars in Common ..
      INTEGER SELdim , SELopt
!     ..
!     .. Common blocks ..
      COMMON /SSLCT / SELopt , SELdim , SELval , SELwr , SELwi
!     ..
!     .. External Functions ..
      REAL SLAMCH , SLANGE
      EXTERNAL SLAMCH , SLANGE
!     ..
!     .. External Subroutines ..
      EXTERNAL SCOPY , SGEESX , SGEMM , SLACPY 
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX , MIN , REAL , SIGN , SQRT
!     ..
!     .. Executable Statements ..
!
!     Check for errors
!
      Info = 0
      IF ( Thresh<ZERO ) THEN
         Info = -3
      ELSEIF ( Nounit<=0 ) THEN
         Info = -5
      ELSEIF ( N<0 ) THEN
         Info = -6
      ELSEIF ( Lda<1 .OR. Lda<N ) THEN
         Info = -8
      ELSEIF ( Ldvs<1 .OR. Ldvs<N ) THEN
         Info = -18
      ELSEIF ( Lwork<3*N ) THEN
         Info = -26
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL xerbla('sget24',-Info)
         RETURN
      ENDIF
!
!     Quick return if nothing to do
!
      DO i = 1 , 17
         Result(i) = -ONE
      ENDDO
!
      IF ( N==0 ) RETURN
!
!     Important constants
!
      smlnum = SLAMCH('Safe minimum')
      ulp = SLAMCH('Precision')
      ulpinv = ONE/ulp
!
!     Perform tests (1)-(13)
!
      SELopt = 0
      liwork = N*N
      DO isort = 0 , 1
         IF ( isort==0 ) THEN
            sort = 'N'
            rsub = 0
         ELSE
            sort = 'S'
            rsub = 6
         ENDIF
!
!        Compute Schur form and Schur vectors, and test them
!
         CALL SLACPY('F',N,N,A,Lda,H,Lda)
         CALL SGEESX('V',sort,sslect,'N',N,H,Lda,sdim,Wr,Wi,Vs,Ldvs,    &
     &               rconde,rcondv,Work,Lwork,Iwork,liwork,Bwork,iinfo)
         IF ( iinfo/=0 .AND. iinfo/=N+2 ) THEN
            Result(1+rsub) = ulpinv
            IF ( Jtype/=22 ) THEN
               WRITE (Nounit,FMT=99002) 'SGEESX1' , iinfo , N , Jtype , &
     &                                  Iseed
            ELSE
               WRITE (Nounit,FMT=99001) 'SGEESX1' , iinfo , N , Iseed(1)
            ENDIF
            Info = ABS(iinfo)
            RETURN
         ENDIF
         IF ( isort==0 ) THEN
            CALL SCOPY(N,Wr,1,Wrtmp,1)
            CALL SCOPY(N,Wi,1,Witmp,1)
         ENDIF
!
!        Do Test (1) or Test (7)
!
         Result(1+rsub) = ZERO
         DO j = 1 , N - 2
            DO i = j + 2 , N
               IF ( H(i,j)/=ZERO ) Result(1+rsub) = ulpinv
            ENDDO
         ENDDO
         DO i = 1 , N - 2
            IF ( H(i+1,i)/=ZERO .AND. H(i+2,i+1)/=ZERO ) Result(1+rsub) &
     &           = ulpinv
         ENDDO
         DO i = 1 , N - 1
            IF ( H(i+1,i)/=ZERO ) THEN
               IF ( H(i,i)/=H(i+1,i+1) .OR. H(i,i+1)==ZERO .OR.         &
     &              SIGN(ONE,H(i+1,i))==SIGN(ONE,H(i,i+1)) )            &
     &              Result(1+rsub) = ulpinv
            ENDIF
         ENDDO
!
!        Test (2) or (8): Compute norm(A - Q*H*Q') / (norm(A) * N * ULP)
!
!        Copy A to VS1, used as workspace
!
         CALL SLACPY(' ',N,N,A,Lda,Vs1,Ldvs)
!
!        Compute Q*H and store in HT.
!
         CALL SGEMM('No transpose','No transpose',N,N,N,ONE,Vs,Ldvs,H,  &
     &              Lda,ZERO,Ht,Lda)
!
!        Compute A - Q*H*Q'
!
         CALL SGEMM('No transpose','Transpose',N,N,N,-ONE,Ht,Lda,Vs,    &
     &              Ldvs,ONE,Vs1,Ldvs)
!
         anorm = MAX(SLANGE('1',N,N,A,Lda,Work),smlnum)
         wnorm = SLANGE('1',N,N,Vs1,Ldvs,Work)
!
         IF ( anorm>wnorm ) THEN
            Result(2+rsub) = (wnorm/anorm)/(N*ulp)
         ELSEIF ( anorm<ONE ) THEN
            Result(2+rsub) = (MIN(wnorm,N*anorm)/anorm)/(N*ulp)
         ELSE
            Result(2+rsub) = MIN(wnorm/anorm,REAL(N))/(N*ulp)
         ENDIF
!
!        Test (3) or (9):  Compute norm( I - Q'*Q ) / ( N * ULP )
!
         CALL sort01('Columns',N,N,Vs,Ldvs,Work,Lwork,Result(3+rsub))
!
!        Do Test (4) or Test (10)
!
         Result(4+rsub) = ZERO
         DO i = 1 , N
            IF ( H(i,i)/=Wr(i) ) Result(4+rsub) = ulpinv
         ENDDO
         IF ( N>1 ) THEN
            IF ( H(2,1)==ZERO .AND. Wi(1)/=ZERO ) Result(4+rsub)        &
     &           = ulpinv
            IF ( H(N,N-1)==ZERO .AND. Wi(N)/=ZERO ) Result(4+rsub)      &
     &           = ulpinv
         ENDIF
         DO i = 1 , N - 1
            IF ( H(i+1,i)/=ZERO ) THEN
               tmp = SQRT(ABS(H(i+1,i)))*SQRT(ABS(H(i,i+1)))
               Result(4+rsub) = MAX(Result(4+rsub),ABS(Wi(i)-tmp)/MAX(  &
     &                          ulp*tmp,smlnum))
               Result(4+rsub) = MAX(Result(4+rsub),ABS(Wi(i+1)+tmp)/MAX(&
     &                          ulp*tmp,smlnum))
            ELSEIF ( i>1 ) THEN
               IF ( H(i+1,i)==ZERO .AND. H(i,i-1)==ZERO .AND.           &
     &              Wi(i)/=ZERO ) Result(4+rsub) = ulpinv
            ENDIF
         ENDDO
!
!        Do Test (5) or Test (11)
!
         CALL SLACPY('F',N,N,A,Lda,Ht,Lda)
         CALL SGEESX('N',sort,sslect,'N',N,Ht,Lda,sdim,Wrt,Wit,Vs,Ldvs, &
     &               rconde,rcondv,Work,Lwork,Iwork,liwork,Bwork,iinfo)
         IF ( iinfo/=0 .AND. iinfo/=N+2 ) THEN
            Result(5+rsub) = ulpinv
            IF ( Jtype/=22 ) THEN
               WRITE (Nounit,FMT=99002) 'SGEESX2' , iinfo , N , Jtype , &
     &                                  Iseed
            ELSE
               WRITE (Nounit,FMT=99001) 'SGEESX2' , iinfo , N , Iseed(1)
            ENDIF
            Info = ABS(iinfo)
            GOTO 100
         ENDIF
!
         Result(5+rsub) = ZERO
         DO j = 1 , N
            DO i = 1 , N
               IF ( H(i,j)/=Ht(i,j) ) Result(5+rsub) = ulpinv
            ENDDO
         ENDDO
!
!        Do Test (6) or Test (12)
!
         Result(6+rsub) = ZERO
         DO i = 1 , N
            IF ( Wr(i)/=Wrt(i) .OR. Wi(i)/=Wit(i) ) Result(6+rsub)      &
     &           = ulpinv
         ENDDO
!
!        Do Test (13)
!
         IF ( isort==1 ) THEN
            Result(13) = ZERO
            knteig = 0
            DO i = 1 , N
               IF ( sslect(Wr(i),Wi(i)) .OR. sslect(Wr(i),-Wi(i)) )     &
     &              knteig = knteig + 1
               IF ( i<N ) THEN
                  IF ( (sslect(Wr(i+1),Wi(i+1)) .OR. sslect(Wr(i+1),-Wi(&
     &                 i+1))) .AND.                                     &
     &                 (.NOT.(sslect(Wr(i),Wi(i)) .OR. sslect(Wr(i),    &
     &                 -Wi(i)))) .AND. iinfo/=N+2 ) Result(13) = ulpinv
               ENDIF
            ENDDO
            IF ( sdim/=knteig ) Result(13) = ulpinv
         ENDIF
!
      ENDDO
!
!     If there is enough workspace, perform tests (14) and (15)
!     as well as (10) through (13)
!
      IF ( Lwork>=N+(N*N)/2 ) THEN
!
!        Compute both RCONDE and RCONDV with VS
!
         sort = 'S'
         Result(14) = ZERO
         Result(15) = ZERO
         CALL SLACPY('F',N,N,A,Lda,Ht,Lda)
         CALL SGEESX('V',sort,sslect,'B',N,Ht,Lda,sdim1,Wrt,Wit,Vs1,    &
     &               Ldvs,rconde,rcondv,Work,Lwork,Iwork,liwork,Bwork,  &
     &               iinfo)
         IF ( iinfo/=0 .AND. iinfo/=N+2 ) THEN
            Result(14) = ulpinv
            Result(15) = ulpinv
            IF ( Jtype/=22 ) THEN
               WRITE (Nounit,FMT=99002) 'SGEESX3' , iinfo , N , Jtype , &
     &                                  Iseed
            ELSE
               WRITE (Nounit,FMT=99001) 'SGEESX3' , iinfo , N , Iseed(1)
            ENDIF
            Info = ABS(iinfo)
            GOTO 100
         ENDIF
!
!        Perform tests (10), (11), (12), and (13)
!
         DO i = 1 , N
            IF ( Wr(i)/=Wrt(i) .OR. Wi(i)/=Wit(i) ) Result(10) = ulpinv
            DO j = 1 , N
               IF ( H(i,j)/=Ht(i,j) ) Result(11) = ulpinv
               IF ( Vs(i,j)/=Vs1(i,j) ) Result(12) = ulpinv
            ENDDO
         ENDDO
         IF ( sdim/=sdim1 ) Result(13) = ulpinv
!
!        Compute both RCONDE and RCONDV without VS, and compare
!
         CALL SLACPY('F',N,N,A,Lda,Ht,Lda)
         CALL SGEESX('N',sort,sslect,'B',N,Ht,Lda,sdim1,Wrt,Wit,Vs1,    &
     &               Ldvs,rcnde1,rcndv1,Work,Lwork,Iwork,liwork,Bwork,  &
     &               iinfo)
         IF ( iinfo/=0 .AND. iinfo/=N+2 ) THEN
            Result(14) = ulpinv
            Result(15) = ulpinv
            IF ( Jtype/=22 ) THEN
               WRITE (Nounit,FMT=99002) 'SGEESX4' , iinfo , N , Jtype , &
     &                                  Iseed
            ELSE
               WRITE (Nounit,FMT=99001) 'SGEESX4' , iinfo , N , Iseed(1)
            ENDIF
            Info = ABS(iinfo)
            GOTO 100
         ENDIF
!
!        Perform tests (14) and (15)
!
         IF ( rcnde1/=rconde ) Result(14) = ulpinv
         IF ( rcndv1/=rcondv ) Result(15) = ulpinv
!
!        Perform tests (10), (11), (12), and (13)
!
         DO i = 1 , N
            IF ( Wr(i)/=Wrt(i) .OR. Wi(i)/=Wit(i) ) Result(10) = ulpinv
            DO j = 1 , N
               IF ( H(i,j)/=Ht(i,j) ) Result(11) = ulpinv
               IF ( Vs(i,j)/=Vs1(i,j) ) Result(12) = ulpinv
            ENDDO
         ENDDO
         IF ( sdim/=sdim1 ) Result(13) = ulpinv
!
!        Compute RCONDE with VS, and compare
!
         CALL SLACPY('F',N,N,A,Lda,Ht,Lda)
         CALL SGEESX('V',sort,sslect,'E',N,Ht,Lda,sdim1,Wrt,Wit,Vs1,    &
     &               Ldvs,rcnde1,rcndv1,Work,Lwork,Iwork,liwork,Bwork,  &
     &               iinfo)
         IF ( iinfo/=0 .AND. iinfo/=N+2 ) THEN
            Result(14) = ulpinv
            IF ( Jtype/=22 ) THEN
               WRITE (Nounit,FMT=99002) 'SGEESX5' , iinfo , N , Jtype , &
     &                                  Iseed
            ELSE
               WRITE (Nounit,FMT=99001) 'SGEESX5' , iinfo , N , Iseed(1)
            ENDIF
            Info = ABS(iinfo)
            GOTO 100
         ENDIF
!
!        Perform test (14)
!
         IF ( rcnde1/=rconde ) Result(14) = ulpinv
!
!        Perform tests (10), (11), (12), and (13)
!
         DO i = 1 , N
            IF ( Wr(i)/=Wrt(i) .OR. Wi(i)/=Wit(i) ) Result(10) = ulpinv
            DO j = 1 , N
               IF ( H(i,j)/=Ht(i,j) ) Result(11) = ulpinv
               IF ( Vs(i,j)/=Vs1(i,j) ) Result(12) = ulpinv
            ENDDO
         ENDDO
         IF ( sdim/=sdim1 ) Result(13) = ulpinv
!
!        Compute RCONDE without VS, and compare
!
         CALL SLACPY('F',N,N,A,Lda,Ht,Lda)
         CALL SGEESX('N',sort,sslect,'E',N,Ht,Lda,sdim1,Wrt,Wit,Vs1,    &
     &               Ldvs,rcnde1,rcndv1,Work,Lwork,Iwork,liwork,Bwork,  &
     &               iinfo)
         IF ( iinfo/=0 .AND. iinfo/=N+2 ) THEN
            Result(14) = ulpinv
            IF ( Jtype/=22 ) THEN
               WRITE (Nounit,FMT=99002) 'SGEESX6' , iinfo , N , Jtype , &
     &                                  Iseed
            ELSE
               WRITE (Nounit,FMT=99001) 'SGEESX6' , iinfo , N , Iseed(1)
            ENDIF
            Info = ABS(iinfo)
            GOTO 100
         ENDIF
!
!        Perform test (14)
!
         IF ( rcnde1/=rconde ) Result(14) = ulpinv
!
!        Perform tests (10), (11), (12), and (13)
!
         DO i = 1 , N
            IF ( Wr(i)/=Wrt(i) .OR. Wi(i)/=Wit(i) ) Result(10) = ulpinv
            DO j = 1 , N
               IF ( H(i,j)/=Ht(i,j) ) Result(11) = ulpinv
               IF ( Vs(i,j)/=Vs1(i,j) ) Result(12) = ulpinv
            ENDDO
         ENDDO
         IF ( sdim/=sdim1 ) Result(13) = ulpinv
!
!        Compute RCONDV with VS, and compare
!
         CALL SLACPY('F',N,N,A,Lda,Ht,Lda)
         CALL SGEESX('V',sort,sslect,'V',N,Ht,Lda,sdim1,Wrt,Wit,Vs1,    &
     &               Ldvs,rcnde1,rcndv1,Work,Lwork,Iwork,liwork,Bwork,  &
     &               iinfo)
         IF ( iinfo/=0 .AND. iinfo/=N+2 ) THEN
            Result(15) = ulpinv
            IF ( Jtype/=22 ) THEN
               WRITE (Nounit,FMT=99002) 'SGEESX7' , iinfo , N , Jtype , &
     &                                  Iseed
            ELSE
               WRITE (Nounit,FMT=99001) 'SGEESX7' , iinfo , N , Iseed(1)
            ENDIF
            Info = ABS(iinfo)
            GOTO 100
         ENDIF
!
!        Perform test (15)
!
         IF ( rcndv1/=rcondv ) Result(15) = ulpinv
!
!        Perform tests (10), (11), (12), and (13)
!
         DO i = 1 , N
            IF ( Wr(i)/=Wrt(i) .OR. Wi(i)/=Wit(i) ) Result(10) = ulpinv
            DO j = 1 , N
               IF ( H(i,j)/=Ht(i,j) ) Result(11) = ulpinv
               IF ( Vs(i,j)/=Vs1(i,j) ) Result(12) = ulpinv
            ENDDO
         ENDDO
         IF ( sdim/=sdim1 ) Result(13) = ulpinv
!
!        Compute RCONDV without VS, and compare
!
         CALL SLACPY('F',N,N,A,Lda,Ht,Lda)
         CALL SGEESX('N',sort,sslect,'V',N,Ht,Lda,sdim1,Wrt,Wit,Vs1,    &
     &               Ldvs,rcnde1,rcndv1,Work,Lwork,Iwork,liwork,Bwork,  &
     &               iinfo)
         IF ( iinfo/=0 .AND. iinfo/=N+2 ) THEN
            Result(15) = ulpinv
            IF ( Jtype/=22 ) THEN
               WRITE (Nounit,FMT=99002) 'SGEESX8' , iinfo , N , Jtype , &
     &                                  Iseed
            ELSE
               WRITE (Nounit,FMT=99001) 'SGEESX8' , iinfo , N , Iseed(1)
            ENDIF
            Info = ABS(iinfo)
            GOTO 100
         ENDIF
!
!        Perform test (15)
!
         IF ( rcndv1/=rcondv ) Result(15) = ulpinv
!
!        Perform tests (10), (11), (12), and (13)
!
         DO i = 1 , N
            IF ( Wr(i)/=Wrt(i) .OR. Wi(i)/=Wit(i) ) Result(10) = ulpinv
            DO j = 1 , N
               IF ( H(i,j)/=Ht(i,j) ) Result(11) = ulpinv
               IF ( Vs(i,j)/=Vs1(i,j) ) Result(12) = ulpinv
            ENDDO
         ENDDO
         IF ( sdim/=sdim1 ) Result(13) = ulpinv
!
      ENDIF
!
!
!     If there are precomputed reciprocal condition numbers, compare
!     computed values with them.
!
 100  IF ( Comp ) THEN
!
!        First set up SELOPT, SELDIM, SELVAL, SELWR, and SELWI so that
!        the logical function sslect selects the eigenvalues specified
!        by NSLCT and ISLCT.
!
         SELdim = N
         SELopt = 1
         eps = MAX(ulp,EPSIN)
         DO i = 1 , N
            ipnt(i) = i
            SELval(i) = .FALSE.
            SELwr(i) = Wrtmp(i)
            SELwi(i) = Witmp(i)
         ENDDO
         DO i = 1 , N - 1
            kmin = i
            vrmin = Wrtmp(i)
            vimin = Witmp(i)
            DO j = i + 1 , N
               IF ( Wrtmp(j)<vrmin ) THEN
                  kmin = j
                  vrmin = Wrtmp(j)
                  vimin = Witmp(j)
               ENDIF
            ENDDO
            Wrtmp(kmin) = Wrtmp(i)
            Witmp(kmin) = Witmp(i)
            Wrtmp(i) = vrmin
            Witmp(i) = vimin
            itmp = ipnt(i)
            ipnt(i) = ipnt(kmin)
            ipnt(kmin) = itmp
         ENDDO
         DO i = 1 , Nslct
            SELval(ipnt(Islct(i))) = .TRUE.
         ENDDO
!
!        Compute condition numbers
!
         CALL SLACPY('F',N,N,A,Lda,Ht,Lda)
         CALL SGEESX('N','S',sslect,'B',N,Ht,Lda,sdim1,Wrt,Wit,Vs1,Ldvs,&
     &               rconde,rcondv,Work,Lwork,Iwork,liwork,Bwork,iinfo)
         IF ( iinfo/=0 .AND. iinfo/=N+2 ) THEN
            Result(16) = ulpinv
            Result(17) = ulpinv
            WRITE (Nounit,FMT=99001) 'SGEESX9' , iinfo , N , Iseed(1)
            Info = ABS(iinfo)
            GOTO 99999
         ENDIF
!
!        Compare condition number for average of selected eigenvalues
!        taking its condition number into account
!
         anorm = SLANGE('1',N,N,A,Lda,Work)
         v = MAX(REAL(N)*eps*anorm,smlnum)
         IF ( anorm==ZERO ) v = ONE
         IF ( v>rcondv ) THEN
            tol = ONE
         ELSE
            tol = v/rcondv
         ENDIF
         IF ( v>Rcdvin ) THEN
            tolin = ONE
         ELSE
            tolin = v/Rcdvin
         ENDIF
         tol = MAX(tol,smlnum/eps)
         tolin = MAX(tolin,smlnum/eps)
         IF ( eps*(Rcdein-tolin)>rconde+tol ) THEN
            Result(16) = ulpinv
         ELSEIF ( Rcdein-tolin>rconde+tol ) THEN
            Result(16) = (Rcdein-tolin)/(rconde+tol)
         ELSEIF ( Rcdein+tolin<eps*(rconde-tol) ) THEN
            Result(16) = ulpinv
         ELSEIF ( Rcdein+tolin<rconde-tol ) THEN
            Result(16) = (rconde-tol)/(Rcdein+tolin)
         ELSE
            Result(16) = ONE
         ENDIF
!
!        Compare condition numbers for right invariant subspace
!        taking its condition number into account
!
         IF ( v>rcondv*rconde ) THEN
            tol = rcondv
         ELSE
            tol = v/rconde
         ENDIF
         IF ( v>Rcdvin*Rcdein ) THEN
            tolin = Rcdvin
         ELSE
            tolin = v/Rcdein
         ENDIF
         tol = MAX(tol,smlnum/eps)
         tolin = MAX(tolin,smlnum/eps)
         IF ( eps*(Rcdvin-tolin)>rcondv+tol ) THEN
            Result(17) = ulpinv
         ELSEIF ( Rcdvin-tolin>rcondv+tol ) THEN
            Result(17) = (Rcdvin-tolin)/(rcondv+tol)
         ELSEIF ( Rcdvin+tolin<eps*(rcondv-tol) ) THEN
            Result(17) = ulpinv
         ELSEIF ( Rcdvin+tolin<rcondv-tol ) THEN
            Result(17) = (rcondv-tol)/(Rcdvin+tolin)
         ELSE
            Result(17) = ONE
         ENDIF
!
!
      ENDIF
!
99001 FORMAT (' sget24: ',A,' returned INFO=',I6,'.',/9X,'N=',I6,       &
     &        ', INPUT EXAMPLE NUMBER = ',I4)
99002 FORMAT (' sget24: ',A,' returned INFO=',I6,'.',/9X,'N=',I6,       &
     &        ', JTYPE=',I6,', ISEED=(',3(I5,','),I5,')')
!
!
!     End of sget24
!
99999 END SUBROUTINE sget24
