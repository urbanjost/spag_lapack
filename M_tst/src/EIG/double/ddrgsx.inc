!*==ddrgsx.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b ddrgsx
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ddrgsx( NSIZE, NCMAX, THRESH, NIN, NOUT, A, LDA, B, AI,
!                          BI, Z, Q, ALPHAR, ALPHAI, BETA, C, LDC, S,
!                          WORK, LWORK, IWORK, LIWORK, BWORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, LDC, LIWORK, LWORK, NCMAX, NIN,
!      $                   NOUT, NSIZE
!       DOUBLE PRECISION   THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            BWORK( * )
!       INTEGER            IWORK( * )
!       DOUBLE PRECISION   A( LDA, * ), AI( LDA, * ), ALPHAI( * ),
!      $                   ALPHAR( * ), B( LDA, * ), BETA( * ),
!      $                   BI( LDA, * ), C( LDC, * ), Q( LDA, * ), S( * ),
!      $                   WORK( * ), Z( LDA, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ddrgsx checks the nonsymmetric generalized eigenvalue (Schur form)
!> problem expert driver DGGESX.
!>
!> DGGESX factors A and B as Q S Z' and Q T Z', where ' means
!> transpose, T is upper triangular, S is in generalized Schur form
!> (block upper triangular, with 1x1 and 2x2 blocks on the diagonal,
!> the 2x2 blocks corresponding to complex conjugate pairs of
!> generalized eigenvalues), and Q and Z are orthogonal.  It also
!> computes the generalized eigenvalues (alpha(1),beta(1)), ...,
!> (alpha(n),beta(n)). Thus, w(j) = alpha(j)/beta(j) is a root of the
!> characteristic equation
!>
!>     det( A - w(j) B ) = 0
!>
!> Optionally it also reorders the eigenvalues so that a selected
!> cluster of eigenvalues appears in the leading diagonal block of the
!> Schur forms; computes a reciprocal condition number for the average
!> of the selected eigenvalues; and computes a reciprocal condition
!> number for the right and left deflating subspaces corresponding to
!> the selected eigenvalues.
!>
!> When ddrgsx is called with NSIZE > 0, five (5) types of built-in
!> matrix pairs are used to test the routine DGGESX.
!>
!> When ddrgsx is called with NSIZE = 0, it reads in test matrix data
!> to test DGGESX.
!>
!> For each matrix pair, the following tests will be performed and
!> compared with the threshold THRESH except for the tests (7) and (9):
!>
!> (1)   | A - Q S Z' | / ( |A| n ulp )
!>
!> (2)   | B - Q T Z' | / ( |B| n ulp )
!>
!> (3)   | I - QQ' | / ( n ulp )
!>
!> (4)   | I - ZZ' | / ( n ulp )
!>
!> (5)   if A is in Schur form (i.e. quasi-triangular form)
!>
!> (6)   maximum over j of D(j)  where:
!>
!>       if alpha(j) is real:
!>                     |alpha(j) - S(j,j)|        |beta(j) - T(j,j)|
!>           D(j) = ------------------------ + -----------------------
!>                  max(|alpha(j)|,|S(j,j)|)   max(|beta(j)|,|T(j,j)|)
!>
!>       if alpha(j) is complex:
!>                                 | det( s S - w T ) |
!>           D(j) = ---------------------------------------------------
!>                  ulp max( s norm(S), |w| norm(T) )*norm( s S - w T )
!>
!>           and S and T are here the 2 x 2 diagonal blocks of S and T
!>           corresponding to the j-th and j+1-th eigenvalues.
!>
!> (7)   if sorting worked and SDIM is the number of eigenvalues
!>       which were selected.
!>
!> (8)   the estimated value DIF does not differ from the true values of
!>       Difu and Difl more than a factor 10*THRESH. If the estimate DIF
!>       equals zero the corresponding true values of Difu and Difl
!>       should be less than EPS*norm(A, B). If the true value of Difu
!>       and Difl equal zero, the estimate DIF should be less than
!>       EPS*norm(A, B).
!>
!> (9)   If INFO = N+3 is returned by DGGESX, the reordering "failed"
!>       and we check that DIF = PL = PR = 0 and that the true value of
!>       Difu and Difl is < EPS*norm(A, B). We count the events when
!>       INFO=N+3.
!>
!> For read-in test matrices, the above tests are run except that the
!> exact value for DIF (and PL) is input data.  Additionally, there is
!> one more test run for read-in test matrices:
!>
!> (10)  the estimated value PL does not differ from the true value of
!>       PLTRU more than a factor THRESH. If the estimate PL equals
!>       zero the corresponding true value of PLTRU should be less than
!>       EPS*norm(A, B). If the true value of PLTRU equal zero, the
!>       estimate PL should be less than EPS*norm(A, B).
!>
!> Note that for the built-in tests, a total of 10*NSIZE*(NSIZE-1)
!> matrix pairs are generated and tested. NSIZE should be kept small.
!>
!> SVD (routine DGESVD) is used for computing the true value of DIF_u
!> and DIF_l when testing the built-in test problems.
!>
!> Built-in Test Matrices
!> ======================
!>
!> All built-in test matrices are the 2 by 2 block of triangular
!> matrices
!>
!>          A = [ A11 A12 ]    and      B = [ B11 B12 ]
!>              [     A22 ]                 [     B22 ]
!>
!> where for different type of A11 and A22 are given as the following.
!> A12 and B12 are chosen so that the generalized Sylvester equation
!>
!>          A11*R - L*A22 = -A12
!>          B11*R - L*B22 = -B12
!>
!> have prescribed solution R and L.
!>
!> Type 1:  A11 = J_m(1,-1) and A_22 = J_k(1-a,1).
!>          B11 = I_m, B22 = I_k
!>          where J_k(a,b) is the k-by-k Jordan block with ``a'' on
!>          diagonal and ``b'' on superdiagonal.
!>
!> Type 2:  A11 = (a_ij) = ( 2(.5-sin(i)) ) and
!>          B11 = (b_ij) = ( 2(.5-sin(ij)) ) for i=1,...,m, j=i,...,m
!>          A22 = (a_ij) = ( 2(.5-sin(i+j)) ) and
!>          B22 = (b_ij) = ( 2(.5-sin(ij)) ) for i=m+1,...,k, j=i,...,k
!>
!> Type 3:  A11, A22 and B11, B22 are chosen as for Type 2, but each
!>          second diagonal block in A_11 and each third diagonal block
!>          in A_22 are made as 2 by 2 blocks.
!>
!> Type 4:  A11 = ( 20(.5 - sin(ij)) ) and B22 = ( 2(.5 - sin(i+j)) )
!>             for i=1,...,m,  j=1,...,m and
!>          A22 = ( 20(.5 - sin(i+j)) ) and B22 = ( 2(.5 - sin(ij)) )
!>             for i=m+1,...,k,  j=m+1,...,k
!>
!> Type 5:  (A,B) and have potentially close or common eigenvalues and
!>          very large departure from block diagonality A_11 is chosen
!>          as the m x m leading submatrix of A_1:
!>                  |  1  b                            |
!>                  | -b  1                            |
!>                  |        1+d  b                    |
!>                  |         -b 1+d                   |
!>           A_1 =  |                  d  1            |
!>                  |                 -1  d            |
!>                  |                        -d  1     |
!>                  |                        -1 -d     |
!>                  |                               1  |
!>          and A_22 is chosen as the k x k leading submatrix of A_2:
!>                  | -1  b                            |
!>                  | -b -1                            |
!>                  |       1-d  b                     |
!>                  |       -b  1-d                    |
!>           A_2 =  |                 d 1+b            |
!>                  |               -1-b d             |
!>                  |                       -d  1+b    |
!>                  |                      -1+b  -d    |
!>                  |                              1-d |
!>          and matrix B are chosen as identity matrices (see dlatm5).
!>
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] NSIZE
!> \verbatim
!>          NSIZE is INTEGER
!>          The maximum size of the matrices to use. NSIZE >= 0.
!>          If NSIZE = 0, no built-in tests matrices are used, but
!>          read-in test matrices are used to test DGGESX.
!> \endverbatim
!>
!> \param[in] NCMAX
!> \verbatim
!>          NCMAX is INTEGER
!>          Maximum allowable NMAX for generating Kroneker matrix
!>          in call to dlakf2
!> \endverbatim
!>
!> \param[in] THRESH
!> \verbatim
!>          THRESH is DOUBLE PRECISION
!>          A test will count as "failed" if the "error", computed as
!>          described above, exceeds THRESH.  Note that the error
!>          is scaled to be O(1), so THRESH should be a reasonably
!>          small multiple of 1, e.g., 10 or 100.  In particular,
!>          it should not depend on the precision (single vs. double)
!>          or the size of the matrix.  THRESH >= 0.
!> \endverbatim
!>
!> \param[in] NIN
!> \verbatim
!>          NIN is INTEGER
!>          The FORTRAN unit number for reading in the data file of
!>          problems to solve.
!> \endverbatim
!>
!> \param[in] NOUT
!> \verbatim
!>          NOUT is INTEGER
!>          The FORTRAN unit number for printing out error messages
!>          (e.g., if a routine returns IINFO not equal to 0.)
!> \endverbatim
!>
!> \param[out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA, NSIZE)
!>          Used to store the matrix whose eigenvalues are to be
!>          computed.  On exit, A contains the last matrix actually used.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of A, B, AI, BI, Z and Q,
!>          LDA >= max( 1, NSIZE ). For the read-in test,
!>          LDA >= max( 1, N ), N is the size of the test matrices.
!> \endverbatim
!>
!> \param[out] B
!> \verbatim
!>          B is DOUBLE PRECISION array, dimension (LDA, NSIZE)
!>          Used to store the matrix whose eigenvalues are to be
!>          computed.  On exit, B contains the last matrix actually used.
!> \endverbatim
!>
!> \param[out] AI
!> \verbatim
!>          AI is DOUBLE PRECISION array, dimension (LDA, NSIZE)
!>          Copy of A, modified by DGGESX.
!> \endverbatim
!>
!> \param[out] BI
!> \verbatim
!>          BI is DOUBLE PRECISION array, dimension (LDA, NSIZE)
!>          Copy of B, modified by DGGESX.
!> \endverbatim
!>
!> \param[out] Z
!> \verbatim
!>          Z is DOUBLE PRECISION array, dimension (LDA, NSIZE)
!>          Z holds the left Schur vectors computed by DGGESX.
!> \endverbatim
!>
!> \param[out] Q
!> \verbatim
!>          Q is DOUBLE PRECISION array, dimension (LDA, NSIZE)
!>          Q holds the right Schur vectors computed by DGGESX.
!> \endverbatim
!>
!> \param[out] ALPHAR
!> \verbatim
!>          ALPHAR is DOUBLE PRECISION array, dimension (NSIZE)
!> \endverbatim
!>
!> \param[out] ALPHAI
!> \verbatim
!>          ALPHAI is DOUBLE PRECISION array, dimension (NSIZE)
!> \endverbatim
!>
!> \param[out] BETA
!> \verbatim
!>          BETA is DOUBLE PRECISION array, dimension (NSIZE)
!>
!>          On exit, (ALPHAR + ALPHAI*i)/BETA are the eigenvalues.
!> \endverbatim
!>
!> \param[out] C
!> \verbatim
!>          C is DOUBLE PRECISION array, dimension (LDC, LDC)
!>          Store the matrix generated by subroutine dlakf2, this is the
!>          matrix formed by Kronecker products used for estimating
!>          DIF.
!> \endverbatim
!>
!> \param[in] LDC
!> \verbatim
!>          LDC is INTEGER
!>          The leading dimension of C. LDC >= max(1, LDA*LDA/2 ).
!> \endverbatim
!>
!> \param[out] S
!> \verbatim
!>          S is DOUBLE PRECISION array, dimension (LDC)
!>          Singular values of C
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (LWORK)
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK.
!>          LWORK >= MAX( 5*NSIZE*NSIZE/2 - 2, 10*(NSIZE+1) )
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (LIWORK)
!> \endverbatim
!>
!> \param[in] LIWORK
!> \verbatim
!>          LIWORK is INTEGER
!>          The dimension of the array IWORK. LIWORK >= NSIZE + 6.
!> \endverbatim
!>
!> \param[out] BWORK
!> \verbatim
!>          BWORK is LOGICAL array, dimension (LDA)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
!>          > 0:  A routine returned an error code.
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
!> \ingroup double_eig
!
!  =====================================================================
      SUBROUTINE ddrgsx(Nsize,Ncmax,Thresh,Nin,Nout,A,Lda,B,Ai,Bi,Z,Q,  &
     &                  Alphar,Alphai,Beta,C,Ldc,S,Work,Lwork,Iwork,    &
     &                  Liwork,Bwork,Info)
      IMPLICIT NONE
!*--DDRGSX363
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2016
!
!     .. Scalar Arguments ..
      INTEGER Info , Lda , Ldc , Liwork , Lwork , Ncmax , Nin , Nout ,  &
     &        Nsize
      DOUBLE PRECISION Thresh
!     ..
!     .. Array Arguments ..
      LOGICAL Bwork(*)
      INTEGER Iwork(*)
      DOUBLE PRECISION A(Lda,*) , Ai(Lda,*) , Alphai(*) , Alphar(*) ,   &
     &                 B(Lda,*) , Beta(*) , Bi(Lda,*) , C(Ldc,*) ,      &
     &                 Q(Lda,*) , S(*) , Work(*) , Z(Lda,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE , TEN
      PARAMETER (ZERO=0.0D+0,ONE=1.0D+0,TEN=1.0D+1)
!     ..
!     .. Local Scalars ..
      LOGICAL ilabad
      CHARACTER sense
      INTEGER bdspac , i , i1 , ifunc , iinfo , j , linfo , maxwrk ,    &
     &        minwrk , mm , mn2 , nerrs , nptknt , ntest , ntestt ,     &
     &        prtype , qba , qbb
      DOUBLE PRECISION abnrm , bignum , diftru , pltru , smlnum ,       &
     &                 temp1 , temp2 , thrsh2 , ulp , ulpinv , weight
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION difest(2) , pl(2) , result(10)
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMCH , DLANGE
      EXTERNAL DLAMCH , DLANGE
!     ..
!     .. External Subroutines ..
      EXTERNAL DGESVD , DGGESX , DLABAD , DLACPY , DLASET 
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX , SQRT
!     ..
!     .. Scalars in Common ..
      LOGICAL FS
      INTEGER K , M , MPLusn , N
!     ..
!     .. Common blocks ..
      COMMON /MN    / M , N , MPLusn , K , FS
!     ..
!     .. Executable Statements ..
!
!     Check for errors
!
      IF ( Nsize<0 ) THEN
         Info = -1
      ELSEIF ( Thresh<ZERO ) THEN
         Info = -2
      ELSEIF ( Nin<=0 ) THEN
         Info = -3
      ELSEIF ( Nout<=0 ) THEN
         Info = -4
      ELSEIF ( Lda<1 .OR. Lda<Nsize ) THEN
         Info = -6
      ELSEIF ( Ldc<1 .OR. Ldc<Nsize*Nsize/2 ) THEN
         Info = -17
      ELSEIF ( Liwork<Nsize+6 ) THEN
         Info = -21
      ENDIF
!
!     Compute workspace
!      (Note: Comments in the code beginning "Workspace:" describe the
!       minimal amount of workspace needed at that point in the code,
!       as well as the preferred amount for good performance.
!       NB refers to the optimal block size for the immediately
!       following subroutine, as returned by ilaenv.)
!
      minwrk = 1
      IF ( Info==0 .AND. Lwork>=1 ) THEN
         minwrk = MAX(10*(Nsize+1),5*Nsize*Nsize/2)
!
!        workspace for sggesx
!
         maxwrk = 9*(Nsize+1) + Nsize*ilaenv(1,'DGEQRF',' ',Nsize,1,    &
     &            Nsize,0)
         maxwrk = MAX(maxwrk,9*(Nsize+1)                                &
     &            +Nsize*ilaenv(1,'DORGQR',' ',Nsize,1,Nsize,-1))
!
!        workspace for dgesvd
!
         bdspac = 5*Nsize*Nsize/2
         maxwrk = MAX(maxwrk,3*Nsize*Nsize/2+                           &
     &            Nsize*Nsize*ilaenv(1,'DGEBRD',' ',Nsize*Nsize/2,      &
     &            Nsize*Nsize/2,-1,-1))
         maxwrk = MAX(maxwrk,bdspac)
!
         maxwrk = MAX(maxwrk,minwrk)
!
         Work(1) = maxwrk
      ENDIF
!
      IF ( Lwork<minwrk ) Info = -19
!
      IF ( Info/=0 ) THEN
         CALL xerbla('ddrgsx',-Info)
         RETURN
      ENDIF
!
!     Important constants
!
      ulp = DLAMCH('P')
      ulpinv = ONE/ulp
      smlnum = DLAMCH('S')/ulp
      bignum = ONE/smlnum
      CALL DLABAD(smlnum,bignum)
      thrsh2 = TEN*Thresh
      ntestt = 0
      nerrs = 0
!
!     Go to the tests for read-in matrix pairs
!
      ifunc = 0
      IF ( Nsize==0 ) THEN
!
!
!     Read in data from file to check accuracy of condition estimation
!     Read input data until N=0
!
         nptknt = 0
         DO
!
            READ (Nin,FMT=*,END=100) MPLusn
            IF ( MPLusn==0 ) EXIT
            READ (Nin,FMT=*,END=100) N
            DO i = 1 , MPLusn
               READ (Nin,FMT=*) (Ai(i,j),j=1,MPLusn)
            ENDDO
            DO i = 1 , MPLusn
               READ (Nin,FMT=*) (Bi(i,j),j=1,MPLusn)
            ENDDO
            READ (Nin,FMT=*) pltru , diftru
!
            nptknt = nptknt + 1
            FS = .TRUE.
            K = 0
            M = MPLusn - N
!
            CALL DLACPY('Full',MPLusn,MPLusn,Ai,Lda,A,Lda)
            CALL DLACPY('Full',MPLusn,MPLusn,Bi,Lda,B,Lda)
!
!     Compute the Schur factorization while swapping the
!     m-by-m (1,1)-blocks with n-by-n (2,2)-blocks.
!
            CALL DGGESX('V','V','S',dlctsx,'B',MPLusn,Ai,Lda,Bi,Lda,mm, &
     &                  Alphar,Alphai,Beta,Q,Lda,Z,Lda,pl,difest,Work,  &
     &                  Lwork,Iwork,Liwork,Bwork,linfo)
!
            IF ( linfo/=0 .AND. linfo/=MPLusn+2 ) THEN
               result(1) = ulpinv
               WRITE (Nout,FMT=99002) 'DGGESX' , linfo , MPLusn , nptknt
               CYCLE
            ENDIF
!
!     Compute the norm(A, B)
!        (should this be norm of (A,B) or (AI,BI)?)
!
            CALL DLACPY('Full',MPLusn,MPLusn,Ai,Lda,Work,MPLusn)
            CALL DLACPY('Full',MPLusn,MPLusn,Bi,Lda,                    &
     &                  Work(MPLusn*MPLusn+1),MPLusn)
            abnrm = DLANGE('Fro',MPLusn,2*MPLusn,Work,MPLusn,Work)
!
!     Do tests (1) to (4)
!
            CALL dget51(1,MPLusn,A,Lda,Ai,Lda,Q,Lda,Z,Lda,Work,result(1)&
     &                  )
            CALL dget51(1,MPLusn,B,Lda,Bi,Lda,Q,Lda,Z,Lda,Work,result(2)&
     &                  )
            CALL dget51(3,MPLusn,B,Lda,Bi,Lda,Q,Lda,Q,Lda,Work,result(3)&
     &                  )
            CALL dget51(3,MPLusn,B,Lda,Bi,Lda,Z,Lda,Z,Lda,Work,result(4)&
     &                  )
!
!     Do tests (5) and (6): check Schur form of A and compare
!     eigenvalues with diagonals.
!
            ntest = 6
            temp1 = ZERO
            result(5) = ZERO
            result(6) = ZERO
!
            DO j = 1 , MPLusn
               ilabad = .FALSE.
               IF ( Alphai(j)==ZERO ) THEN
                  temp2 = (ABS(Alphar(j)-Ai(j,j))                       &
     &                    /MAX(smlnum,ABS(Alphar(j)),ABS(Ai(j,j)))      &
     &                    +ABS(Beta(j)-Bi(j,j))                         &
     &                    /MAX(smlnum,ABS(Beta(j)),ABS(Bi(j,j))))/ulp
                  IF ( j<MPLusn ) THEN
                     IF ( Ai(j+1,j)/=ZERO ) THEN
                        ilabad = .TRUE.
                        result(5) = ulpinv
                     ENDIF
                  ENDIF
                  IF ( j>1 ) THEN
                     IF ( Ai(j,j-1)/=ZERO ) THEN
                        ilabad = .TRUE.
                        result(5) = ulpinv
                     ENDIF
                  ENDIF
               ELSE
                  IF ( Alphai(j)>ZERO ) THEN
                     i1 = j
                  ELSE
                     i1 = j - 1
                  ENDIF
                  IF ( i1<=0 .OR. i1>=MPLusn ) THEN
                     ilabad = .TRUE.
                  ELSEIF ( i1<MPLusn-1 ) THEN
                     IF ( Ai(i1+2,i1+1)/=ZERO ) THEN
                        ilabad = .TRUE.
                        result(5) = ulpinv
                     ENDIF
                  ELSEIF ( i1>1 ) THEN
                     IF ( Ai(i1,i1-1)/=ZERO ) THEN
                        ilabad = .TRUE.
                        result(5) = ulpinv
                     ENDIF
                  ENDIF
                  IF ( .NOT.ilabad ) THEN
                     CALL dget53(Ai(i1,i1),Lda,Bi(i1,i1),Lda,Beta(j),   &
     &                           Alphar(j),Alphai(j),temp2,iinfo)
                     IF ( iinfo>=3 ) THEN
                        WRITE (Nout,FMT=99003) iinfo , j , MPLusn ,     &
     &                         nptknt
                        Info = ABS(iinfo)
                     ENDIF
                  ELSE
                     temp2 = ulpinv
                  ENDIF
               ENDIF
               temp1 = MAX(temp1,temp2)
               IF ( ilabad ) WRITE (Nout,FMT=99004) j , MPLusn , nptknt
            ENDDO
            result(6) = temp1
!
!     Test (7) (if sorting worked)  <--------- need to be checked.
!
            ntest = 7
            result(7) = ZERO
            IF ( linfo==MPLusn+3 ) result(7) = ulpinv
!
!     Test (8): compare the estimated value of DIF and its true value.
!
            ntest = 8
            result(8) = ZERO
            IF ( difest(2)==ZERO ) THEN
               IF ( diftru>abnrm*ulp ) result(8) = ulpinv
            ELSEIF ( diftru==ZERO ) THEN
               IF ( difest(2)>abnrm*ulp ) result(8) = ulpinv
            ELSEIF ( (diftru>thrsh2*difest(2)) .OR.                     &
     &               (diftru*thrsh2<difest(2)) ) THEN
               result(8) = MAX(diftru/difest(2),difest(2)/diftru)
            ENDIF
!
!     Test (9)
!
            ntest = 9
            result(9) = ZERO
            IF ( linfo==(MPLusn+2) ) THEN
               IF ( diftru>abnrm*ulp ) result(9) = ulpinv
               IF ( (ifunc>1) .AND. (difest(2)/=ZERO) ) result(9)       &
     &              = ulpinv
               IF ( (ifunc==1) .AND. (pl(1)/=ZERO) ) result(9) = ulpinv
            ENDIF
!
!     Test (10): compare the estimated value of PL and it true value.
!
            ntest = 10
            result(10) = ZERO
            IF ( pl(1)==ZERO ) THEN
               IF ( pltru>abnrm*ulp ) result(10) = ulpinv
            ELSEIF ( pltru==ZERO ) THEN
               IF ( pl(1)>abnrm*ulp ) result(10) = ulpinv
            ELSEIF ( (pltru>Thresh*pl(1)) .OR. (pltru*Thresh<pl(1)) )   &
     &               THEN
               result(10) = ulpinv
            ENDIF
!
            ntestt = ntestt + ntest
!
!     Print out tests which fail.
!
            DO j = 1 , ntest
               IF ( result(j)>=Thresh ) THEN
!
!           If this is the first test to fail,
!           print a header to the data file.
!
                  IF ( nerrs==0 ) THEN
                     WRITE (Nout,FMT=99005) 'DGX'
!
!              Matrix types
!
                     WRITE (Nout,FMT=99006)
!
!              Tests performed
!
                     WRITE (Nout,FMT=99008) 'orthogonal' , '''' ,       &
     &                      'transpose' , ('''',i=1,4)
!
                  ENDIF
                  nerrs = nerrs + 1
                  IF ( result(j)<10000.0D0 ) THEN
                     WRITE (Nout,FMT=99011) nptknt , MPLusn , j ,       &
     &                      result(j)
                  ELSE
                     WRITE (Nout,FMT=99012) nptknt , MPLusn , j ,       &
     &                      result(j)
                  ENDIF
               ENDIF
!
!
            ENDDO
         ENDDO
      ELSE
!
!     Test the built-in matrix pairs.
!     Loop over different functions (IFUNC) of DGGESX, types (PRTYPE)
!     of test matrices, different size (M+N)
!
         prtype = 0
         qba = 3
         qbb = 4
         weight = SQRT(ulp)
!
         DO ifunc = 0 , 3
            DO prtype = 1 , 5
               DO M = 1 , Nsize - 1
                  DO N = 1 , Nsize - M
!
                     weight = ONE/weight
                     MPLusn = M + N
!
!                 Generate test matrices
!
                     FS = .TRUE.
                     K = 0
!
                     CALL DLASET('Full',MPLusn,MPLusn,ZERO,ZERO,Ai,Lda)
                     CALL DLASET('Full',MPLusn,MPLusn,ZERO,ZERO,Bi,Lda)
!
                     CALL dlatm5(prtype,M,N,Ai,Lda,Ai(M+1,M+1),Lda,     &
     &                           Ai(1,M+1),Lda,Bi,Lda,Bi(M+1,M+1),Lda,  &
     &                           Bi(1,M+1),Lda,Q,Lda,Z,Lda,weight,qba,  &
     &                           qbb)
!
!                 Compute the Schur factorization and swapping the
!                 m-by-m (1,1)-blocks with n-by-n (2,2)-blocks.
!                 Swapping is accomplished via the function dlctsx
!                 which is supplied below.
!
                     IF ( ifunc==0 ) THEN
                        sense = 'N'
                     ELSEIF ( ifunc==1 ) THEN
                        sense = 'E'
                     ELSEIF ( ifunc==2 ) THEN
                        sense = 'V'
                     ELSEIF ( ifunc==3 ) THEN
                        sense = 'B'
                     ENDIF
!
                     CALL DLACPY('Full',MPLusn,MPLusn,Ai,Lda,A,Lda)
                     CALL DLACPY('Full',MPLusn,MPLusn,Bi,Lda,B,Lda)
!
                     CALL DGGESX('V','V','S',dlctsx,sense,MPLusn,Ai,Lda,&
     &                           Bi,Lda,mm,Alphar,Alphai,Beta,Q,Lda,Z,  &
     &                           Lda,pl,difest,Work,Lwork,Iwork,Liwork, &
     &                           Bwork,linfo)
!
                     IF ( linfo/=0 .AND. linfo/=MPLusn+2 ) THEN
                        result(1) = ulpinv
                        WRITE (Nout,FMT=99001) 'DGGESX' , linfo ,       &
     &                         MPLusn , prtype
                        Info = linfo
                        CYCLE
                     ENDIF
!
!                 Compute the norm(A, B)
!
                     CALL DLACPY('Full',MPLusn,MPLusn,Ai,Lda,Work,      &
     &                           MPLusn)
                     CALL DLACPY('Full',MPLusn,MPLusn,Bi,Lda,           &
     &                           Work(MPLusn*MPLusn+1),MPLusn)
                     abnrm = DLANGE('Fro',MPLusn,2*MPLusn,Work,MPLusn,  &
     &                       Work)
!
!                 Do tests (1) to (4)
!
                     CALL dget51(1,MPLusn,A,Lda,Ai,Lda,Q,Lda,Z,Lda,Work,&
     &                           result(1))
                     CALL dget51(1,MPLusn,B,Lda,Bi,Lda,Q,Lda,Z,Lda,Work,&
     &                           result(2))
                     CALL dget51(3,MPLusn,B,Lda,Bi,Lda,Q,Lda,Q,Lda,Work,&
     &                           result(3))
                     CALL dget51(3,MPLusn,B,Lda,Bi,Lda,Z,Lda,Z,Lda,Work,&
     &                           result(4))
                     ntest = 4
!
!                 Do tests (5) and (6): check Schur form of A and
!                 compare eigenvalues with diagonals.
!
                     temp1 = ZERO
                     result(5) = ZERO
                     result(6) = ZERO
!
                     DO j = 1 , MPLusn
                        ilabad = .FALSE.
                        IF ( Alphai(j)==ZERO ) THEN
                           temp2 = (ABS(Alphar(j)-Ai(j,j))              &
     &                             /MAX(smlnum,ABS(Alphar(j)),          &
     &                             ABS(Ai(j,j)))+ABS(Beta(j)-Bi(j,j))   &
     &                             /MAX(smlnum,ABS(Beta(j)),ABS(Bi(j,j))&
     &                             ))/ulp
                           IF ( j<MPLusn ) THEN
                              IF ( Ai(j+1,j)/=ZERO ) THEN
                                 ilabad = .TRUE.
                                 result(5) = ulpinv
                              ENDIF
                           ENDIF
                           IF ( j>1 ) THEN
                              IF ( Ai(j,j-1)/=ZERO ) THEN
                                 ilabad = .TRUE.
                                 result(5) = ulpinv
                              ENDIF
                           ENDIF
                        ELSE
                           IF ( Alphai(j)>ZERO ) THEN
                              i1 = j
                           ELSE
                              i1 = j - 1
                           ENDIF
                           IF ( i1<=0 .OR. i1>=MPLusn ) THEN
                              ilabad = .TRUE.
                           ELSEIF ( i1<MPLusn-1 ) THEN
                              IF ( Ai(i1+2,i1+1)/=ZERO ) THEN
                                 ilabad = .TRUE.
                                 result(5) = ulpinv
                              ENDIF
                           ELSEIF ( i1>1 ) THEN
                              IF ( Ai(i1,i1-1)/=ZERO ) THEN
                                 ilabad = .TRUE.
                                 result(5) = ulpinv
                              ENDIF
                           ENDIF
                           IF ( .NOT.ilabad ) THEN
                              CALL dget53(Ai(i1,i1),Lda,Bi(i1,i1),Lda,  &
     &                           Beta(j),Alphar(j),Alphai(j),temp2,     &
     &                           iinfo)
                              IF ( iinfo>=3 ) THEN
                                 WRITE (Nout,FMT=99003) iinfo , j ,     &
     &                                  MPLusn , prtype
                                 Info = ABS(iinfo)
                              ENDIF
                           ELSE
                              temp2 = ulpinv
                           ENDIF
                        ENDIF
                        temp1 = MAX(temp1,temp2)
                        IF ( ilabad ) WRITE (Nout,FMT=99004) j ,        &
     &                       MPLusn , prtype
                     ENDDO
                     result(6) = temp1
                     ntest = ntest + 2
!
!                 Test (7) (if sorting worked)
!
                     result(7) = ZERO
                     IF ( linfo==MPLusn+3 ) THEN
                        result(7) = ulpinv
                     ELSEIF ( mm/=N ) THEN
                        result(7) = ulpinv
                     ENDIF
                     ntest = ntest + 1
!
!                 Test (8): compare the estimated value DIF and its
!                 value. first, compute the exact DIF.
!
                     result(8) = ZERO
                     mn2 = mm*(MPLusn-mm)*2
                     IF ( ifunc>=2 .AND. mn2<=Ncmax*Ncmax ) THEN
!
!                    Note: for either following two causes, there are
!                    almost same number of test cases fail the test.
!
                        CALL dlakf2(mm,MPLusn-mm,Ai,Lda,Ai(mm+1,mm+1),  &
     &                              Bi,Bi(mm+1,mm+1),C,Ldc)
!
                        CALL DGESVD('N','N',mn2,mn2,C,Ldc,S,Work,1,     &
     &                              Work(2),1,Work(3),Lwork-2,Info)
                        diftru = S(mn2)
!
                        IF ( difest(2)==ZERO ) THEN
                           IF ( diftru>abnrm*ulp ) result(8) = ulpinv
                        ELSEIF ( diftru==ZERO ) THEN
                           IF ( difest(2)>abnrm*ulp ) result(8) = ulpinv
                        ELSEIF ( (diftru>thrsh2*difest(2)) .OR.         &
     &                           (diftru*thrsh2<difest(2)) ) THEN
                           result(8) = MAX(diftru/difest(2),difest(2)/  &
     &                                 diftru)
                        ENDIF
                        ntest = ntest + 1
                     ENDIF
!
!                 Test (9)
!
                     result(9) = ZERO
                     IF ( linfo==(MPLusn+2) ) THEN
                        IF ( diftru>abnrm*ulp ) result(9) = ulpinv
                        IF ( (ifunc>1) .AND. (difest(2)/=ZERO) )        &
     &                       result(9) = ulpinv
                        IF ( (ifunc==1) .AND. (pl(1)/=ZERO) ) result(9) &
     &                       = ulpinv
                        ntest = ntest + 1
                     ENDIF
!
                     ntestt = ntestt + ntest
!
!                 Print out tests which fail.
!
                     DO j = 1 , 9
                        IF ( result(j)>=Thresh ) THEN
!
!                       If this is the first test to fail,
!                       print a header to the data file.
!
                           IF ( nerrs==0 ) THEN
                              WRITE (Nout,FMT=99005) 'DGX'
!
!                          Matrix types
!
                              WRITE (Nout,FMT=99007)
!
!                          Tests performed
!
                              WRITE (Nout,FMT=99008) 'orthogonal' ,     &
     &                               '''' , 'transpose' , ('''',i=1,4)
!
                           ENDIF
                           nerrs = nerrs + 1
                           IF ( result(j)<10000.0D0 ) THEN
                              WRITE (Nout,FMT=99009) MPLusn , prtype ,  &
     &                               weight , M , j , result(j)
                           ELSE
                              WRITE (Nout,FMT=99010) MPLusn , prtype ,  &
     &                               weight , M , j , result(j)
                           ENDIF
                        ENDIF
                     ENDDO
!
                  ENDDO
               ENDDO
            ENDDO
!
         ENDDO
      ENDIF
!
!
!     Summary
!
 100  CALL alasvm('DGX',Nout,nerrs,ntestt,0)
!
      Work(1) = maxwrk
!
      RETURN
!
99001 FORMAT (' ddrgsx: ',A,' returned INFO=',I6,'.',/9X,'N=',I6,       &
     &        ', JTYPE=',I6,')')
!
99002 FORMAT (' ddrgsx: ',A,' returned INFO=',I6,'.',/9X,'N=',I6,       &
     &        ', Input Example #',I2,')')
!
99003 FORMAT (' ddrgsx: dget53 returned INFO=',I1,' for eigenvalue ',I6,&
     &        '.',/9X,'N=',I6,', JTYPE=',I6,')')
!
99004 FORMAT (' ddrgsx: S not in Schur form at eigenvalue ',I6,'.',/9X, &
     &        'N=',I6,', JTYPE=',I6,')')
!
99005 FORMAT (/1X,A3,' -- Real Expert Generalized Schur form',          &
     &        ' problem driver')
!
99006 FORMAT ('Input Example')
!
99007 FORMAT (' Matrix types: ',                                        &
     &        /'  1:  A is a block diagonal matrix of Jordan blocks ',  &
     &        'and B is the identity ',/'      matrix, ',               &
     &        /'  2:  A and B are upper triangular matrices, ',         &
     &        /'  3:  A and B are as type 2, but each second diagonal ',&
     &        'block in A_11 and ',                                     &
     &        /'      each third diaongal block in A_22 are 2x2 blocks,'&
     &        ,/'  4:  A and B are block diagonal matrices, ',          &
     &        /'  5:  (A,B) has potentially close or common ',          &
     &        'eigenvalues.',/)
!
99008 FORMAT (/' Tests performed:  (S is Schur, T is triangular, ',     &
     &        'Q and Z are ',A,',',/19X,' a is alpha, b is beta, and ', &
     &        A,' means ',A,'.)',/'  1 = | A - Q S Z',A,                &
     &        ' | / ( |A| n ulp )      2 = | B - Q T Z',A,              &
     &        ' | / ( |B| n ulp )',/'  3 = | I - QQ',A,                 &
     &        ' | / ( n ulp )             4 = | I - ZZ',A,              &
     &        ' | / ( n ulp )',/'  5 = 1/ULP  if A is not in ',         &
     &        'Schur form S',/'  6 = difference between (alpha,beta)',  &
     &        ' and diagonals of (S,T)',                                &
     &        /'  7 = 1/ULP  if SDIM is not the correct number of ',    &
     &        'selected eigenvalues',                                   &
     &        /'  8 = 1/ULP  if DIFEST/DIFTRU > 10*THRESH or ',         &
     &        'DIFTRU/DIFEST > 10*THRESH',                              &
     &        /'  9 = 1/ULP  if DIFEST <> 0 or DIFTRU > ULP*norm(A,B) ',&
     &        'when reordering fails',                                  &
     &        /' 10 = 1/ULP  if PLEST/PLTRU > THRESH or ',              &
     &        'PLTRU/PLEST > THRESH',                                   &
     &        /'    ( Test 10 is only for input examples )',/)
99009 FORMAT (' Matrix order=',I2,', type=',I2,', a=',D10.3,            &
     &        ', order(A_11)=',I2,', result ',I2,' is ',0P,F8.2)
99010 FORMAT (' Matrix order=',I2,', type=',I2,', a=',D10.3,            &
     &        ', order(A_11)=',I2,', result ',I2,' is ',0P,D10.3)
99011 FORMAT (' Input example #',I2,', matrix order=',I4,',',' result ',&
     &        I2,' is',0P,F8.2)
99012 FORMAT (' Input example #',I2,', matrix order=',I4,',',' result ',&
     &        I2,' is',1P,D10.3)
!
!     End of ddrgsx
!
      END SUBROUTINE ddrgsx
