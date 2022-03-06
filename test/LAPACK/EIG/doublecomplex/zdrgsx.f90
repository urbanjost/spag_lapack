!*==zdrgsx.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b ZDRGSX
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZDRGSX( NSIZE, NCMAX, THRESH, NIN, NOUT, A, LDA, B, AI,
!                          BI, Z, Q, ALPHA, BETA, C, LDC, S, WORK, LWORK,
!                          RWORK, IWORK, LIWORK, BWORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, LDC, LIWORK, LWORK, NCMAX, NIN,
!      $                   NOUT, NSIZE
!       DOUBLE PRECISION   THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            BWORK( * )
!       INTEGER            IWORK( * )
!       DOUBLE PRECISION   RWORK( * ), S( * )
!       COMPLEX*16         A( LDA, * ), AI( LDA, * ), ALPHA( * ),
!      $                   B( LDA, * ), BETA( * ), BI( LDA, * ),
!      $                   C( LDC, * ), Q( LDA, * ), WORK( * ),
!      $                   Z( LDA, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZDRGSX checks the nonsymmetric generalized eigenvalue (Schur form)
!> problem expert driver ZGGESX.
!>
!> ZGGES factors A and B as Q*S*Z'  and Q*T*Z' , where ' means conjugate
!> transpose, S and T are  upper triangular (i.e., in generalized Schur
!> form), and Q and Z are unitary. It also computes the generalized
!> eigenvalues (alpha(j),beta(j)), j=1,...,n.  Thus,
!> w(j) = alpha(j)/beta(j) is a root of the characteristic equation
!>
!>                 det( A - w(j) B ) = 0
!>
!> Optionally it also reorders the eigenvalues so that a selected
!> cluster of eigenvalues appears in the leading diagonal block of the
!> Schur forms; computes a reciprocal condition number for the average
!> of the selected eigenvalues; and computes a reciprocal condition
!> number for the right and left deflating subspaces corresponding to
!> the selected eigenvalues.
!>
!> When ZDRGSX is called with NSIZE > 0, five (5) types of built-in
!> matrix pairs are used to test the routine ZGGESX.
!>
!> When ZDRGSX is called with NSIZE = 0, it reads in test matrix data
!> to test ZGGESX.
!> (need more details on what kind of read-in data are needed).
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
!> (5)   if A is in Schur form (i.e. triangular form)
!>
!> (6)   maximum over j of D(j)  where:
!>
!>                     |alpha(j) - S(j,j)|        |beta(j) - T(j,j)|
!>           D(j) = ------------------------ + -----------------------
!>                  max(|alpha(j)|,|S(j,j)|)   max(|beta(j)|,|T(j,j)|)
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
!> (9)   If INFO = N+3 is returned by ZGGESX, the reordering "failed"
!>       and we check that DIF = PL = PR = 0 and that the true value of
!>       Difu and Difl is < EPS*norm(A, B). We count the events when
!>       INFO=N+3.
!>
!> For read-in test matrices, the same tests are run except that the
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
!> SVD (routine ZGESVD) is used for computing the true value of DIF_u
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
!>          and matrix B are chosen as identity matrices (see DLATM5).
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
!>          in call to ZLAKF2
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
!>          (e.g., if a routine returns INFO not equal to 0.)
!> \endverbatim
!>
!> \param[out] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA, NSIZE)
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
!>          B is COMPLEX*16 array, dimension (LDA, NSIZE)
!>          Used to store the matrix whose eigenvalues are to be
!>          computed.  On exit, B contains the last matrix actually used.
!> \endverbatim
!>
!> \param[out] AI
!> \verbatim
!>          AI is COMPLEX*16 array, dimension (LDA, NSIZE)
!>          Copy of A, modified by ZGGESX.
!> \endverbatim
!>
!> \param[out] BI
!> \verbatim
!>          BI is COMPLEX*16 array, dimension (LDA, NSIZE)
!>          Copy of B, modified by ZGGESX.
!> \endverbatim
!>
!> \param[out] Z
!> \verbatim
!>          Z is COMPLEX*16 array, dimension (LDA, NSIZE)
!>          Z holds the left Schur vectors computed by ZGGESX.
!> \endverbatim
!>
!> \param[out] Q
!> \verbatim
!>          Q is COMPLEX*16 array, dimension (LDA, NSIZE)
!>          Q holds the right Schur vectors computed by ZGGESX.
!> \endverbatim
!>
!> \param[out] ALPHA
!> \verbatim
!>          ALPHA is COMPLEX*16 array, dimension (NSIZE)
!> \endverbatim
!>
!> \param[out] BETA
!> \verbatim
!>          BETA is COMPLEX*16 array, dimension (NSIZE)
!>
!>          On exit, ALPHA/BETA are the eigenvalues.
!> \endverbatim
!>
!> \param[out] C
!> \verbatim
!>          C is COMPLEX*16 array, dimension (LDC, LDC)
!>          Store the matrix generated by subroutine ZLAKF2, this is the
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
!>          WORK is COMPLEX*16 array, dimension (LWORK)
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK.  LWORK >= 3*NSIZE*NSIZE/2
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array,
!>                                 dimension (5*NSIZE*NSIZE/2 - 4)
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
!>          The dimension of the array IWORK. LIWORK >= NSIZE + 2.
!> \endverbatim
!>
!> \param[out] BWORK
!> \verbatim
!>          BWORK is LOGICAL array, dimension (NSIZE)
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
!> \ingroup complex16_eig
!
!  =====================================================================
      SUBROUTINE ZDRGSX(Nsize,Ncmax,Thresh,Nin,Nout,A,Lda,B,Ai,Bi,Z,Q,  &
     &                  Alpha,Beta,C,Ldc,S,Work,Lwork,Rwork,Iwork,      &
     &                  Liwork,Bwork,Info)
      IMPLICIT NONE
!*--ZDRGSX353
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
      DOUBLE PRECISION Rwork(*) , S(*)
      COMPLEX*16 A(Lda,*) , Ai(Lda,*) , Alpha(*) , B(Lda,*) , Beta(*) , &
     &           Bi(Lda,*) , C(Ldc,*) , Q(Lda,*) , Work(*) , Z(Lda,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE , TEN
      PARAMETER (ZERO=0.0D+0,ONE=1.0D+0,TEN=1.0D+1)
      COMPLEX*16 CZERO
      PARAMETER (CZERO=(0.0D+0,0.0D+0))
!     ..
!     .. Local Scalars ..
      LOGICAL ilabad
      CHARACTER sense
      INTEGER bdspac , i , ifunc , j , linfo , maxwrk , minwrk , mm ,   &
     &        mn2 , nerrs , nptknt , ntest , ntestt , prtype , qba , qbb
      DOUBLE PRECISION abnrm , bignum , diftru , pltru , smlnum ,       &
     &                 temp1 , temp2 , thrsh2 , ulp , ulpinv , weight
      COMPLEX*16 x
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION difest(2) , pl(2) , result(10)
!     ..
!     .. External Functions ..
      LOGICAL ZLCTSX
      INTEGER ILAENV
      DOUBLE PRECISION DLAMCH , ZLANGE
      EXTERNAL ZLCTSX , ILAENV , DLAMCH , ZLANGE
!     ..
!     .. External Subroutines ..
      EXTERNAL ALASVM , DLABAD , XERBLA , ZGESVD , ZGET51 , ZGGESX ,    &
     &         ZLACPY , ZLAKF2 , ZLASET , ZLATM5
!     ..
!     .. Scalars in Common ..
      LOGICAL FS
      INTEGER K , M , MPLusn , N
!     ..
!     .. Common blocks ..
      COMMON /MN    / M , N , MPLusn , K , FS
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , DBLE , DIMAG , MAX , SQRT
!     ..
!     .. Statement Functions ..
      DOUBLE PRECISION ABS1
!     ..
!     .. Statement Function definitions ..
      ABS1(x) = ABS(DBLE(x)) + ABS(DIMAG(x))
!     ..
!     .. Executable Statements ..
!
!     Check for errors
!
      Info = 0
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
         Info = -15
      ELSEIF ( Liwork<Nsize+2 ) THEN
         Info = -21
      ENDIF
!
!     Compute workspace
!      (Note: Comments in the code beginning "Workspace:" describe the
!       minimal amount of workspace needed at that point in the code,
!       as well as the preferred amount for good performance.
!       NB refers to the optimal block size for the immediately
!       following subroutine, as returned by ILAENV.)
!
      minwrk = 1
      IF ( Info==0 .AND. Lwork>=1 ) THEN
         minwrk = 3*Nsize*Nsize/2
!
!        workspace for cggesx
!
         maxwrk = Nsize*(1+ILAENV(1,'ZGEQRF',' ',Nsize,1,Nsize,0))
         maxwrk = MAX(maxwrk,                                           &
     &            Nsize*(1+ILAENV(1,'ZUNGQR',' ',Nsize,1,Nsize,-1)))
!
!        workspace for zgesvd
!
         bdspac = 3*Nsize*Nsize/2
         maxwrk = MAX(maxwrk,                                           &
     &            Nsize*Nsize*(1+ILAENV(1,'ZGEBRD',' ',Nsize*Nsize/2,   &
     &            Nsize*Nsize/2,-1,-1)))
         maxwrk = MAX(maxwrk,bdspac)
!
         maxwrk = MAX(maxwrk,minwrk)
!
         Work(1) = maxwrk
      ENDIF
!
      IF ( Lwork<minwrk ) Info = -18
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('ZDRGSX',-Info)
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
            CALL ZLACPY('Full',MPLusn,MPLusn,Ai,Lda,A,Lda)
            CALL ZLACPY('Full',MPLusn,MPLusn,Bi,Lda,B,Lda)
!
!     Compute the Schur factorization while swapping the
!     m-by-m (1,1)-blocks with n-by-n (2,2)-blocks.
!
            CALL ZGGESX('V','V','S',ZLCTSX,'B',MPLusn,Ai,Lda,Bi,Lda,mm, &
     &                  Alpha,Beta,Q,Lda,Z,Lda,pl,difest,Work,Lwork,    &
     &                  Rwork,Iwork,Liwork,Bwork,linfo)
!
            IF ( linfo/=0 .AND. linfo/=MPLusn+2 ) THEN
               result(1) = ulpinv
               WRITE (Nout,FMT=99002) 'ZGGESX' , linfo , MPLusn , nptknt
               CYCLE
            ENDIF
!
!     Compute the norm(A, B)
!        (should this be norm of (A,B) or (AI,BI)?)
!
            CALL ZLACPY('Full',MPLusn,MPLusn,Ai,Lda,Work,MPLusn)
            CALL ZLACPY('Full',MPLusn,MPLusn,Bi,Lda,                    &
     &                  Work(MPLusn*MPLusn+1),MPLusn)
            abnrm = ZLANGE('Fro',MPLusn,2*MPLusn,Work,MPLusn,Rwork)
!
!     Do tests (1) to (4)
!
            CALL ZGET51(1,MPLusn,A,Lda,Ai,Lda,Q,Lda,Z,Lda,Work,Rwork,   &
     &                  result(1))
            CALL ZGET51(1,MPLusn,B,Lda,Bi,Lda,Q,Lda,Z,Lda,Work,Rwork,   &
     &                  result(2))
            CALL ZGET51(3,MPLusn,B,Lda,Bi,Lda,Q,Lda,Q,Lda,Work,Rwork,   &
     &                  result(3))
            CALL ZGET51(3,MPLusn,B,Lda,Bi,Lda,Z,Lda,Z,Lda,Work,Rwork,   &
     &                  result(4))
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
               temp2 = (ABS1(Alpha(j)-Ai(j,j))                          &
     &                 /MAX(smlnum,ABS1(Alpha(j)),ABS1(Ai(j,j)))        &
     &                 +ABS1(Beta(j)-Bi(j,j))                           &
     &                 /MAX(smlnum,ABS1(Beta(j)),ABS1(Bi(j,j))))/ulp
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
               temp1 = MAX(temp1,temp2)
               IF ( ilabad ) WRITE (Nout,FMT=99003) j , MPLusn , nptknt
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
                     WRITE (Nout,FMT=99004) 'ZGX'
!
!              Matrix types
!
                     WRITE (Nout,FMT=99005)
!
!              Tests performed
!
                     WRITE (Nout,FMT=99007) 'unitary' , '''' ,          &
     &                      'transpose' , ('''',i=1,4)
!
                  ENDIF
                  nerrs = nerrs + 1
                  IF ( result(j)<10000.0D0 ) THEN
                     WRITE (Nout,FMT=99010) nptknt , MPLusn , j ,       &
     &                      result(j)
                  ELSE
                     WRITE (Nout,FMT=99011) nptknt , MPLusn , j ,       &
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
!     Loop over different functions (IFUNC) of ZGGESX, types (PRTYPE)
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
                     CALL ZLASET('Full',MPLusn,MPLusn,CZERO,CZERO,Ai,   &
     &                           Lda)
                     CALL ZLASET('Full',MPLusn,MPLusn,CZERO,CZERO,Bi,   &
     &                           Lda)
!
                     CALL ZLATM5(prtype,M,N,Ai,Lda,Ai(M+1,M+1),Lda,     &
     &                           Ai(1,M+1),Lda,Bi,Lda,Bi(M+1,M+1),Lda,  &
     &                           Bi(1,M+1),Lda,Q,Lda,Z,Lda,weight,qba,  &
     &                           qbb)
!
!                 Compute the Schur factorization and swapping the
!                 m-by-m (1,1)-blocks with n-by-n (2,2)-blocks.
!                 Swapping is accomplished via the function ZLCTSX
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
                     CALL ZLACPY('Full',MPLusn,MPLusn,Ai,Lda,A,Lda)
                     CALL ZLACPY('Full',MPLusn,MPLusn,Bi,Lda,B,Lda)
!
                     CALL ZGGESX('V','V','S',ZLCTSX,sense,MPLusn,Ai,Lda,&
     &                           Bi,Lda,mm,Alpha,Beta,Q,Lda,Z,Lda,pl,   &
     &                           difest,Work,Lwork,Rwork,Iwork,Liwork,  &
     &                           Bwork,linfo)
!
                     IF ( linfo/=0 .AND. linfo/=MPLusn+2 ) THEN
                        result(1) = ulpinv
                        WRITE (Nout,FMT=99001) 'ZGGESX' , linfo ,       &
     &                         MPLusn , prtype
                        Info = linfo
                        CYCLE
                     ENDIF
!
!                 Compute the norm(A, B)
!
                     CALL ZLACPY('Full',MPLusn,MPLusn,Ai,Lda,Work,      &
     &                           MPLusn)
                     CALL ZLACPY('Full',MPLusn,MPLusn,Bi,Lda,           &
     &                           Work(MPLusn*MPLusn+1),MPLusn)
                     abnrm = ZLANGE('Fro',MPLusn,2*MPLusn,Work,MPLusn,  &
     &                       Rwork)
!
!                 Do tests (1) to (4)
!
                     result(2) = ZERO
                     CALL ZGET51(1,MPLusn,A,Lda,Ai,Lda,Q,Lda,Z,Lda,Work,&
     &                           Rwork,result(1))
                     CALL ZGET51(1,MPLusn,B,Lda,Bi,Lda,Q,Lda,Z,Lda,Work,&
     &                           Rwork,result(2))
                     CALL ZGET51(3,MPLusn,B,Lda,Bi,Lda,Q,Lda,Q,Lda,Work,&
     &                           Rwork,result(3))
                     CALL ZGET51(3,MPLusn,B,Lda,Bi,Lda,Z,Lda,Z,Lda,Work,&
     &                           Rwork,result(4))
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
                        temp2 = (ABS1(Alpha(j)-Ai(j,j))                 &
     &                          /MAX(smlnum,ABS1(Alpha(j)),ABS1(Ai(j,j))&
     &                          )+ABS1(Beta(j)-Bi(j,j))                 &
     &                          /MAX(smlnum,ABS1(Beta(j)),ABS1(Bi(j,j)))&
     &                          )/ulp
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
                        temp1 = MAX(temp1,temp2)
                        IF ( ilabad ) WRITE (Nout,FMT=99003) j ,        &
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
!                    Note: for either following two cases, there are
!                    almost same number of test cases fail the test.
!
                        CALL ZLAKF2(mm,MPLusn-mm,Ai,Lda,Ai(mm+1,mm+1),  &
     &                              Bi,Bi(mm+1,mm+1),C,Ldc)
!
                        CALL ZGESVD('N','N',mn2,mn2,C,Ldc,S,Work,1,     &
     &                              Work(2),1,Work(3),Lwork-2,Rwork,    &
     &                              Info)
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
                              WRITE (Nout,FMT=99004) 'ZGX'
!
!                          Matrix types
!
                              WRITE (Nout,FMT=99006)
!
!                          Tests performed
!
                              WRITE (Nout,FMT=99007) 'unitary' , '''' , &
     &                               'transpose' , ('''',i=1,4)
!
                           ENDIF
                           nerrs = nerrs + 1
                           IF ( result(j)<10000.0D0 ) THEN
                              WRITE (Nout,FMT=99008) MPLusn , prtype ,  &
     &                               weight , M , j , result(j)
                           ELSE
                              WRITE (Nout,FMT=99009) MPLusn , prtype ,  &
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
 100  CALL ALASVM('ZGX',Nout,nerrs,ntestt,0)
!
      Work(1) = maxwrk
!
      RETURN
!
99001 FORMAT (' ZDRGSX: ',A,' returned INFO=',I6,'.',/9X,'N=',I6,       &
     &        ', JTYPE=',I6,')')
!
99002 FORMAT (' ZDRGSX: ',A,' returned INFO=',I6,'.',/9X,'N=',I6,       &
     &        ', Input Example #',I2,')')
!
99003 FORMAT (' ZDRGSX: S not in Schur form at eigenvalue ',I6,'.',/9X, &
     &        'N=',I6,', JTYPE=',I6,')')
!
99004 FORMAT (/1X,A3,' -- Complex Expert Generalized Schur form',       &
     &        ' problem driver')
!
99005 FORMAT ('Input Example')
!
99006 FORMAT (' Matrix types: ',                                        &
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
99007 FORMAT (/' Tests performed:  (S is Schur, T is triangular, ',     &
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
99008 FORMAT (' Matrix order=',I2,', type=',I2,', a=',D10.3,            &
     &        ', order(A_11)=',I2,', result ',I2,' is ',0P,F8.2)
99009 FORMAT (' Matrix order=',I2,', type=',I2,', a=',D10.3,            &
     &        ', order(A_11)=',I2,', result ',I2,' is ',0P,D10.3)
99010 FORMAT (' Input example #',I2,', matrix order=',I4,',',' result ',&
     &        I2,' is',0P,F8.2)
99011 FORMAT (' Input example #',I2,', matrix order=',I4,',',' result ',&
     &        I2,' is',1P,D10.3)
!
!     End of ZDRGSX
!
      END SUBROUTINE ZDRGSX
