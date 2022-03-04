!*==cchkhb.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b CCHKHB
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CCHKHB( NSIZES, NN, NWDTHS, KK, NTYPES, DOTYPE, ISEED,
!                          THRESH, NOUNIT, A, LDA, SD, SE, U, LDU, WORK,
!                          LWORK, RWORK, RESULT, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, LDU, LWORK, NOUNIT, NSIZES, NTYPES,
!      $                   NWDTHS
!       REAL               THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            DOTYPE( * )
!       INTEGER            ISEED( 4 ), KK( * ), NN( * )
!       REAL               RESULT( * ), RWORK( * ), SD( * ), SE( * )
!       COMPLEX            A( LDA, * ), U( LDU, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CCHKHB tests the reduction of a Hermitian band matrix to tridiagonal
!> from, used with the Hermitian eigenvalue problem.
!>
!> CHBTRD factors a Hermitian band matrix A as  U S U* , where * means
!> conjugate transpose, S is symmetric tridiagonal, and U is unitary.
!> CHBTRD can use either just the lower or just the upper triangle
!> of A; CCHKHB checks both cases.
!>
!> When CCHKHB is called, a number of matrix "sizes" ("n's"), a number
!> of bandwidths ("k's"), and a number of matrix "types" are
!> specified.  For each size ("n"), each bandwidth ("k") less than or
!> equal to "n", and each type of matrix, one matrix will be generated
!> and used to test the hermitian banded reduction routine.  For each
!> matrix, a number of tests will be performed:
!>
!> (1)     | A - V S V* | / ( |A| n ulp )  computed by CHBTRD with
!>                                         UPLO='U'
!>
!> (2)     | I - UU* | / ( n ulp )
!>
!> (3)     | A - V S V* | / ( |A| n ulp )  computed by CHBTRD with
!>                                         UPLO='L'
!>
!> (4)     | I - UU* | / ( n ulp )
!>
!> The "sizes" are specified by an array NN(1:NSIZES); the value of
!> each element NN(j) specifies one size.
!> The "types" are specified by a logical array DOTYPE( 1:NTYPES );
!> if DOTYPE(j) is .TRUE., then matrix type "j" will be generated.
!> Currently, the list of possible types is:
!>
!> (1)  The zero matrix.
!> (2)  The identity matrix.
!>
!> (3)  A diagonal matrix with evenly spaced entries
!>      1, ..., ULP  and random signs.
!>      (ULP = (first number larger than 1) - 1 )
!> (4)  A diagonal matrix with geometrically spaced entries
!>      1, ..., ULP  and random signs.
!> (5)  A diagonal matrix with "clustered" entries 1, ULP, ..., ULP
!>      and random signs.
!>
!> (6)  Same as (4), but multiplied by SQRT( overflow threshold )
!> (7)  Same as (4), but multiplied by SQRT( underflow threshold )
!>
!> (8)  A matrix of the form  U* D U, where U is unitary and
!>      D has evenly spaced entries 1, ..., ULP with random signs
!>      on the diagonal.
!>
!> (9)  A matrix of the form  U* D U, where U is unitary and
!>      D has geometrically spaced entries 1, ..., ULP with random
!>      signs on the diagonal.
!>
!> (10) A matrix of the form  U* D U, where U is unitary and
!>      D has "clustered" entries 1, ULP,..., ULP with random
!>      signs on the diagonal.
!>
!> (11) Same as (8), but multiplied by SQRT( overflow threshold )
!> (12) Same as (8), but multiplied by SQRT( underflow threshold )
!>
!> (13) Hermitian matrix with random entries chosen from (-1,1).
!> (14) Same as (13), but multiplied by SQRT( overflow threshold )
!> (15) Same as (13), but multiplied by SQRT( underflow threshold )
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] NSIZES
!> \verbatim
!>          NSIZES is INTEGER
!>          The number of sizes of matrices to use.  If it is zero,
!>          CCHKHB does nothing.  It must be at least zero.
!> \endverbatim
!>
!> \param[in] NN
!> \verbatim
!>          NN is INTEGER array, dimension (NSIZES)
!>          An array containing the sizes to be used for the matrices.
!>          Zero values will be skipped.  The values must be at least
!>          zero.
!> \endverbatim
!>
!> \param[in] NWDTHS
!> \verbatim
!>          NWDTHS is INTEGER
!>          The number of bandwidths to use.  If it is zero,
!>          CCHKHB does nothing.  It must be at least zero.
!> \endverbatim
!>
!> \param[in] KK
!> \verbatim
!>          KK is INTEGER array, dimension (NWDTHS)
!>          An array containing the bandwidths to be used for the band
!>          matrices.  The values must be at least zero.
!> \endverbatim
!>
!> \param[in] NTYPES
!> \verbatim
!>          NTYPES is INTEGER
!>          The number of elements in DOTYPE.   If it is zero, CCHKHB
!>          does nothing.  It must be at least zero.  If it is MAXTYP+1
!>          and NSIZES is 1, then an additional type, MAXTYP+1 is
!>          defined, which is to use whatever matrix is in A.  This
!>          is only useful if DOTYPE(1:MAXTYP) is .FALSE. and
!>          DOTYPE(MAXTYP+1) is .TRUE. .
!> \endverbatim
!>
!> \param[in] DOTYPE
!> \verbatim
!>          DOTYPE is LOGICAL array, dimension (NTYPES)
!>          If DOTYPE(j) is .TRUE., then for each size in NN a
!>          matrix of that size and of type j will be generated.
!>          If NTYPES is smaller than the maximum number of types
!>          defined (PARAMETER MAXTYP), then types NTYPES+1 through
!>          MAXTYP will not be generated.  If NTYPES is larger
!>          than MAXTYP, DOTYPE(MAXTYP+1) through DOTYPE(NTYPES)
!>          will be ignored.
!> \endverbatim
!>
!> \param[in,out] ISEED
!> \verbatim
!>          ISEED is INTEGER array, dimension (4)
!>          On entry ISEED specifies the seed of the random number
!>          generator. The array elements should be between 0 and 4095;
!>          if not they will be reduced mod 4096.  Also, ISEED(4) must
!>          be odd.  The random number generator uses a linear
!>          congruential sequence limited to small integers, and so
!>          should produce machine independent random numbers. The
!>          values of ISEED are changed on exit, and can be used in the
!>          next call to CCHKHB to continue the same random number
!>          sequence.
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
!>          (e.g., if a routine returns IINFO not equal to 0.)
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX array, dimension
!>                            (LDA, max(NN))
!>          Used to hold the matrix whose eigenvalues are to be
!>          computed.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of A.  It must be at least 2 (not 1!)
!>          and at least max( KK )+1.
!> \endverbatim
!>
!> \param[out] SD
!> \verbatim
!>          SD is REAL array, dimension (max(NN))
!>          Used to hold the diagonal of the tridiagonal matrix computed
!>          by CHBTRD.
!> \endverbatim
!>
!> \param[out] SE
!> \verbatim
!>          SE is REAL array, dimension (max(NN))
!>          Used to hold the off-diagonal of the tridiagonal matrix
!>          computed by CHBTRD.
!> \endverbatim
!>
!> \param[out] U
!> \verbatim
!>          U is COMPLEX array, dimension (LDU, max(NN))
!>          Used to hold the unitary matrix computed by CHBTRD.
!> \endverbatim
!>
!> \param[in] LDU
!> \verbatim
!>          LDU is INTEGER
!>          The leading dimension of U.  It must be at least 1
!>          and at least max( NN ).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension (LWORK)
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The number of entries in WORK.  This must be at least
!>          max( LDA+1, max(NN)+1 )*max(NN).
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array
!> \endverbatim
!>
!> \param[out] RESULT
!> \verbatim
!>          RESULT is REAL array, dimension (4)
!>          The values computed by the tests described above.
!>          The values are currently limited to 1/ulp, to avoid
!>          overflow.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          If 0, then everything ran OK.
!>
!>-----------------------------------------------------------------------
!>
!>       Some Local Variables and Parameters:
!>       ---- ----- --------- --- ----------
!>       ZERO, ONE       Real 0 and 1.
!>       MAXTYP          The number of types defined.
!>       NTEST           The number of tests performed, or which can
!>                       be performed so far, for the current matrix.
!>       NTESTT          The total number of tests performed so far.
!>       NMAX            Largest value in NN.
!>       NMATS           The number of matrices generated so far.
!>       NERRS           The number of tests which have exceeded THRESH
!>                       so far.
!>       COND, IMODE     Values to be passed to the matrix generators.
!>       ANORM           Norm of A; passed to matrix generators.
!>
!>       OVFL, UNFL      Overflow and underflow thresholds.
!>       ULP, ULPINV     Finest relative precision and its inverse.
!>       RTOVFL, RTUNFL  Square roots of the previous 2 values.
!>               The following four arrays decode JTYPE:
!>       KTYPE(j)        The general type (1-10) for type "j".
!>       KMODE(j)        The MODE value to be passed to the matrix
!>                       generator for type "j".
!>       KMAGN(j)        The order of magnitude ( O(1),
!>                       O(overflow^(1/2) ), O(underflow^(1/2) )
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
!> \ingroup complex_eig
!
!  =====================================================================
      SUBROUTINE CCHKHB(Nsizes,Nn,Nwdths,Kk,Ntypes,Dotype,Iseed,Thresh, &
     &                  Nounit,A,Lda,Sd,Se,U,Ldu,Work,Lwork,Rwork,      &
     &                  Result,Info)
      IMPLICIT NONE
!*--CCHKHB302
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Info , Lda , Ldu , Lwork , Nounit , Nsizes , Ntypes ,     &
     &        Nwdths
      REAL Thresh
!     ..
!     .. Array Arguments ..
      LOGICAL Dotype(*)
      INTEGER Iseed(4) , Kk(*) , Nn(*)
      REAL Result(*) , Rwork(*) , Sd(*) , Se(*)
      COMPLEX A(Lda,*) , U(Ldu,*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX CZERO , CONE
      PARAMETER (CZERO=(0.0E+0,0.0E+0),CONE=(1.0E+0,0.0E+0))
      REAL ZERO , ONE , TWO , TEN
      PARAMETER (ZERO=0.0E+0,ONE=1.0E+0,TWO=2.0E+0,TEN=10.0E+0)
      REAL HALF
      PARAMETER (HALF=ONE/TWO)
      INTEGER MAXTYP
      PARAMETER (MAXTYP=15)
!     ..
!     .. Local Scalars ..
      LOGICAL badnn , badnnb
      INTEGER i , iinfo , imode , itype , j , jc , jcol , jr , jsize ,  &
     &        jtype , jwidth , k , kmax , mtypes , n , nerrs , nmats ,  &
     &        nmax , ntest , ntestt
      REAL aninv , anorm , cond , ovfl , rtovfl , rtunfl , temp1 , ulp ,&
     &     ulpinv , unfl
!     ..
!     .. Local Arrays ..
      INTEGER idumma(1) , ioldsd(4) , kmagn(MAXTYP) , kmode(MAXTYP) ,   &
     &        ktype(MAXTYP)
!     ..
!     .. External Functions ..
      REAL SLAMCH
      EXTERNAL SLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL CHBT21 , CHBTRD , CLACPY , CLATMR , CLATMS , CLASET ,    &
     &         SLASUM , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , CONJG , MAX , MIN , REAL , SQRT
!     ..
!     .. Data statements ..
      DATA ktype/1 , 2 , 5*4 , 5*5 , 3*8/
      DATA kmagn/2*1 , 1 , 1 , 1 , 2 , 3 , 1 , 1 , 1 , 2 , 3 , 1 , 2 ,  &
     &     3/
      DATA kmode/2*0 , 4 , 3 , 1 , 4 , 4 , 4 , 3 , 1 , 4 , 4 , 0 , 0 ,  &
     &     0/
!     ..
!     .. Executable Statements ..
!
!     Check for errors
!
      ntestt = 0
      Info = 0
!
!     Important constants
!
      badnn = .FALSE.
      nmax = 1
      DO j = 1 , Nsizes
         nmax = MAX(nmax,Nn(j))
         IF ( Nn(j)<0 ) badnn = .TRUE.
      ENDDO
!
      badnnb = .FALSE.
      kmax = 0
      DO j = 1 , Nsizes
         kmax = MAX(kmax,Kk(j))
         IF ( Kk(j)<0 ) badnnb = .TRUE.
      ENDDO
      kmax = MIN(nmax-1,kmax)
!
!     Check for errors
!
      IF ( Nsizes<0 ) THEN
         Info = -1
      ELSEIF ( badnn ) THEN
         Info = -2
      ELSEIF ( Nwdths<0 ) THEN
         Info = -3
      ELSEIF ( badnnb ) THEN
         Info = -4
      ELSEIF ( Ntypes<0 ) THEN
         Info = -5
      ELSEIF ( Lda<kmax+1 ) THEN
         Info = -11
      ELSEIF ( Ldu<nmax ) THEN
         Info = -15
      ELSEIF ( (MAX(Lda,nmax)+1)*nmax>Lwork ) THEN
         Info = -17
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('CCHKHB',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( Nsizes==0 .OR. Ntypes==0 .OR. Nwdths==0 ) RETURN
!
!     More Important constants
!
      unfl = SLAMCH('Safe minimum')
      ovfl = ONE/unfl
      ulp = SLAMCH('Epsilon')*SLAMCH('Base')
      ulpinv = ONE/ulp
      rtunfl = SQRT(unfl)
      rtovfl = SQRT(ovfl)
!
!     Loop over sizes, types
!
      nerrs = 0
      nmats = 0
!
      DO jsize = 1 , Nsizes
         n = Nn(jsize)
         aninv = ONE/REAL(MAX(1,n))
!
         DO jwidth = 1 , Nwdths
            k = Kk(jwidth)
            IF ( k<=n ) THEN
               k = MAX(0,MIN(n-1,k))
!
               IF ( Nsizes/=1 ) THEN
                  mtypes = MIN(MAXTYP,Ntypes)
               ELSE
                  mtypes = MIN(MAXTYP+1,Ntypes)
               ENDIF
!
               DO jtype = 1 , mtypes
                  IF ( Dotype(jtype) ) THEN
                     nmats = nmats + 1
                     ntest = 0
!
                     DO j = 1 , 4
                        ioldsd(j) = Iseed(j)
                     ENDDO
!
!              Compute "A".
!              Store as "Upper"; later, we will copy to other format.
!
!              Control parameters:
!
!                  KMAGN  KMODE        KTYPE
!              =1  O(1)   clustered 1  zero
!              =2  large  clustered 2  identity
!              =3  small  exponential  (none)
!              =4         arithmetic   diagonal, (w/ eigenvalues)
!              =5         random log   hermitian, w/ eigenvalues
!              =6         random       (none)
!              =7                      random diagonal
!              =8                      random hermitian
!              =9                      positive definite
!              =10                     diagonally dominant tridiagonal
!
                     IF ( mtypes<=MAXTYP ) THEN
!
                        itype = ktype(jtype)
                        imode = kmode(jtype)
!
!              Compute norm
!
                        IF ( kmagn(jtype)==2 ) THEN
!
                           anorm = (rtovfl*ulp)*aninv
                        ELSEIF ( kmagn(jtype)==3 ) THEN
!
                           anorm = rtunfl*n*ulpinv
                        ELSE
!
                           anorm = ONE
                        ENDIF
!
!
                        CALL CLASET('Full',Lda,n,CZERO,CZERO,A,Lda)
                        iinfo = 0
                        IF ( jtype<=15 ) THEN
                           cond = ulpinv
                        ELSE
                           cond = ulpinv*aninv/TEN
                        ENDIF
!
!              Special Matrices -- Identity & Jordan block
!
!                 Zero
!
                        IF ( itype==1 ) THEN
                           iinfo = 0
!
                        ELSEIF ( itype==2 ) THEN
!
!                 Identity
!
                           DO jcol = 1 , n
                              A(k+1,jcol) = anorm
                           ENDDO
!
                        ELSEIF ( itype==4 ) THEN
!
!                 Diagonal Matrix, [Eigen]values Specified
!
                           CALL CLATMS(n,n,'S',Iseed,'H',Rwork,imode,   &
     &                                 cond,anorm,0,0,'Q',A(k+1,1),Lda, &
     &                                 Work,iinfo)
!
                        ELSEIF ( itype==5 ) THEN
!
!                 Hermitian, eigenvalues specified
!
                           CALL CLATMS(n,n,'S',Iseed,'H',Rwork,imode,   &
     &                                 cond,anorm,k,k,'Q',A,Lda,Work,   &
     &                                 iinfo)
!
                        ELSEIF ( itype==7 ) THEN
!
!                 Diagonal, random eigenvalues
!
                           CALL CLATMR(n,n,'S',Iseed,'H',Work,6,ONE,    &
     &                                 CONE,'T','N',Work(n+1),1,ONE,    &
     &                                 Work(2*n+1),1,ONE,'N',idumma,0,0,&
     &                                 ZERO,anorm,'Q',A(k+1,1),Lda,     &
     &                                 idumma,iinfo)
!
                        ELSEIF ( itype==8 ) THEN
!
!                 Hermitian, random eigenvalues
!
                           CALL CLATMR(n,n,'S',Iseed,'H',Work,6,ONE,    &
     &                                 CONE,'T','N',Work(n+1),1,ONE,    &
     &                                 Work(2*n+1),1,ONE,'N',idumma,k,k,&
     &                                 ZERO,anorm,'Q',A,Lda,idumma,     &
     &                                 iinfo)
!
                        ELSEIF ( itype==9 ) THEN
!
!                 Positive definite, eigenvalues specified.
!
                           CALL CLATMS(n,n,'S',Iseed,'P',Rwork,imode,   &
     &                                 cond,anorm,k,k,'Q',A,Lda,        &
     &                                 Work(n+1),iinfo)
!
                        ELSEIF ( itype==10 ) THEN
!
!                 Positive definite tridiagonal, eigenvalues specified.
!
                           IF ( n>1 ) k = MAX(1,k)
                           CALL CLATMS(n,n,'S',Iseed,'P',Rwork,imode,   &
     &                                 cond,anorm,1,1,'Q',A(k,1),Lda,   &
     &                                 Work,iinfo)
                           DO i = 2 , n
                              temp1 = ABS(A(k,i))                       &
     &                                /SQRT(ABS(A(k+1,i-1)*A(k+1,i)))
                              IF ( temp1>HALF ) A(k,i)                  &
     &                             = HALF*SQRT(ABS(A(k+1,i-1)*A(k+1,i)))
                           ENDDO
!
                        ELSE
!
                           iinfo = 1
                        ENDIF
!
                        IF ( iinfo/=0 ) THEN
                           WRITE (Nounit,FMT=99001) 'Generator' ,       &
     &                            iinfo , n , jtype , ioldsd
                           Info = ABS(iinfo)
                           RETURN
                        ENDIF
                     ENDIF
!
!
!              Call CHBTRD to compute S and U from upper triangle.
!
                     CALL CLACPY(' ',k+1,n,A,Lda,Work,Lda)
!
                     ntest = 1
                     CALL CHBTRD('V','U',n,k,Work,Lda,Sd,Se,U,Ldu,      &
     &                           Work(Lda*n+1),iinfo)
!
                     IF ( iinfo/=0 ) THEN
                        WRITE (Nounit,FMT=99001) 'CHBTRD(U)' , iinfo ,  &
     &                         n , jtype , ioldsd
                        Info = ABS(iinfo)
                        IF ( iinfo<0 ) THEN
                           RETURN
                        ELSE
                           Result(1) = ulpinv
                           GOTO 2
                        ENDIF
                     ENDIF
!
!              Do tests 1 and 2
!
                     CALL CHBT21('Upper',n,k,1,A,Lda,Sd,Se,U,Ldu,Work,  &
     &                           Rwork,Result(1))
!
!              Convert A from Upper-Triangle-Only storage to
!              Lower-Triangle-Only storage.
!
                     DO jc = 1 , n
                        DO jr = 0 , MIN(k,n-jc)
                           A(jr+1,jc) = CONJG(A(k+1-jr,jc+jr))
                        ENDDO
                     ENDDO
                     DO jc = n + 1 - k , n
                        DO jr = MIN(k,n-jc) + 1 , k
                           A(jr+1,jc) = ZERO
                        ENDDO
                     ENDDO
!
!              Call CHBTRD to compute S and U from lower triangle
!
                     CALL CLACPY(' ',k+1,n,A,Lda,Work,Lda)
!
                     ntest = 3
                     CALL CHBTRD('V','L',n,k,Work,Lda,Sd,Se,U,Ldu,      &
     &                           Work(Lda*n+1),iinfo)
!
                     IF ( iinfo/=0 ) THEN
                        WRITE (Nounit,FMT=99001) 'CHBTRD(L)' , iinfo ,  &
     &                         n , jtype , ioldsd
                        Info = ABS(iinfo)
                        IF ( iinfo<0 ) THEN
                           RETURN
                        ELSE
                           Result(3) = ulpinv
                           GOTO 2
                        ENDIF
                     ENDIF
                     ntest = 4
!
!              Do tests 3 and 4
!
                     CALL CHBT21('Lower',n,k,1,A,Lda,Sd,Se,U,Ldu,Work,  &
     &                           Rwork,Result(3))
!
!              End of Loop -- Check for RESULT(j) > THRESH
!
 2                   ntestt = ntestt + ntest
!
!              Print out tests which fail.
!
                     DO jr = 1 , ntest
                        IF ( Result(jr)>=Thresh ) THEN
!
!                    If this is the first test to fail,
!                    print a header to the data file.
!
                           IF ( nerrs==0 ) THEN
                              WRITE (Nounit,FMT=99002) 'CHB'
                              WRITE (Nounit,FMT=99003)
                              WRITE (Nounit,FMT=99004)
                              WRITE (Nounit,FMT=99005) 'Hermitian'
                              WRITE (Nounit,FMT=99006) 'unitary' , '*' ,&
     &                               'conjugate transpose' , ('*',j=1,4)
                           ENDIF
                           nerrs = nerrs + 1
                           WRITE (Nounit,FMT=99007) n , k , ioldsd ,    &
     &                            jtype , jr , Result(jr)
                        ENDIF
                     ENDDO
                  ENDIF
!
               ENDDO
            ENDIF
         ENDDO
      ENDDO
!
!     Summary
!
      CALL SLASUM('CHB',Nounit,nerrs,ntestt)
      RETURN
!
99001 FORMAT (' CCHKHB: ',A,' returned INFO=',I6,'.',/9X,'N=',I6,       &
     &        ', JTYPE=',I6,', ISEED=(',3(I5,','),I5,')')
99002 FORMAT (/1X,A3,                                                   &
     &     ' -- Complex Hermitian Banded Tridiagonal Reduction Routines'&
     &     )
99003 FORMAT (' Matrix types (see SCHK23 for details): ')
!
99004 FORMAT (/' Special Matrices:',                                    &
     &        /'  1=Zero matrix.                        ',              &
     &        '  5=Diagonal: clustered entries.',                       &
     &        /'  2=Identity matrix.                    ',              &
     &        '  6=Diagonal: large, evenly spaced.',                    &
     &        /'  3=Diagonal: evenly spaced entries.    ',              &
     &        '  7=Diagonal: small, evenly spaced.',                    &
     &        /'  4=Diagonal: geometr. spaced entries.')
99005 FORMAT (' Dense ',A,' Banded Matrices:',                          &
     &        /'  8=Evenly spaced eigenvals.            ',              &
     &        ' 12=Small, evenly spaced eigenvals.',                    &
     &        /'  9=Geometrically spaced eigenvals.     ',              &
     &        ' 13=Matrix with random O(1) entries.',                   &
     &        /' 10=Clustered eigenvalues.              ',              &
     &        ' 14=Matrix with large random entries.',                  &
     &        /' 11=Large, evenly spaced eigenvals.     ',              &
     &        ' 15=Matrix with small random entries.')
!
99006 FORMAT (/' Tests performed:   (S is Tridiag,  U is ',A,',',/20X,A,&
     &        ' means ',A,'.',/' UPLO=''U'':',/'  1= | A - U S U',A1,   &
     &        ' | / ( |A| n ulp )     ','  2= | I - U U',A1,            &
     &        ' | / ( n ulp )',/' UPLO=''L'':',/'  3= | A - U S U',A1,  &
     &        ' | / ( |A| n ulp )     ','  4= | I - U U',A1,            &
     &        ' | / ( n ulp )')
99007 FORMAT (' N=',I5,', K=',I4,', seed=',4(I4,','),' type ',I2,       &
     &        ', test(',I2,')=',G10.3)
!
!     End of CCHKHB
!
      END SUBROUTINE CCHKHB
