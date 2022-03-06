!*==zebchvxx.f90  processed by SPAG 7.51RB at 17:35 on  4 Mar 2022
!> \brief \b ZEBCHVXX
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!      SUBROUTINE ZEBCHVXX( THRESH, PATH )
!
!     .. Scalar Arguments ..
!      DOUBLE PRECISION  THRESH
!      CHARACTER*3       PATH
!       ..
!
!  Purpose
!  ======
!
!> \details \b Purpose:
!> \verbatim
!>
!>  ZEBCHVXX will run Z**SVXX on a series of Hilbert matrices and then
!>  compare the error bounds returned by Z**SVXX to see if the returned
!>  answer indeed falls within those bounds.
!>
!>  Eight test ratios will be computed.  The tests will pass if they are .LT.
!>  THRESH.  There are two cases that are determined by 1 / (SQRT( N ) * EPS).
!>  If that value is .LE. to the component wise reciprocal condition number,
!>  it uses the guaranteed case, other wise it uses the unguaranteed case.
!>
!>  Test ratios:
!>     Let Xc be X_computed and Xt be X_truth.
!>     The norm used is the infinity norm.
!>
!>     Let A be the guaranteed case and B be the unguaranteed case.
!>
!>       1. Normwise guaranteed forward error bound.
!>       A: norm ( abs( Xc - Xt ) / norm ( Xt ) .LE. ERRBND( *, nwise_i, bnd_i ) and
!>          ERRBND( *, nwise_i, bnd_i ) .LE. MAX(SQRT(N),10) * EPS.
!>          If these conditions are met, the test ratio is set to be
!>          ERRBND( *, nwise_i, bnd_i ) / MAX(SQRT(N), 10).  Otherwise it is 1/EPS.
!>       B: For this case, CGESVXX should just return 1.  If it is less than
!>          one, treat it the same as in 1A.  Otherwise it fails. (Set test
!>          ratio to ERRBND( *, nwise_i, bnd_i ) * THRESH?)
!>
!>       2. Componentwise guaranteed forward error bound.
!>       A: norm ( abs( Xc(j) - Xt(j) ) ) / norm (Xt(j)) .LE. ERRBND( *, cwise_i, bnd_i )
!>          for all j .AND. ERRBND( *, cwise_i, bnd_i ) .LE. MAX(SQRT(N), 10) * EPS.
!>          If these conditions are met, the test ratio is set to be
!>          ERRBND( *, cwise_i, bnd_i ) / MAX(SQRT(N), 10).  Otherwise it is 1/EPS.
!>       B: Same as normwise test ratio.
!>
!>       3. Backwards error.
!>       A: The test ratio is set to BERR/EPS.
!>       B: Same test ratio.
!>
!>       4. Reciprocal condition number.
!>       A: A condition number is computed with Xt and compared with the one
!>          returned from CGESVXX.  Let RCONDc be the RCOND returned by CGESVXX
!>          and RCONDt be the RCOND from the truth value.  Test ratio is set to
!>          MAX(RCONDc/RCONDt, RCONDt/RCONDc).
!>       B: Test ratio is set to 1 / (EPS * RCONDc).
!>
!>       5. Reciprocal normwise condition number.
!>       A: The test ratio is set to
!>          MAX(ERRBND( *, nwise_i, cond_i ) / NCOND, NCOND / ERRBND( *, nwise_i, cond_i )).
!>       B: Test ratio is set to 1 / (EPS * ERRBND( *, nwise_i, cond_i )).
!>
!>       6. Reciprocal componentwise condition number.
!>       A: Test ratio is set to
!>          MAX(ERRBND( *, cwise_i, cond_i ) / CCOND, CCOND / ERRBND( *, cwise_i, cond_i )).
!>       B: Test ratio is set to 1 / (EPS * ERRBND( *, cwise_i, cond_i )).
!>
!>     .. Parameters ..
!>     NMAX is determined by the largest number in the inverse of the hilbert
!>     matrix.  Precision is exhausted when the largest entry in it is greater
!>     than 2 to the power of the number of bits in the fraction of the data
!>     type used plus one, which is 24 for single precision.
!>     NMAX should be 6 for single and 11 for double.
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
!> \ingroup complex16_lin
!
!  =====================================================================
      SUBROUTINE ZEBCHVXX(Thresh,Path)
      IMPLICIT NONE
!*--ZEBCHVXX100
!     .. Scalar Arguments ..
      DOUBLE PRECISION Thresh
      CHARACTER*3 Path
 
      INTEGER NMAX , NPARAMS , NERRBND , NTESTS , kl , ku
      PARAMETER (NMAX=10,NPARAMS=2,NERRBND=3,NTESTS=6)
 
!     .. Local Scalars ..
      INTEGER n , nrhs , info , i , j , k , nfail , lda , n_aux_tests , &
     &        ldab , ldafb
      CHARACTER fact , trans , uplo , equed
      CHARACTER*2 c2
      CHARACTER(3) nguar , cguar
      LOGICAL printed_guide
      DOUBLE PRECISION ncond , ccond , m , normdif , normt , rcond ,    &
     &                 rnorm , rinorm , sumr , sumri , eps , berr(NMAX) &
     &                 , rpvgrw , orcond , cwise_err , nwise_err ,      &
     &                 cwise_bnd , nwise_bnd , cwise_rcond ,            &
     &                 nwise_rcond , condthresh , errthresh
      COMPLEX*16 zdum
 
!     .. Local Arrays ..
      DOUBLE PRECISION tstrat(NTESTS) , rinv(NMAX) , params(NPARAMS) ,  &
     &                 s(NMAX) , r(NMAX) , c(NMAX) , rwork(3*NMAX) ,    &
     &                 diff(NMAX,NMAX) , errbnd_n(NMAX*3) ,             &
     &                 errbnd_c(NMAX*3)
      INTEGER ipiv(NMAX)
      COMPLEX*16 a(NMAX,NMAX) , invhilb(NMAX,NMAX) , x(NMAX,NMAX) ,     &
     &           work(NMAX*3*5) , af(NMAX,NMAX) , b(NMAX,NMAX) ,        &
     &           acopy(NMAX,NMAX) , ab((NMAX-1)+(NMAX-1)+1,NMAX) ,      &
     &           abcopy((NMAX-1)+(NMAX-1)+1,NMAX) ,                     &
     &           afb(2*(NMAX-1)+(NMAX-1)+1,NMAX)
 
!     .. External Functions ..
      DOUBLE PRECISION DLAMCH
 
!     .. External Subroutines ..
EXTERNAL    ZGESVXX,   ZPOSVXX,   ZSYSVXX,   ZGBSVXX
EXTERNAL    ZLACPY,   LSAMEN
      LOGICAL LSAMEN
 
!     .. Intrinsic Functions ..
      INTRINSIC SQRT , MAX , ABS , DBLE , DIMAG
 
!     .. Statement Functions ..
      DOUBLE PRECISION CABS1
 
!     .. Statement Function Definitions ..
      CABS1(zdum) = ABS(DBLE(zdum)) + ABS(DIMAG(zdum))
 
!     .. Parameters ..
      INTEGER BND_I , COND_I
      PARAMETER (BND_I=2,COND_I=3)
 
!  Create the loop to test out the Hilbert matrices
 
      fact = 'E'
      uplo = 'U'
      trans = 'N'
      equed = 'N'
      eps = DLAMCH('Epsilon')
      nfail = 0
      n_aux_tests = 0
      lda = NMAX
      ldab = (NMAX-1) + (NMAX-1) + 1
      ldafb = 2*(NMAX-1) + (NMAX-1) + 1
      c2 = Path(2:3)
 
!     Main loop to test the different Hilbert Matrices.
 
      printed_guide = .FALSE.
 
      DO n = 1 , NMAX
         params(1) = -1
         params(2) = -1
 
         kl = n - 1
         ku = n - 1
         nrhs = n
         m = MAX(SQRT(DBLE(n)),10.0D+0)
 
!        Generate the Hilbert matrix, its inverse, and the
!        right hand side, all scaled by the LCM(1,..,2N-1).
         CALL zlahilb(n,n,a,lda,invhilb,lda,b,lda,work,info,Path)
 
!        Copy A into ACOPY.
         CALL ZLACPY('ALL',n,n,a,NMAX,acopy,NMAX)
 
!        Store A in band format for GB tests
         DO j = 1 , n
            DO i = 1 , kl + ku + 1
               ab(i,j) = (0.0D+0,0.0D+0)
            ENDDO
         ENDDO
         DO j = 1 , n
            DO i = MAX(1,j-ku) , MIN(n,j+kl)
               ab(ku+1+i-j,j) = a(i,j)
            ENDDO
         ENDDO
 
!        Copy AB into ABCOPY.
         DO j = 1 , n
            DO i = 1 , kl + ku + 1
               abcopy(i,j) = (0.0D+0,0.0D+0)
            ENDDO
         ENDDO
         CALL ZLACPY('ALL',kl+ku+1,n,ab,ldab,abcopy,ldab)
 
!        Call Z**SVXX with default PARAMS and N_ERR_BND = 3.
         IF ( LSAMEN(2,c2,'SY') ) THEN
            CALL ZSYSVXX(fact,uplo,n,nrhs,acopy,lda,af,lda,ipiv,equed,s,&
     &                   b,lda,x,lda,orcond,rpvgrw,berr,NERRBND,        &
     &                   errbnd_n,errbnd_c,NPARAMS,params,work,rwork,   &
     &                   info)
         ELSEIF ( LSAMEN(2,c2,'PO') ) THEN
            CALL ZPOSVXX(fact,uplo,n,nrhs,acopy,lda,af,lda,equed,s,b,   &
     &                   lda,x,lda,orcond,rpvgrw,berr,NERRBND,errbnd_n, &
     &                   errbnd_c,NPARAMS,params,work,rwork,info)
         ELSEIF ( LSAMEN(2,c2,'HE') ) THEN
            CALL ZHESVXX(fact,uplo,n,nrhs,acopy,lda,af,lda,ipiv,equed,s,&
     &                   b,lda,x,lda,orcond,rpvgrw,berr,NERRBND,        &
     &                   errbnd_n,errbnd_c,NPARAMS,params,work,rwork,   &
     &                   info)
         ELSEIF ( LSAMEN(2,c2,'GB') ) THEN
            CALL ZGBSVXX(fact,trans,n,kl,ku,nrhs,abcopy,ldab,afb,ldafb, &
     &                   ipiv,equed,r,c,b,lda,x,lda,orcond,rpvgrw,berr, &
     &                   NERRBND,errbnd_n,errbnd_c,NPARAMS,params,work, &
     &                   rwork,info)
         ELSE
            CALL ZGESVXX(fact,trans,n,nrhs,acopy,lda,af,lda,ipiv,equed, &
     &                   r,c,b,lda,x,lda,orcond,rpvgrw,berr,NERRBND,    &
     &                   errbnd_n,errbnd_c,NPARAMS,params,work,rwork,   &
     &                   info)
         ENDIF
 
         n_aux_tests = n_aux_tests + 1
         IF ( orcond>=eps ) THEN
!        Either factorization failed or the matrix is flagged, and 1 <=
!        INFO <= N+1. We don't decide based on rcond anymore.
!            IF (INFO .EQ. 0 .OR. INFO .GT. N+1) THEN
!               NFAIL = NFAIL + 1
!               WRITE (*, FMT=8000) N, INFO, ORCOND, RCOND
!            END IF
!        Either everything succeeded (INFO == 0) or some solution failed
!        to converge (INFO > N+1).
            IF ( info>0 .AND. info<=n+1 ) THEN
               nfail = nfail + 1
               WRITE (*,FMT=99012) c2 , n , info , orcond , rcond
            ENDIF
         ENDIF
 
!        Calculating the difference between Z**SVXX's X and the true X.
         DO i = 1 , n
            DO j = 1 , nrhs
               diff(i,j) = x(i,j) - invhilb(i,j)
            ENDDO
         ENDDO
 
!        Calculating the RCOND
         rnorm = 0
         rinorm = 0
         IF ( LSAMEN(2,c2,'PO') .OR. LSAMEN(2,c2,'SY') .OR.             &
     &        LSAMEN(2,c2,'HE') ) THEN
            DO i = 1 , n
               sumr = 0
               sumri = 0
               DO j = 1 , n
                  sumr = sumr + s(i)*CABS1(a(i,j))*s(j)
                  sumri = sumri + CABS1(invhilb(i,j))/(s(j)*s(i))
               ENDDO
               rnorm = MAX(rnorm,sumr)
               rinorm = MAX(rinorm,sumri)
            ENDDO
         ELSEIF ( LSAMEN(2,c2,'GE') .OR. LSAMEN(2,c2,'GB') ) THEN
            DO i = 1 , n
               sumr = 0
               sumri = 0
               DO j = 1 , n
                  sumr = sumr + r(i)*CABS1(a(i,j))*c(j)
                  sumri = sumri + CABS1(invhilb(i,j))/(r(j)*c(i))
               ENDDO
               rnorm = MAX(rnorm,sumr)
               rinorm = MAX(rinorm,sumri)
            ENDDO
         ENDIF
 
         rnorm = rnorm/CABS1(a(1,1))
         rcond = 1.0D+0/(rnorm*rinorm)
 
!        Calculating the R for normwise rcond.
         DO i = 1 , n
            rinv(i) = 0.0D+0
         ENDDO
         DO j = 1 , n
            DO i = 1 , n
               rinv(i) = rinv(i) + CABS1(a(i,j))
            ENDDO
         ENDDO
 
!        Calculating the Normwise rcond.
         rinorm = 0.0D+0
         DO i = 1 , n
            sumri = 0.0D+0
            DO j = 1 , n
               sumri = sumri + CABS1(invhilb(i,j)*rinv(j))
            ENDDO
            rinorm = MAX(rinorm,sumri)
         ENDDO
 
!        invhilb is the inverse *unscaled* Hilbert matrix, so scale its norm
!        by 1/A(1,1) to make the scaling match A (the scaled Hilbert matrix)
         ncond = CABS1(a(1,1))/rinorm
 
         condthresh = m*eps
         errthresh = m*eps
 
         DO k = 1 , nrhs
            normt = 0.0D+0
            normdif = 0.0D+0
            cwise_err = 0.0D+0
            DO i = 1 , n
               normt = MAX(CABS1(invhilb(i,k)),normt)
               normdif = MAX(CABS1(x(i,k)-invhilb(i,k)),normdif)
               IF ( invhilb(i,k)/=0.0D+0 ) THEN
                  cwise_err = MAX(CABS1(x(i,k)-invhilb(i,k))            &
     &                        /CABS1(invhilb(i,k)),cwise_err)
               ELSEIF ( x(i,k)/=0.0D+0 ) THEN
                  cwise_err = DLAMCH('OVERFLOW')
               ENDIF
            ENDDO
            IF ( normt/=0.0D+0 ) THEN
               nwise_err = normdif/normt
            ELSEIF ( normdif/=0.0D+0 ) THEN
               nwise_err = DLAMCH('OVERFLOW')
            ELSE
               nwise_err = 0.0D+0
            ENDIF
 
            DO i = 1 , n
               rinv(i) = 0.0D+0
            ENDDO
            DO j = 1 , n
               DO i = 1 , n
                  rinv(i) = rinv(i) + CABS1(a(i,j)*invhilb(j,k))
               ENDDO
            ENDDO
            rinorm = 0.0D+0
            DO i = 1 , n
               sumri = 0.0D+0
               DO j = 1 , n
                  sumri = sumri +                                       &
     &                    CABS1(invhilb(i,j)*rinv(j)/invhilb(i,k))
               ENDDO
               rinorm = MAX(rinorm,sumri)
            ENDDO
!        invhilb is the inverse *unscaled* Hilbert matrix, so scale its norm
!        by 1/A(1,1) to make the scaling match A (the scaled Hilbert matrix)
            ccond = CABS1(a(1,1))/rinorm
 
!        Forward error bound tests
            nwise_bnd = errbnd_n(k+(BND_I-1)*nrhs)
            cwise_bnd = errbnd_c(k+(BND_I-1)*nrhs)
            nwise_rcond = errbnd_n(k+(COND_I-1)*nrhs)
            cwise_rcond = errbnd_c(k+(COND_I-1)*nrhs)
!            write (*,*) 'nwise : ', n, k, ncond, nwise_rcond,
!     $           condthresh, ncond.ge.condthresh
!            write (*,*) 'nwise2: ', k, nwise_bnd, nwise_err, errthresh
            IF ( ncond>=condthresh ) THEN
               nguar = 'YES'
               IF ( nwise_bnd>errthresh ) THEN
                  tstrat(1) = 1/(2.0D+0*eps)
               ELSE
                  IF ( nwise_bnd/=0.0D+0 ) THEN
                     tstrat(1) = nwise_err/nwise_bnd
                  ELSEIF ( nwise_err/=0.0D+0 ) THEN
                     tstrat(1) = 1/(16.0*eps)
                  ELSE
                     tstrat(1) = 0.0D+0
                  ENDIF
                  IF ( tstrat(1)>1.0D+0 ) tstrat(1) = 1/(4.0D+0*eps)
               ENDIF
            ELSE
               nguar = 'NO'
               IF ( nwise_bnd<1.0D+0 ) THEN
                  tstrat(1) = 1/(8.0D+0*eps)
               ELSE
                  tstrat(1) = 1.0D+0
               ENDIF
            ENDIF
!            write (*,*) 'cwise : ', n, k, ccond, cwise_rcond,
!     $           condthresh, ccond.ge.condthresh
!            write (*,*) 'cwise2: ', k, cwise_bnd, cwise_err, errthresh
            IF ( ccond>=condthresh ) THEN
               cguar = 'YES'
               IF ( cwise_bnd>errthresh ) THEN
                  tstrat(2) = 1/(2.0D+0*eps)
               ELSE
                  IF ( cwise_bnd/=0.0D+0 ) THEN
                     tstrat(2) = cwise_err/cwise_bnd
                  ELSEIF ( cwise_err/=0.0D+0 ) THEN
                     tstrat(2) = 1/(16.0D+0*eps)
                  ELSE
                     tstrat(2) = 0.0D+0
                  ENDIF
                  IF ( tstrat(2)>1.0D+0 ) tstrat(2) = 1/(4.0D+0*eps)
               ENDIF
            ELSE
               cguar = 'NO'
               IF ( cwise_bnd<1.0D+0 ) THEN
                  tstrat(2) = 1/(8.0D+0*eps)
               ELSE
                  tstrat(2) = 1.0D+0
               ENDIF
            ENDIF
 
!     Backwards error test
            tstrat(3) = berr(k)/eps
 
!     Condition number tests
            tstrat(4) = rcond/orcond
            IF ( rcond>=condthresh .AND. tstrat(4)<1.0D+0 ) tstrat(4)   &
     &           = 1.0D+0/tstrat(4)
 
            tstrat(5) = ncond/nwise_rcond
            IF ( ncond>=condthresh .AND. tstrat(5)<1.0D+0 ) tstrat(5)   &
     &           = 1.0D+0/tstrat(5)
 
            tstrat(6) = ccond/nwise_rcond
            IF ( ccond>=condthresh .AND. tstrat(6)<1.0D+0 ) tstrat(6)   &
     &           = 1.0D+0/tstrat(6)
 
            DO i = 1 , NTESTS
               IF ( tstrat(i)>Thresh ) THEN
                  IF ( .NOT.printed_guide ) THEN
                     WRITE (*,*)
                     WRITE (*,99001) 1
!     Test ratios.
99001                FORMAT (3X,I2,                                     &
     &                       ': Normwise guaranteed forward error',/5X, &
     &                       'Guaranteed case: if norm ( abs( Xc - Xt )'&
     &                       ,                                          &
     &           ' / norm ( Xt ) .LE. ERRBND( *, nwise_i, bnd_i ), then'&
     &           ,/5X,                                                  &
     &         'ERRBND( *, nwise_i, bnd_i ) .LE. MAX(SQRT(N), 10) * EPS'&
     &         )
                     WRITE (*,99002) 2
99002                FORMAT (3X,I2,                                     &
     &                       ': Componentwise guaranteed forward error')
                     WRITE (*,99003) 3
99003                FORMAT (3X,I2,': Backwards error')
                     WRITE (*,99004) 4
99004                FORMAT (3X,I2,': Reciprocal condition number')
                     WRITE (*,99005) 5
99005                FORMAT (3X,I2,                                     &
     &                       ': Reciprocal normwise condition number')
                     WRITE (*,99006) 6
99006                FORMAT (3X,I2,': Raw normwise error estimate')
                     WRITE (*,99007) 7
99007                FORMAT (3X,I2,                                     &
     &                     ': Reciprocal componentwise condition number'&
     &                     )
                     WRITE (*,99008) 8
99008                FORMAT (3X,I2,': Raw componentwise error estimate')
                     WRITE (*,*)
                     printed_guide = .TRUE.
                  ENDIF
                  WRITE (*,99009) c2 , n , k , nguar , cguar , i ,      &
     &                            tstrat(i)
99009             FORMAT (' Z',A2,'SVXX: N =',I2,', RHS = ',I2,         &
     &                    ', NWISE GUAR. = ',A,', CWISE GUAR. = ',A,    &
     &                    ' test(',I1,') =',G12.5)
                  nfail = nfail + 1
               ENDIF
            ENDDO
         ENDDO
 
!$$$         WRITE(*,*)
!$$$         WRITE(*,*) 'Normwise Error Bounds'
!$$$         WRITE(*,*) 'Guaranteed error bound: ',ERRBND(NRHS,nwise_i,bnd_i)
!$$$         WRITE(*,*) 'Reciprocal condition number: ',ERRBND(NRHS,nwise_i,cond_i)
!$$$         WRITE(*,*) 'Raw error estimate: ',ERRBND(NRHS,nwise_i,rawbnd_i)
!$$$         WRITE(*,*)
!$$$         WRITE(*,*) 'Componentwise Error Bounds'
!$$$         WRITE(*,*) 'Guaranteed error bound: ',ERRBND(NRHS,cwise_i,bnd_i)
!$$$         WRITE(*,*) 'Reciprocal condition number: ',ERRBND(NRHS,cwise_i,cond_i)
!$$$         WRITE(*,*) 'Raw error estimate: ',ERRBND(NRHS,cwise_i,rawbnd_i)
!$$$         print *, 'Info: ', info
!$$$         WRITE(*,*)
!         WRITE(*,*) 'TSTRAT: ',TSTRAT
 
      ENDDO
 
      WRITE (*,*)
      IF ( nfail>0 ) THEN
         WRITE (*,99010) c2 , nfail , NTESTS*n + n_aux_tests
99010    FORMAT (' Z',A2,'SVXX: ',I6,' out of ',I6,                     &
     &           ' tests failed to pass the threshold')
      ELSE
         WRITE (*,99011) c2
99011    FORMAT (' Z',A2,'SVXX passed the tests of error bounds')
      ENDIF
 
99012 FORMAT (' Z',A2,'SVXX: N =',I2,', INFO = ',I3,', ORCOND = ',G12.5,&
     &        ', real RCOND = ',G12.5)
 
      END SUBROUTINE ZEBCHVXX
