!*==cchktz.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b cchktz
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE cchktz( DOTYPE, NM, MVAL, NN, NVAL, THRESH, TSTERR, A,
!                          COPYA, S, TAU, WORK, RWORK, NOUT )
!
!       .. Scalar Arguments ..
!       LOGICAL            TSTERR
!       INTEGER            NM, NN, NOUT
!       REAL               THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            DOTYPE( * )
!       INTEGER            MVAL( * ), NVAL( * )
!       REAL               S( * ), RWORK( * )
!       COMPLEX            A( * ), COPYA( * ), TAU( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> cchktz tests CTZRZF.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] DOTYPE
!> \verbatim
!>          DOTYPE is LOGICAL array, dimension (NTYPES)
!>          The matrix types to be used for testing.  Matrices of type j
!>          (for 1 <= j <= NTYPES) are used for testing if DOTYPE(j) =
!>          .TRUE.; if DOTYPE(j) = .FALSE., then type j is not used.
!> \endverbatim
!>
!> \param[in] NM
!> \verbatim
!>          NM is INTEGER
!>          The number of values of M contained in the vector MVAL.
!> \endverbatim
!>
!> \param[in] MVAL
!> \verbatim
!>          MVAL is INTEGER array, dimension (NM)
!>          The values of the matrix row dimension M.
!> \endverbatim
!>
!> \param[in] NN
!> \verbatim
!>          NN is INTEGER
!>          The number of values of N contained in the vector NVAL.
!> \endverbatim
!>
!> \param[in] NVAL
!> \verbatim
!>          NVAL is INTEGER array, dimension (NN)
!>          The values of the matrix column dimension N.
!> \endverbatim
!>
!> \param[in] THRESH
!> \verbatim
!>          THRESH is REAL
!>          The threshold value for the test ratios.  A result is
!>          included in the output file if RESULT >= THRESH.  To have
!>          every test ratio printed, use THRESH = 0.
!> \endverbatim
!>
!> \param[in] TSTERR
!> \verbatim
!>          TSTERR is LOGICAL
!>          Flag that indicates whether error exits are to be tested.
!> \endverbatim
!>
!> \param[out] A
!> \verbatim
!>          A is COMPLEX array, dimension (MMAX*NMAX)
!>          where MMAX is the maximum value of M in MVAL and NMAX is the
!>          maximum value of N in NVAL.
!> \endverbatim
!>
!> \param[out] COPYA
!> \verbatim
!>          COPYA is COMPLEX array, dimension (MMAX*NMAX)
!> \endverbatim
!>
!> \param[out] S
!> \verbatim
!>          S is REAL array, dimension
!>                      (min(MMAX,NMAX))
!> \endverbatim
!>
!> \param[out] TAU
!> \verbatim
!>          TAU is COMPLEX array, dimension (MMAX)
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension
!>                      (MMAX*NMAX + 4*NMAX + MMAX)
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension (2*NMAX)
!> \endverbatim
!>
!> \param[in] NOUT
!> \verbatim
!>          NOUT is INTEGER
!>          The unit number for output.
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
!> \ingroup complex_lin
!
!  =====================================================================
      SUBROUTINE cchktz(Dotype,Nm,Mval,Nn,Nval,Thresh,Tsterr,A,Copya,S, &
     &                  Tau,Work,Rwork,Nout)
      IMPLICIT NONE
!*--CCHKTZ141
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      LOGICAL Tsterr
      INTEGER Nm , Nn , Nout
      REAL Thresh
!     ..
!     .. Array Arguments ..
      LOGICAL Dotype(*)
      INTEGER Mval(*) , Nval(*)
      REAL S(*) , Rwork(*)
      COMPLEX A(*) , Copya(*) , Tau(*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER NTYPES
      PARAMETER (NTYPES=3)
      INTEGER NTESTS
      PARAMETER (NTESTS=3)
      REAL ONE , ZERO
      PARAMETER (ONE=1.0E0,ZERO=0.0E0)
!     ..
!     .. Local Scalars ..
      CHARACTER*3 path
      INTEGER i , im , imode , in , info , k , lda , lwork , m , mnmin ,&
     &        mode , n , nerrs , nfail , nrun
      REAL eps
!     ..
!     .. Local Arrays ..
      INTEGER iseed(4) , iseedy(4)
      REAL result(NTESTS)
!     ..
!     .. External Functions ..
      REAL SLAMCH
EXTERNAL    SLAMCH
!     ..
!     .. External Subroutines ..
EXTERNAL    CGEQR2,   CLACPY,   CLASET
EXTERNAL    CTZRZF
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC CMPLX , MAX , MIN
!     ..
!     .. Scalars in Common ..
      LOGICAL LERr , OK
      CHARACTER*32 SRNamt
      INTEGER INFot , IOUnit
!     ..
!     .. Common blocks ..
      COMMON /INFOC / INFot , IOUnit , OK , LERr
      COMMON /SRNAMC/ SRNamt
!     ..
!     .. Data statements ..
      DATA iseedy/1988 , 1989 , 1990 , 1991/
!     ..
!     .. Executable Statements ..
!
!     Initialize constants and the random number seed.
!
      path(1:1) = 'Complex precision'
      path(2:3) = 'TZ'
      nrun = 0
      nfail = 0
      nerrs = 0
      DO i = 1 , 4
         iseed(i) = iseedy(i)
      ENDDO
      eps = SLAMCH('Epsilon')
!
!     Test the error exits
!
      IF ( Tsterr ) CALL cerrtz(path,Nout)
      INFot = 0
!
      DO im = 1 , Nm
!
!        Do for each value of M in MVAL.
!
         m = Mval(im)
         lda = MAX(1,m)
!
         DO in = 1 , Nn
!
!           Do for each value of N in NVAL for which M .LE. N.
!
            n = Nval(in)
            mnmin = MIN(m,n)
            lwork = MAX(1,n*n+4*m+n)
!
            IF ( m<=n ) THEN
               DO imode = 1 , NTYPES
                  IF ( Dotype(imode) ) THEN
!
!                 Do for each type of singular value distribution.
!                    0:  zero matrix
!                    1:  one small singular value
!                    2:  exponential distribution
!
                     mode = imode - 1
!
!                 Test CTZRZF
!
!                 Generate test matrix of size m by n using
!                 singular value distribution indicated by `mode'.
!
                     IF ( mode==0 ) THEN
                        CALL CLASET('Full',m,n,CMPLX(ZERO),CMPLX(ZERO), &
     &                              A,lda)
                        DO i = 1 , mnmin
                           S(i) = ZERO
                        ENDDO
                     ELSE
                        CALL CLATMS(m,n,'Uniform',iseed,'Nonsymmetric', &
     &                              S,imode,ONE/eps,ONE,m,n,            &
     &                              'No packing',A,lda,Work,info)
                        CALL CGEQR2(m,n,A,lda,Work,Work(mnmin+1),info)
                        CALL CLASET('Lower',m-1,n,CMPLX(ZERO),          &
     &                              CMPLX(ZERO),A(2),lda)
                        CALL slaord('Decreasing',mnmin,S,1)
                     ENDIF
!
!                 Save A and its singular values
!
                     CALL CLACPY('All',m,n,A,lda,Copya,lda)
!
!                 Call CTZRZF to reduce the upper trapezoidal matrix to
!                 upper triangular form.
!
                     SRNamt = 'CTZRZF'
                     CALL CTZRZF(m,n,A,lda,Tau,Work,lwork,info)
!
!                 Compute norm(svd(a) - svd(r))
!
                     result(1) = cqrt12(m,m,A,lda,S,Work,lwork,Rwork)
!
!                 Compute norm( A - R*Q )
!
                     result(2) = crzt01(m,n,Copya,A,lda,Tau,Work,lwork)
!
!                 Compute norm(Q'*Q - I).
!
                     result(3) = crzt02(m,n,A,lda,Tau,Work,lwork)
!
!                 Print information about the tests that did not pass
!                 the threshold.
!
                     DO k = 1 , NTESTS
                        IF ( result(k)>=Thresh ) THEN
                           IF ( nfail==0 .AND. nerrs==0 )               &
     &                          CALL alahd(Nout,path)
                           WRITE (Nout,FMT=99001) m , n , imode , k ,   &
     &                            result(k)
                           nfail = nfail + 1
                        ENDIF
                     ENDDO
                     nrun = nrun + 3
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
      ENDDO
!
!     Print a summary of the results.
!
      CALL alasum(path,Nout,nfail,nrun,nerrs)
!
99001 FORMAT (' M =',I5,', N =',I5,', type ',I2,', test ',I2,           &
     &        ', ratio =',G12.5)
!
!     End if cchktz
!
      END SUBROUTINE cchktz
