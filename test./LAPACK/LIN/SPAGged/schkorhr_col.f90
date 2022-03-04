!*==schkorhr_col.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b SCHKORHR_COL
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SCHKORHR_COL( THRESH, TSTERR, NM, MVAL, NN, NVAL, NNB,
!                                NBVAL, NOUT )
!
!       .. Scalar Arguments ..
!       LOGICAL            TSTERR
!       INTEGER            NM, NN, NNB, NOUT
!       REAL               THRESH
!       ..
!       .. Array Arguments ..
!       INTEGER            MVAL( * ), NBVAL( * ), NVAL( * )
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SCHKORHR_COL tests:
!>   1) SORGTSQR and SORHR_COL using SLATSQR, SGEMQRT,
!>   2) SORGTSQR_ROW and SORHR_COL inside DGETSQRHRT
!>      (which calls SLATSQR, SORGTSQR_ROW and SORHR_COL) using SGEMQRT.
!> Therefore, SLATSQR (part of SGEQR), SGEMQRT (part of SGEMQR)
!> have to be tested before this test.
!>
!> \endverbatim
!
!  Arguments:
!  ==========
!
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
!> \param[in] NNB
!> \verbatim
!>          NNB is INTEGER
!>          The number of values of NB contained in the vector NBVAL.
!> \endverbatim
!>
!> \param[in] NBVAL
!> \verbatim
!>          NBVAL is INTEGER array, dimension (NBVAL)
!>          The values of the blocksize NB.
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
!> \date November 2020
!
!> \ingroup single_lin
!
!  =====================================================================
      SUBROUTINE SCHKORHR_COL(Thresh,Tsterr,Nm,Mval,Nn,Nval,Nnb,Nbval,  &
     &                        Nout)
      IMPLICIT NONE
!*--SCHKORHR_COL112
!
!  -- LAPACK test routine (version 3.10.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2020
!
!     .. Scalar Arguments ..
      LOGICAL Tsterr
      INTEGER Nm , Nn , Nnb , Nout
      REAL Thresh
!     ..
!     .. Array Arguments ..
      INTEGER Mval(*) , Nbval(*) , Nval(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER NTESTS
      PARAMETER (NTESTS=6)
!     ..
!     .. Local Scalars ..
      CHARACTER(LEN=3) path
      INTEGER i , imb1 , inb1 , inb2 , j , t , m , n , mb1 , nb1 , nb2 ,&
     &        nfail , nerrs , nrun
!
!     .. Local Arrays ..
      REAL result(NTESTS)
!     ..
!     .. External Subroutines ..
      EXTERNAL ALAHD , ALASUM , SERRORHR_COL , SORHR_COL01 , SORHR_COL02
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX , MIN
!     ..
!     .. Scalars in Common ..
      LOGICAL LERr , OK
      CHARACTER(LEN=32) SRNamt
      INTEGER INFot , NUNit
!     ..
!     .. Common blocks ..
      COMMON /INFOC / INFot , NUNit , OK , LERr
      COMMON /SRNAMC/ SRNamt
!     ..
!     .. Executable Statements ..
!
!     Initialize constants
!
      path(1:1) = 'S'
      path(2:3) = 'HH'
      nrun = 0
      nfail = 0
      nerrs = 0
!
!     Test the error exits
!
      IF ( Tsterr ) CALL SERRORHR_COL(path,Nout)
      INFot = 0
!
!     Do for each value of M in MVAL.
!
      DO i = 1 , Nm
         m = Mval(i)
!
!        Do for each value of N in NVAL.
!
         DO j = 1 , Nn
            n = Nval(j)
!
!           Only for M >= N
!
            IF ( MIN(m,n)>0 .AND. m>=n ) THEN
!
!              Do for each possible value of MB1
!
               DO imb1 = 1 , Nnb
                  mb1 = Nbval(imb1)
!
!                 Only for MB1 > N
!
                  IF ( mb1>n ) THEN
!
!                    Do for each possible value of NB1
!
                     DO inb1 = 1 , Nnb
                        nb1 = Nbval(inb1)
!
!                       Do for each possible value of NB2
!
                        DO inb2 = 1 , Nnb
                           nb2 = Nbval(inb2)
!
                           IF ( nb1>0 .AND. nb2>0 ) THEN
!
!                             Test SORHR_COL
!
                              CALL SORHR_COL01(m,n,mb1,nb1,nb2,result)
!
!                             Print information about the tests that did
!                             not pass the threshold.
!
                              DO t = 1 , NTESTS
                                 IF ( result(t)>=Thresh ) THEN
                                    IF ( nfail==0 .AND. nerrs==0 )      &
     &                                 CALL ALAHD(Nout,path)
                                    WRITE (Nout,FMT=99001) m , n , mb1 ,&
     &                                 nb1 , nb2 , t , result(t)
                                    nfail = nfail + 1
                                 ENDIF
                              ENDDO
                              nrun = nrun + NTESTS
                           ENDIF
                        ENDDO
                     ENDDO
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
      ENDDO
!
!     Do for each value of M in MVAL.
!
      DO i = 1 , Nm
         m = Mval(i)
!
!        Do for each value of N in NVAL.
!
         DO j = 1 , Nn
            n = Nval(j)
!
!           Only for M >= N
!
            IF ( MIN(m,n)>0 .AND. m>=n ) THEN
!
!              Do for each possible value of MB1
!
               DO imb1 = 1 , Nnb
                  mb1 = Nbval(imb1)
!
!                 Only for MB1 > N
!
                  IF ( mb1>n ) THEN
!
!                    Do for each possible value of NB1
!
                     DO inb1 = 1 , Nnb
                        nb1 = Nbval(inb1)
!
!                       Do for each possible value of NB2
!
                        DO inb2 = 1 , Nnb
                           nb2 = Nbval(inb2)
!
                           IF ( nb1>0 .AND. nb2>0 ) THEN
!
!                             Test SORHR_COL
!
                              CALL SORHR_COL02(m,n,mb1,nb1,nb2,result)
!
!                             Print information about the tests that did
!                             not pass the threshold.
!
                              DO t = 1 , NTESTS
                                 IF ( result(t)>=Thresh ) THEN
                                    IF ( nfail==0 .AND. nerrs==0 )      &
     &                                 CALL ALAHD(Nout,path)
                                    WRITE (Nout,FMT=99002) m , n , mb1 ,&
     &                                 nb1 , nb2 , t , result(t)
                                    nfail = nfail + 1
                                 ENDIF
                              ENDDO
                              nrun = nrun + NTESTS
                           ENDIF
                        ENDDO
                     ENDDO
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
      ENDDO
!
!     Print a summary of the results.
!
      CALL ALASUM(path,Nout,nfail,nrun,nerrs)
!
99001 FORMAT ('SORGTSQR and SORHR_COL: M=',I5,', N=',I5,', MB1=',I5,    &
     &        ', NB1=',I5,', NB2=',I5,' test(',I2,')=',G12.5)
99002 FORMAT ('SORGTSQR_ROW and SORHR_COL: M=',I5,', N=',I5,', MB1=',I5,&
     &        ', NB1=',I5,', NB2=',I5,' test(',I2,')=',G12.5)
!
!     End of SCHKORHR_COL
!
      END SUBROUTINE SCHKORHR_COL
