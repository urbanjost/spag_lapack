!*==alaerh.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b ALAERH
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ALAERH( PATH, SUBNAM, INFO, INFOE, OPTS, M, N, KL, KU,
!                          N5, IMAT, NFAIL, NERRS, NOUT )
!
!       .. Scalar Arguments ..
!       CHARACTER*3        PATH
!       CHARACTER*( * )    SUBNAM
!       CHARACTER*( * )    OPTS
!       INTEGER            IMAT, INFO, INFOE, KL, KU, M, N, N5, NERRS,
!      $                   NFAIL, NOUT
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ALAERH is an error handler for the LAPACK routines.  It prints the
!> header if this is the first error message and prints the error code
!> and form of recovery, if any.  The character evaluations in this
!> routine may make it slow, but it should not be called once the LAPACK
!> routines are fully debugged.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] PATH
!> \verbatim
!>          PATH is CHARACTER*3
!>          The LAPACK path name of subroutine SUBNAM.
!> \endverbatim
!>
!> \param[in] SUBNAM
!> \verbatim
!>          SUBNAM is CHARACTER*(*)
!>          The name of the subroutine that returned an error code.
!> \endverbatim
!>
!> \param[in] INFO
!> \verbatim
!>          INFO is INTEGER
!>          The error code returned from routine SUBNAM.
!> \endverbatim
!>
!> \param[in] INFOE
!> \verbatim
!>          INFOE is INTEGER
!>          The expected error code from routine SUBNAM, if SUBNAM were
!>          error-free.  If INFOE = 0, an error message is printed, but
!>          if INFOE.NE.0, we assume only the return code INFO is wrong.
!> \endverbatim
!>
!> \param[in] OPTS
!> \verbatim
!>          OPTS is CHARACTER*(*)
!>          The character options to the subroutine SUBNAM, concatenated
!>          into a single character string.  For example, UPLO = 'U',
!>          TRANS = 'T', and DIAG = 'N' for a triangular routine would
!>          be specified as OPTS = 'UTN'.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The matrix row dimension.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The matrix column dimension.  Accessed only if PATH = xGE or
!>          xGB.
!> \endverbatim
!>
!> \param[in] KL
!> \verbatim
!>          KL is INTEGER
!>          The number of sub-diagonals of the matrix.  Accessed only if
!>          PATH = xGB, xPB, or xTB.  Also used for NRHS for PATH = xLS.
!> \endverbatim
!>
!> \param[in] KU
!> \verbatim
!>          KU is INTEGER
!>          The number of super-diagonals of the matrix.  Accessed only
!>          if PATH = xGB.
!> \endverbatim
!>
!> \param[in] N5
!> \verbatim
!>          N5 is INTEGER
!>          A fifth integer parameter, may be the blocksize NB or the
!>          number of right hand sides NRHS.
!> \endverbatim
!>
!> \param[in] IMAT
!> \verbatim
!>          IMAT is INTEGER
!>          The matrix type.
!> \endverbatim
!>
!> \param[in] NFAIL
!> \verbatim
!>          NFAIL is INTEGER
!>          The number of prior tests that did not pass the threshold;
!>          used to determine if the header should be printed.
!> \endverbatim
!>
!> \param[in,out] NERRS
!> \verbatim
!>          NERRS is INTEGER
!>          On entry, the number of errors already detected; used to
!>          determine if the header should be printed.
!>          On exit, NERRS is increased by 1.
!> \endverbatim
!>
!> \param[in] NOUT
!> \verbatim
!>          NOUT is INTEGER
!>          The unit number on which results are to be printed.
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
!> \ingroup aux_lin
!
!  =====================================================================
      SUBROUTINE ALAERH(Path,Subnam,Info,Infoe,Opts,M,N,Kl,Ku,N5,Imat,  &
     &                  Nfail,Nerrs,Nout)
      IMPLICIT NONE
!*--ALAERH151
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER*3 Path
      CHARACTER*(*) Subnam
      CHARACTER*(*) Opts
      INTEGER Imat , Info , Infoe , Kl , Ku , M , N , N5 , Nerrs ,      &
     &        Nfail , Nout
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      CHARACTER uplo
      CHARACTER*2 p2
      CHARACTER*3 c3
!     ..
!     .. External Functions ..
      LOGICAL LSAME , LSAMEN
      EXTERNAL LSAME , LSAMEN
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC LEN_TRIM
!     ..
!     .. External Subroutines ..
      EXTERNAL ALADHD , ALAHD
!     ..
!     .. Executable Statements ..
!
      IF ( Info==0 ) RETURN
      p2 = Path(2:3)
      c3 = Subnam(4:6)
!
!     Print the header if this is the first error message.
!
      IF ( Nfail==0 .AND. Nerrs==0 ) THEN
         IF ( LSAMEN(3,c3,'SV ') .OR. LSAMEN(3,c3,'SVX') ) THEN
            CALL ALADHD(Nout,Path)
         ELSE
            CALL ALAHD(Nout,Path)
         ENDIF
      ENDIF
      Nerrs = Nerrs + 1
!
!     Print the message detailing the error and form of recovery,
!     if any.
!
      IF ( LSAMEN(2,p2,'GE') ) THEN
!
!        xGE:  General matrices
!
         IF ( LSAMEN(3,c3,'TRF') ) THEN
            IF ( Info/=Infoe .AND. Infoe/=0 ) THEN
               WRITE (Nout,FMT=99012) Subnam(1:LEN_TRIM(Subnam)) ,      &
     &                                Info , Infoe , M , N , N5 , Imat
            ELSE
               WRITE (Nout,FMT=99025) Subnam(1:LEN_TRIM(Subnam)) ,      &
     &                                Info , M , N , N5 , Imat
            ENDIF
            IF ( Info/=0 ) WRITE (Nout,FMT=99051)
!
         ELSEIF ( LSAMEN(3,c3,'SV ') ) THEN
!
            IF ( Info/=Infoe .AND. Infoe/=0 ) THEN
               WRITE (Nout,FMT=99016) Subnam(1:LEN_TRIM(Subnam)) ,      &
     &                                Info , Infoe , N , N5 , Imat
            ELSE
               WRITE (Nout,FMT=99030) Subnam(1:LEN_TRIM(Subnam)) ,      &
     &                                Info , N , N5 , Imat
            ENDIF
!
         ELSEIF ( LSAMEN(3,c3,'SVX') ) THEN
!
            IF ( Info/=Infoe .AND. Infoe/=0 ) THEN
               WRITE (Nout,FMT=99008) Subnam(1:LEN_TRIM(Subnam)) ,      &
     &                                Info , Infoe , Opts(1:1) ,        &
     &                                Opts(2:2) , N , N5 , Imat
            ELSE
               WRITE (Nout,FMT=99003) Subnam(1:LEN_TRIM(Subnam)) ,      &
     &                                Info , Opts(1:1) , Opts(2:2) , N ,&
     &                                N5 , Imat
            ENDIF
!
         ELSEIF ( LSAMEN(3,c3,'TRI') ) THEN
!
            WRITE (Nout,FMT=99029) Subnam(1:LEN_TRIM(Subnam)) , Info ,  &
     &                             N , N5 , Imat
!
         ELSEIF ( LSAMEN(5,Subnam(2:6),'LATMS') ) THEN
!
            WRITE (Nout,FMT=99022) Subnam(1:LEN_TRIM(Subnam)) , Info ,  &
     &                             M , N , Imat
!
         ELSEIF ( LSAMEN(3,c3,'CON') ) THEN
!
            WRITE (Nout,FMT=99031) Subnam(1:LEN_TRIM(Subnam)) , Info ,  &
     &                             Opts(1:1) , M , Imat
!
         ELSEIF ( LSAMEN(3,c3,'LS ') ) THEN
!
            WRITE (Nout,FMT=99035) Subnam(1:LEN_TRIM(Subnam)) , Info ,  &
     &                             Opts(1:1) , M , N , Kl , N5 , Imat
!
         ELSEIF ( LSAMEN(3,c3,'LSX') .OR. LSAMEN(3,c3,'LSS') ) THEN
!
            WRITE (Nout,FMT=99026) Subnam(1:LEN_TRIM(Subnam)) , Info ,  &
     &                             M , N , Kl , N5 , Imat
!
         ELSE
!
            WRITE (Nout,FMT=99037) Subnam(1:LEN_TRIM(Subnam)) , Info ,  &
     &                             Opts(1:1) , M , N5 , Imat
         ENDIF
!
      ELSEIF ( LSAMEN(2,p2,'GB') ) THEN
!
!        xGB:  General band matrices
!
         IF ( LSAMEN(3,c3,'TRF') ) THEN
            IF ( Info/=Infoe .AND. Infoe/=0 ) THEN
               WRITE (Nout,FMT=99011) Subnam(1:LEN_TRIM(Subnam)) ,      &
     &                                Info , Infoe , M , N , Kl , Ku ,  &
     &                                N5 , Imat
            ELSE
               WRITE (Nout,FMT=99024) Subnam(1:LEN_TRIM(Subnam)) ,      &
     &                                Info , M , N , Kl , Ku , N5 , Imat
            ENDIF
            IF ( Info/=0 ) WRITE (Nout,FMT=99051)
!
         ELSEIF ( LSAMEN(3,c3,'SV ') ) THEN
!
            IF ( Info/=Infoe .AND. Infoe/=0 ) THEN
               WRITE (Nout,FMT=99014) Subnam(1:LEN_TRIM(Subnam)) ,      &
     &                                Info , Infoe , N , Kl , Ku , N5 , &
     &                                Imat
            ELSE
               WRITE (Nout,FMT=99028) Subnam(1:LEN_TRIM(Subnam)) ,      &
     &                                Info , N , Kl , Ku , N5 , Imat
            ENDIF
!
         ELSEIF ( LSAMEN(3,c3,'SVX') ) THEN
!
            IF ( Info/=Infoe .AND. Infoe/=0 ) THEN
               WRITE (Nout,FMT=99007) Subnam(1:LEN_TRIM(Subnam)) ,      &
     &                                Info , Infoe , Opts(1:1) ,        &
     &                                Opts(2:2) , N , Kl , Ku , N5 ,    &
     &                                Imat
            ELSE
               WRITE (Nout,FMT=99002) Subnam(1:LEN_TRIM(Subnam)) ,      &
     &                                Info , Opts(1:1) , Opts(2:2) , N ,&
     &                                Kl , Ku , N5 , Imat
            ENDIF
!
         ELSEIF ( LSAMEN(5,Subnam(2:6),'LATMS') ) THEN
!
            WRITE (Nout,FMT=99023) Subnam(1:LEN_TRIM(Subnam)) , Info ,  &
     &                             M , N , Kl , Ku , Imat
!
         ELSEIF ( LSAMEN(3,c3,'CON') ) THEN
!
            WRITE (Nout,FMT=99032) Subnam(1:LEN_TRIM(Subnam)) , Info ,  &
     &                             Opts(1:1) , M , Kl , Ku , Imat
!
         ELSE
!
            WRITE (Nout,FMT=99036) Subnam(1:LEN_TRIM(Subnam)) , Info ,  &
     &                             Opts(1:1) , M , Kl , Ku , N5 , Imat
         ENDIF
!
      ELSEIF ( LSAMEN(2,p2,'GT') ) THEN
!
!        xGT:  General tridiagonal matrices
!
         IF ( LSAMEN(3,c3,'TRF') ) THEN
            IF ( Info/=Infoe .AND. Infoe/=0 ) THEN
               WRITE (Nout,FMT=99013) Subnam(1:LEN_TRIM(Subnam)) ,      &
     &                                Info , Infoe , N , Imat
            ELSE
               WRITE (Nout,FMT=99027) Subnam(1:LEN_TRIM(Subnam)) ,      &
     &                                Info , N , Imat
            ENDIF
            IF ( Info/=0 ) WRITE (Nout,FMT=99051)
!
         ELSEIF ( LSAMEN(3,c3,'SV ') ) THEN
!
            IF ( Info/=Infoe .AND. Infoe/=0 ) THEN
               WRITE (Nout,FMT=99016) Subnam(1:LEN_TRIM(Subnam)) ,      &
     &                                Info , Infoe , N , N5 , Imat
            ELSE
               WRITE (Nout,FMT=99030) Subnam(1:LEN_TRIM(Subnam)) ,      &
     &                                Info , N , N5 , Imat
            ENDIF
!
         ELSEIF ( LSAMEN(3,c3,'SVX') ) THEN
!
            IF ( Info/=Infoe .AND. Infoe/=0 ) THEN
               WRITE (Nout,FMT=99008) Subnam(1:LEN_TRIM(Subnam)) ,      &
     &                                Info , Infoe , Opts(1:1) ,        &
     &                                Opts(2:2) , N , N5 , Imat
            ELSE
               WRITE (Nout,FMT=99003) Subnam(1:LEN_TRIM(Subnam)) ,      &
     &                                Info , Opts(1:1) , Opts(2:2) , N ,&
     &                                N5 , Imat
            ENDIF
!
         ELSEIF ( LSAMEN(3,c3,'CON') ) THEN
!
            WRITE (Nout,FMT=99031) Subnam(1:LEN_TRIM(Subnam)) , Info ,  &
     &                             Opts(1:1) , M , Imat
!
         ELSE
!
            WRITE (Nout,FMT=99037) Subnam(1:LEN_TRIM(Subnam)) , Info ,  &
     &                             Opts(1:1) , M , N5 , Imat
         ENDIF
!
      ELSEIF ( LSAMEN(2,p2,'PO') ) THEN
!
!        xPO:  Symmetric or Hermitian positive definite matrices
!
         uplo = Opts(1:1)
         IF ( LSAMEN(3,c3,'TRF') ) THEN
            IF ( Info/=Infoe .AND. Infoe/=0 ) THEN
               WRITE (Nout,FMT=99020) Subnam(1:LEN_TRIM(Subnam)) ,      &
     &                                Info , Infoe , uplo , M , N5 ,    &
     &                                Imat
            ELSE
               WRITE (Nout,FMT=99044) Subnam(1:LEN_TRIM(Subnam)) ,      &
     &                                Info , uplo , M , N5 , Imat
            ENDIF
            IF ( Info/=0 ) WRITE (Nout,FMT=99051)
!
         ELSEIF ( LSAMEN(3,c3,'SV ') ) THEN
!
            IF ( Info/=Infoe .AND. Infoe/=0 ) THEN
               WRITE (Nout,FMT=99021) Subnam(1:LEN_TRIM(Subnam)) ,      &
     &                                Info , Infoe , uplo , N , N5 ,    &
     &                                Imat
            ELSE
               WRITE (Nout,FMT=99045) Subnam(1:LEN_TRIM(Subnam)) ,      &
     &                                Info , uplo , N , N5 , Imat
            ENDIF
!
         ELSEIF ( LSAMEN(3,c3,'SVX') ) THEN
!
            IF ( Info/=Infoe .AND. Infoe/=0 ) THEN
               WRITE (Nout,FMT=99010) Subnam(1:LEN_TRIM(Subnam)) ,      &
     &                                Info , Infoe , Opts(1:1) ,        &
     &                                Opts(2:2) , N , N5 , Imat
            ELSE
               WRITE (Nout,FMT=99005) Subnam(1:LEN_TRIM(Subnam)) ,      &
     &                                Info , Opts(1:1) , Opts(2:2) , N ,&
     &                                N5 , Imat
            ENDIF
!
         ELSEIF ( LSAMEN(3,c3,'TRI') ) THEN
!
            WRITE (Nout,FMT=99044) Subnam(1:LEN_TRIM(Subnam)) , Info ,  &
     &                             uplo , M , N5 , Imat
!
         ELSEIF ( LSAMEN(5,Subnam(2:6),'LATMS') .OR. LSAMEN(3,c3,'CON') &
     &            ) THEN
!
            WRITE (Nout,FMT=99040) Subnam(1:LEN_TRIM(Subnam)) , Info ,  &
     &                             uplo , M , Imat
!
         ELSE
!
            WRITE (Nout,FMT=99045) Subnam(1:LEN_TRIM(Subnam)) , Info ,  &
     &                             uplo , M , N5 , Imat
         ENDIF
!
      ELSEIF ( LSAMEN(2,p2,'PS') ) THEN
!
!        xPS:  Symmetric or Hermitian positive semi-definite matrices
!
         uplo = Opts(1:1)
         IF ( LSAMEN(3,c3,'TRF') ) THEN
            IF ( Info/=Infoe .AND. Infoe/=0 ) THEN
               WRITE (Nout,FMT=99020) Subnam , Info , Infoe , uplo , M ,&
     &                                N5 , Imat
            ELSE
               WRITE (Nout,FMT=99044) Subnam , Info , uplo , M , N5 ,   &
     &                                Imat
            ENDIF
            IF ( Info/=0 ) WRITE (Nout,FMT=99051)
!
         ELSEIF ( LSAMEN(3,c3,'SV ') ) THEN
!
            IF ( Info/=Infoe .AND. Infoe/=0 ) THEN
               WRITE (Nout,FMT=99021) Subnam , Info , Infoe , uplo , N ,&
     &                                N5 , Imat
            ELSE
               WRITE (Nout,FMT=99045) Subnam , Info , uplo , N , N5 ,   &
     &                                Imat
            ENDIF
!
         ELSEIF ( LSAMEN(3,c3,'SVX') ) THEN
!
            IF ( Info/=Infoe .AND. Infoe/=0 ) THEN
               WRITE (Nout,FMT=99010) Subnam , Info , Infoe , Opts(1:1) &
     &                                , Opts(2:2) , N , N5 , Imat
            ELSE
               WRITE (Nout,FMT=99005) Subnam , Info , Opts(1:1) ,       &
     &                                Opts(2:2) , N , N5 , Imat
            ENDIF
!
         ELSEIF ( LSAMEN(3,c3,'TRI') ) THEN
!
            WRITE (Nout,FMT=99044) Subnam , Info , uplo , M , N5 , Imat
!
         ELSEIF ( LSAMEN(5,Subnam(2:6),'LATMT') .OR. LSAMEN(3,c3,'CON') &
     &            ) THEN
!
            WRITE (Nout,FMT=99040) Subnam , Info , uplo , M , Imat
!
         ELSE
!
            WRITE (Nout,FMT=99045) Subnam , Info , uplo , M , N5 , Imat
         ENDIF
!
      ELSEIF ( LSAMEN(2,p2,'SY') .OR. LSAMEN(2,p2,'SR') .OR.            &
     &         LSAMEN(2,p2,'SK') .OR. LSAMEN(2,p2,'HE') .OR.            &
     &         LSAMEN(2,p2,'HR') .OR. LSAMEN(2,p2,'HK') .OR.            &
     &         LSAMEN(2,p2,'HA') ) THEN
!
!        xSY: symmetric indefinite matrices
!             with partial (Bunch-Kaufman) pivoting;
!        xSR: symmetric indefinite matrices
!             with rook (bounded Bunch-Kaufman) pivoting;
!        xSK: symmetric indefinite matrices
!             with rook (bounded Bunch-Kaufman) pivoting,
!             new storage format;
!        xHE: Hermitian indefinite matrices
!             with partial (Bunch-Kaufman) pivoting.
!        xHR: Hermitian indefinite matrices
!             with rook (bounded Bunch-Kaufman) pivoting;
!        xHK: Hermitian indefinite matrices
!             with rook (bounded Bunch-Kaufman) pivoting,
!             new storage format;
!        xHA: Hermitian matrices
!             Aasen Algorithm
!
         uplo = Opts(1:1)
         IF ( LSAMEN(3,c3,'TRF') ) THEN
            IF ( Info/=Infoe .AND. Infoe/=0 ) THEN
               WRITE (Nout,FMT=99020) Subnam(1:LEN_TRIM(Subnam)) ,      &
     &                                Info , Infoe , uplo , M , N5 ,    &
     &                                Imat
            ELSE
               WRITE (Nout,FMT=99044) Subnam(1:LEN_TRIM(Subnam)) ,      &
     &                                Info , uplo , M , N5 , Imat
            ENDIF
            IF ( Info/=0 ) WRITE (Nout,FMT=99051)
!
         ELSEIF ( LSAMEN(2,c3,'SV') ) THEN
!
            IF ( Info/=Infoe .AND. Infoe/=0 ) THEN
               WRITE (Nout,FMT=99021) Subnam(1:LEN_TRIM(Subnam)) ,      &
     &                                Info , Infoe , uplo , N , N5 ,    &
     &                                Imat
            ELSE
               WRITE (Nout,FMT=99045) Subnam(1:LEN_TRIM(Subnam)) ,      &
     &                                Info , uplo , N , N5 , Imat
            ENDIF
!
         ELSEIF ( LSAMEN(3,c3,'SVX') ) THEN
!
            IF ( Info/=Infoe .AND. Infoe/=0 ) THEN
               WRITE (Nout,FMT=99010) Subnam(1:LEN_TRIM(Subnam)) ,      &
     &                                Info , Infoe , Opts(1:1) ,        &
     &                                Opts(2:2) , N , N5 , Imat
            ELSE
               WRITE (Nout,FMT=99005) Subnam(1:LEN_TRIM(Subnam)) ,      &
     &                                Info , Opts(1:1) , Opts(2:2) , N ,&
     &                                N5 , Imat
            ENDIF
!
         ELSEIF ( LSAMEN(5,Subnam(2:6),'LATMS') .OR. LSAMEN(3,c3,'TRI') &
     &            .OR. LSAMEN(3,c3,'CON') ) THEN
!
            WRITE (Nout,FMT=99040) Subnam(1:LEN_TRIM(Subnam)) , Info ,  &
     &                             uplo , M , Imat
!
         ELSE
!
            WRITE (Nout,FMT=99045) Subnam(1:LEN_TRIM(Subnam)) , Info ,  &
     &                             uplo , M , N5 , Imat
         ENDIF
!
      ELSEIF ( LSAMEN(2,p2,'PP') .OR. LSAMEN(2,p2,'SP') .OR.            &
     &         LSAMEN(2,p2,'HP') ) THEN
!
!        xPP, xHP, or xSP:  Symmetric or Hermitian packed matrices
!
         uplo = Opts(1:1)
         IF ( LSAMEN(3,c3,'TRF') ) THEN
            IF ( Info/=Infoe .AND. Infoe/=0 ) THEN
               WRITE (Nout,FMT=99017) Subnam(1:LEN_TRIM(Subnam)) ,      &
     &                                Info , Infoe , uplo , M , Imat
            ELSE
               WRITE (Nout,FMT=99040) Subnam(1:LEN_TRIM(Subnam)) ,      &
     &                                Info , uplo , M , Imat
            ENDIF
            IF ( Info/=0 ) WRITE (Nout,FMT=99051)
!
         ELSEIF ( LSAMEN(3,c3,'SV ') ) THEN
!
            IF ( Info/=Infoe .AND. Infoe/=0 ) THEN
               WRITE (Nout,FMT=99021) Subnam(1:LEN_TRIM(Subnam)) ,      &
     &                                Info , Infoe , uplo , N , N5 ,    &
     &                                Imat
            ELSE
               WRITE (Nout,FMT=99045) Subnam(1:LEN_TRIM(Subnam)) ,      &
     &                                Info , uplo , N , N5 , Imat
            ENDIF
!
         ELSEIF ( LSAMEN(3,c3,'SVX') ) THEN
!
            IF ( Info/=Infoe .AND. Infoe/=0 ) THEN
               WRITE (Nout,FMT=99010) Subnam(1:LEN_TRIM(Subnam)) ,      &
     &                                Info , Infoe , Opts(1:1) ,        &
     &                                Opts(2:2) , N , N5 , Imat
            ELSE
               WRITE (Nout,FMT=99005) Subnam(1:LEN_TRIM(Subnam)) ,      &
     &                                Info , Opts(1:1) , Opts(2:2) , N ,&
     &                                N5 , Imat
            ENDIF
!
         ELSEIF ( LSAMEN(5,Subnam(2:6),'LATMS') .OR. LSAMEN(3,c3,'TRI') &
     &            .OR. LSAMEN(3,c3,'CON') ) THEN
!
            WRITE (Nout,FMT=99040) Subnam(1:LEN_TRIM(Subnam)) , Info ,  &
     &                             uplo , M , Imat
!
         ELSE
!
            WRITE (Nout,FMT=99045) Subnam(1:LEN_TRIM(Subnam)) , Info ,  &
     &                             uplo , M , N5 , Imat
         ENDIF
!
      ELSEIF ( LSAMEN(2,p2,'PB') ) THEN
!
!        xPB:  Symmetric (Hermitian) positive definite band matrix
!
         uplo = Opts(1:1)
         IF ( LSAMEN(3,c3,'TRF') ) THEN
            IF ( Info/=Infoe .AND. Infoe/=0 ) THEN
               WRITE (Nout,FMT=99018) Subnam(1:LEN_TRIM(Subnam)) ,      &
     &                                Info , Infoe , uplo , M , Kl ,    &
     &                                N5 , Imat
            ELSE
               WRITE (Nout,FMT=99042) Subnam(1:LEN_TRIM(Subnam)) ,      &
     &                                Info , uplo , M , Kl , N5 , Imat
            ENDIF
            IF ( Info/=0 ) WRITE (Nout,FMT=99051)
!
         ELSEIF ( LSAMEN(3,c3,'SV ') ) THEN
!
            IF ( Info/=Infoe .AND. Infoe/=0 ) THEN
               WRITE (Nout,FMT=99019) Subnam(1:LEN_TRIM(Subnam)) ,      &
     &                                Info , Infoe , uplo , N , Kl ,    &
     &                                N5 , Imat
            ELSE
               WRITE (Nout,FMT=99043) Subnam(1:LEN_TRIM(Subnam)) ,      &
     &                                Info , uplo , N , Kl , N5 , Imat
            ENDIF
!
         ELSEIF ( LSAMEN(3,c3,'SVX') ) THEN
!
            IF ( Info/=Infoe .AND. Infoe/=0 ) THEN
               WRITE (Nout,FMT=99009) Subnam(1:LEN_TRIM(Subnam)) ,      &
     &                                Info , Infoe , Opts(1:1) ,        &
     &                                Opts(2:2) , N , Kl , N5 , Imat
            ELSE
               WRITE (Nout,FMT=99004) Subnam(1:LEN_TRIM(Subnam)) ,      &
     &                                Info , Opts(1:1) , Opts(2:2) , N ,&
     &                                Kl , N5 , Imat
            ENDIF
!
         ELSEIF ( LSAMEN(5,Subnam(2:6),'LATMS') .OR. LSAMEN(3,c3,'CON') &
     &            ) THEN
!
            WRITE (Nout,FMT=99041) Subnam(1:LEN_TRIM(Subnam)) , Info ,  &
     &                             uplo , M , Kl , Imat
!
         ELSE
!
            WRITE (Nout,FMT=99043) Subnam(1:LEN_TRIM(Subnam)) , Info ,  &
     &                             uplo , M , Kl , N5 , Imat
         ENDIF
!
      ELSEIF ( LSAMEN(2,p2,'PT') ) THEN
!
!        xPT:  Positive definite tridiagonal matrices
!
         IF ( LSAMEN(3,c3,'TRF') ) THEN
            IF ( Info/=Infoe .AND. Infoe/=0 ) THEN
               WRITE (Nout,FMT=99013) Subnam(1:LEN_TRIM(Subnam)) ,      &
     &                                Info , Infoe , N , Imat
            ELSE
               WRITE (Nout,FMT=99027) Subnam(1:LEN_TRIM(Subnam)) ,      &
     &                                Info , N , Imat
            ENDIF
            IF ( Info/=0 ) WRITE (Nout,FMT=99051)
!
         ELSEIF ( LSAMEN(3,c3,'SV ') ) THEN
!
            IF ( Info/=Infoe .AND. Infoe/=0 ) THEN
               WRITE (Nout,FMT=99016) Subnam(1:LEN_TRIM(Subnam)) ,      &
     &                                Info , Infoe , N , N5 , Imat
            ELSE
               WRITE (Nout,FMT=99030) Subnam(1:LEN_TRIM(Subnam)) ,      &
     &                                Info , N , N5 , Imat
            ENDIF
!
         ELSEIF ( LSAMEN(3,c3,'SVX') ) THEN
!
            IF ( Info/=Infoe .AND. Infoe/=0 ) THEN
               WRITE (Nout,FMT=99006) Subnam(1:LEN_TRIM(Subnam)) ,      &
     &                                Info , Infoe , Opts(1:1) , N ,    &
     &                                N5 , Imat
            ELSE
               WRITE (Nout,FMT=99001) Subnam(1:LEN_TRIM(Subnam)) ,      &
     &                                Info , Opts(1:1) , N , N5 , Imat
            ENDIF
!
         ELSEIF ( LSAMEN(3,c3,'CON') ) THEN
!
            IF ( LSAME(Subnam(1:1),'S') .OR. LSAME(Subnam(1:1),'D') )   &
     &           THEN
               WRITE (Nout,FMT=99027) Subnam(1:LEN_TRIM(Subnam)) ,      &
     &                                Info , M , Imat
            ELSE
               WRITE (Nout,FMT=99031) Subnam(1:LEN_TRIM(Subnam)) ,      &
     &                                Info , Opts(1:1) , M , Imat
            ENDIF
!
         ELSE
!
            WRITE (Nout,FMT=99037) Subnam(1:LEN_TRIM(Subnam)) , Info ,  &
     &                             Opts(1:1) , M , N5 , Imat
         ENDIF
!
      ELSEIF ( LSAMEN(2,p2,'TR') ) THEN
!
!        xTR:  Triangular matrix
!
         IF ( LSAMEN(3,c3,'TRI') ) THEN
            WRITE (Nout,FMT=99039) Subnam(1:LEN_TRIM(Subnam)) , Info ,  &
     &                             Opts(1:1) , Opts(2:2) , M , N5 , Imat
         ELSEIF ( LSAMEN(3,c3,'CON') ) THEN
            WRITE (Nout,FMT=99033) Subnam(1:LEN_TRIM(Subnam)) , Info ,  &
     &                             Opts(1:1) , Opts(2:2) , Opts(3:3) ,  &
     &                             M , Imat
         ELSEIF ( LSAMEN(5,Subnam(2:6),'LATRS') ) THEN
            WRITE (Nout,FMT=99048) Subnam(1:LEN_TRIM(Subnam)) , Info ,  &
     &                             Opts(1:1) , Opts(2:2) , Opts(3:3) ,  &
     &                             Opts(4:4) , M , Imat
         ELSE
            WRITE (Nout,FMT=99047) Subnam(1:LEN_TRIM(Subnam)) , Info ,  &
     &                             Opts(1:1) , Opts(2:2) , Opts(3:3) ,  &
     &                             M , N5 , Imat
         ENDIF
!
      ELSEIF ( LSAMEN(2,p2,'TP') ) THEN
!
!        xTP:  Triangular packed matrix
!
         IF ( LSAMEN(3,c3,'TRI') ) THEN
            WRITE (Nout,FMT=99038) Subnam(1:LEN_TRIM(Subnam)) , Info ,  &
     &                             Opts(1:1) , Opts(2:2) , M , Imat
         ELSEIF ( LSAMEN(3,c3,'CON') ) THEN
            WRITE (Nout,FMT=99033) Subnam(1:LEN_TRIM(Subnam)) , Info ,  &
     &                             Opts(1:1) , Opts(2:2) , Opts(3:3) ,  &
     &                             M , Imat
         ELSEIF ( LSAMEN(5,Subnam(2:6),'LATPS') ) THEN
            WRITE (Nout,FMT=99048) Subnam(1:LEN_TRIM(Subnam)) , Info ,  &
     &                             Opts(1:1) , Opts(2:2) , Opts(3:3) ,  &
     &                             Opts(4:4) , M , Imat
         ELSE
            WRITE (Nout,FMT=99047) Subnam(1:LEN_TRIM(Subnam)) , Info ,  &
     &                             Opts(1:1) , Opts(2:2) , Opts(3:3) ,  &
     &                             M , N5 , Imat
         ENDIF
!
      ELSEIF ( LSAMEN(2,p2,'TB') ) THEN
!
!        xTB:  Triangular band matrix
!
         IF ( LSAMEN(3,c3,'CON') ) THEN
            WRITE (Nout,FMT=99034) Subnam(1:LEN_TRIM(Subnam)) , Info ,  &
     &                             Opts(1:1) , Opts(2:2) , Opts(3:3) ,  &
     &                             M , Kl , Imat
         ELSEIF ( LSAMEN(5,Subnam(2:6),'LATBS') ) THEN
            WRITE (Nout,FMT=99049) Subnam(1:LEN_TRIM(Subnam)) , Info ,  &
     &                             Opts(1:1) , Opts(2:2) , Opts(3:3) ,  &
     &                             Opts(4:4) , M , Kl , Imat
         ELSE
            WRITE (Nout,FMT=99046) Subnam(1:LEN_TRIM(Subnam)) , Info ,  &
     &                             Opts(1:1) , Opts(2:2) , Opts(3:3) ,  &
     &                             M , Kl , N5 , Imat
         ENDIF
!
      ELSEIF ( LSAMEN(2,p2,'QR') ) THEN
!
!        xQR:  QR factorization
!
         IF ( LSAMEN(3,c3,'QRS') ) THEN
            WRITE (Nout,FMT=99026) Subnam(1:LEN_TRIM(Subnam)) , Info ,  &
     &                             M , N , Kl , N5 , Imat
         ELSEIF ( LSAMEN(5,Subnam(2:6),'LATMS') ) THEN
            WRITE (Nout,FMT=99022) Subnam(1:LEN_TRIM(Subnam)) , Info ,  &
     &                             M , N , Imat
         ENDIF
!
      ELSEIF ( LSAMEN(2,p2,'LQ') ) THEN
!
!        xLQ:  LQ factorization
!
         IF ( LSAMEN(3,c3,'LQS') ) THEN
            WRITE (Nout,FMT=99026) Subnam(1:LEN_TRIM(Subnam)) , Info ,  &
     &                             M , N , Kl , N5 , Imat
         ELSEIF ( LSAMEN(5,Subnam(2:6),'LATMS') ) THEN
            WRITE (Nout,FMT=99022) Subnam(1:LEN_TRIM(Subnam)) , Info ,  &
     &                             M , N , Imat
         ENDIF
!
      ELSEIF ( LSAMEN(2,p2,'QL') ) THEN
!
!        xQL:  QL factorization
!
         IF ( LSAMEN(3,c3,'QLS') ) THEN
            WRITE (Nout,FMT=99026) Subnam(1:LEN_TRIM(Subnam)) , Info ,  &
     &                             M , N , Kl , N5 , Imat
         ELSEIF ( LSAMEN(5,Subnam(2:6),'LATMS') ) THEN
            WRITE (Nout,FMT=99022) Subnam(1:LEN_TRIM(Subnam)) , Info ,  &
     &                             M , N , Imat
         ENDIF
!
      ELSEIF ( LSAMEN(2,p2,'RQ') ) THEN
!
!        xRQ:  RQ factorization
!
         IF ( LSAMEN(3,c3,'RQS') ) THEN
            WRITE (Nout,FMT=99026) Subnam(1:LEN_TRIM(Subnam)) , Info ,  &
     &                             M , N , Kl , N5 , Imat
         ELSEIF ( LSAMEN(5,Subnam(2:6),'LATMS') ) THEN
            WRITE (Nout,FMT=99022) Subnam(1:LEN_TRIM(Subnam)) , Info ,  &
     &                             M , N , Imat
         ENDIF
!
      ELSEIF ( LSAMEN(2,p2,'LU') ) THEN
!
         IF ( Info/=Infoe .AND. Infoe/=0 ) THEN
            WRITE (Nout,FMT=99012) Subnam(1:LEN_TRIM(Subnam)) , Info ,  &
     &                             Infoe , M , N , N5 , Imat
         ELSE
            WRITE (Nout,FMT=99025) Subnam(1:LEN_TRIM(Subnam)) , Info ,  &
     &                             M , N , N5 , Imat
         ENDIF
!
      ELSEIF ( LSAMEN(2,p2,'CH') ) THEN
!
         IF ( Info/=Infoe .AND. Infoe/=0 ) THEN
            WRITE (Nout,FMT=99015) Subnam(1:LEN_TRIM(Subnam)) , Info ,  &
     &                             Infoe , M , N5 , Imat
         ELSE
            WRITE (Nout,FMT=99029) Subnam(1:LEN_TRIM(Subnam)) , Info ,  &
     &                             M , N5 , Imat
         ENDIF
!
      ELSE
!
!        Print a generic message if the path is unknown.
!
         WRITE (Nout,FMT=99050) Subnam(1:LEN_TRIM(Subnam)) , Info
      ENDIF
!
!     Description of error message (alphabetical, left to right)
!
!     SUBNAM, INFO, FACT, N, NRHS, IMAT
!
99001 FORMAT (' *** Error code from ',A,'=',I5,', FACT=''',A1,''', N=', &
     &        I5,', NRHS=',I4,', type ',I2)
!
!     SUBNAM, INFO, FACT, TRANS, N, KL, KU, NRHS, IMAT
!
99002 FORMAT (' *** Error code from ',A,' =',I5,/' ==> FACT=''',A1,     &
     &        ''', TRANS=''',A1,''', N=',I5,', KL=',I5,', KU=',I5,      &
     &        ', NRHS=',I4,', type ',I1)
!
!     SUBNAM, INFO, FACT, TRANS, N, NRHS, IMAT
!
99003 FORMAT (' *** Error code from ',A,' =',I5,/' ==> FACT=''',A1,     &
     &        ''', TRANS=''',A1,''', N =',I5,', NRHS =',I4,', type ',I2)
!
!     SUBNAM, INFO, FACT, UPLO, N, KD, NRHS, IMAT
!
99004 FORMAT (' *** Error code from ',A,' =',I5,/' ==> FACT=''',A1,     &
     &        ''', UPLO=''',A1,''', N=',I5,', KD=',I5,', NRHS=',I4,     &
     &        ', type ',I2)
!
!     SUBNAM, INFO, FACT, UPLO, N, NRHS, IMAT
!
99005 FORMAT (' *** Error code from ',A,' =',I5,/' ==> FACT=''',A1,     &
     &        ''', UPLO=''',A1,''', N =',I5,', NRHS =',I4,', type ',I2)
!
!     SUBNAM, INFO, INFOE, FACT, N, NRHS, IMAT
!
99006 FORMAT (' *** ',A,' returned with INFO =',I5,' instead of ',I2,   &
     &        /' ==> FACT=''',A1,''', N =',I5,', NRHS =',I4,', type ',  &
     &        I2)
!
!     SUBNAM, INFO, INFOE, FACT, TRANS, N, KL, KU, NRHS, IMAT
!
99007 FORMAT (' *** ',A,' returned with INFO =',I5,' instead of ',I2,   &
     &        /' ==> FACT=''',A1,''', TRANS=''',A1,''', N=',I5,', KL=', &
     &        I5,', KU=',I5,', NRHS=',I4,', type ',I1)
!
!     SUBNAM, INFO, INFOE, FACT, TRANS, N, NRHS, IMAT
!
99008 FORMAT (' *** ',A,' returned with INFO =',I5,' instead of ',I2,   &
     &        /' ==> FACT=''',A1,''', TRANS=''',A1,''', N =',I5,        &
     &        ', NRHS =',I4,', type ',I2)
!
!     SUBNAM, INFO, INFOE, FACT, UPLO, N, KD, NRHS, IMAT
!
99009 FORMAT (' *** ',A,' returned with INFO =',I5,' instead of ',I2,   &
     &        /' ==> FACT=''',A1,''', UPLO=''',A1,''', N=',I5,', KD=',  &
     &        I5,', NRHS=',I4,', type ',I2)
!
!     SUBNAM, INFO, INFOE, FACT, UPLO, N, NRHS, IMAT
!
99010 FORMAT (' *** ',A,' returned with INFO =',I5,' instead of ',I2,   &
     &        /' ==> FACT=''',A1,''', UPLO=''',A1,''', N =',I5,         &
     &        ', NRHS =',I4,', type ',I2)
!
!     SUBNAM, INFO, INFOE, M, N, KL, KU, NB, IMAT
!
99011 FORMAT (' *** ',A,' returned with INFO =',I5,' instead of ',I2,   &
     &        /' ==> M = ',I5,', N =',I5,', KL =',I5,', KU =',I5,       &
     &        ', NB =',I4,', type ',I2)
!
!     SUBNAM, INFO, INFOE, M, N, NB, IMAT
!
99012 FORMAT (' *** ',A,' returned with INFO =',I5,' instead of ',I2,   &
     &        /' ==> M =',I5,', N =',I5,', NB =',I4,', type ',I2)
!
!     SUBNAM, INFO, INFOE, N, IMAT
!
99013 FORMAT (' *** ',A,' returned with INFO =',I5,' instead of ',I2,   &
     &        ' for N=',I5,', type ',I2)
!
!     SUBNAM, INFO, INFOE, N, KL, KU, NRHS, IMAT
!
99014 FORMAT (' *** ',A,' returned with INFO =',I5,' instead of ',I2,   &
     &        /' ==> N =',I5,', KL =',I5,', KU =',I5,', NRHS =',I4,     &
     &        ', type ',I2)
!
!     SUBNAM, INFO, INFOE, N, NB, IMAT
!
99015 FORMAT (' *** ',A,' returned with INFO =',I5,' instead of ',I2,   &
     &        /' ==> N =',I5,', NB =',I4,', type ',I2)
!
!     SUBNAM, INFO, INFOE, N, NRHS, IMAT
!
99016 FORMAT (' *** ',A,' returned with INFO =',I5,' instead of ',I2,   &
     &        /' ==> N =',I5,', NRHS =',I4,', type ',I2)
!
!     SUBNAM, INFO, INFOE, UPLO, N, IMAT
!
99017 FORMAT (' *** ',A,' returned with INFO =',I5,' instead of ',I2,   &
     &        /' ==> UPLO = ''',A1,''', N =',I5,', type ',I2)
!
!     SUBNAM, INFO, INFOE, UPLO, N, KD, NB, IMAT
!
99018 FORMAT (' *** ',A,' returned with INFO =',I5,' instead of ',I2,   &
     &        /' ==> UPLO = ''',A1,''', N =',I5,', KD =',I5,', NB =',I4,&
     &        ', type ',I2)
!
!     SUBNAM, INFO, INFOE, UPLO, N, KD, NRHS, IMAT
!
99019 FORMAT (' *** ',A,' returned with INFO =',I5,' instead of ',I2,   &
     &        /' ==> UPLO=''',A1,''', N =',I5,', KD =',I5,', NRHS =',I4,&
     &        ', type ',I2)
!
!     SUBNAM, INFO, INFOE, UPLO, N, NB, IMAT
!
99020 FORMAT (' *** ',A,' returned with INFO =',I5,' instead of ',I2,   &
     &        /' ==> UPLO = ''',A1,''', N =',I5,', NB =',I4,', type ',  &
     &        I2)
!
!     SUBNAM, INFO, INFOE, UPLO, N, NRHS, IMAT
!
99021 FORMAT (' *** ',A,' returned with INFO =',I5,' instead of ',I2,   &
     &        /' ==> UPLO = ''',A1,''', N =',I5,', NRHS =',I4,', type ',&
     &        I2)
!
!     SUBNAM, INFO, M, N, IMAT
!
99022 FORMAT (' *** Error code from ',A,' =',I5,' for M =',I5,', N =',  &
     &        I5,', type ',I2)
!
!     SUBNAM, INFO, M, N, KL, KU, IMAT
!
99023 FORMAT (' *** Error code from ',A,' =',I5,/' ==> M = ',I5,', N =',&
     &        I5,', KL =',I5,', KU =',I5,', type ',I2)
!
!     SUBNAM, INFO, M, N, KL, KU, NB, IMAT
!
99024 FORMAT (' *** Error code from ',A,' =',I5,/' ==> M = ',I5,', N =',&
     &        I5,', KL =',I5,', KU =',I5,', NB =',I4,', type ',I2)
!
!     SUBNAM, INFO, M, N, NB, IMAT
!
99025 FORMAT (' *** Error code from ',A,'=',I5,' for M=',I5,', N=',I5,  &
     &        ', NB=',I4,', type ',I2)
!
!     SUBNAM, INFO, M, N, NRHS, NB, IMAT
!
99026 FORMAT (' *** Error code from ',A,'=',I5,/' ==> M =',I5,', N =',  &
     &        I5,', NRHS =',I4,', NB =',I4,', type ',I2)
!
!     SUBNAM, INFO, N, IMAT
!
99027 FORMAT (' *** Error code from ',A,' =',I5,' for N =',I5,', type ',&
     &        I2)
!
!     SUBNAM, INFO, N, KL, KU, NRHS, IMAT
!
99028 FORMAT (' *** Error code from ',A,' =',I5,/' ==> N =',I5,', KL =',&
     &        I5,', KU =',I5,', NRHS =',I4,', type ',I2)
!
!     SUBNAM, INFO, N, NB, IMAT
!
99029 FORMAT (' *** Error code from ',A,'=',I5,' for N=',I5,', NB=',I4, &
     &        ', type ',I2)
!
!     SUBNAM, INFO, N, NRHS, IMAT
!
99030 FORMAT (' *** Error code from ',A,' =',I5,' for N =',I5,          &
     &        ', NRHS =',I4,', type ',I2)
!
!     SUBNAM, INFO, NORM, N, IMAT
!
99031 FORMAT (' *** Error code from ',A,' =',I5,' for NORM = ''',A1,    &
     &        ''', N =',I5,', type ',I2)
!
!     SUBNAM, INFO, NORM, N, KL, KU, IMAT
!
99032 FORMAT (' *** Error code from ',A,' =',I5,/' ==> NORM =''',A1,    &
     &        ''', N =',I5,', KL =',I5,', KU =',I5,', type ',I2)
!
!     SUBNAM, INFO, NORM, UPLO, DIAG, N, IMAT
!
99033 FORMAT (' *** Error code from ',A,' =',I5,/' ==> NORM=''',A1,     &
     &        ''', UPLO =''',A1,''', DIAG=''',A1,''', N =',I5,', type ',&
     &        I2)
!
!     SUBNAM, INFO, NORM, UPLO, DIAG, N, KD, IMAT
!
99034 FORMAT (' *** Error code from ',A,' =',I5,/' ==> NORM=''',A1,     &
     &        ''', UPLO =''',A1,''', DIAG=''',A1,''', N=',I5,', KD=',I5,&
     &        ', type ',I2)
!
!     SUBNAM, INFO, TRANS, M, N, NRHS, NB, IMAT
!
99035 FORMAT (' *** Error code from ',A,' =',I5,/' ==> TRANS = ''',A1,  &
     &        ''', M =',I5,', N =',I5,', NRHS =',I4,', NB =',I4,        &
     &        ', type ',I2)
!
!     SUBNAM, INFO, TRANS, N, KL, KU, NRHS, IMAT
!
99036 FORMAT (' *** Error code from ',A,'=',I5,/' ==> TRANS=''',A1,     &
     &        ''', N =',I5,', KL =',I5,', KU =',I5,', NRHS =',I4,       &
     &        ', type ',I2)
!
!     SUBNAM, INFO, TRANS, N, NRHS, IMAT
!
99037 FORMAT (' *** Error code from ',A,' =',I5,/' ==> TRANS = ''',A1,  &
     &        ''', N =',I5,', NRHS =',I4,', type ',I2)
!
!     SUBNAM, INFO, UPLO, DIAG, N, IMAT
!
99038 FORMAT (' *** Error code from ',A,' =',I5,/' ==> UPLO=''',A1,     &
     &        ''', DIAG =''',A1,''', N =',I5,', type ',I2)
!
!     SUBNAM, INFO, UPLO, DIAG, N, NB, IMAT
!
99039 FORMAT (' *** Error code from ',A,' =',I5,/' ==> UPLO=''',A1,     &
     &        ''', DIAG =''',A1,''', N =',I5,', NB =',I4,', type ',I2)
!
!     SUBNAM, INFO, UPLO, N, IMAT
!
99040 FORMAT (' *** Error code from ',A,' =',I5,' for UPLO = ''',A1,    &
     &        ''', N =',I5,', type ',I2)
!
!     SUBNAM, INFO, UPLO, N, KD, IMAT
!
99041 FORMAT (' *** Error code from ',A,' =',I5,/' ==> UPLO = ''',A1,   &
     &        ''', N =',I5,', KD =',I5,', type ',I2)
!
!     SUBNAM, INFO, UPLO, N, KD, NB, IMAT
!
99042 FORMAT (' *** Error code from ',A,' =',I5,/' ==> UPLO = ''',A1,   &
     &        ''', N =',I5,', KD =',I5,', NB =',I4,', type ',I2)
!
!     SUBNAM, INFO, UPLO, N, KD, NRHS, IMAT
!
99043 FORMAT (' *** Error code from ',A,'=',I5,/' ==> UPLO = ''',A1,    &
     &        ''', N =',I5,', KD =',I5,', NRHS =',I4,', type ',I2)
!
!     SUBNAM, INFO, UPLO, N, NB, IMAT
!
99044 FORMAT (' *** Error code from ',A,' =',I5,/' ==> UPLO = ''',A1,   &
     &        ''', N =',I5,', NB =',I4,', type ',I2)
!
!     SUBNAM, INFO, UPLO, N, NRHS, IMAT
!
99045 FORMAT (' *** Error code from ',A,' =',I5,/' ==> UPLO = ''',A1,   &
     &        ''', N =',I5,', NRHS =',I4,', type ',I2)
!
!     SUBNAM, INFO, UPLO, TRANS, DIAG, N, KD, NRHS, IMAT
!
99046 FORMAT (' *** Error code from ',A,' =',I5,/' ==> UPLO=''',A1,     &
     &        ''', TRANS=''',A1,''', DIAG=''',A1,''', N=',I5,', KD=',I5,&
     &        ', NRHS=',I4,', type ',I2)
!
!     SUBNAM, INFO, UPLO, TRANS, DIAG, N, NRHS, IMAT
!
99047 FORMAT (' *** Error code from ',A,' =',I5,/' ==> UPLO=''',A1,     &
     &        ''', TRANS=''',A1,''', DIAG=''',A1,''', N =',I5,          &
     &        ', NRHS =',I4,', type ',I2)
!
!     SUBNAM, INFO, UPLO, TRANS, DIAG, NORMIN, N, IMAT
!
99048 FORMAT (' *** Error code from ',A,' =',I5,/' ==> UPLO=''',A1,     &
     &        ''', TRANS=''',A1,''', DIAG=''',A1,''', NORMIN=''',A1,    &
     &        ''', N =',I5,', type ',I2)
!
!     SUBNAM, INFO, UPLO, TRANS, DIAG, NORMIN, N, KD, IMAT
!
99049 FORMAT (' *** Error code from ',A,' =',I5,/' ==> UPLO=''',A1,     &
     &        ''', TRANS=''',A1,''', DIAG=''',A1,''', NORMIN=''',A1,    &
     &        ''', N=',I5,', KD=',I5,', type ',I2)
!
!     Unknown type
!
99050 FORMAT (' *** Error code from ',A,' =',I5)
!
!     What we do next
!
99051 FORMAT (' ==> Doing only the condition estimate for this case')
!
!
!     End of ALAERH
!
      END SUBROUTINE ALAERH
