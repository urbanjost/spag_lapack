!*==ilaenv.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b ILAENV
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ILAENV + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ilaenv.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ilaenv.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ilaenv.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       INTEGER FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
!
!       .. Scalar Arguments ..
!       CHARACTER*( * )    NAME, OPTS
!       INTEGER            ISPEC, N1, N2, N3, N4
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ILAENV is called from the LAPACK routines to choose problem-dependent
!> parameters for the local environment.  See ISPEC for a description of
!> the parameters.
!>
!> ILAENV returns an INTEGER
!> if ILAENV >= 0: ILAENV returns the value of the parameter specified by ISPEC
!> if ILAENV < 0:  if ILAENV = -k, the k-th argument had an illegal value.
!>
!> This version provides a set of parameters which should give good,
!> but not optimal, performance on many of the currently available
!> computers.  Users are encouraged to modify this subroutine to set
!> the tuning parameters for their particular machine using the option
!> and problem size information in the arguments.
!>
!> This routine will not function correctly if it is converted to all
!> lower case.  Converting it to all upper case is allowed.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] ISPEC
!> \verbatim
!>          ISPEC is INTEGER
!>          Specifies the parameter to be returned as the value of
!>          ILAENV.
!>          = 1: the optimal blocksize; if this value is 1, an unblocked
!>               algorithm will give the best performance.
!>          = 2: the minimum block size for which the block routine
!>               should be used; if the usable block size is less than
!>               this value, an unblocked routine should be used.
!>          = 3: the crossover point (in a block routine, for N less
!>               than this value, an unblocked routine should be used)
!>          = 4: the number of shifts, used in the nonsymmetric
!>               eigenvalue routines (DEPRECATED)
!>          = 5: the minimum column dimension for blocking to be used;
!>               rectangular blocks must have dimension at least k by m,
!>               where k is given by ILAENV(2,...) and m by ILAENV(5,...)
!>          = 6: the crossover point for the SVD (when reducing an m by n
!>               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds
!>               this value, a QR factorization is used first to reduce
!>               the matrix to a triangular form.)
!>          = 7: the number of processors
!>          = 8: the crossover point for the multishift QR method
!>               for nonsymmetric eigenvalue problems (DEPRECATED)
!>          = 9: maximum size of the subproblems at the bottom of the
!>               computation tree in the divide-and-conquer algorithm
!>               (used by xGELSD and xGESDD)
!>          =10: ieee infinity and NaN arithmetic can be trusted not to trap
!>          =11: infinity arithmetic can be trusted not to trap
!>          12 <= ISPEC <= 16:
!>               xHSEQR or related subroutines,
!>               see IPARMQ for detailed explanation
!> \endverbatim
!>
!> \param[in] NAME
!> \verbatim
!>          NAME is CHARACTER*(*)
!>          The name of the calling subroutine, in either upper case or
!>          lower case.
!> \endverbatim
!>
!> \param[in] OPTS
!> \verbatim
!>          OPTS is CHARACTER*(*)
!>          The character options to the subroutine NAME, concatenated
!>          into a single character string.  For example, UPLO = 'U',
!>          TRANS = 'T', and DIAG = 'N' for a triangular routine would
!>          be specified as OPTS = 'UTN'.
!> \endverbatim
!>
!> \param[in] N1
!> \verbatim
!>          N1 is INTEGER
!> \endverbatim
!>
!> \param[in] N2
!> \verbatim
!>          N2 is INTEGER
!> \endverbatim
!>
!> \param[in] N3
!> \verbatim
!>          N3 is INTEGER
!> \endverbatim
!>
!> \param[in] N4
!> \verbatim
!>          N4 is INTEGER
!>          Problem dimensions for the subroutine NAME; these may not all
!>          be required.
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
!> \date November 2019
!
!> \ingroup OTHERauxiliary
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The following conventions have been used when calling ILAENV from the
!>  LAPACK routines:
!>  1)  OPTS is a concatenation of all of the character options to
!>      subroutine NAME, in the same order that they appear in the
!>      argument list for NAME, even if they are not used in determining
!>      the value of the parameter specified by ISPEC.
!>  2)  The problem dimensions N1, N2, N3, N4 are specified in the order
!>      that they appear in the argument list for NAME.  N1 is used
!>      first, N2 second, and so on, and unused problem dimensions are
!>      passed a value of -1.
!>  3)  The parameter value returned by ILAENV is checked for validity in
!>      the calling subroutine.  For example, ILAENV is used to retrieve
!>      the optimal blocksize for STRTRI as follows:
!>
!>      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )
!>      IF( NB.LE.1 ) NB = MAX( 1, N )
!> \endverbatim
!>
!  =====================================================================
      INTEGER FUNCTION ILAENV(Ispec,Name,Opts,N1,N2,N3,N4)
      IMPLICIT NONE
!*--ILAENV166
!
!  -- LAPACK auxiliary routine (version 3.9.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2019
!
!     .. Scalar Arguments ..
      CHARACTER*(*) Name , Opts
      INTEGER Ispec , N1 , N2 , N3 , N4
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      INTEGER i , ic , iz , nb , nbmin , nx
      LOGICAL cname , sname , twostage
      CHARACTER c1*1 , c2*2 , c4*2 , c3*3 , subnam*16
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC CHAR , ICHAR , INT , MIN , REAL
!     ..
!     .. External Functions ..
      INTEGER IEEECK , IPARMQ , IPARAM2STAGE
      EXTERNAL IEEECK , IPARMQ , IPARAM2STAGE
!     ..
!     .. Executable Statements ..
!
      IF ( Ispec==1 .OR. Ispec==2 .OR. Ispec==3 ) THEN
!
!
!     Convert NAME to upper case if the first character is lower case.
!
         ILAENV = 1
         subnam = Name
         ic = ICHAR(subnam(1:1))
         iz = ICHAR('Z')
         IF ( iz==90 .OR. iz==122 ) THEN
!
!        ASCII character set
!
            IF ( ic>=97 .AND. ic<=122 ) THEN
               subnam(1:1) = CHAR(ic-32)
               DO i = 2 , 6
                  ic = ICHAR(subnam(i:i))
                  IF ( ic>=97 .AND. ic<=122 ) subnam(i:i) = CHAR(ic-32)
               ENDDO
            ENDIF
!
         ELSEIF ( iz==233 .OR. iz==169 ) THEN
!
!        EBCDIC character set
!
            IF ( (ic>=129 .AND. ic<=137) .OR. (ic>=145 .AND. ic<=153)   &
     &           .OR. (ic>=162 .AND. ic<=169) ) THEN
               subnam(1:1) = CHAR(ic+64)
               DO i = 2 , 6
                  ic = ICHAR(subnam(i:i))
                  IF ( (ic>=129 .AND. ic<=137) .OR.                     &
     &                 (ic>=145 .AND. ic<=153) .OR.                     &
     &                 (ic>=162 .AND. ic<=169) ) subnam(i:i)            &
     &                 = CHAR(ic+64)
               ENDDO
            ENDIF
!
         ELSEIF ( iz==218 .OR. iz==250 ) THEN
!
!        Prime machines:  ASCII+128
!
            IF ( ic>=225 .AND. ic<=250 ) THEN
               subnam(1:1) = CHAR(ic-32)
               DO i = 2 , 6
                  ic = ICHAR(subnam(i:i))
                  IF ( ic>=225 .AND. ic<=250 ) subnam(i:i) = CHAR(ic-32)
               ENDDO
            ENDIF
         ENDIF
!
         c1 = subnam(1:1)
         sname = c1=='S' .OR. c1=='D'
         cname = c1=='C' .OR. c1=='Z'
         IF ( .NOT.(cname .OR. sname) ) RETURN
         c2 = subnam(2:3)
         c3 = subnam(4:6)
         c4 = c3(2:3)
         twostage = LEN(subnam)>=11 .AND. subnam(11:11)=='2'
!
         IF ( Ispec==2 ) THEN
!
!
!     ISPEC = 2:  minimum block size
!
            nbmin = 2
            IF ( c2=='GE' ) THEN
               IF ( c3=='QRF' .OR. c3=='RQF' .OR. c3=='LQF' .OR.        &
     &              c3=='QLF' ) THEN
                  IF ( sname ) THEN
                     nbmin = 2
                  ELSE
                     nbmin = 2
                  ENDIF
               ELSEIF ( c3=='HRD' ) THEN
                  IF ( sname ) THEN
                     nbmin = 2
                  ELSE
                     nbmin = 2
                  ENDIF
               ELSEIF ( c3=='BRD' ) THEN
                  IF ( sname ) THEN
                     nbmin = 2
                  ELSE
                     nbmin = 2
                  ENDIF
               ELSEIF ( c3=='TRI' ) THEN
                  IF ( sname ) THEN
                     nbmin = 2
                  ELSE
                     nbmin = 2
                  ENDIF
               ENDIF
            ELSEIF ( c2=='SY' ) THEN
               IF ( c3=='TRF' ) THEN
                  IF ( sname ) THEN
                     nbmin = 8
                  ELSE
                     nbmin = 8
                  ENDIF
               ELSEIF ( sname .AND. c3=='TRD' ) THEN
                  nbmin = 2
               ENDIF
            ELSEIF ( cname .AND. c2=='HE' ) THEN
               IF ( c3=='TRD' ) nbmin = 2
            ELSEIF ( sname .AND. c2=='OR' ) THEN
               IF ( c3(1:1)=='G' ) THEN
                  IF ( c4=='QR' .OR. c4=='RQ' .OR. c4=='LQ' .OR.        &
     &                 c4=='QL' .OR. c4=='HR' .OR. c4=='TR' .OR.        &
     &                 c4=='BR' ) nbmin = 2
               ELSEIF ( c3(1:1)=='M' ) THEN
                  IF ( c4=='QR' .OR. c4=='RQ' .OR. c4=='LQ' .OR.        &
     &                 c4=='QL' .OR. c4=='HR' .OR. c4=='TR' .OR.        &
     &                 c4=='BR' ) nbmin = 2
               ENDIF
            ELSEIF ( cname .AND. c2=='UN' ) THEN
               IF ( c3(1:1)=='G' ) THEN
                  IF ( c4=='QR' .OR. c4=='RQ' .OR. c4=='LQ' .OR.        &
     &                 c4=='QL' .OR. c4=='HR' .OR. c4=='TR' .OR.        &
     &                 c4=='BR' ) nbmin = 2
               ELSEIF ( c3(1:1)=='M' ) THEN
                  IF ( c4=='QR' .OR. c4=='RQ' .OR. c4=='LQ' .OR.        &
     &                 c4=='QL' .OR. c4=='HR' .OR. c4=='TR' .OR.        &
     &                 c4=='BR' ) nbmin = 2
               ENDIF
            ELSEIF ( c2=='GG' ) THEN
               nbmin = 2
               IF ( c3=='HD3' ) nbmin = 2
            ENDIF
            ILAENV = nbmin
            RETURN
         ELSEIF ( Ispec==3 ) THEN
!
!
!     ISPEC = 3:  crossover point
!
            nx = 0
            IF ( c2=='GE' ) THEN
               IF ( c3=='QRF' .OR. c3=='RQF' .OR. c3=='LQF' .OR.        &
     &              c3=='QLF' ) THEN
                  IF ( sname ) THEN
                     nx = 128
                  ELSE
                     nx = 128
                  ENDIF
               ELSEIF ( c3=='HRD' ) THEN
                  IF ( sname ) THEN
                     nx = 128
                  ELSE
                     nx = 128
                  ENDIF
               ELSEIF ( c3=='BRD' ) THEN
                  IF ( sname ) THEN
                     nx = 128
                  ELSE
                     nx = 128
                  ENDIF
               ENDIF
            ELSEIF ( c2=='SY' ) THEN
               IF ( sname .AND. c3=='TRD' ) nx = 32
            ELSEIF ( cname .AND. c2=='HE' ) THEN
               IF ( c3=='TRD' ) nx = 32
            ELSEIF ( sname .AND. c2=='OR' ) THEN
               IF ( c3(1:1)=='G' ) THEN
                  IF ( c4=='QR' .OR. c4=='RQ' .OR. c4=='LQ' .OR.        &
     &                 c4=='QL' .OR. c4=='HR' .OR. c4=='TR' .OR.        &
     &                 c4=='BR' ) nx = 128
               ENDIF
            ELSEIF ( cname .AND. c2=='UN' ) THEN
               IF ( c3(1:1)=='G' ) THEN
                  IF ( c4=='QR' .OR. c4=='RQ' .OR. c4=='LQ' .OR.        &
     &                 c4=='QL' .OR. c4=='HR' .OR. c4=='TR' .OR.        &
     &                 c4=='BR' ) nx = 128
               ENDIF
            ELSEIF ( c2=='GG' ) THEN
               nx = 128
               IF ( c3=='HD3' ) nx = 128
            ENDIF
            ILAENV = nx
            RETURN
         ELSE
!
!
!     ISPEC = 1:  block size
!
!     In these examples, separate code is provided for setting NB for
!     real and complex.  We assume that NB will take the same value in
!     single or double precision.
!
            nb = 1
!
            IF ( subnam(2:6)=='LAORH' ) THEN
!
!        This is for *LAORHR_GETRFNP routine
!
               IF ( sname ) THEN
                  nb = 32
               ELSE
                  nb = 32
               ENDIF
            ELSEIF ( c2=='GE' ) THEN
               IF ( c3=='TRF' ) THEN
                  IF ( sname ) THEN
                     nb = 64
                  ELSE
                     nb = 64
                  ENDIF
               ELSEIF ( c3=='QRF' .OR. c3=='RQF' .OR. c3=='LQF' .OR.    &
     &                  c3=='QLF' ) THEN
                  IF ( sname ) THEN
                     nb = 32
                  ELSE
                     nb = 32
                  ENDIF
               ELSEIF ( c3=='QR ' ) THEN
                  IF ( N3==1 ) THEN
                     IF ( sname ) THEN
!     M*N
                        IF ( (N1*N2<=131072) .OR. (N1<=8192) ) THEN
                           nb = N1
                        ELSE
                           nb = 32768/N2
                        ENDIF
                     ELSEIF ( (N1*N2<=131072) .OR. (N1<=8192) ) THEN
                        nb = N1
                     ELSE
                        nb = 32768/N2
                     ENDIF
                  ELSEIF ( sname ) THEN
                     nb = 1
                  ELSE
                     nb = 1
                  ENDIF
               ELSEIF ( c3=='LQ ' ) THEN
                  IF ( N3==2 ) THEN
                     IF ( sname ) THEN
!     M*N
                        IF ( (N1*N2<=131072) .OR. (N1<=8192) ) THEN
                           nb = N1
                        ELSE
                           nb = 32768/N2
                        ENDIF
                     ELSEIF ( (N1*N2<=131072) .OR. (N1<=8192) ) THEN
                        nb = N1
                     ELSE
                        nb = 32768/N2
                     ENDIF
                  ELSEIF ( sname ) THEN
                     nb = 1
                  ELSE
                     nb = 1
                  ENDIF
               ELSEIF ( c3=='HRD' ) THEN
                  IF ( sname ) THEN
                     nb = 32
                  ELSE
                     nb = 32
                  ENDIF
               ELSEIF ( c3=='BRD' ) THEN
                  IF ( sname ) THEN
                     nb = 32
                  ELSE
                     nb = 32
                  ENDIF
               ELSEIF ( c3=='TRI' ) THEN
                  IF ( sname ) THEN
                     nb = 64
                  ELSE
                     nb = 64
                  ENDIF
               ENDIF
            ELSEIF ( c2=='PO' ) THEN
               IF ( c3=='TRF' ) THEN
                  IF ( sname ) THEN
                     nb = 64
                  ELSE
                     nb = 64
                  ENDIF
               ENDIF
            ELSEIF ( c2=='SY' ) THEN
               IF ( c3=='TRF' ) THEN
                  IF ( sname ) THEN
                     IF ( twostage ) THEN
                        nb = 192
                     ELSE
                        nb = 64
                     ENDIF
                  ELSEIF ( twostage ) THEN
                     nb = 192
                  ELSE
                     nb = 64
                  ENDIF
               ELSEIF ( sname .AND. c3=='TRD' ) THEN
                  nb = 32
               ELSEIF ( sname .AND. c3=='GST' ) THEN
                  nb = 64
               ENDIF
            ELSEIF ( cname .AND. c2=='HE' ) THEN
               IF ( c3=='TRF' ) THEN
                  IF ( twostage ) THEN
                     nb = 192
                  ELSE
                     nb = 64
                  ENDIF
               ELSEIF ( c3=='TRD' ) THEN
                  nb = 32
               ELSEIF ( c3=='GST' ) THEN
                  nb = 64
               ENDIF
            ELSEIF ( sname .AND. c2=='OR' ) THEN
               IF ( c3(1:1)=='G' ) THEN
                  IF ( c4=='QR' .OR. c4=='RQ' .OR. c4=='LQ' .OR.        &
     &                 c4=='QL' .OR. c4=='HR' .OR. c4=='TR' .OR.        &
     &                 c4=='BR' ) nb = 32
               ELSEIF ( c3(1:1)=='M' ) THEN
                  IF ( c4=='QR' .OR. c4=='RQ' .OR. c4=='LQ' .OR.        &
     &                 c4=='QL' .OR. c4=='HR' .OR. c4=='TR' .OR.        &
     &                 c4=='BR' ) nb = 32
               ENDIF
            ELSEIF ( cname .AND. c2=='UN' ) THEN
               IF ( c3(1:1)=='G' ) THEN
                  IF ( c4=='QR' .OR. c4=='RQ' .OR. c4=='LQ' .OR.        &
     &                 c4=='QL' .OR. c4=='HR' .OR. c4=='TR' .OR.        &
     &                 c4=='BR' ) nb = 32
               ELSEIF ( c3(1:1)=='M' ) THEN
                  IF ( c4=='QR' .OR. c4=='RQ' .OR. c4=='LQ' .OR.        &
     &                 c4=='QL' .OR. c4=='HR' .OR. c4=='TR' .OR.        &
     &                 c4=='BR' ) nb = 32
               ENDIF
            ELSEIF ( c2=='GB' ) THEN
               IF ( c3=='TRF' ) THEN
                  IF ( sname ) THEN
                     IF ( N4<=64 ) THEN
                        nb = 1
                     ELSE
                        nb = 32
                     ENDIF
                  ELSEIF ( N4<=64 ) THEN
                     nb = 1
                  ELSE
                     nb = 32
                  ENDIF
               ENDIF
            ELSEIF ( c2=='PB' ) THEN
               IF ( c3=='TRF' ) THEN
                  IF ( sname ) THEN
                     IF ( N2<=64 ) THEN
                        nb = 1
                     ELSE
                        nb = 32
                     ENDIF
                  ELSEIF ( N2<=64 ) THEN
                     nb = 1
                  ELSE
                     nb = 32
                  ENDIF
               ENDIF
            ELSEIF ( c2=='TR' ) THEN
               IF ( c3=='TRI' ) THEN
                  IF ( sname ) THEN
                     nb = 64
                  ELSE
                     nb = 64
                  ENDIF
               ELSEIF ( c3=='EVC' ) THEN
                  IF ( sname ) THEN
                     nb = 64
                  ELSE
                     nb = 64
                  ENDIF
               ENDIF
            ELSEIF ( c2=='LA' ) THEN
               IF ( c3=='UUM' ) THEN
                  IF ( sname ) THEN
                     nb = 64
                  ELSE
                     nb = 64
                  ENDIF
               ENDIF
            ELSEIF ( sname .AND. c2=='ST' ) THEN
               IF ( c3=='EBZ' ) nb = 1
            ELSEIF ( c2=='GG' ) THEN
               nb = 32
               IF ( c3=='HD3' ) THEN
                  IF ( sname ) THEN
                     nb = 32
                  ELSE
                     nb = 32
                  ENDIF
               ENDIF
            ENDIF
            ILAENV = nb
            RETURN
         ENDIF
      ELSEIF ( Ispec==4 ) THEN
!
!
!     ISPEC = 4:  number of shifts (used by xHSEQR)
!
         ILAENV = 6
         RETURN
      ELSEIF ( Ispec==5 ) THEN
!
!
!     ISPEC = 5:  minimum column dimension (not used)
!
         ILAENV = 2
         RETURN
      ELSEIF ( Ispec==6 ) THEN
!
!
!     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD)
!
         ILAENV = INT(REAL(MIN(N1,N2))*1.6E0)
         RETURN
      ELSEIF ( Ispec==7 ) THEN
!
!
!     ISPEC = 7:  number of processors (not used)
!
         ILAENV = 1
         RETURN
      ELSEIF ( Ispec==8 ) THEN
!
!
!     ISPEC = 8:  crossover point for multishift (used by xHSEQR)
!
         ILAENV = 50
         RETURN
      ELSEIF ( Ispec==9 ) THEN
!
!
!     ISPEC = 9:  maximum size of the subproblems at the bottom of the
!                 computation tree in the divide-and-conquer algorithm
!                 (used by xGELSD and xGESDD)
!
         ILAENV = 25
         RETURN
      ELSEIF ( Ispec==10 ) THEN
!
!
!     ISPEC = 10: ieee and infinity NaN arithmetic can be trusted not to trap
!
!     ILAENV = 0
         ILAENV = 1
         IF ( ILAENV==1 ) ILAENV = IEEECK(1,0.0,1.0)
         RETURN
      ELSEIF ( Ispec==11 ) THEN
!
!
!     ISPEC = 11: ieee infinity arithmetic can be trusted not to trap
!
!     ILAENV = 0
         ILAENV = 1
         IF ( ILAENV==1 ) ILAENV = IEEECK(0,0.0,1.0)
         RETURN
      ELSEIF ( Ispec==12 .OR. Ispec==13 .OR. Ispec==14 .OR.             &
     &         Ispec==15 .OR. Ispec==16 ) THEN
!
!
!     12 <= ISPEC <= 16: xHSEQR or related subroutines.
!
         ILAENV = IPARMQ(Ispec,Name,Opts,N1,N2,N3,N4)
         GOTO 99999
      ENDIF
!
!     Invalid value for ISPEC
!
      ILAENV = -1
      RETURN
!
!     End of ILAENV
!
99999 END FUNCTION ILAENV
