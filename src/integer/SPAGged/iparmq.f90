!*==iparmq.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b IPARMQ
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download IPARMQ + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/iparmq.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/iparmq.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/iparmq.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       INTEGER FUNCTION IPARMQ( ISPEC, NAME, OPTS, N, ILO, IHI, LWORK )
!
!       .. Scalar Arguments ..
!       INTEGER            IHI, ILO, ISPEC, LWORK, N
!       CHARACTER          NAME*( * ), OPTS*( * )
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>      This program sets problem and machine dependent parameters
!>      useful for xHSEQR and related subroutines for eigenvalue
!>      problems. It is called whenever
!>      IPARMQ is called with 12 <= ISPEC <= 16
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] ISPEC
!> \verbatim
!>          ISPEC is INTEGER
!>              ISPEC specifies which tunable parameter IPARMQ should
!>              return.
!>
!>              ISPEC=12: (INMIN)  Matrices of order nmin or less
!>                        are sent directly to xLAHQR, the implicit
!>                        double shift QR algorithm.  NMIN must be
!>                        at least 11.
!>
!>              ISPEC=13: (INWIN)  Size of the deflation window.
!>                        This is best set greater than or equal to
!>                        the number of simultaneous shifts NS.
!>                        Larger matrices benefit from larger deflation
!>                        windows.
!>
!>              ISPEC=14: (INIBL) Determines when to stop nibbling and
!>                        invest in an (expensive) multi-shift QR sweep.
!>                        If the aggressive early deflation subroutine
!>                        finds LD converged eigenvalues from an order
!>                        NW deflation window and LD > (NW*NIBBLE)/100,
!>                        then the next QR sweep is skipped and early
!>                        deflation is applied immediately to the
!>                        remaining active diagonal block.  Setting
!>                        IPARMQ(ISPEC=14) = 0 causes TTQRE to skip a
!>                        multi-shift QR sweep whenever early deflation
!>                        finds a converged eigenvalue.  Setting
!>                        IPARMQ(ISPEC=14) greater than or equal to 100
!>                        prevents TTQRE from skipping a multi-shift
!>                        QR sweep.
!>
!>              ISPEC=15: (NSHFTS) The number of simultaneous shifts in
!>                        a multi-shift QR iteration.
!>
!>              ISPEC=16: (IACC22) IPARMQ is set to 0, 1 or 2 with the
!>                        following meanings.
!>                        0:  During the multi-shift QR/QZ sweep,
!>                            blocked eigenvalue reordering, blocked
!>                            Hessenberg-triangular reduction,
!>                            reflections and/or rotations are not
!>                            accumulated when updating the
!>                            far-from-diagonal matrix entries.
!>                        1:  During the multi-shift QR/QZ sweep,
!>                            blocked eigenvalue reordering, blocked
!>                            Hessenberg-triangular reduction,
!>                            reflections and/or rotations are
!>                            accumulated, and matrix-matrix
!>                            multiplication is used to update the
!>                            far-from-diagonal matrix entries.
!>                        2:  During the multi-shift QR/QZ sweep,
!>                            blocked eigenvalue reordering, blocked
!>                            Hessenberg-triangular reduction,
!>                            reflections and/or rotations are
!>                            accumulated, and 2-by-2 block structure
!>                            is exploited during matrix-matrix
!>                            multiplies.
!>                        (If xTRMM is slower than xGEMM, then
!>                        IPARMQ(ISPEC=16)=1 may be more efficient than
!>                        IPARMQ(ISPEC=16)=2 despite the greater level of
!>                        arithmetic work implied by the latter choice.)
!> \endverbatim
!>
!> \param[in] NAME
!> \verbatim
!>          NAME is CHARACTER string
!>               Name of the calling subroutine
!> \endverbatim
!>
!> \param[in] OPTS
!> \verbatim
!>          OPTS is CHARACTER string
!>               This is a concatenation of the string arguments to
!>               TTQRE.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>               N is the order of the Hessenberg matrix H.
!> \endverbatim
!>
!> \param[in] ILO
!> \verbatim
!>          ILO is INTEGER
!> \endverbatim
!>
!> \param[in] IHI
!> \verbatim
!>          IHI is INTEGER
!>               It is assumed that H is already upper triangular
!>               in rows and columns 1:ILO-1 and IHI+1:N.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>               The amount of workspace available.
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
!> \date June 2017
!
!> \ingroup OTHERauxiliary
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>       Little is known about how best to choose these parameters.
!>       It is possible to use different values of the parameters
!>       for each of CHSEQR, DHSEQR, SHSEQR and ZHSEQR.
!>
!>       It is probably best to choose different parameters for
!>       different matrices and different parameters at different
!>       times during the iteration, but this has not been
!>       implemented --- yet.
!>
!>
!>       The best choices of most of the parameters depend
!>       in an ill-understood way on the relative execution
!>       rate of xLAQR3 and xLAQR5 and on the nature of each
!>       particular eigenvalue problem.  Experiment may be the
!>       only practical way to determine which choices are most
!>       effective.
!>
!>       Following is a list of default values supplied by IPARMQ.
!>       These defaults may be adjusted in order to attain better
!>       performance in any particular computational environment.
!>
!>       IPARMQ(ISPEC=12) The xLAHQR vs xLAQR0 crossover point.
!>                        Default: 75. (Must be at least 11.)
!>
!>       IPARMQ(ISPEC=13) Recommended deflation window size.
!>                        This depends on ILO, IHI and NS, the
!>                        number of simultaneous shifts returned
!>                        by IPARMQ(ISPEC=15).  The default for
!>                        (IHI-ILO+1) <= 500 is NS.  The default
!>                        for (IHI-ILO+1) > 500 is 3*NS/2.
!>
!>       IPARMQ(ISPEC=14) Nibble crossover point.  Default: 14.
!>
!>       IPARMQ(ISPEC=15) Number of simultaneous shifts, NS.
!>                        a multi-shift QR iteration.
!>
!>                        If IHI-ILO+1 is ...
!>
!>                        greater than      ...but less    ... the
!>                        or equal to ...      than        default is
!>
!>                                0               30       NS =   2+
!>                               30               60       NS =   4+
!>                               60              150       NS =  10
!>                              150              590       NS =  **
!>                              590             3000       NS =  64
!>                             3000             6000       NS = 128
!>                             6000             infinity   NS = 256
!>
!>                    (+)  By default matrices of this order are
!>                         passed to the implicit double shift routine
!>                         xLAHQR.  See IPARMQ(ISPEC=12) above.   These
!>                         values of NS are used only in case of a rare
!>                         xLAHQR failure.
!>
!>                    (**) The asterisks (**) indicate an ad-hoc
!>                         function increasing from 10 to 64.
!>
!>       IPARMQ(ISPEC=16) Select structured matrix multiply.
!>                        (See ISPEC=16 above for details.)
!>                        Default: 3.
!> \endverbatim
!>
!  =====================================================================
      FUNCTION IPARMQ(Ispec,Name,Opts,N,Ilo,Ihi,Lwork)
      IMPLICIT NONE
!*--IPARMQ226
!
! PARAMETER definitions rewritten by SPAG
!
      INTEGER , PARAMETER  ::  INMIN = 12 , INWIN = 13 , INIBL = 14 ,   &
     &                         ISHFTS = 15 , IACC22 = 16 , NMIN = 75 ,  &
     &                         K22MIN = 14 , KACMIN = 14 , NIBBLE = 14 ,&
     &                         KNWSWP = 500
      REAL , PARAMETER  ::  TWO = 2.0
      INTEGER :: IPARMQ
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: Ispec
      CHARACTER(*) , INTENT(IN) :: Name
      CHARACTER(*) :: Opts
      INTEGER :: N
      INTEGER , INTENT(IN) :: Ilo
      INTEGER , INTENT(IN) :: Ihi
      INTEGER :: Lwork
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: i , ic , iz , nh , ns
      CHARACTER(6) :: subnam
!
! End of declarations rewritten by SPAG
!
!
!  ================================================================
!     .. Parameters ..
!     ..
!     .. Local Scalars ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
      IF ( (Ispec==ISHFTS) .OR. (Ispec==INWIN) .OR. (Ispec==IACC22) )   &
     &     THEN
!
!        ==== Set the number simultaneous shifts ====
!
         nh = Ihi - Ilo + 1
         ns = 2
         IF ( nh>=30 ) ns = 4
         IF ( nh>=60 ) ns = 10
         IF ( nh>=150 ) ns = MAX(10,nh/NINT(LOG(REAL(nh))/LOG(TWO)))
         IF ( nh>=590 ) ns = 64
         IF ( nh>=3000 ) ns = 128
         IF ( nh>=6000 ) ns = 256
         ns = MAX(2,ns-MOD(ns,2))
      ENDIF
!
      IF ( Ispec==INMIN ) THEN
!
!
!        ===== Matrices of order smaller than NMIN get sent
!        .     to xLAHQR, the classic double shift algorithm.
!        .     This must be at least 11. ====
!
         IPARMQ = NMIN
!
      ELSEIF ( Ispec==INIBL ) THEN
!
!        ==== INIBL: skip a multi-shift qr iteration and
!        .    whenever aggressive early deflation finds
!        .    at least (NIBBLE*(window size)/100) deflations. ====
!
         IPARMQ = NIBBLE
!
      ELSEIF ( Ispec==ISHFTS ) THEN
!
!        ==== NSHFTS: The number of simultaneous shifts =====
!
         IPARMQ = ns
!
      ELSEIF ( Ispec==INWIN ) THEN
!
!        ==== NW: deflation window size.  ====
!
         IF ( nh<=KNWSWP ) THEN
            IPARMQ = ns
         ELSE
            IPARMQ = 3*ns/2
         ENDIF
!
      ELSEIF ( Ispec==IACC22 ) THEN
!
!        ==== IACC22: Whether to accumulate reflections
!        .     before updating the far-from-diagonal elements
!        .     and whether to use 2-by-2 block structure while
!        .     doing it.  A small amount of work could be saved
!        .     by making this choice dependent also upon the
!        .     NH=IHI-ILO+1.
!
!
!        Convert NAME to upper case if the first character is lower case.
!
         IPARMQ = 0
         subnam = Name
         ic = ICHAR(subnam(1:1))
         iz = ICHAR('Z')
         IF ( iz==90 .OR. iz==122 ) THEN
!
!           ASCII character set
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
!           EBCDIC character set
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
!           Prime machines:  ASCII+128
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
         IF ( subnam(2:6)=='GGHRD' .OR. subnam(2:6)=='GGHD3' ) THEN
            IPARMQ = 1
            IF ( nh>=K22MIN ) IPARMQ = 2
         ELSEIF ( subnam(4:6)=='EXC' ) THEN
            IF ( nh>=KACMIN ) IPARMQ = 1
            IF ( nh>=K22MIN ) IPARMQ = 2
         ELSEIF ( subnam(2:6)=='HSEQR' .OR. subnam(2:5)=='LAQR' ) THEN
            IF ( ns>=KACMIN ) IPARMQ = 1
            IF ( ns>=K22MIN ) IPARMQ = 2
         ENDIF
!
      ELSE
!        ===== invalid value of ispec =====
         IPARMQ = -1
!
      ENDIF
!
!     ==== End of IPARMQ ====
!
      END FUNCTION IPARMQ
