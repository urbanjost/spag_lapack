!*==chseqr.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b CHSEQR
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CHSEQR + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chseqr.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chseqr.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chseqr.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CHSEQR( JOB, COMPZ, N, ILO, IHI, H, LDH, W, Z, LDZ,
!                          WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            IHI, ILO, INFO, LDH, LDZ, LWORK, N
!       CHARACTER          COMPZ, JOB
!       ..
!       .. Array Arguments ..
!       COMPLEX            H( LDH, * ), W( * ), WORK( * ), Z( LDZ, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    CHSEQR computes the eigenvalues of a Hessenberg matrix H
!>    and, optionally, the matrices T and Z from the Schur decomposition
!>    H = Z T Z**H, where T is an upper triangular matrix (the
!>    Schur form), and Z is the unitary matrix of Schur vectors.
!>
!>    Optionally Z may be postmultiplied into an input unitary
!>    matrix Q so that this routine can give the Schur factorization
!>    of a matrix A which has been reduced to the Hessenberg form H
!>    by the unitary matrix Q:  A = Q*H*Q**H = (QZ)*T*(QZ)**H.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] JOB
!> \verbatim
!>          JOB is CHARACTER*1
!>           = 'E':  compute eigenvalues only;
!>           = 'S':  compute eigenvalues and the Schur form T.
!> \endverbatim
!>
!> \param[in] COMPZ
!> \verbatim
!>          COMPZ is CHARACTER*1
!>           = 'N':  no Schur vectors are computed;
!>           = 'I':  Z is initialized to the unit matrix and the matrix Z
!>                   of Schur vectors of H is returned;
!>           = 'V':  Z must contain an unitary matrix Q on entry, and
!>                   the product Q*Z is returned.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>           The order of the matrix H.  N >= 0.
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
!>
!>           It is assumed that H is already upper triangular in rows
!>           and columns 1:ILO-1 and IHI+1:N. ILO and IHI are normally
!>           set by a previous call to CGEBAL, and then passed to ZGEHRD
!>           when the matrix output by CGEBAL is reduced to Hessenberg
!>           form. Otherwise ILO and IHI should be set to 1 and N
!>           respectively.  If N > 0, then 1 <= ILO <= IHI <= N.
!>           If N = 0, then ILO = 1 and IHI = 0.
!> \endverbatim
!>
!> \param[in,out] H
!> \verbatim
!>          H is COMPLEX array, dimension (LDH,N)
!>           On entry, the upper Hessenberg matrix H.
!>           On exit, if INFO = 0 and JOB = 'S', H contains the upper
!>           triangular matrix T from the Schur decomposition (the
!>           Schur form). If INFO = 0 and JOB = 'E', the contents of
!>           H are unspecified on exit.  (The output value of H when
!>           INFO > 0 is given under the description of INFO below.)
!>
!>           Unlike earlier versions of CHSEQR, this subroutine may
!>           explicitly H(i,j) = 0 for i > j and j = 1, 2, ... ILO-1
!>           or j = IHI+1, IHI+2, ... N.
!> \endverbatim
!>
!> \param[in] LDH
!> \verbatim
!>          LDH is INTEGER
!>           The leading dimension of the array H. LDH >= max(1,N).
!> \endverbatim
!>
!> \param[out] W
!> \verbatim
!>          W is COMPLEX array, dimension (N)
!>           The computed eigenvalues. If JOB = 'S', the eigenvalues are
!>           stored in the same order as on the diagonal of the Schur
!>           form returned in H, with W(i) = H(i,i).
!> \endverbatim
!>
!> \param[in,out] Z
!> \verbatim
!>          Z is COMPLEX array, dimension (LDZ,N)
!>           If COMPZ = 'N', Z is not referenced.
!>           If COMPZ = 'I', on entry Z need not be set and on exit,
!>           if INFO = 0, Z contains the unitary matrix Z of the Schur
!>           vectors of H.  If COMPZ = 'V', on entry Z must contain an
!>           N-by-N matrix Q, which is assumed to be equal to the unit
!>           matrix except for the submatrix Z(ILO:IHI,ILO:IHI). On exit,
!>           if INFO = 0, Z contains Q*Z.
!>           Normally Q is the unitary matrix generated by CUNGHR
!>           after the call to CGEHRD which formed the Hessenberg matrix
!>           H. (The output value of Z when INFO > 0 is given under
!>           the description of INFO below.)
!> \endverbatim
!>
!> \param[in] LDZ
!> \verbatim
!>          LDZ is INTEGER
!>           The leading dimension of the array Z.  if COMPZ = 'I' or
!>           COMPZ = 'V', then LDZ >= MAX(1,N).  Otherwise, LDZ >= 1.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension (LWORK)
!>           On exit, if INFO = 0, WORK(1) returns an estimate of
!>           the optimal value for LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>           The dimension of the array WORK.  LWORK >= max(1,N)
!>           is sufficient and delivers very good and sometimes
!>           optimal performance.  However, LWORK as large as 11*N
!>           may be required for optimal performance.  A workspace
!>           query is recommended to determine the optimal workspace
!>           size.
!>
!>           If LWORK = -1, then CHSEQR does a workspace query.
!>           In this case, CHSEQR checks the input parameters and
!>           estimates the optimal workspace size for the given
!>           values of N, ILO and IHI.  The estimate is returned
!>           in WORK(1).  No error message related to LWORK is
!>           issued by XERBLA.  Neither H nor Z are accessed.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>             = 0:  successful exit
!>             < 0:  if INFO = -i, the i-th argument had an illegal
!>                    value
!>             > 0:  if INFO = i, CHSEQR failed to compute all of
!>                the eigenvalues.  Elements 1:ilo-1 and i+1:n of W
!>                contain those eigenvalues which have been
!>                successfully computed.  (Failures are rare.)
!>
!>                If INFO > 0 and JOB = 'E', then on exit, the
!>                remaining unconverged eigenvalues are the eigen-
!>                values of the upper Hessenberg matrix rows and
!>                columns ILO through INFO of the final, output
!>                value of H.
!>
!>                If INFO > 0 and JOB   = 'S', then on exit
!>
!>           (*)  (initial value of H)*U  = U*(final value of H)
!>
!>                where U is a unitary matrix.  The final
!>                value of  H is upper Hessenberg and triangular in
!>                rows and columns INFO+1 through IHI.
!>
!>                If INFO > 0 and COMPZ = 'V', then on exit
!>
!>                  (final value of Z)  =  (initial value of Z)*U
!>
!>                where U is the unitary matrix in (*) (regard-
!>                less of the value of JOB.)
!>
!>                If INFO > 0 and COMPZ = 'I', then on exit
!>                      (final value of Z)  = U
!>                where U is the unitary matrix in (*) (regard-
!>                less of the value of JOB.)
!>
!>                If INFO > 0 and COMPZ = 'N', then Z is not
!>                accessed.
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
!> \ingroup complexOTHERcomputational
!
!> \par Contributors:
!  ==================
!>
!>       Karen Braman and Ralph Byers, Department of Mathematics,
!>       University of Kansas, USA
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>             Default values supplied by
!>             ILAENV(ISPEC,'CHSEQR',JOB(:1)//COMPZ(:1),N,ILO,IHI,LWORK).
!>             It is suggested that these defaults be adjusted in order
!>             to attain best performance in each particular
!>             computational environment.
!>
!>            ISPEC=12: The CLAHQR vs CLAQR0 crossover point.
!>                      Default: 75. (Must be at least 11.)
!>
!>            ISPEC=13: Recommended deflation window size.
!>                      This depends on ILO, IHI and NS.  NS is the
!>                      number of simultaneous shifts returned
!>                      by ILAENV(ISPEC=15).  (See ISPEC=15 below.)
!>                      The default for (IHI-ILO+1) <= 500 is NS.
!>                      The default for (IHI-ILO+1) >  500 is 3*NS/2.
!>
!>            ISPEC=14: Nibble crossover point. (See IPARMQ for
!>                      details.)  Default: 14% of deflation window
!>                      size.
!>
!>            ISPEC=15: Number of simultaneous shifts in a multishift
!>                      QR iteration.
!>
!>                      If IHI-ILO+1 is ...
!>
!>                      greater than      ...but less    ... the
!>                      or equal to ...      than        default is
!>
!>                           1               30          NS =   2(+)
!>                          30               60          NS =   4(+)
!>                          60              150          NS =  10(+)
!>                         150              590          NS =  **
!>                         590             3000          NS =  64
!>                        3000             6000          NS = 128
!>                        6000             infinity      NS = 256
!>
!>                  (+)  By default some or all matrices of this order
!>                       are passed to the implicit double shift routine
!>                       CLAHQR and this parameter is ignored.  See
!>                       ISPEC=12 above and comments in IPARMQ for
!>                       details.
!>
!>                 (**)  The asterisks (**) indicate an ad-hoc
!>                       function of N increasing from 10 to 64.
!>
!>            ISPEC=16: Select structured matrix multiply.
!>                      If the number of simultaneous shifts (specified
!>                      by ISPEC=15) is less than 14, then the default
!>                      for ISPEC=16 is 0.  Otherwise the default for
!>                      ISPEC=16 is 2.
!> \endverbatim
!
!> \par References:
!  ================
!>
!>       K. Braman, R. Byers and R. Mathias, The Multi-Shift QR
!>       Algorithm Part I: Maintaining Well Focused Shifts, and Level 3
!>       Performance, SIAM Journal of Matrix Analysis, volume 23, pages
!>       929--947, 2002.
!> \n
!>       K. Braman, R. Byers and R. Mathias, The Multi-Shift QR
!>       Algorithm Part II: Aggressive Early Deflation, SIAM Journal
!>       of Matrix Analysis, volume 23, pages 948--973, 2002.
!
!  =====================================================================
      SUBROUTINE CHSEQR(Job,Compz,N,Ilo,Ihi,H,Ldh,W,Z,Ldz,Work,Lwork,   &
     &                  Info)
      USE S_CCOPY
      USE S_CLACPY
      USE S_CLAHQR
      USE S_CLAQR0
      USE S_CLASET
      USE S_ILAENV
      USE S_LSAME
      USE S_XERBLA
      IMPLICIT NONE
!*--CHSEQR311
!
! PARAMETER definitions rewritten by SPAG
!
      INTEGER , PARAMETER  ::  NTINY = 15 , NL = 49
      COMPLEX , PARAMETER  ::  ZERO = (0.0E0,0.0E0) ,                   &
     &                         ONE = (1.0E0,0.0E0)
      REAL , PARAMETER  ::  RZERO = 0.0E0
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Job
      CHARACTER :: Compz
      INTEGER :: N
      INTEGER :: Ilo
      INTEGER :: Ihi
      COMPLEX , DIMENSION(Ldh,*) :: H
      INTEGER :: Ldh
      COMPLEX , DIMENSION(*) :: W
      COMPLEX , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      COMPLEX , DIMENSION(NL,NL) :: hl
      LOGICAL :: initz , lquery , wantt , wantz
      INTEGER :: kbot , nmin
      COMPLEX , DIMENSION(NL) :: workl
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
!
!     ==== Matrices of order NTINY or smaller must be processed by
!     .    CLAHQR because of insufficient subdiagonal scratch space.
!     .    (This is a hard limit.) ====
!
!     ==== NL allocates some local workspace to help small matrices
!     .    through a rare CLAHQR failure.  NL > NTINY = 15 is
!     .    required and NL <= NMIN = ILAENV(ISPEC=12,...) is recom-
!     .    mended.  (The default value of NMIN is 75.)  Using NL = 49
!     .    allows up to six simultaneous shifts and a 16-by-16
!     .    deflation window.  ====
!     ..
!     .. Local Arrays ..
!     ..
!     .. Local Scalars ..
!     ..
!     .. External Functions ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
!     ==== Decode and check the input parameters. ====
!
      wantt = LSAME(Job,'S')
      initz = LSAME(Compz,'I')
      wantz = initz .OR. LSAME(Compz,'V')
      Work(1) = CMPLX(REAL(MAX(1,N)),RZERO)
      lquery = Lwork== - 1
!
      Info = 0
      IF ( .NOT.LSAME(Job,'E') .AND. .NOT.wantt ) THEN
         Info = -1
      ELSEIF ( .NOT.LSAME(Compz,'N') .AND. .NOT.wantz ) THEN
         Info = -2
      ELSEIF ( N<0 ) THEN
         Info = -3
      ELSEIF ( Ilo<1 .OR. Ilo>MAX(1,N) ) THEN
         Info = -4
      ELSEIF ( Ihi<MIN(Ilo,N) .OR. Ihi>N ) THEN
         Info = -5
      ELSEIF ( Ldh<MAX(1,N) ) THEN
         Info = -7
      ELSEIF ( Ldz<1 .OR. (wantz .AND. Ldz<MAX(1,N)) ) THEN
         Info = -10
      ELSEIF ( Lwork<MAX(1,N) .AND. .NOT.lquery ) THEN
         Info = -12
      ENDIF
!
      IF ( Info/=0 ) THEN
!
!        ==== Quick return in case of invalid argument. ====
!
         CALL XERBLA('CHSEQR',-Info)
         RETURN
!
      ELSEIF ( N==0 ) THEN
!
!        ==== Quick return in case N = 0; nothing to do. ====
!
         RETURN
!
      ELSEIF ( lquery ) THEN
!
!        ==== Quick return in case of a workspace query ====
!
         CALL CLAQR0(wantt,wantz,N,Ilo,Ihi,H,Ldh,W,Ilo,Ihi,Z,Ldz,Work,  &
     &               Lwork,Info)
!        ==== Ensure reported workspace size is backward-compatible with
!        .    previous LAPACK versions. ====
         Work(1) = CMPLX(MAX(REAL(Work(1)),REAL(MAX(1,N))),RZERO)
         RETURN
!
      ELSE
!
!        ==== copy eigenvalues isolated by CGEBAL ====
!
         IF ( Ilo>1 ) CALL CCOPY(Ilo-1,H,Ldh+1,W,1)
         IF ( Ihi<N ) CALL CCOPY(N-Ihi,H(Ihi+1,Ihi+1),Ldh+1,W(Ihi+1),1)
!
!        ==== Initialize Z, if requested ====
!
         IF ( initz ) CALL CLASET('A',N,N,ZERO,ONE,Z,Ldz)
!
!        ==== Quick return if possible ====
!
         IF ( Ilo==Ihi ) THEN
            W(Ilo) = H(Ilo,Ilo)
            RETURN
         ENDIF
!
!        ==== CLAHQR/CLAQR0 crossover point ====
!
         nmin = ILAENV(12,'CHSEQR',Job(:1)//Compz(:1),N,Ilo,Ihi,Lwork)
         nmin = MAX(NTINY,nmin)
!
!        ==== CLAQR0 for big matrices; CLAHQR for small ones ====
!
         IF ( N>nmin ) THEN
            CALL CLAQR0(wantt,wantz,N,Ilo,Ihi,H,Ldh,W,Ilo,Ihi,Z,Ldz,    &
     &                  Work,Lwork,Info)
         ELSE
!
!           ==== Small matrix ====
!
            CALL CLAHQR(wantt,wantz,N,Ilo,Ihi,H,Ldh,W,Ilo,Ihi,Z,Ldz,    &
     &                  Info)
!
            IF ( Info>0 ) THEN
!
!              ==== A rare CLAHQR failure!  CLAQR0 sometimes succeeds
!              .    when CLAHQR fails. ====
!
               kbot = Info
!
               IF ( N>=NL ) THEN
!
!                 ==== Larger matrices have enough subdiagonal scratch
!                 .    space to call CLAQR0 directly. ====
!
                  CALL CLAQR0(wantt,wantz,N,Ilo,kbot,H,Ldh,W,Ilo,Ihi,Z, &
     &                        Ldz,Work,Lwork,Info)
!
               ELSE
!
!                 ==== Tiny matrices don't have enough subdiagonal
!                 .    scratch space to benefit from CLAQR0.  Hence,
!                 .    tiny matrices must be copied into a larger
!                 .    array before calling CLAQR0. ====
!
                  CALL CLACPY('A',N,N,H,Ldh,hl,NL)
                  hl(N+1,N) = ZERO
                  CALL CLASET('A',NL,NL-N,ZERO,ZERO,hl(1,N+1),NL)
                  CALL CLAQR0(wantt,wantz,NL,Ilo,kbot,hl,NL,W,Ilo,Ihi,Z,&
     &                        Ldz,workl,NL,Info)
                  IF ( wantt .OR. Info/=0 )                             &
     &                 CALL CLACPY('A',N,N,hl,NL,H,Ldh)
               ENDIF
            ENDIF
         ENDIF
!
!        ==== Clear out the trash, if necessary. ====
!
         IF ( (wantt .OR. Info/=0) .AND. N>2 )                          &
     &        CALL CLASET('L',N-2,N-2,ZERO,ZERO,H(3,1),Ldh)
!
!        ==== Ensure reported workspace size is backward-compatible with
!        .    previous LAPACK versions. ====
!
         Work(1) = CMPLX(MAX(REAL(MAX(1,N)),REAL(Work(1))),RZERO)
      ENDIF
!
!     ==== End of CHSEQR ====
!
      END SUBROUTINE CHSEQR