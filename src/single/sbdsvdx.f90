!*==sbdsvdx.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b SBDSVDX
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SBDSVDX + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sbdsvdx.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sbdsvdx.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sbdsvdx.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!     SUBROUTINE SBDSVDX( UPLO, JOBZ, RANGE, N, D, E, VL, VU, IL, IU,
!    $                    NS, S, Z, LDZ, WORK, IWORK, INFO )
!
!     .. Scalar Arguments ..
!      CHARACTER          JOBZ, RANGE, UPLO
!      INTEGER            IL, INFO, IU, LDZ, N, NS
!      REAL               VL, VU
!     ..
!     .. Array Arguments ..
!      INTEGER            IWORK( * )
!      REAL               D( * ), E( * ), S( * ), WORK( * ),
!                         Z( LDZ, * )
!       ..
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>  SBDSVDX computes the singular value decomposition (SVD) of a real
!>  N-by-N (upper or lower) bidiagonal matrix B, B = U * S * VT,
!>  where S is a diagonal matrix with non-negative diagonal elements
!>  (the singular values of B), and U and VT are orthogonal matrices
!>  of left and right singular vectors, respectively.
!>
!>  Given an upper bidiagonal B with diagonal D = [ d_1 d_2 ... d_N ]
!>  and superdiagonal E = [ e_1 e_2 ... e_N-1 ], SBDSVDX computes the
!>  singular value decompositon of B through the eigenvalues and
!>  eigenvectors of the N*2-by-N*2 tridiagonal matrix
!>
!>        |  0  d_1                |
!>        | d_1  0  e_1            |
!>  TGK = |     e_1  0  d_2        |
!>        |         d_2  .   .     |
!>        |              .   .   . |
!>
!>  If (s,u,v) is a singular triplet of B with ||u|| = ||v|| = 1, then
!>  (+/-s,q), ||q|| = 1, are eigenpairs of TGK, with q = P * ( u' +/-v' ) /
!>  sqrt(2) = ( v_1 u_1 v_2 u_2 ... v_n u_n ) / sqrt(2), and
!>  P = [ e_{n+1} e_{1} e_{n+2} e_{2} ... ].
!>
!>  Given a TGK matrix, one can either a) compute -s,-v and change signs
!>  so that the singular values (and corresponding vectors) are already in
!>  descending order (as in SGESVD/SGESDD) or b) compute s,v and reorder
!>  the values (and corresponding vectors). SBDSVDX implements a) by
!>  calling SSTEVX (bisection plus inverse iteration, to be replaced
!>  with a version of the Multiple Relative Robust Representation
!>  algorithm. (See P. Willems and B. Lang, A framework for the MR^3
!>  algorithm: theory and implementation, SIAM J. Sci. Comput.,
!>  35:740-766, 2013.)
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          = 'U':  B is upper bidiagonal;
!>          = 'L':  B is lower bidiagonal.
!> \endverbatim
!>
!> \param[in] JOBZ
!> \verbatim
!>          JOBZ is CHARACTER*1
!>          = 'N':  Compute singular values only;
!>          = 'V':  Compute singular values and singular vectors.
!> \endverbatim
!>
!> \param[in] RANGE
!> \verbatim
!>          RANGE is CHARACTER*1
!>          = 'A': all singular values will be found.
!>          = 'V': all singular values in the half-open interval [VL,VU)
!>                 will be found.
!>          = 'I': the IL-th through IU-th singular values will be found.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the bidiagonal matrix.  N >= 0.
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is REAL array, dimension (N)
!>          The n diagonal elements of the bidiagonal matrix B.
!> \endverbatim
!>
!> \param[in] E
!> \verbatim
!>          E is REAL array, dimension (max(1,N-1))
!>          The (n-1) superdiagonal elements of the bidiagonal matrix
!>          B in elements 1 to N-1.
!> \endverbatim
!>
!> \param[in] VL
!> \verbatim
!>         VL is REAL
!>          If RANGE='V', the lower bound of the interval to
!>          be searched for singular values. VU > VL.
!>          Not referenced if RANGE = 'A' or 'I'.
!> \endverbatim
!>
!> \param[in] VU
!> \verbatim
!>         VU is REAL
!>          If RANGE='V', the upper bound of the interval to
!>          be searched for singular values. VU > VL.
!>          Not referenced if RANGE = 'A' or 'I'.
!> \endverbatim
!>
!> \param[in] IL
!> \verbatim
!>          IL is INTEGER
!>          If RANGE='I', the index of the
!>          smallest singular value to be returned.
!>          1 <= IL <= IU <= min(M,N), if min(M,N) > 0.
!>          Not referenced if RANGE = 'A' or 'V'.
!> \endverbatim
!>
!> \param[in] IU
!> \verbatim
!>          IU is INTEGER
!>          If RANGE='I', the index of the
!>          largest singular value to be returned.
!>          1 <= IL <= IU <= min(M,N), if min(M,N) > 0.
!>          Not referenced if RANGE = 'A' or 'V'.
!> \endverbatim
!>
!> \param[out] NS
!> \verbatim
!>          NS is INTEGER
!>          The total number of singular values found.  0 <= NS <= N.
!>          If RANGE = 'A', NS = N, and if RANGE = 'I', NS = IU-IL+1.
!> \endverbatim
!>
!> \param[out] S
!> \verbatim
!>          S is REAL array, dimension (N)
!>          The first NS elements contain the selected singular values in
!>          ascending order.
!> \endverbatim
!>
!> \param[out] Z
!> \verbatim
!>          Z is REAL array, dimension (2*N,K)
!>          If JOBZ = 'V', then if INFO = 0 the first NS columns of Z
!>          contain the singular vectors of the matrix B corresponding to
!>          the selected singular values, with U in rows 1 to N and V
!>          in rows N+1 to N*2, i.e.
!>          Z = [ U ]
!>              [ V ]
!>          If JOBZ = 'N', then Z is not referenced.
!>          Note: The user must ensure that at least K = NS+1 columns are
!>          supplied in the array Z; if RANGE = 'V', the exact value of
!>          NS is not known in advance and an upper bound must be used.
!> \endverbatim
!>
!> \param[in] LDZ
!> \verbatim
!>          LDZ is INTEGER
!>          The leading dimension of the array Z. LDZ >= 1, and if
!>          JOBZ = 'V', LDZ >= max(2,N*2).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (14*N)
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (12*N)
!>          If JOBZ = 'V', then if INFO = 0, the first NS elements of
!>          IWORK are zero. If INFO > 0, then IWORK contains the indices
!>          of the eigenvectors that failed to converge in DSTEVX.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!>          > 0:  if INFO = i, then i eigenvectors failed to converge
!>                   in SSTEVX. The indices of the eigenvectors
!>                   (as returned by SSTEVX) are stored in the
!>                   array IWORK.
!>                if INFO = N*2 + 1, an internal error occurred.
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
!> \ingroup realOTHEReigen
!
!  =====================================================================
      SUBROUTINE SBDSVDX(Uplo,Jobz,Range,N,D,E,Vl,Vu,Il,Iu,Ns,S,Z,Ldz,  &
     &                   Work,Iwork,Info)
      IMPLICIT NONE
!*--SBDSVDX230
!
!  -- LAPACK driver routine (version 3.8.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2017
!
!     .. Scalar Arguments ..
      CHARACTER Jobz , Range , Uplo
      INTEGER Il , Info , Iu , Ldz , N , Ns
      REAL Vl , Vu
!     ..
!     .. Array Arguments ..
      INTEGER Iwork(*)
      REAL D(*) , E(*) , S(*) , Work(*) , Z(Ldz,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE , TEN , HNDRD , MEIGTH
      PARAMETER (ZERO=0.0E0,ONE=1.0E0,TEN=10.0E0,HNDRD=100.0E0,         &
     &           MEIGTH=-0.1250E0)
      REAL FUDGE
      PARAMETER (FUDGE=2.0E0)
!     ..
!     .. Local Scalars ..
      CHARACTER rngvx
      LOGICAL allsv , indsv , lower , split , sveq0 , valsv , wantz
      INTEGER i , icolz , idbeg , idend , idtgk , idptr , ieptr ,       &
     &        ietgk , iifail , iiwork , iltgk , irowu , irowv , irowz , &
     &        isbeg , isplt , itemp , iutgk , j , k , ntgk , nru , nrv ,&
     &        nsl
      REAL abstol , eps , emin , mu , nrmu , nrmv , ortol , smax ,      &
     &     smin , sqrt2 , thresh , tol , ulp , vltgk , vutgk , zjtji
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      INTEGER ISAMAX
      REAL SDOT , SLAMCH , SNRM2
      EXTERNAL ISAMAX , LSAME , SAXPY , SDOT , SLAMCH , SNRM2
!     ..
!     .. External Subroutines ..
      EXTERNAL SCOPY , SLASET , SSCAL , SSWAP , SSTEVX , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , REAL , SIGN , SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      allsv = LSAME(Range,'A')
      valsv = LSAME(Range,'V')
      indsv = LSAME(Range,'I')
      wantz = LSAME(Jobz,'V')
      lower = LSAME(Uplo,'L')
!
      Info = 0
      IF ( .NOT.LSAME(Uplo,'U') .AND. .NOT.lower ) THEN
         Info = -1
      ELSEIF ( .NOT.(wantz .OR. LSAME(Jobz,'N')) ) THEN
         Info = -2
      ELSEIF ( .NOT.(allsv .OR. valsv .OR. indsv) ) THEN
         Info = -3
      ELSEIF ( N<0 ) THEN
         Info = -4
      ELSEIF ( N>0 ) THEN
         IF ( valsv ) THEN
            IF ( Vl<ZERO ) THEN
               Info = -7
            ELSEIF ( Vu<=Vl ) THEN
               Info = -8
            ENDIF
         ELSEIF ( indsv ) THEN
            IF ( Il<1 .OR. Il>MAX(1,N) ) THEN
               Info = -9
            ELSEIF ( Iu<MIN(N,Il) .OR. Iu>N ) THEN
               Info = -10
            ENDIF
         ENDIF
      ENDIF
      IF ( Info==0 ) THEN
         IF ( Ldz<1 .OR. (wantz .AND. Ldz<N*2) ) Info = -14
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('SBDSVDX',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible (N.LE.1)
!
      Ns = 0
      IF ( N==0 ) RETURN
!
      IF ( N==1 ) THEN
         IF ( allsv .OR. indsv ) THEN
            Ns = 1
            S(1) = ABS(D(1))
         ELSEIF ( Vl<ABS(D(1)) .AND. Vu>=ABS(D(1)) ) THEN
            Ns = 1
            S(1) = ABS(D(1))
         ENDIF
         IF ( wantz ) THEN
            Z(1,1) = SIGN(ONE,D(1))
            Z(2,1) = ONE
         ENDIF
         RETURN
      ENDIF
!
      abstol = 2*SLAMCH('Safe Minimum')
      ulp = SLAMCH('Precision')
      eps = SLAMCH('Epsilon')
      sqrt2 = SQRT(2.0E0)
      ortol = SQRT(ulp)
!
!     Criterion for splitting is taken from SBDSQR when singular
!     values are computed to relative accuracy TOL. (See J. Demmel and
!     W. Kahan, Accurate singular values of bidiagonal matrices, SIAM
!     J. Sci. and Stat. Comput., 11:873–912, 1990.)
!
      tol = MAX(TEN,MIN(HNDRD,eps**MEIGTH))*eps
!
!     Compute approximate maximum, minimum singular values.
!
      i = ISAMAX(N,D,1)
      smax = ABS(D(i))
      i = ISAMAX(N-1,E,1)
      smax = MAX(smax,ABS(E(i)))
!
!     Compute threshold for neglecting D's and E's.
!
      smin = ABS(D(1))
      IF ( smin/=ZERO ) THEN
         mu = smin
         DO i = 2 , N
            mu = ABS(D(i))*(mu/(mu+ABS(E(i-1))))
            smin = MIN(smin,mu)
            IF ( smin==ZERO ) EXIT
         ENDDO
      ENDIF
      smin = smin/SQRT(REAL(N))
      thresh = tol*smin
!
!     Check for zeros in D and E (splits), i.e. submatrices.
!
      DO i = 1 , N - 1
         IF ( ABS(D(i))<=thresh ) D(i) = ZERO
         IF ( ABS(E(i))<=thresh ) E(i) = ZERO
      ENDDO
      IF ( ABS(D(N))<=thresh ) D(N) = ZERO
!
!     Pointers for arrays used by SSTEVX.
!
      idtgk = 1
      ietgk = idtgk + N*2
      itemp = ietgk + N*2
      iifail = 1
      iiwork = iifail + N*2
!
!     Set RNGVX, which corresponds to RANGE for SSTEVX in TGK mode.
!     VL,VU or IL,IU are redefined to conform to implementation a)
!     described in the leading comments.
!
      iltgk = 0
      iutgk = 0
      vltgk = ZERO
      vutgk = ZERO
!
      IF ( allsv ) THEN
!
!        All singular values will be found. We aim at -s (see
!        leading comments) with RNGVX = 'I'. IL and IU are set
!        later (as ILTGK and IUTGK) according to the dimension
!        of the active submatrix.
!
         rngvx = 'I'
         IF ( wantz ) CALL SLASET('F',N*2,N+1,ZERO,ZERO,Z,Ldz)
      ELSEIF ( valsv ) THEN
!
!        Find singular values in a half-open interval. We aim
!        at -s (see leading comments) and we swap VL and VU
!        (as VUTGK and VLTGK), changing their signs.
!
         rngvx = 'V'
         vltgk = -Vu
         vutgk = -Vl
         Work(idtgk:idtgk+2*N-1) = ZERO
         CALL SCOPY(N,D,1,Work(ietgk),2)
         CALL SCOPY(N-1,E,1,Work(ietgk+1),2)
         CALL SSTEVX('N','V',N*2,Work(idtgk),Work(ietgk),vltgk,vutgk,   &
     &               iltgk,iltgk,abstol,Ns,S,Z,Ldz,Work(itemp),         &
     &               Iwork(iiwork),Iwork(iifail),Info)
         IF ( Ns==0 ) THEN
            RETURN
         ELSE
            IF ( wantz ) CALL SLASET('F',N*2,Ns,ZERO,ZERO,Z,Ldz)
         ENDIF
      ELSEIF ( indsv ) THEN
!
!        Find the IL-th through the IU-th singular values. We aim
!        at -s (see leading comments) and indices are mapped into
!        values, therefore mimicking SSTEBZ, where
!
!        GL = GL - FUDGE*TNORM*ULP*N - FUDGE*TWO*PIVMIN
!        GU = GU + FUDGE*TNORM*ULP*N + FUDGE*PIVMIN
!
         iltgk = Il
         iutgk = Iu
         rngvx = 'V'
         Work(idtgk:idtgk+2*N-1) = ZERO
         CALL SCOPY(N,D,1,Work(ietgk),2)
         CALL SCOPY(N-1,E,1,Work(ietgk+1),2)
         CALL SSTEVX('N','I',N*2,Work(idtgk),Work(ietgk),vltgk,vltgk,   &
     &               iltgk,iltgk,abstol,Ns,S,Z,Ldz,Work(itemp),         &
     &               Iwork(iiwork),Iwork(iifail),Info)
         vltgk = S(1) - FUDGE*smax*ulp*N
         Work(idtgk:idtgk+2*N-1) = ZERO
         CALL SCOPY(N,D,1,Work(ietgk),2)
         CALL SCOPY(N-1,E,1,Work(ietgk+1),2)
         CALL SSTEVX('N','I',N*2,Work(idtgk),Work(ietgk),vutgk,vutgk,   &
     &               iutgk,iutgk,abstol,Ns,S,Z,Ldz,Work(itemp),         &
     &               Iwork(iiwork),Iwork(iifail),Info)
         vutgk = S(1) + FUDGE*smax*ulp*N
         vutgk = MIN(vutgk,ZERO)
!
!        If VLTGK=VUTGK, SSTEVX returns an error message,
!        so if needed we change VUTGK slightly.
!
         IF ( vltgk==vutgk ) vltgk = vltgk - tol
!
         IF ( wantz ) CALL SLASET('F',N*2,Iu-Il+1,ZERO,ZERO,Z,Ldz)
      ENDIF
!
!     Initialize variables and pointers for S, Z, and WORK.
!
!     NRU, NRV: number of rows in U and V for the active submatrix
!     IDBEG, ISBEG: offsets for the entries of D and S
!     IROWZ, ICOLZ: offsets for the rows and columns of Z
!     IROWU, IROWV: offsets for the rows of U and V
!
      Ns = 0
      nru = 0
      nrv = 0
      idbeg = 1
      isbeg = 1
      irowz = 1
      icolz = 1
      irowu = 2
      irowv = 1
      split = .FALSE.
      sveq0 = .FALSE.
!
!     Form the tridiagonal TGK matrix.
!
      S(1:N) = ZERO
      Work(ietgk+2*N-1) = ZERO
      Work(idtgk:idtgk+2*N-1) = ZERO
      CALL SCOPY(N,D,1,Work(ietgk),2)
      CALL SCOPY(N-1,E,1,Work(ietgk+1),2)
!
!
!     Check for splits in two levels, outer level
!     in E and inner level in D.
!
      DO ieptr = 2 , N*2 , 2
         IF ( Work(ietgk+ieptr-1)==ZERO ) THEN
!
!           Split in E (this piece of B is square) or bottom
!           of the (input bidiagonal) matrix.
!
            isplt = idbeg
            idend = ieptr - 1
            DO idptr = idbeg , idend , 2
               IF ( Work(ietgk+idptr-1)==ZERO ) THEN
!
!                 Split in D (rectangular submatrix). Set the number
!                 of rows in U and V (NRU and NRV) accordingly.
!
                  IF ( idptr==idbeg ) THEN
!
!                    D=0 at the top.
!
                     sveq0 = .TRUE.
                     IF ( idbeg==idend ) THEN
                        nru = 1
                        nrv = 1
                     ENDIF
                  ELSEIF ( idptr==idend ) THEN
!
!                    D=0 at the bottom.
!
                     sveq0 = .TRUE.
                     nru = (idend-isplt)/2 + 1
                     nrv = nru
                     IF ( isplt/=idbeg ) nru = nru + 1
                  ELSEIF ( isplt==idbeg ) THEN
!
!                       Split: top rectangular submatrix.
!
                     nru = (idptr-idbeg)/2
                     nrv = nru + 1
                  ELSE
!
!                       Split: middle square submatrix.
!
                     nru = (idptr-isplt)/2 + 1
                     nrv = nru
                  ENDIF
               ELSEIF ( idptr==idend ) THEN
!
!                 Last entry of D in the active submatrix.
!
                  IF ( isplt==idbeg ) THEN
!
!                    No split (trivial case).
!
                     nru = (idend-idbeg)/2 + 1
                     nrv = nru
                  ELSE
!
!                    Split: bottom rectangular submatrix.
!
                     nrv = (idend-isplt)/2 + 1
                     nru = nrv + 1
                  ENDIF
               ENDIF
!
               ntgk = nru + nrv
!
               IF ( ntgk>0 ) THEN
!
!                 Compute eigenvalues/vectors of the active
!                 submatrix according to RANGE:
!                 if RANGE='A' (ALLSV) then RNGVX = 'I'
!                 if RANGE='V' (VALSV) then RNGVX = 'V'
!                 if RANGE='I' (INDSV) then RNGVX = 'V'
!
                  iltgk = 1
                  iutgk = ntgk/2
                  IF ( allsv .OR. vutgk==ZERO ) THEN
!                        Special case: eigenvalue equal to zero or very
!                        small, additional eigenvector is needed.
                     IF ( sveq0 .OR. smin<eps .OR. MOD(ntgk,2)>0 )      &
     &                    iutgk = iutgk + 1
                  ENDIF
!
!                 Workspace needed by SSTEVX:
!                 WORK( ITEMP: ): 2*5*NTGK
!                 IWORK( 1: ): 2*6*NTGK
!
                  CALL SSTEVX(Jobz,rngvx,ntgk,Work(idtgk+isplt-1),      &
     &                        Work(ietgk+isplt-1),vltgk,vutgk,iltgk,    &
     &                        iutgk,abstol,nsl,S(isbeg),Z(irowz,icolz), &
     &                        Ldz,Work(itemp),Iwork(iiwork),            &
     &                        Iwork(iifail),Info)
!                    Exit with the error code from SSTEVX.
                  IF ( Info/=0 ) RETURN
                  emin = ABS(MAXVAL(S(isbeg:isbeg+nsl-1)))
!
                  IF ( nsl>0 .AND. wantz ) THEN
!
!                    Normalize u=Z([2,4,...],:) and v=Z([1,3,...],:),
!                    changing the sign of v as discussed in the leading
!                    comments. The norms of u and v may be (slightly)
!                    different from 1/sqrt(2) if the corresponding
!                    eigenvalues are very small or too close. We check
!                    those norms and, if needed, reorthogonalize the
!                    vectors.
!
                     IF ( nsl>1 .AND. vutgk==ZERO .AND. MOD(ntgk,2)     &
     &                    ==0 .AND. emin==0 .AND. .NOT.split ) THEN
!
!                       D=0 at the top or bottom of the active submatrix:
!                       one eigenvalue is equal to zero; concatenate the
!                       eigenvectors corresponding to the two smallest
!                       eigenvalues.
!
                        Z(irowz:irowz+ntgk-1,icolz+nsl-2)               &
     &                     = Z(irowz:irowz+ntgk-1,icolz+nsl-2)          &
     &                     + Z(irowz:irowz+ntgk-1,icolz+nsl-1)
                        Z(irowz:irowz+ntgk-1,icolz+nsl-1) = ZERO
!                       IF( IUTGK*2.GT.NTGK ) THEN
!                          Eigenvalue equal to zero or very small.
!                          NSL = NSL - 1
!                       END IF
                     ENDIF
!
                     DO i = 0 , MIN(nsl-1,nru-1)
                        nrmu = SNRM2(nru,Z(irowu,icolz+i),2)
                        IF ( nrmu==ZERO ) THEN
                           Info = N*2 + 1
                           RETURN
                        ENDIF
                        CALL SSCAL(nru,ONE/nrmu,Z(irowu,icolz+i),2)
                        IF ( nrmu/=ONE .AND. ABS(nrmu-ortol)*sqrt2>ONE )&
     &                       THEN
                           DO j = 0 , i - 1
                              zjtji = -SDOT(nru,Z(irowu,icolz+j),2,     &
     &                                Z(irowu,icolz+i),2)
                              CALL SAXPY(nru,zjtji,Z(irowu,icolz+j),2,  &
     &                           Z(irowu,icolz+i),2)
                           ENDDO
                           nrmu = SNRM2(nru,Z(irowu,icolz+i),2)
                           CALL SSCAL(nru,ONE/nrmu,Z(irowu,icolz+i),2)
                        ENDIF
                     ENDDO
                     DO i = 0 , MIN(nsl-1,nrv-1)
                        nrmv = SNRM2(nrv,Z(irowv,icolz+i),2)
                        IF ( nrmv==ZERO ) THEN
                           Info = N*2 + 1
                           RETURN
                        ENDIF
                        CALL SSCAL(nrv,-ONE/nrmv,Z(irowv,icolz+i),2)
                        IF ( nrmv/=ONE .AND. ABS(nrmv-ortol)*sqrt2>ONE )&
     &                       THEN
                           DO j = 0 , i - 1
                              zjtji = -SDOT(nrv,Z(irowv,icolz+j),2,     &
     &                                Z(irowv,icolz+i),2)
                              CALL SAXPY(nru,zjtji,Z(irowv,icolz+j),2,  &
     &                           Z(irowv,icolz+i),2)
                           ENDDO
                           nrmv = SNRM2(nrv,Z(irowv,icolz+i),2)
                           CALL SSCAL(nrv,ONE/nrmv,Z(irowv,icolz+i),2)
                        ENDIF
                     ENDDO
                     IF ( vutgk==ZERO .AND. idptr<idend .AND.           &
     &                    MOD(ntgk,2)>0 ) THEN
!
!                       D=0 in the middle of the active submatrix (one
!                       eigenvalue is equal to zero): save the corresponding
!                       eigenvector for later use (when bottom of the
!                       active submatrix is reached).
!
                        split = .TRUE.
                        Z(irowz:irowz+ntgk-1,N+1)                       &
     &                     = Z(irowz:irowz+ntgk-1,Ns+nsl)
                        Z(irowz:irowz+ntgk-1,Ns+nsl) = ZERO
                     ENDIF
                  ENDIF  !** WANTZ **!
!
                  nsl = MIN(nsl,nru)
                  sveq0 = .FALSE.
!
!                 Absolute values of the eigenvalues of TGK.
!
                  DO i = 0 , nsl - 1
                     S(isbeg+i) = ABS(S(isbeg+i))
                  ENDDO
!
!                 Update pointers for TGK, S and Z.
!
                  isbeg = isbeg + nsl
                  irowz = irowz + ntgk
                  icolz = icolz + nsl
                  irowu = irowz
                  irowv = irowz + 1
                  isplt = idptr + 1
                  Ns = Ns + nsl
                  nru = 0
                  nrv = 0
               ENDIF  !** NTGK.GT.0 **!
               IF ( irowz<N*2 .AND. wantz ) Z(1:irowz-1,icolz) = ZERO
            ENDDO  !** IDPTR loop **!
            IF ( split .AND. wantz ) THEN
!
!              Bring back eigenvector corresponding
!              to eigenvalue equal to zero.
!
               Z(idbeg:idend-ntgk+1,isbeg-1)                            &
     &            = Z(idbeg:idend-ntgk+1,isbeg-1)                       &
     &            + Z(idbeg:idend-ntgk+1,N+1)
               Z(idbeg:idend-ntgk+1,N+1) = 0
            ENDIF
            irowv = irowv - 1
            irowu = irowu + 1
            idbeg = ieptr + 1
            sveq0 = .FALSE.
            split = .FALSE.
         ENDIF  !** Check for split in E **!
      ENDDO  !** IEPTR loop **!
!
!     Sort the singular values into decreasing order (insertion sort on
!     singular values, but only one transposition per singular vector)
!
      DO i = 1 , Ns - 1
         k = 1
         smin = S(1)
         DO j = 2 , Ns + 1 - i
            IF ( S(j)<=smin ) THEN
               k = j
               smin = S(j)
            ENDIF
         ENDDO
         IF ( k/=Ns+1-i ) THEN
            S(k) = S(Ns+1-i)
            S(Ns+1-i) = smin
            IF ( wantz ) CALL SSWAP(N*2,Z(1,k),1,Z(1,Ns+1-i),1)
         ENDIF
      ENDDO
!
!     If RANGE=I, check for singular values/vectors to be discarded.
!
      IF ( indsv ) THEN
         k = Iu - Il + 1
         IF ( k<Ns ) THEN
            S(k+1:Ns) = ZERO
            IF ( wantz ) Z(1:N*2,k+1:Ns) = ZERO
            Ns = k
         ENDIF
      ENDIF
!
!     Reorder Z: U = Z( 1:N,1:NS ), V = Z( N+1:N*2,1:NS ).
!     If B is a lower diagonal, swap U and V.
!
      IF ( wantz ) THEN
         DO i = 1 , Ns
            CALL SCOPY(N*2,Z(1,i),1,Work,1)
            IF ( lower ) THEN
               CALL SCOPY(N,Work(2),2,Z(N+1,i),1)
               CALL SCOPY(N,Work(1),2,Z(1,i),1)
            ELSE
               CALL SCOPY(N,Work(2),2,Z(1,i),1)
               CALL SCOPY(N,Work(1),2,Z(N+1,i),1)
            ENDIF
         ENDDO
      ENDIF
!
!
!     End of SBDSVDX
!
      END SUBROUTINE SBDSVDX
