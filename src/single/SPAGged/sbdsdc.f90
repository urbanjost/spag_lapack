!*==sbdsdc.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
 
!> \brief \b SBDSDC
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SBDSDC + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sbdsdc.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sbdsdc.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sbdsdc.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SBDSDC( UPLO, COMPQ, N, D, E, U, LDU, VT, LDVT, Q, IQ,
!                          WORK, IWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          COMPQ, UPLO
!       INTEGER            INFO, LDU, LDVT, N
!       ..
!       .. Array Arguments ..
!       INTEGER            IQ( * ), IWORK( * )
!       REAL               D( * ), E( * ), Q( * ), U( LDU, * ),
!      $                   VT( LDVT, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SBDSDC computes the singular value decomposition (SVD) of a real
!> N-by-N (upper or lower) bidiagonal matrix B:  B = U * S * VT,
!> using a divide and conquer method, where S is a diagonal matrix
!> with non-negative diagonal elements (the singular values of B), and
!> U and VT are orthogonal matrices of left and right singular vectors,
!> respectively. SBDSDC can be used to compute all singular values,
!> and optionally, singular vectors or singular vectors in compact form.
!>
!> This code makes very mild assumptions about floating point
!> arithmetic. It will work on machines with a guard digit in
!> add/subtract, or on those binary machines without guard digits
!> which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or Cray-2.
!> It could conceivably fail on hexadecimal or decimal machines
!> without guard digits, but we know of none.  See SLASD3 for details.
!>
!> The code currently calls SLASDQ if singular values only are desired.
!> However, it can be slightly modified to compute singular values
!> using the divide and conquer method.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          = 'U':  B is upper bidiagonal.
!>          = 'L':  B is lower bidiagonal.
!> \endverbatim
!>
!> \param[in] COMPQ
!> \verbatim
!>          COMPQ is CHARACTER*1
!>          Specifies whether singular vectors are to be computed
!>          as follows:
!>          = 'N':  Compute singular values only;
!>          = 'P':  Compute singular values and compute singular
!>                  vectors in compact form;
!>          = 'I':  Compute singular values and singular vectors.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix B.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] D
!> \verbatim
!>          D is REAL array, dimension (N)
!>          On entry, the n diagonal elements of the bidiagonal matrix B.
!>          On exit, if INFO=0, the singular values of B.
!> \endverbatim
!>
!> \param[in,out] E
!> \verbatim
!>          E is REAL array, dimension (N-1)
!>          On entry, the elements of E contain the offdiagonal
!>          elements of the bidiagonal matrix whose SVD is desired.
!>          On exit, E has been destroyed.
!> \endverbatim
!>
!> \param[out] U
!> \verbatim
!>          U is REAL array, dimension (LDU,N)
!>          If  COMPQ = 'I', then:
!>             On exit, if INFO = 0, U contains the left singular vectors
!>             of the bidiagonal matrix.
!>          For other values of COMPQ, U is not referenced.
!> \endverbatim
!>
!> \param[in] LDU
!> \verbatim
!>          LDU is INTEGER
!>          The leading dimension of the array U.  LDU >= 1.
!>          If singular vectors are desired, then LDU >= max( 1, N ).
!> \endverbatim
!>
!> \param[out] VT
!> \verbatim
!>          VT is REAL array, dimension (LDVT,N)
!>          If  COMPQ = 'I', then:
!>             On exit, if INFO = 0, VT**T contains the right singular
!>             vectors of the bidiagonal matrix.
!>          For other values of COMPQ, VT is not referenced.
!> \endverbatim
!>
!> \param[in] LDVT
!> \verbatim
!>          LDVT is INTEGER
!>          The leading dimension of the array VT.  LDVT >= 1.
!>          If singular vectors are desired, then LDVT >= max( 1, N ).
!> \endverbatim
!>
!> \param[out] Q
!> \verbatim
!>          Q is REAL array, dimension (LDQ)
!>          If  COMPQ = 'P', then:
!>             On exit, if INFO = 0, Q and IQ contain the left
!>             and right singular vectors in a compact form,
!>             requiring O(N log N) space instead of 2*N**2.
!>             In particular, Q contains all the REAL data in
!>             LDQ >= N*(11 + 2*SMLSIZ + 8*INT(LOG_2(N/(SMLSIZ+1))))
!>             words of memory, where SMLSIZ is returned by ILAENV and
!>             is equal to the maximum size of the subproblems at the
!>             bottom of the computation tree (usually about 25).
!>          For other values of COMPQ, Q is not referenced.
!> \endverbatim
!>
!> \param[out] IQ
!> \verbatim
!>          IQ is INTEGER array, dimension (LDIQ)
!>          If  COMPQ = 'P', then:
!>             On exit, if INFO = 0, Q and IQ contain the left
!>             and right singular vectors in a compact form,
!>             requiring O(N log N) space instead of 2*N**2.
!>             In particular, IQ contains all INTEGER data in
!>             LDIQ >= N*(3 + 3*INT(LOG_2(N/(SMLSIZ+1))))
!>             words of memory, where SMLSIZ is returned by ILAENV and
!>             is equal to the maximum size of the subproblems at the
!>             bottom of the computation tree (usually about 25).
!>          For other values of COMPQ, IQ is not referenced.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (MAX(1,LWORK))
!>          If COMPQ = 'N' then LWORK >= (4 * N).
!>          If COMPQ = 'P' then LWORK >= (6 * N).
!>          If COMPQ = 'I' then LWORK >= (3 * N**2 + 4 * N).
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (8*N)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit.
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
!>          > 0:  The algorithm failed to compute a singular value.
!>                The update process of divide and conquer failed.
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
!> \ingroup auxOTHERcomputational
!
!> \par Contributors:
!  ==================
!>
!>     Ming Gu and Huan Ren, Computer Science Division, University of
!>     California at Berkeley, USA
!>
!  =====================================================================
      SUBROUTINE SBDSDC(Uplo,Compq,N,D,E,U,Ldu,Vt,Ldvt,Q,Iq,Work,Iwork, &
     &                  Info)
      USE S_ILAENV
      USE S_LSAME
      USE S_SCOPY
      USE S_SLAMCH
      USE S_SLANST
      USE S_SLARTG
      USE S_SLASCL
      USE S_SLASD0
      USE S_SLASDA
      USE S_SLASDQ
      USE S_SLASET
      USE S_SLASR
      USE S_SSWAP
      USE S_XERBLA
      IMPLICIT NONE
!*--SBDSDC224
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0 , TWO = 2.0E+0
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Uplo
      CHARACTER :: Compq
      INTEGER :: N
      REAL , INTENT(INOUT) , DIMENSION(*) :: D
      REAL , INTENT(INOUT) , DIMENSION(*) :: E
      REAL , DIMENSION(Ldu,*) :: U
      INTEGER :: Ldu
      REAL , DIMENSION(Ldvt,*) :: Vt
      INTEGER :: Ldvt
      REAL , DIMENSION(*) :: Q
      INTEGER , DIMENSION(*) :: Iq
      REAL , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL :: cs , eps , orgnrm , p , r , sn
      INTEGER :: difl , difr , givcol , givnum , givptr , i , ic ,      &
     &           icompq , ierr , ii , is , iu , iuplo , ivt , j , k ,   &
     &           kk , mlvl , nm1 , nsize , perm , poles , qstart ,      &
     &           smlsiz , smlszp , sqre , start , wstart , z
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
!  =====================================================================
!  Changed dimension statement in comment describing E from (N) to
!  (N-1).  Sven, 17 Feb 05.
!  =====================================================================
!
!     .. Parameters ..
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
!     Test the input parameters.
!
      Info = 0
!
      iuplo = 0
      IF ( LSAME(Uplo,'U') ) iuplo = 1
      IF ( LSAME(Uplo,'L') ) iuplo = 2
      IF ( LSAME(Compq,'N') ) THEN
         icompq = 0
      ELSEIF ( LSAME(Compq,'P') ) THEN
         icompq = 1
      ELSEIF ( LSAME(Compq,'I') ) THEN
         icompq = 2
      ELSE
         icompq = -1
      ENDIF
      IF ( iuplo==0 ) THEN
         Info = -1
      ELSEIF ( icompq<0 ) THEN
         Info = -2
      ELSEIF ( N<0 ) THEN
         Info = -3
      ELSEIF ( (Ldu<1) .OR. ((icompq==2) .AND. (Ldu<N)) ) THEN
         Info = -7
      ELSEIF ( (Ldvt<1) .OR. ((icompq==2) .AND. (Ldvt<N)) ) THEN
         Info = -9
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('SBDSDC',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 ) RETURN
      smlsiz = ILAENV(9,'SBDSDC',' ',0,0,0,0)
      IF ( N==1 ) THEN
         IF ( icompq==1 ) THEN
            Q(1) = SIGN(ONE,D(1))
            Q(1+smlsiz*N) = ONE
         ELSEIF ( icompq==2 ) THEN
            U(1,1) = SIGN(ONE,D(1))
            Vt(1,1) = ONE
         ENDIF
         D(1) = ABS(D(1))
         RETURN
      ENDIF
      nm1 = N - 1
!
!     If matrix lower bidiagonal, rotate to be upper bidiagonal
!     by applying Givens rotations on the left
!
      wstart = 1
      qstart = 3
      IF ( icompq==1 ) THEN
         CALL SCOPY(N,D,1,Q(1),1)
         CALL SCOPY(N-1,E,1,Q(N+1),1)
      ENDIF
      IF ( iuplo==2 ) THEN
         qstart = 5
         IF ( icompq==2 ) wstart = 2*N - 1
         DO i = 1 , N - 1
            CALL SLARTG(D(i),E(i),cs,sn,r)
            D(i) = r
            E(i) = sn*D(i+1)
            D(i+1) = cs*D(i+1)
            IF ( icompq==1 ) THEN
               Q(i+2*N) = cs
               Q(i+3*N) = sn
            ELSEIF ( icompq==2 ) THEN
               Work(i) = cs
               Work(nm1+i) = -sn
            ENDIF
         ENDDO
      ENDIF
!
!     If ICOMPQ = 0, use SLASDQ to compute the singular values.
!
      IF ( icompq==0 ) THEN
!        Ignore WSTART, instead using WORK( 1 ), since the two vectors
!        for CS and -SN above are added only if ICOMPQ == 2,
!        and adding them exceeds documented WORK size of 4*n.
         CALL SLASDQ('U',0,N,0,0,0,D,E,Vt,Ldvt,U,Ldu,U,Ldu,Work(1),Info)
         GOTO 100
      ENDIF
!
!     If N is smaller than the minimum divide size SMLSIZ, then solve
!     the problem with another solver.
!
      IF ( N<=smlsiz ) THEN
         IF ( icompq==2 ) THEN
            CALL SLASET('A',N,N,ZERO,ONE,U,Ldu)
            CALL SLASET('A',N,N,ZERO,ONE,Vt,Ldvt)
            CALL SLASDQ('U',0,N,N,N,0,D,E,Vt,Ldvt,U,Ldu,U,Ldu,          &
     &                  Work(wstart),Info)
         ELSEIF ( icompq==1 ) THEN
            iu = 1
            ivt = iu + N
            CALL SLASET('A',N,N,ZERO,ONE,Q(iu+(qstart-1)*N),N)
            CALL SLASET('A',N,N,ZERO,ONE,Q(ivt+(qstart-1)*N),N)
            CALL SLASDQ('U',0,N,N,N,0,D,E,Q(ivt+(qstart-1)*N),N,        &
     &                  Q(iu+(qstart-1)*N),N,Q(iu+(qstart-1)*N),N,      &
     &                  Work(wstart),Info)
         ENDIF
         GOTO 100
      ENDIF
!
      IF ( icompq==2 ) THEN
         CALL SLASET('A',N,N,ZERO,ONE,U,Ldu)
         CALL SLASET('A',N,N,ZERO,ONE,Vt,Ldvt)
      ENDIF
!
!     Scale.
!
      orgnrm = SLANST('M',N,D,E)
      IF ( orgnrm==ZERO ) RETURN
      CALL SLASCL('G',0,0,orgnrm,ONE,N,1,D,N,ierr)
      CALL SLASCL('G',0,0,orgnrm,ONE,nm1,1,E,nm1,ierr)
!
      eps = SLAMCH('Epsilon')
!
      mlvl = INT(LOG(REAL(N)/REAL(smlsiz+1))/LOG(TWO)) + 1
      smlszp = smlsiz + 1
!
      IF ( icompq==1 ) THEN
         iu = 1
         ivt = 1 + smlsiz
         difl = ivt + smlszp
         difr = difl + mlvl
         z = difr + mlvl*2
         ic = z + mlvl
         is = ic + 1
         poles = is + 1
         givnum = poles + 2*mlvl
!
         k = 1
         givptr = 2
         perm = 3
         givcol = perm + mlvl
      ENDIF
!
      DO i = 1 , N
         IF ( ABS(D(i))<eps ) D(i) = SIGN(eps,D(i))
      ENDDO
!
      start = 1
      sqre = 0
!
      DO i = 1 , nm1
         IF ( (ABS(E(i))<eps) .OR. (i==nm1) ) THEN
!
!        Subproblem found. First determine its size and then
!        apply divide and conquer on it.
!
            IF ( i<nm1 ) THEN
!
!        A subproblem with E(I) small for I < NM1.
!
               nsize = i - start + 1
            ELSEIF ( ABS(E(i))>=eps ) THEN
!
!        A subproblem with E(NM1) not too small but I = NM1.
!
               nsize = N - start + 1
            ELSE
!
!        A subproblem with E(NM1) small. This implies an
!        1-by-1 subproblem at D(N). Solve this 1-by-1 problem
!        first.
!
               nsize = i - start + 1
               IF ( icompq==2 ) THEN
                  U(N,N) = SIGN(ONE,D(N))
                  Vt(N,N) = ONE
               ELSEIF ( icompq==1 ) THEN
                  Q(N+(qstart-1)*N) = SIGN(ONE,D(N))
                  Q(N+(smlsiz+qstart-1)*N) = ONE
               ENDIF
               D(N) = ABS(D(N))
            ENDIF
            IF ( icompq==2 ) THEN
               CALL SLASD0(nsize,sqre,D(start),E(start),U(start,start), &
     &                     Ldu,Vt(start,start),Ldvt,smlsiz,Iwork,       &
     &                     Work(wstart),Info)
            ELSE
               CALL SLASDA(icompq,smlsiz,nsize,sqre,D(start),E(start),  &
     &                     Q(start+(iu+qstart-2)*N),N,                  &
     &                     Q(start+(ivt+qstart-2)*N),Iq(start+k*N),     &
     &                     Q(start+(difl+qstart-2)*N),                  &
     &                     Q(start+(difr+qstart-2)*N),                  &
     &                     Q(start+(z+qstart-2)*N),                     &
     &                     Q(start+(poles+qstart-2)*N),                 &
     &                     Iq(start+givptr*N),Iq(start+givcol*N),N,     &
     &                     Iq(start+perm*N),Q(start+(givnum+qstart-2)*N)&
     &                     ,Q(start+(ic+qstart-2)*N),                   &
     &                     Q(start+(is+qstart-2)*N),Work(wstart),Iwork, &
     &                     Info)
            ENDIF
            IF ( Info/=0 ) RETURN
            start = i + 1
         ENDIF
      ENDDO
!
!     Unscale
!
      CALL SLASCL('G',0,0,ONE,orgnrm,N,1,D,N,ierr)
!
!     Use Selection Sort to minimize swaps of singular vectors
!
 100  DO ii = 2 , N
         i = ii - 1
         kk = i
         p = D(i)
         DO j = ii , N
            IF ( D(j)>p ) THEN
               kk = j
               p = D(j)
            ENDIF
         ENDDO
         IF ( kk/=i ) THEN
            D(kk) = D(i)
            D(i) = p
            IF ( icompq==1 ) THEN
               Iq(i) = kk
            ELSEIF ( icompq==2 ) THEN
               CALL SSWAP(N,U(1,i),1,U(1,kk),1)
               CALL SSWAP(N,Vt(i,1),Ldvt,Vt(kk,1),Ldvt)
            ENDIF
         ELSEIF ( icompq==1 ) THEN
            Iq(i) = i
         ENDIF
      ENDDO
!
!     If ICOMPQ = 1, use IQ(N,1) as the indicator for UPLO
!
      IF ( icompq==1 ) THEN
         IF ( iuplo==1 ) THEN
            Iq(N) = 1
         ELSE
            Iq(N) = 0
         ENDIF
      ENDIF
!
!     If B is lower bidiagonal, update U by those Givens rotations
!     which rotated B to be upper bidiagonal
!
      IF ( (iuplo==2) .AND. (icompq==2) )                               &
     &     CALL SLASR('L','V','B',N,N,Work(1),Work(N),U,Ldu)
!
!
!     End of SBDSDC
!
      END SUBROUTINE SBDSDC
