!*==slalsa.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b SLALSA computes the SVD of the coefficient matrix in compact form. Used by sgelsd.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SLALSA + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slalsa.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slalsa.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slalsa.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLALSA( ICOMPQ, SMLSIZ, N, NRHS, B, LDB, BX, LDBX, U,
!                          LDU, VT, K, DIFL, DIFR, Z, POLES, GIVPTR,
!                          GIVCOL, LDGCOL, PERM, GIVNUM, C, S, WORK,
!                          IWORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            ICOMPQ, INFO, LDB, LDBX, LDGCOL, LDU, N, NRHS,
!      $                   SMLSIZ
!       ..
!       .. Array Arguments ..
!       INTEGER            GIVCOL( LDGCOL, * ), GIVPTR( * ), IWORK( * ),
!      $                   K( * ), PERM( LDGCOL, * )
!       REAL               B( LDB, * ), BX( LDBX, * ), C( * ),
!      $                   DIFL( LDU, * ), DIFR( LDU, * ),
!      $                   GIVNUM( LDU, * ), POLES( LDU, * ), S( * ),
!      $                   U( LDU, * ), VT( LDU, * ), WORK( * ),
!      $                   Z( LDU, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLALSA is an itermediate step in solving the least squares problem
!> by computing the SVD of the coefficient matrix in compact form (The
!> singular vectors are computed as products of simple orthorgonal
!> matrices.).
!>
!> If ICOMPQ = 0, SLALSA applies the inverse of the left singular vector
!> matrix of an upper bidiagonal matrix to the right hand side; and if
!> ICOMPQ = 1, SLALSA applies the right singular vector matrix to the
!> right hand side. The singular vector matrices were generated in
!> compact form by SLALSA.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] ICOMPQ
!> \verbatim
!>          ICOMPQ is INTEGER
!>         Specifies whether the left or the right singular vector
!>         matrix is involved.
!>         = 0: Left singular vector matrix
!>         = 1: Right singular vector matrix
!> \endverbatim
!>
!> \param[in] SMLSIZ
!> \verbatim
!>          SMLSIZ is INTEGER
!>         The maximum size of the subproblems at the bottom of the
!>         computation tree.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>         The row and column dimensions of the upper bidiagonal matrix.
!> \endverbatim
!>
!> \param[in] NRHS
!> \verbatim
!>          NRHS is INTEGER
!>         The number of columns of B and BX. NRHS must be at least 1.
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is REAL array, dimension ( LDB, NRHS )
!>         On input, B contains the right hand sides of the least
!>         squares problem in rows 1 through M.
!>         On output, B contains the solution X in rows 1 through N.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>         The leading dimension of B in the calling subprogram.
!>         LDB must be at least max(1,MAX( M, N ) ).
!> \endverbatim
!>
!> \param[out] BX
!> \verbatim
!>          BX is REAL array, dimension ( LDBX, NRHS )
!>         On exit, the result of applying the left or right singular
!>         vector matrix to B.
!> \endverbatim
!>
!> \param[in] LDBX
!> \verbatim
!>          LDBX is INTEGER
!>         The leading dimension of BX.
!> \endverbatim
!>
!> \param[in] U
!> \verbatim
!>          U is REAL array, dimension ( LDU, SMLSIZ ).
!>         On entry, U contains the left singular vector matrices of all
!>         subproblems at the bottom level.
!> \endverbatim
!>
!> \param[in] LDU
!> \verbatim
!>          LDU is INTEGER, LDU = > N.
!>         The leading dimension of arrays U, VT, DIFL, DIFR,
!>         POLES, GIVNUM, and Z.
!> \endverbatim
!>
!> \param[in] VT
!> \verbatim
!>          VT is REAL array, dimension ( LDU, SMLSIZ+1 ).
!>         On entry, VT**T contains the right singular vector matrices of
!>         all subproblems at the bottom level.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is INTEGER array, dimension ( N ).
!> \endverbatim
!>
!> \param[in] DIFL
!> \verbatim
!>          DIFL is REAL array, dimension ( LDU, NLVL ).
!>         where NLVL = INT(log_2 (N/(SMLSIZ+1))) + 1.
!> \endverbatim
!>
!> \param[in] DIFR
!> \verbatim
!>          DIFR is REAL array, dimension ( LDU, 2 * NLVL ).
!>         On entry, DIFL(*, I) and DIFR(*, 2 * I -1) record
!>         distances between singular values on the I-th level and
!>         singular values on the (I -1)-th level, and DIFR(*, 2 * I)
!>         record the normalizing factors of the right singular vectors
!>         matrices of subproblems on I-th level.
!> \endverbatim
!>
!> \param[in] Z
!> \verbatim
!>          Z is REAL array, dimension ( LDU, NLVL ).
!>         On entry, Z(1, I) contains the components of the deflation-
!>         adjusted updating row vector for subproblems on the I-th
!>         level.
!> \endverbatim
!>
!> \param[in] POLES
!> \verbatim
!>          POLES is REAL array, dimension ( LDU, 2 * NLVL ).
!>         On entry, POLES(*, 2 * I -1: 2 * I) contains the new and old
!>         singular values involved in the secular equations on the I-th
!>         level.
!> \endverbatim
!>
!> \param[in] GIVPTR
!> \verbatim
!>          GIVPTR is INTEGER array, dimension ( N ).
!>         On entry, GIVPTR( I ) records the number of Givens
!>         rotations performed on the I-th problem on the computation
!>         tree.
!> \endverbatim
!>
!> \param[in] GIVCOL
!> \verbatim
!>          GIVCOL is INTEGER array, dimension ( LDGCOL, 2 * NLVL ).
!>         On entry, for each I, GIVCOL(*, 2 * I - 1: 2 * I) records the
!>         locations of Givens rotations performed on the I-th level on
!>         the computation tree.
!> \endverbatim
!>
!> \param[in] LDGCOL
!> \verbatim
!>          LDGCOL is INTEGER, LDGCOL = > N.
!>         The leading dimension of arrays GIVCOL and PERM.
!> \endverbatim
!>
!> \param[in] PERM
!> \verbatim
!>          PERM is INTEGER array, dimension ( LDGCOL, NLVL ).
!>         On entry, PERM(*, I) records permutations done on the I-th
!>         level of the computation tree.
!> \endverbatim
!>
!> \param[in] GIVNUM
!> \verbatim
!>          GIVNUM is REAL array, dimension ( LDU, 2 * NLVL ).
!>         On entry, GIVNUM(*, 2 *I -1 : 2 * I) records the C- and S-
!>         values of Givens rotations performed on the I-th level on the
!>         computation tree.
!> \endverbatim
!>
!> \param[in] C
!> \verbatim
!>          C is REAL array, dimension ( N ).
!>         On entry, if the I-th subproblem is not square,
!>         C( I ) contains the C-value of a Givens rotation related to
!>         the right null space of the I-th subproblem.
!> \endverbatim
!>
!> \param[in] S
!> \verbatim
!>          S is REAL array, dimension ( N ).
!>         On entry, if the I-th subproblem is not square,
!>         S( I ) contains the S-value of a Givens rotation related to
!>         the right null space of the I-th subproblem.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (N)
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (3*N)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit.
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
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
!> \ingroup realOTHERcomputational
!
!> \par Contributors:
!  ==================
!>
!>     Ming Gu and Ren-Cang Li, Computer Science Division, University of
!>       California at Berkeley, USA \n
!>     Osni Marques, LBNL/NERSC, USA \n
!
!  =====================================================================
      SUBROUTINE SLALSA(Icompq,Smlsiz,N,Nrhs,B,Ldb,Bx,Ldbx,U,Ldu,Vt,K,  &
     &                  Difl,Difr,Z,Poles,Givptr,Givcol,Ldgcol,Perm,    &
     &                  Givnum,C,S,Work,Iwork,Info)
      USE S_SCOPY
      USE S_SGEMM
      USE S_SLALS0
      USE S_SLASDT
      USE S_XERBLA
      IMPLICIT NONE
!*--SLALSA275
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER :: Icompq
      INTEGER :: Smlsiz
      INTEGER :: N
      INTEGER :: Nrhs
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , DIMENSION(Ldbx,*) :: Bx
      INTEGER :: Ldbx
      REAL , DIMENSION(Ldu,*) :: U
      INTEGER :: Ldu
      REAL , DIMENSION(Ldu,*) :: Vt
      INTEGER , DIMENSION(*) :: K
      REAL , DIMENSION(Ldu,*) :: Difl
      REAL , DIMENSION(Ldu,*) :: Difr
      REAL , DIMENSION(Ldu,*) :: Z
      REAL , DIMENSION(Ldu,*) :: Poles
      INTEGER , DIMENSION(*) :: Givptr
      INTEGER , DIMENSION(Ldgcol,*) :: Givcol
      INTEGER :: Ldgcol
      INTEGER , DIMENSION(Ldgcol,*) :: Perm
      REAL , DIMENSION(Ldu,*) :: Givnum
      REAL , DIMENSION(*) :: C
      REAL , DIMENSION(*) :: S
      REAL , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: i , i1 , ic , im1 , inode , j , lf , ll , lvl , lvl2 , &
     &           nd , ndb1 , ndiml , ndimr , nl , nlf , nlp1 , nlvl ,   &
     &           nr , nrf , nrp1 , sqre
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
!     ..
!     .. Local Scalars ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      Info = 0
!
      IF ( (Icompq<0) .OR. (Icompq>1) ) THEN
         Info = -1
      ELSEIF ( Smlsiz<3 ) THEN
         Info = -2
      ELSEIF ( N<Smlsiz ) THEN
         Info = -3
      ELSEIF ( Nrhs<1 ) THEN
         Info = -4
      ELSEIF ( Ldb<N ) THEN
         Info = -6
      ELSEIF ( Ldbx<N ) THEN
         Info = -8
      ELSEIF ( Ldu<N ) THEN
         Info = -10
      ELSEIF ( Ldgcol<N ) THEN
         Info = -19
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('SLALSA',-Info)
         RETURN
      ENDIF
!
!     Book-keeping and  setting up the computation tree.
!
      inode = 1
      ndiml = inode + N
      ndimr = ndiml + N
!
      CALL SLASDT(N,nlvl,nd,Iwork(inode),Iwork(ndiml),Iwork(ndimr),     &
     &            Smlsiz)
!
!     The following code applies back the left singular vector factors.
!     For applying back the right singular vector factors, go to 50.
!
      IF ( Icompq==1 ) THEN
!
!     ICOMPQ = 1: applying back the right singular vector factors.
!
!
!     First now go through the right singular vector matrices of all
!     the tree nodes top-down.
!
         j = 0
         DO lvl = 1 , nlvl
            lvl2 = 2*lvl - 1
!
!        Find the first node LF and last node LL on
!        the current level LVL.
!
            IF ( lvl==1 ) THEN
               lf = 1
               ll = 1
            ELSE
               lf = 2**(lvl-1)
               ll = 2*lf - 1
            ENDIF
            DO i = ll , lf , -1
               im1 = i - 1
               ic = Iwork(inode+im1)
               nl = Iwork(ndiml+im1)
               nr = Iwork(ndimr+im1)
               nlf = ic - nl
               nrf = ic + 1
               IF ( i==ll ) THEN
                  sqre = 0
               ELSE
                  sqre = 1
               ENDIF
               j = j + 1
               CALL SLALS0(Icompq,nl,nr,sqre,Nrhs,B(nlf,1),Ldb,Bx(nlf,1)&
     &                     ,Ldbx,Perm(nlf,lvl),Givptr(j),               &
     &                     Givcol(nlf,lvl2),Ldgcol,Givnum(nlf,lvl2),Ldu,&
     &                     Poles(nlf,lvl2),Difl(nlf,lvl),Difr(nlf,lvl2),&
     &                     Z(nlf,lvl),K(j),C(j),S(j),Work,Info)
            ENDDO
         ENDDO
!
!     The nodes on the bottom level of the tree were solved
!     by SLASDQ. The corresponding right singular vector
!     matrices are in explicit form. Apply them back.
!
         ndb1 = (nd+1)/2
         DO i = ndb1 , nd
            i1 = i - 1
            ic = Iwork(inode+i1)
            nl = Iwork(ndiml+i1)
            nr = Iwork(ndimr+i1)
            nlp1 = nl + 1
            IF ( i==nd ) THEN
               nrp1 = nr
            ELSE
               nrp1 = nr + 1
            ENDIF
            nlf = ic - nl
            nrf = ic + 1
            CALL SGEMM('T','N',nlp1,Nrhs,nlp1,ONE,Vt(nlf,1),Ldu,B(nlf,1)&
     &                 ,Ldb,ZERO,Bx(nlf,1),Ldbx)
            CALL SGEMM('T','N',nrp1,Nrhs,nrp1,ONE,Vt(nrf,1),Ldu,B(nrf,1)&
     &                 ,Ldb,ZERO,Bx(nrf,1),Ldbx)
         ENDDO
         GOTO 99999
      ENDIF
!
!     The nodes on the bottom level of the tree were solved
!     by SLASDQ. The corresponding left and right singular vector
!     matrices are in explicit form. First apply back the left
!     singular vector matrices.
!
      ndb1 = (nd+1)/2
      DO i = ndb1 , nd
!
!        IC : center row of each node
!        NL : number of rows of left  subproblem
!        NR : number of rows of right subproblem
!        NLF: starting row of the left   subproblem
!        NRF: starting row of the right  subproblem
!
         i1 = i - 1
         ic = Iwork(inode+i1)
         nl = Iwork(ndiml+i1)
         nr = Iwork(ndimr+i1)
         nlf = ic - nl
         nrf = ic + 1
         CALL SGEMM('T','N',nl,Nrhs,nl,ONE,U(nlf,1),Ldu,B(nlf,1),Ldb,   &
     &              ZERO,Bx(nlf,1),Ldbx)
         CALL SGEMM('T','N',nr,Nrhs,nr,ONE,U(nrf,1),Ldu,B(nrf,1),Ldb,   &
     &              ZERO,Bx(nrf,1),Ldbx)
      ENDDO
!
!     Next copy the rows of B that correspond to unchanged rows
!     in the bidiagonal matrix to BX.
!
      DO i = 1 , nd
         ic = Iwork(inode+i-1)
         CALL SCOPY(Nrhs,B(ic,1),Ldb,Bx(ic,1),Ldbx)
      ENDDO
!
!     Finally go through the left singular vector matrices of all
!     the other subproblems bottom-up on the tree.
!
      j = 2**nlvl
      sqre = 0
!
      DO lvl = nlvl , 1 , -1
         lvl2 = 2*lvl - 1
!
!        find the first node LF and last node LL on
!        the current level LVL
!
         IF ( lvl==1 ) THEN
            lf = 1
            ll = 1
         ELSE
            lf = 2**(lvl-1)
            ll = 2*lf - 1
         ENDIF
         DO i = lf , ll
            im1 = i - 1
            ic = Iwork(inode+im1)
            nl = Iwork(ndiml+im1)
            nr = Iwork(ndimr+im1)
            nlf = ic - nl
            nrf = ic + 1
            j = j - 1
            CALL SLALS0(Icompq,nl,nr,sqre,Nrhs,Bx(nlf,1),Ldbx,B(nlf,1), &
     &                  Ldb,Perm(nlf,lvl),Givptr(j),Givcol(nlf,lvl2),   &
     &                  Ldgcol,Givnum(nlf,lvl2),Ldu,Poles(nlf,lvl2),    &
     &                  Difl(nlf,lvl),Difr(nlf,lvl2),Z(nlf,lvl),K(j),   &
     &                  C(j),S(j),Work,Info)
         ENDDO
      ENDDO
!
!
!
!     End of SLALSA
!
99999 END SUBROUTINE SLALSA
