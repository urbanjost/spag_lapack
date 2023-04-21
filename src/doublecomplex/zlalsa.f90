!*==zlalsa.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b ZLALSA computes the SVD of the coefficient matrix in compact form. Used by sgelsd.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZLALSA + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlalsa.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlalsa.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlalsa.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZLALSA( ICOMPQ, SMLSIZ, N, NRHS, B, LDB, BX, LDBX, U,
!                          LDU, VT, K, DIFL, DIFR, Z, POLES, GIVPTR,
!                          GIVCOL, LDGCOL, PERM, GIVNUM, C, S, RWORK,
!                          IWORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            ICOMPQ, INFO, LDB, LDBX, LDGCOL, LDU, N, NRHS,
!      $                   SMLSIZ
!       ..
!       .. Array Arguments ..
!       INTEGER            GIVCOL( LDGCOL, * ), GIVPTR( * ), IWORK( * ),
!      $                   K( * ), PERM( LDGCOL, * )
!       DOUBLE PRECISION   C( * ), DIFL( LDU, * ), DIFR( LDU, * ),
!      $                   GIVNUM( LDU, * ), POLES( LDU, * ), RWORK( * ),
!      $                   S( * ), U( LDU, * ), VT( LDU, * ), Z( LDU, * )
!       COMPLEX*16         B( LDB, * ), BX( LDBX, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZLALSA is an itermediate step in solving the least squares problem
!> by computing the SVD of the coefficient matrix in compact form (The
!> singular vectors are computed as products of simple orthorgonal
!> matrices.).
!>
!> If ICOMPQ = 0, ZLALSA applies the inverse of the left singular vector
!> matrix of an upper bidiagonal matrix to the right hand side; and if
!> ICOMPQ = 1, ZLALSA applies the right singular vector matrix to the
!> right hand side. The singular vector matrices were generated in
!> compact form by ZLALSA.
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
!>          B is COMPLEX*16 array, dimension ( LDB, NRHS )
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
!>          BX is COMPLEX*16 array, dimension ( LDBX, NRHS )
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
!>          U is DOUBLE PRECISION array, dimension ( LDU, SMLSIZ ).
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
!>          VT is DOUBLE PRECISION array, dimension ( LDU, SMLSIZ+1 ).
!>         On entry, VT**H contains the right singular vector matrices of
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
!>          DIFL is DOUBLE PRECISION array, dimension ( LDU, NLVL ).
!>         where NLVL = INT(log_2 (N/(SMLSIZ+1))) + 1.
!> \endverbatim
!>
!> \param[in] DIFR
!> \verbatim
!>          DIFR is DOUBLE PRECISION array, dimension ( LDU, 2 * NLVL ).
!>         On entry, DIFL(*, I) and DIFR(*, 2 * I -1) record
!>         distances between singular values on the I-th level and
!>         singular values on the (I -1)-th level, and DIFR(*, 2 * I)
!>         record the normalizing factors of the right singular vectors
!>         matrices of subproblems on I-th level.
!> \endverbatim
!>
!> \param[in] Z
!> \verbatim
!>          Z is DOUBLE PRECISION array, dimension ( LDU, NLVL ).
!>         On entry, Z(1, I) contains the components of the deflation-
!>         adjusted updating row vector for subproblems on the I-th
!>         level.
!> \endverbatim
!>
!> \param[in] POLES
!> \verbatim
!>          POLES is DOUBLE PRECISION array, dimension ( LDU, 2 * NLVL ).
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
!>          GIVNUM is DOUBLE PRECISION array, dimension ( LDU, 2 * NLVL ).
!>         On entry, GIVNUM(*, 2 *I -1 : 2 * I) records the C- and S-
!>         values of Givens rotations performed on the I-th level on the
!>         computation tree.
!> \endverbatim
!>
!> \param[in] C
!> \verbatim
!>          C is DOUBLE PRECISION array, dimension ( N ).
!>         On entry, if the I-th subproblem is not square,
!>         C( I ) contains the C-value of a Givens rotation related to
!>         the right null space of the I-th subproblem.
!> \endverbatim
!>
!> \param[in] S
!> \verbatim
!>          S is DOUBLE PRECISION array, dimension ( N ).
!>         On entry, if the I-th subproblem is not square,
!>         S( I ) contains the S-value of a Givens rotation related to
!>         the right null space of the I-th subproblem.
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array, dimension at least
!>         MAX( (SMLSZ+1)*NRHS*3, N*(1+NRHS) + 2*NRHS ).
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
!> \ingroup complex16OTHERcomputational
!
!> \par Contributors:
!  ==================
!>
!>     Ming Gu and Ren-Cang Li, Computer Science Division, University of
!>       California at Berkeley, USA \n
!>     Osni Marques, LBNL/NERSC, USA \n
!
!  =====================================================================
      SUBROUTINE ZLALSA(Icompq,Smlsiz,N,Nrhs,B,Ldb,Bx,Ldbx,U,Ldu,Vt,K,  &
     &                  Difl,Difr,Z,Poles,Givptr,Givcol,Ldgcol,Perm,    &
     &                  Givnum,C,S,Rwork,Iwork,Info)
      IMPLICIT NONE
!*--ZLALSA270
!
!  -- LAPACK computational routine (version 3.7.1) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2017
!
!     .. Scalar Arguments ..
      INTEGER Icompq , Info , Ldb , Ldbx , Ldgcol , Ldu , N , Nrhs ,    &
     &        Smlsiz
!     ..
!     .. Array Arguments ..
      INTEGER Givcol(Ldgcol,*) , Givptr(*) , Iwork(*) , K(*) ,          &
     &        Perm(Ldgcol,*)
      DOUBLE PRECISION C(*) , Difl(Ldu,*) , Difr(Ldu,*) , Givnum(Ldu,*) &
     &                 , Poles(Ldu,*) , Rwork(*) , S(*) , U(Ldu,*) ,    &
     &                 Vt(Ldu,*) , Z(Ldu,*)
      COMPLEX*16 B(Ldb,*) , Bx(Ldbx,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D0,ONE=1.0D0)
!     ..
!     .. Local Scalars ..
      INTEGER i , i1 , ic , im1 , inode , j , jcol , jimag , jreal ,    &
     &        jrow , lf , ll , lvl , lvl2 , nd , ndb1 , ndiml , ndimr , &
     &        nl , nlf , nlp1 , nlvl , nr , nrf , nrp1 , sqre
!     ..
!     .. External Subroutines ..
      EXTERNAL DGEMM , DLASDT , XERBLA , ZCOPY , ZLALS0
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DBLE , DCMPLX , DIMAG
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
         CALL XERBLA('ZLALSA',-Info)
         RETURN
      ENDIF
!
!     Book-keeping and  setting up the computation tree.
!
      inode = 1
      ndiml = inode + N
      ndimr = ndiml + N
!
      CALL DLASDT(N,nlvl,nd,Iwork(inode),Iwork(ndiml),Iwork(ndimr),     &
     &            Smlsiz)
!
!     The following code applies back the left singular vector factors.
!     For applying back the right singular vector factors, go to 170.
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
               CALL ZLALS0(Icompq,nl,nr,sqre,Nrhs,B(nlf,1),Ldb,Bx(nlf,1)&
     &                     ,Ldbx,Perm(nlf,lvl),Givptr(j),               &
     &                     Givcol(nlf,lvl2),Ldgcol,Givnum(nlf,lvl2),Ldu,&
     &                     Poles(nlf,lvl2),Difl(nlf,lvl),Difr(nlf,lvl2),&
     &                     Z(nlf,lvl),K(j),C(j),S(j),Rwork,Info)
            ENDDO
         ENDDO
!
!     The nodes on the bottom level of the tree were solved
!     by DLASDQ. The corresponding right singular vector
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
!
!        Since B and BX are complex, the following call to DGEMM is
!        performed in two steps (real and imaginary parts).
!
!        CALL DGEMM( 'T', 'N', NLP1, NRHS, NLP1, ONE, VT( NLF, 1 ), LDU,
!    $               B( NLF, 1 ), LDB, ZERO, BX( NLF, 1 ), LDBX )
!
            j = nlp1*Nrhs*2
            DO jcol = 1 , Nrhs
               DO jrow = nlf , nlf + nlp1 - 1
                  j = j + 1
                  Rwork(j) = DBLE(B(jrow,jcol))
               ENDDO
            ENDDO
            CALL DGEMM('T','N',nlp1,Nrhs,nlp1,ONE,Vt(nlf,1),Ldu,        &
     &                 Rwork(1+nlp1*Nrhs*2),nlp1,ZERO,Rwork(1),nlp1)
            j = nlp1*Nrhs*2
            DO jcol = 1 , Nrhs
               DO jrow = nlf , nlf + nlp1 - 1
                  j = j + 1
                  Rwork(j) = DIMAG(B(jrow,jcol))
               ENDDO
            ENDDO
            CALL DGEMM('T','N',nlp1,Nrhs,nlp1,ONE,Vt(nlf,1),Ldu,        &
     &                 Rwork(1+nlp1*Nrhs*2),nlp1,ZERO,Rwork(1+nlp1*Nrhs)&
     &                 ,nlp1)
            jreal = 0
            jimag = nlp1*Nrhs
            DO jcol = 1 , Nrhs
               DO jrow = nlf , nlf + nlp1 - 1
                  jreal = jreal + 1
                  jimag = jimag + 1
                  Bx(jrow,jcol) = DCMPLX(Rwork(jreal),Rwork(jimag))
               ENDDO
            ENDDO
!
!        Since B and BX are complex, the following call to DGEMM is
!        performed in two steps (real and imaginary parts).
!
!        CALL DGEMM( 'T', 'N', NRP1, NRHS, NRP1, ONE, VT( NRF, 1 ), LDU,
!    $               B( NRF, 1 ), LDB, ZERO, BX( NRF, 1 ), LDBX )
!
            j = nrp1*Nrhs*2
            DO jcol = 1 , Nrhs
               DO jrow = nrf , nrf + nrp1 - 1
                  j = j + 1
                  Rwork(j) = DBLE(B(jrow,jcol))
               ENDDO
            ENDDO
            CALL DGEMM('T','N',nrp1,Nrhs,nrp1,ONE,Vt(nrf,1),Ldu,        &
     &                 Rwork(1+nrp1*Nrhs*2),nrp1,ZERO,Rwork(1),nrp1)
            j = nrp1*Nrhs*2
            DO jcol = 1 , Nrhs
               DO jrow = nrf , nrf + nrp1 - 1
                  j = j + 1
                  Rwork(j) = DIMAG(B(jrow,jcol))
               ENDDO
            ENDDO
            CALL DGEMM('T','N',nrp1,Nrhs,nrp1,ONE,Vt(nrf,1),Ldu,        &
     &                 Rwork(1+nrp1*Nrhs*2),nrp1,ZERO,Rwork(1+nrp1*Nrhs)&
     &                 ,nrp1)
            jreal = 0
            jimag = nrp1*Nrhs
            DO jcol = 1 , Nrhs
               DO jrow = nrf , nrf + nrp1 - 1
                  jreal = jreal + 1
                  jimag = jimag + 1
                  Bx(jrow,jcol) = DCMPLX(Rwork(jreal),Rwork(jimag))
               ENDDO
            ENDDO
!
         ENDDO
         GOTO 99999
      ENDIF
!
!     The nodes on the bottom level of the tree were solved
!     by DLASDQ. The corresponding left and right singular vector
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
!
!        Since B and BX are complex, the following call to DGEMM
!        is performed in two steps (real and imaginary parts).
!
!        CALL DGEMM( 'T', 'N', NL, NRHS, NL, ONE, U( NLF, 1 ), LDU,
!     $               B( NLF, 1 ), LDB, ZERO, BX( NLF, 1 ), LDBX )
!
         j = nl*Nrhs*2
         DO jcol = 1 , Nrhs
            DO jrow = nlf , nlf + nl - 1
               j = j + 1
               Rwork(j) = DBLE(B(jrow,jcol))
            ENDDO
         ENDDO
         CALL DGEMM('T','N',nl,Nrhs,nl,ONE,U(nlf,1),Ldu,                &
     &              Rwork(1+nl*Nrhs*2),nl,ZERO,Rwork(1),nl)
         j = nl*Nrhs*2
         DO jcol = 1 , Nrhs
            DO jrow = nlf , nlf + nl - 1
               j = j + 1
               Rwork(j) = DIMAG(B(jrow,jcol))
            ENDDO
         ENDDO
         CALL DGEMM('T','N',nl,Nrhs,nl,ONE,U(nlf,1),Ldu,                &
     &              Rwork(1+nl*Nrhs*2),nl,ZERO,Rwork(1+nl*Nrhs),nl)
         jreal = 0
         jimag = nl*Nrhs
         DO jcol = 1 , Nrhs
            DO jrow = nlf , nlf + nl - 1
               jreal = jreal + 1
               jimag = jimag + 1
               Bx(jrow,jcol) = DCMPLX(Rwork(jreal),Rwork(jimag))
            ENDDO
         ENDDO
!
!        Since B and BX are complex, the following call to DGEMM
!        is performed in two steps (real and imaginary parts).
!
!        CALL DGEMM( 'T', 'N', NR, NRHS, NR, ONE, U( NRF, 1 ), LDU,
!    $               B( NRF, 1 ), LDB, ZERO, BX( NRF, 1 ), LDBX )
!
         j = nr*Nrhs*2
         DO jcol = 1 , Nrhs
            DO jrow = nrf , nrf + nr - 1
               j = j + 1
               Rwork(j) = DBLE(B(jrow,jcol))
            ENDDO
         ENDDO
         CALL DGEMM('T','N',nr,Nrhs,nr,ONE,U(nrf,1),Ldu,                &
     &              Rwork(1+nr*Nrhs*2),nr,ZERO,Rwork(1),nr)
         j = nr*Nrhs*2
         DO jcol = 1 , Nrhs
            DO jrow = nrf , nrf + nr - 1
               j = j + 1
               Rwork(j) = DIMAG(B(jrow,jcol))
            ENDDO
         ENDDO
         CALL DGEMM('T','N',nr,Nrhs,nr,ONE,U(nrf,1),Ldu,                &
     &              Rwork(1+nr*Nrhs*2),nr,ZERO,Rwork(1+nr*Nrhs),nr)
         jreal = 0
         jimag = nr*Nrhs
         DO jcol = 1 , Nrhs
            DO jrow = nrf , nrf + nr - 1
               jreal = jreal + 1
               jimag = jimag + 1
               Bx(jrow,jcol) = DCMPLX(Rwork(jreal),Rwork(jimag))
            ENDDO
         ENDDO
!
      ENDDO
!
!     Next copy the rows of B that correspond to unchanged rows
!     in the bidiagonal matrix to BX.
!
      DO i = 1 , nd
         ic = Iwork(inode+i-1)
         CALL ZCOPY(Nrhs,B(ic,1),Ldb,Bx(ic,1),Ldbx)
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
            CALL ZLALS0(Icompq,nl,nr,sqre,Nrhs,Bx(nlf,1),Ldbx,B(nlf,1), &
     &                  Ldb,Perm(nlf,lvl),Givptr(j),Givcol(nlf,lvl2),   &
     &                  Ldgcol,Givnum(nlf,lvl2),Ldu,Poles(nlf,lvl2),    &
     &                  Difl(nlf,lvl),Difr(nlf,lvl2),Z(nlf,lvl),K(j),   &
     &                  C(j),S(j),Rwork,Info)
         ENDDO
      ENDDO
!
!
!
!     End of ZLALSA
!
99999 END SUBROUTINE ZLALSA
