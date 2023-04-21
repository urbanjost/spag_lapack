!*==slasda.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
 
!> \brief \b SLASDA computes the singular value decomposition (SVD) of a real upper bidiagonal matrix with diagonal d and off-diagonal e. Used by sbdsdc.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SLASDA + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slasda.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slasda.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slasda.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLASDA( ICOMPQ, SMLSIZ, N, SQRE, D, E, U, LDU, VT, K,
!                          DIFL, DIFR, Z, POLES, GIVPTR, GIVCOL, LDGCOL,
!                          PERM, GIVNUM, C, S, WORK, IWORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            ICOMPQ, INFO, LDGCOL, LDU, N, SMLSIZ, SQRE
!       ..
!       .. Array Arguments ..
!       INTEGER            GIVCOL( LDGCOL, * ), GIVPTR( * ), IWORK( * ),
!      $                   K( * ), PERM( LDGCOL, * )
!       REAL               C( * ), D( * ), DIFL( LDU, * ), DIFR( LDU, * ),
!      $                   E( * ), GIVNUM( LDU, * ), POLES( LDU, * ),
!      $                   S( * ), U( LDU, * ), VT( LDU, * ), WORK( * ),
!      $                   Z( LDU, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> Using a divide and conquer approach, SLASDA computes the singular
!> value decomposition (SVD) of a real upper bidiagonal N-by-M matrix
!> B with diagonal D and offdiagonal E, where M = N + SQRE. The
!> algorithm computes the singular values in the SVD B = U * S * VT.
!> The orthogonal matrices U and VT are optionally computed in
!> compact form.
!>
!> A related subroutine, SLASD0, computes the singular values and
!> the singular vectors in explicit form.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] ICOMPQ
!> \verbatim
!>          ICOMPQ is INTEGER
!>         Specifies whether singular vectors are to be computed
!>         in compact form, as follows
!>         = 0: Compute singular values only.
!>         = 1: Compute singular vectors of upper bidiagonal
!>              matrix in compact form.
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
!>         The row dimension of the upper bidiagonal matrix. This is
!>         also the dimension of the main diagonal array D.
!> \endverbatim
!>
!> \param[in] SQRE
!> \verbatim
!>          SQRE is INTEGER
!>         Specifies the column dimension of the bidiagonal matrix.
!>         = 0: The bidiagonal matrix has column dimension M = N;
!>         = 1: The bidiagonal matrix has column dimension M = N + 1.
!> \endverbatim
!>
!> \param[in,out] D
!> \verbatim
!>          D is REAL array, dimension ( N )
!>         On entry D contains the main diagonal of the bidiagonal
!>         matrix. On exit D, if INFO = 0, contains its singular values.
!> \endverbatim
!>
!> \param[in] E
!> \verbatim
!>          E is REAL array, dimension ( M-1 )
!>         Contains the subdiagonal entries of the bidiagonal matrix.
!>         On exit, E has been destroyed.
!> \endverbatim
!>
!> \param[out] U
!> \verbatim
!>          U is REAL array,
!>         dimension ( LDU, SMLSIZ ) if ICOMPQ = 1, and not referenced
!>         if ICOMPQ = 0. If ICOMPQ = 1, on exit, U contains the left
!>         singular vector matrices of all subproblems at the bottom
!>         level.
!> \endverbatim
!>
!> \param[in] LDU
!> \verbatim
!>          LDU is INTEGER, LDU = > N.
!>         The leading dimension of arrays U, VT, DIFL, DIFR, POLES,
!>         GIVNUM, and Z.
!> \endverbatim
!>
!> \param[out] VT
!> \verbatim
!>          VT is REAL array,
!>         dimension ( LDU, SMLSIZ+1 ) if ICOMPQ = 1, and not referenced
!>         if ICOMPQ = 0. If ICOMPQ = 1, on exit, VT**T contains the right
!>         singular vector matrices of all subproblems at the bottom
!>         level.
!> \endverbatim
!>
!> \param[out] K
!> \verbatim
!>          K is INTEGER array, dimension ( N )
!>         if ICOMPQ = 1 and dimension 1 if ICOMPQ = 0.
!>         If ICOMPQ = 1, on exit, K(I) is the dimension of the I-th
!>         secular equation on the computation tree.
!> \endverbatim
!>
!> \param[out] DIFL
!> \verbatim
!>          DIFL is REAL array, dimension ( LDU, NLVL ),
!>         where NLVL = floor(log_2 (N/SMLSIZ))).
!> \endverbatim
!>
!> \param[out] DIFR
!> \verbatim
!>          DIFR is REAL array,
!>                  dimension ( LDU, 2 * NLVL ) if ICOMPQ = 1 and
!>                  dimension ( N ) if ICOMPQ = 0.
!>         If ICOMPQ = 1, on exit, DIFL(1:N, I) and DIFR(1:N, 2 * I - 1)
!>         record distances between singular values on the I-th
!>         level and singular values on the (I -1)-th level, and
!>         DIFR(1:N, 2 * I ) contains the normalizing factors for
!>         the right singular vector matrix. See SLASD8 for details.
!> \endverbatim
!>
!> \param[out] Z
!> \verbatim
!>          Z is REAL array,
!>                  dimension ( LDU, NLVL ) if ICOMPQ = 1 and
!>                  dimension ( N ) if ICOMPQ = 0.
!>         The first K elements of Z(1, I) contain the components of
!>         the deflation-adjusted updating row vector for subproblems
!>         on the I-th level.
!> \endverbatim
!>
!> \param[out] POLES
!> \verbatim
!>          POLES is REAL array,
!>         dimension ( LDU, 2 * NLVL ) if ICOMPQ = 1, and not referenced
!>         if ICOMPQ = 0. If ICOMPQ = 1, on exit, POLES(1, 2*I - 1) and
!>         POLES(1, 2*I) contain  the new and old singular values
!>         involved in the secular equations on the I-th level.
!> \endverbatim
!>
!> \param[out] GIVPTR
!> \verbatim
!>          GIVPTR is INTEGER array,
!>         dimension ( N ) if ICOMPQ = 1, and not referenced if
!>         ICOMPQ = 0. If ICOMPQ = 1, on exit, GIVPTR( I ) records
!>         the number of Givens rotations performed on the I-th
!>         problem on the computation tree.
!> \endverbatim
!>
!> \param[out] GIVCOL
!> \verbatim
!>          GIVCOL is INTEGER array,
!>         dimension ( LDGCOL, 2 * NLVL ) if ICOMPQ = 1, and not
!>         referenced if ICOMPQ = 0. If ICOMPQ = 1, on exit, for each I,
!>         GIVCOL(1, 2 *I - 1) and GIVCOL(1, 2 *I) record the locations
!>         of Givens rotations performed on the I-th level on the
!>         computation tree.
!> \endverbatim
!>
!> \param[in] LDGCOL
!> \verbatim
!>          LDGCOL is INTEGER, LDGCOL = > N.
!>         The leading dimension of arrays GIVCOL and PERM.
!> \endverbatim
!>
!> \param[out] PERM
!> \verbatim
!>          PERM is INTEGER array, dimension ( LDGCOL, NLVL )
!>         if ICOMPQ = 1, and not referenced
!>         if ICOMPQ = 0. If ICOMPQ = 1, on exit, PERM(1, I) records
!>         permutations done on the I-th level of the computation tree.
!> \endverbatim
!>
!> \param[out] GIVNUM
!> \verbatim
!>          GIVNUM is REAL array,
!>         dimension ( LDU,  2 * NLVL ) if ICOMPQ = 1, and not
!>         referenced if ICOMPQ = 0. If ICOMPQ = 1, on exit, for each I,
!>         GIVNUM(1, 2 *I - 1) and GIVNUM(1, 2 *I) record the C- and S-
!>         values of Givens rotations performed on the I-th level on
!>         the computation tree.
!> \endverbatim
!>
!> \param[out] C
!> \verbatim
!>          C is REAL array,
!>         dimension ( N ) if ICOMPQ = 1, and dimension 1 if ICOMPQ = 0.
!>         If ICOMPQ = 1 and the I-th subproblem is not square, on exit,
!>         C( I ) contains the C-value of a Givens rotation related to
!>         the right null space of the I-th subproblem.
!> \endverbatim
!>
!> \param[out] S
!> \verbatim
!>          S is REAL array, dimension ( N ) if
!>         ICOMPQ = 1, and dimension 1 if ICOMPQ = 0. If ICOMPQ = 1
!>         and the I-th subproblem is not square, on exit, S( I )
!>         contains the S-value of a Givens rotation related to
!>         the right null space of the I-th subproblem.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension
!>         (6 * N + (SMLSIZ + 1)*(SMLSIZ + 1)).
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (7*N).
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit.
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
!>          > 0:  if INFO = 1, a singular value did not converge
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
!> \ingroup OTHERauxiliary
!
!> \par Contributors:
!  ==================
!>
!>     Ming Gu and Huan Ren, Computer Science Division, University of
!>     California at Berkeley, USA
!>
!  =====================================================================
      SUBROUTINE SLASDA(Icompq,Smlsiz,N,Sqre,D,E,U,Ldu,Vt,K,Difl,Difr,Z,&
     &                  Poles,Givptr,Givcol,Ldgcol,Perm,Givnum,C,S,Work,&
     &                  Iwork,Info)
      IMPLICIT NONE
!*--SLASDA278
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Icompq , Info , Ldgcol , Ldu , N , Smlsiz , Sqre
!     ..
!     .. Array Arguments ..
      INTEGER Givcol(Ldgcol,*) , Givptr(*) , Iwork(*) , K(*) ,          &
     &        Perm(Ldgcol,*)
      REAL C(*) , D(*) , Difl(Ldu,*) , Difr(Ldu,*) , E(*) ,             &
     &     Givnum(Ldu,*) , Poles(Ldu,*) , S(*) , U(Ldu,*) , Vt(Ldu,*) , &
     &     Work(*) , Z(Ldu,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE
      PARAMETER (ZERO=0.0E+0,ONE=1.0E+0)
!     ..
!     .. Local Scalars ..
      INTEGER i , i1 , ic , idxq , idxqi , im1 , inode , itemp , iwk ,  &
     &        j , lf , ll , lvl , lvl2 , m , ncc , nd , ndb1 , ndiml ,  &
     &        ndimr , nl , nlf , nlp1 , nlvl , nr , nrf , nrp1 , nru ,  &
     &        nwork1 , nwork2 , smlszp , sqrei , vf , vfi , vl , vli
      REAL alpha , beta
!     ..
!     .. External Subroutines ..
      EXTERNAL SCOPY , SLASD6 , SLASDQ , SLASDT , SLASET , XERBLA
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
      ELSEIF ( N<0 ) THEN
         Info = -3
      ELSEIF ( (Sqre<0) .OR. (Sqre>1) ) THEN
         Info = -4
      ELSEIF ( Ldu<(N+Sqre) ) THEN
         Info = -8
      ELSEIF ( Ldgcol<N ) THEN
         Info = -17
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('SLASDA',-Info)
         RETURN
      ENDIF
!
      m = N + Sqre
!
!     If the input matrix is too small, call SLASDQ to find the SVD.
!
      IF ( N<=Smlsiz ) THEN
         IF ( Icompq==0 ) THEN
            CALL SLASDQ('U',Sqre,N,0,0,0,D,E,Vt,Ldu,U,Ldu,U,Ldu,Work,   &
     &                  Info)
         ELSE
            CALL SLASDQ('U',Sqre,N,m,N,0,D,E,Vt,Ldu,U,Ldu,U,Ldu,Work,   &
     &                  Info)
         ENDIF
         RETURN
      ENDIF
!
!     Book-keeping and  set up the computation tree.
!
      inode = 1
      ndiml = inode + N
      ndimr = ndiml + N
      idxq = ndimr + N
      iwk = idxq + N
!
      ncc = 0
      nru = 0
!
      smlszp = Smlsiz + 1
      vf = 1
      vl = vf + m
      nwork1 = vl + m
      nwork2 = nwork1 + smlszp*smlszp
!
      CALL SLASDT(N,nlvl,nd,Iwork(inode),Iwork(ndiml),Iwork(ndimr),     &
     &            Smlsiz)
!
!     for the nodes on bottom level of the tree, solve
!     their subproblems by SLASDQ.
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
         nlp1 = nl + 1
         nr = Iwork(ndimr+i1)
         nlf = ic - nl
         nrf = ic + 1
         idxqi = idxq + nlf - 2
         vfi = vf + nlf - 1
         vli = vl + nlf - 1
         sqrei = 1
         IF ( Icompq==0 ) THEN
            CALL SLASET('A',nlp1,nlp1,ZERO,ONE,Work(nwork1),smlszp)
            CALL SLASDQ('U',sqrei,nl,nlp1,nru,ncc,D(nlf),E(nlf),        &
     &                  Work(nwork1),smlszp,Work(nwork2),nl,Work(nwork2)&
     &                  ,nl,Work(nwork2),Info)
            itemp = nwork1 + nl*smlszp
            CALL SCOPY(nlp1,Work(nwork1),1,Work(vfi),1)
            CALL SCOPY(nlp1,Work(itemp),1,Work(vli),1)
         ELSE
            CALL SLASET('A',nl,nl,ZERO,ONE,U(nlf,1),Ldu)
            CALL SLASET('A',nlp1,nlp1,ZERO,ONE,Vt(nlf,1),Ldu)
            CALL SLASDQ('U',sqrei,nl,nlp1,nl,ncc,D(nlf),E(nlf),Vt(nlf,1)&
     &                  ,Ldu,U(nlf,1),Ldu,U(nlf,1),Ldu,Work(nwork1),    &
     &                  Info)
            CALL SCOPY(nlp1,Vt(nlf,1),1,Work(vfi),1)
            CALL SCOPY(nlp1,Vt(nlf,nlp1),1,Work(vli),1)
         ENDIF
         IF ( Info/=0 ) RETURN
         DO j = 1 , nl
            Iwork(idxqi+j) = j
         ENDDO
         IF ( (i==nd) .AND. (Sqre==0) ) THEN
            sqrei = 0
         ELSE
            sqrei = 1
         ENDIF
         idxqi = idxqi + nlp1
         vfi = vfi + nlp1
         vli = vli + nlp1
         nrp1 = nr + sqrei
         IF ( Icompq==0 ) THEN
            CALL SLASET('A',nrp1,nrp1,ZERO,ONE,Work(nwork1),smlszp)
            CALL SLASDQ('U',sqrei,nr,nrp1,nru,ncc,D(nrf),E(nrf),        &
     &                  Work(nwork1),smlszp,Work(nwork2),nr,Work(nwork2)&
     &                  ,nr,Work(nwork2),Info)
            itemp = nwork1 + (nrp1-1)*smlszp
            CALL SCOPY(nrp1,Work(nwork1),1,Work(vfi),1)
            CALL SCOPY(nrp1,Work(itemp),1,Work(vli),1)
         ELSE
            CALL SLASET('A',nr,nr,ZERO,ONE,U(nrf,1),Ldu)
            CALL SLASET('A',nrp1,nrp1,ZERO,ONE,Vt(nrf,1),Ldu)
            CALL SLASDQ('U',sqrei,nr,nrp1,nr,ncc,D(nrf),E(nrf),Vt(nrf,1)&
     &                  ,Ldu,U(nrf,1),Ldu,U(nrf,1),Ldu,Work(nwork1),    &
     &                  Info)
            CALL SCOPY(nrp1,Vt(nrf,1),1,Work(vfi),1)
            CALL SCOPY(nrp1,Vt(nrf,nrp1),1,Work(vli),1)
         ENDIF
         IF ( Info/=0 ) RETURN
         DO j = 1 , nr
            Iwork(idxqi+j) = j
         ENDDO
      ENDDO
!
!     Now conquer each subproblem bottom-up.
!
      j = 2**nlvl
      DO lvl = nlvl , 1 , -1
         lvl2 = lvl*2 - 1
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
         DO i = lf , ll
            im1 = i - 1
            ic = Iwork(inode+im1)
            nl = Iwork(ndiml+im1)
            nr = Iwork(ndimr+im1)
            nlf = ic - nl
            nrf = ic + 1
            IF ( i==ll ) THEN
               sqrei = Sqre
            ELSE
               sqrei = 1
            ENDIF
            vfi = vf + nlf - 1
            vli = vl + nlf - 1
            idxqi = idxq + nlf - 1
            alpha = D(ic)
            beta = E(ic)
            IF ( Icompq==0 ) THEN
               CALL SLASD6(Icompq,nl,nr,sqrei,D(nlf),Work(vfi),Work(vli)&
     &                     ,alpha,beta,Iwork(idxqi),Perm,Givptr(1),     &
     &                     Givcol,Ldgcol,Givnum,Ldu,Poles,Difl,Difr,Z,  &
     &                     K(1),C(1),S(1),Work(nwork1),Iwork(iwk),Info)
            ELSE
               j = j - 1
               CALL SLASD6(Icompq,nl,nr,sqrei,D(nlf),Work(vfi),Work(vli)&
     &                     ,alpha,beta,Iwork(idxqi),Perm(nlf,lvl),      &
     &                     Givptr(j),Givcol(nlf,lvl2),Ldgcol,           &
     &                     Givnum(nlf,lvl2),Ldu,Poles(nlf,lvl2),        &
     &                     Difl(nlf,lvl),Difr(nlf,lvl2),Z(nlf,lvl),K(j),&
     &                     C(j),S(j),Work(nwork1),Iwork(iwk),Info)
            ENDIF
            IF ( Info/=0 ) RETURN
         ENDDO
      ENDDO
!
!
!     End of SLASDA
!
      END SUBROUTINE SLASDA
