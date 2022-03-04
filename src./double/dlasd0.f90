!*==dlasd0.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DLASD0 computes the singular values of a real upper bidiagonal n-by-m matrix B with diagonal d and off-diagonal e. Used by sbdsdc.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLASD0 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlasd0.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlasd0.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlasd0.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLASD0( N, SQRE, D, E, U, LDU, VT, LDVT, SMLSIZ, IWORK,
!                          WORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDU, LDVT, N, SMLSIZ, SQRE
!       ..
!       .. Array Arguments ..
!       INTEGER            IWORK( * )
!       DOUBLE PRECISION   D( * ), E( * ), U( LDU, * ), VT( LDVT, * ),
!      $                   WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> Using a divide and conquer approach, DLASD0 computes the singular
!> value decomposition (SVD) of a real upper bidiagonal N-by-M
!> matrix B with diagonal D and offdiagonal E, where M = N + SQRE.
!> The algorithm computes orthogonal matrices U and VT such that
!> B = U * S * VT. The singular values S are overwritten on D.
!>
!> A related subroutine, DLASDA, computes only the singular values,
!> and optionally, the singular vectors in compact form.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>         On entry, the row dimension of the upper bidiagonal matrix.
!>         This is also the dimension of the main diagonal array D.
!> \endverbatim
!>
!> \param[in] SQRE
!> \verbatim
!>          SQRE is INTEGER
!>         Specifies the column dimension of the bidiagonal matrix.
!>         = 0: The bidiagonal matrix has column dimension M = N;
!>         = 1: The bidiagonal matrix has column dimension M = N+1;
!> \endverbatim
!>
!> \param[in,out] D
!> \verbatim
!>          D is DOUBLE PRECISION array, dimension (N)
!>         On entry D contains the main diagonal of the bidiagonal
!>         matrix.
!>         On exit D, if INFO = 0, contains its singular values.
!> \endverbatim
!>
!> \param[in,out] E
!> \verbatim
!>          E is DOUBLE PRECISION array, dimension (M-1)
!>         Contains the subdiagonal entries of the bidiagonal matrix.
!>         On exit, E has been destroyed.
!> \endverbatim
!>
!> \param[out] U
!> \verbatim
!>          U is DOUBLE PRECISION array, dimension (LDU, N)
!>         On exit, U contains the left singular vectors.
!> \endverbatim
!>
!> \param[in] LDU
!> \verbatim
!>          LDU is INTEGER
!>         On entry, leading dimension of U.
!> \endverbatim
!>
!> \param[out] VT
!> \verbatim
!>          VT is DOUBLE PRECISION array, dimension (LDVT, M)
!>         On exit, VT**T contains the right singular vectors.
!> \endverbatim
!>
!> \param[in] LDVT
!> \verbatim
!>          LDVT is INTEGER
!>         On entry, leading dimension of VT.
!> \endverbatim
!>
!> \param[in] SMLSIZ
!> \verbatim
!>          SMLSIZ is INTEGER
!>         On entry, maximum size of the subproblems at the
!>         bottom of the computation tree.
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (8*N)
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (3*M**2+2*M)
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
!> \date June 2017
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
      SUBROUTINE DLASD0(N,Sqre,D,E,U,Ldu,Vt,Ldvt,Smlsiz,Iwork,Work,Info)
      USE F77KINDS                        
      USE S_DLASD1
      USE S_DLASDQ
      USE S_DLASDT
      USE S_XERBLA
      IMPLICIT NONE
!*--DLASD0158
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER :: N
      INTEGER :: Sqre
      REAL(R8KIND) , DIMENSION(*) :: D
      REAL(R8KIND) , DIMENSION(*) :: E
      REAL(R8KIND) , DIMENSION(Ldu,*) :: U
      INTEGER :: Ldu
      REAL(R8KIND) , DIMENSION(Ldvt,*) :: Vt
      INTEGER :: Ldvt
      INTEGER :: Smlsiz
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: alpha , beta
      INTEGER :: i , i1 , ic , idxq , idxqc , im1 , inode , itemp ,     &
     &           iwk , j , lf , ll , lvl , m , ncc , nd , ndb1 , ndiml ,&
     &           ndimr , nl , nlf , nlp1 , nlvl , nr , nrf , nrp1 ,     &
     &           sqrei
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
!  =====================================================================
!
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
      IF ( N<0 ) THEN
         Info = -1
      ELSEIF ( (Sqre<0) .OR. (Sqre>1) ) THEN
         Info = -2
      ENDIF
!
      m = N + Sqre
!
      IF ( Ldu<N ) THEN
         Info = -6
      ELSEIF ( Ldvt<m ) THEN
         Info = -8
      ELSEIF ( Smlsiz<3 ) THEN
         Info = -9
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('DLASD0',-Info)
         RETURN
      ENDIF
!
!     If the input matrix is too small, call DLASDQ to find the SVD.
!
      IF ( N<=Smlsiz ) THEN
         CALL DLASDQ('U',Sqre,N,m,N,0,D,E,Vt,Ldvt,U,Ldu,U,Ldu,Work,Info)
         RETURN
      ENDIF
!
!     Set up the computation tree.
!
      inode = 1
      ndiml = inode + N
      ndimr = ndiml + N
      idxq = ndimr + N
      iwk = idxq + N
      CALL DLASDT(N,nlvl,nd,Iwork(inode),Iwork(ndiml),Iwork(ndimr),     &
     &            Smlsiz)
!
!     For the nodes on bottom level of the tree, solve
!     their subproblems by DLASDQ.
!
      ndb1 = (nd+1)/2
      ncc = 0
      DO i = ndb1 , nd
!
!     IC : center row of each node
!     NL : number of rows of left  subproblem
!     NR : number of rows of right subproblem
!     NLF: starting row of the left   subproblem
!     NRF: starting row of the right  subproblem
!
         i1 = i - 1
         ic = Iwork(inode+i1)
         nl = Iwork(ndiml+i1)
         nlp1 = nl + 1
         nr = Iwork(ndimr+i1)
         nrp1 = nr + 1
         nlf = ic - nl
         nrf = ic + 1
         sqrei = 1
         CALL DLASDQ('U',sqrei,nl,nlp1,nl,ncc,D(nlf),E(nlf),Vt(nlf,nlf),&
     &               Ldvt,U(nlf,nlf),Ldu,U(nlf,nlf),Ldu,Work,Info)
         IF ( Info/=0 ) RETURN
         itemp = idxq + nlf - 2
         DO j = 1 , nl
            Iwork(itemp+j) = j
         ENDDO
         IF ( i==nd ) THEN
            sqrei = Sqre
         ELSE
            sqrei = 1
         ENDIF
         nrp1 = nr + sqrei
         CALL DLASDQ('U',sqrei,nr,nrp1,nr,ncc,D(nrf),E(nrf),Vt(nrf,nrf),&
     &               Ldvt,U(nrf,nrf),Ldu,U(nrf,nrf),Ldu,Work,Info)
         IF ( Info/=0 ) RETURN
         itemp = idxq + ic
         DO j = 1 , nr
            Iwork(itemp+j-1) = j
         ENDDO
      ENDDO
!
!     Now conquer each subproblem bottom-up.
!
      DO lvl = nlvl , 1 , -1
!
!        Find the first node LF and last node LL on the
!        current level LVL.
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
            IF ( (Sqre==0) .AND. (i==ll) ) THEN
               sqrei = Sqre
            ELSE
               sqrei = 1
            ENDIF
            idxqc = idxq + nlf - 1
            alpha = D(ic)
            beta = E(ic)
            CALL DLASD1(nl,nr,sqrei,D(nlf),alpha,beta,U(nlf,nlf),Ldu,   &
     &                  Vt(nlf,nlf),Ldvt,Iwork(idxqc),Iwork(iwk),Work,  &
     &                  Info)
!
!        Report the possible convergence failure.
!
            IF ( Info/=0 ) RETURN
         ENDDO
      ENDDO
!
!
!     End of DLASD0
!
      END SUBROUTINE DLASD0
