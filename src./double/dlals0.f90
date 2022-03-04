!*==dlals0.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DLALS0 applies back multiplying factors in solving the least squares problem using divide and conquer SVD approach. Used by sgelsd.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLALS0 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlals0.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlals0.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlals0.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLALS0( ICOMPQ, NL, NR, SQRE, NRHS, B, LDB, BX, LDBX,
!                          PERM, GIVPTR, GIVCOL, LDGCOL, GIVNUM, LDGNUM,
!                          POLES, DIFL, DIFR, Z, K, C, S, WORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            GIVPTR, ICOMPQ, INFO, K, LDB, LDBX, LDGCOL,
!      $                   LDGNUM, NL, NR, NRHS, SQRE
!       DOUBLE PRECISION   C, S
!       ..
!       .. Array Arguments ..
!       INTEGER            GIVCOL( LDGCOL, * ), PERM( * )
!       DOUBLE PRECISION   B( LDB, * ), BX( LDBX, * ), DIFL( * ),
!      $                   DIFR( LDGNUM, * ), GIVNUM( LDGNUM, * ),
!      $                   POLES( LDGNUM, * ), WORK( * ), Z( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLALS0 applies back the multiplying factors of either the left or the
!> right singular vector matrix of a diagonal matrix appended by a row
!> to the right hand side matrix B in solving the least squares problem
!> using the divide-and-conquer SVD approach.
!>
!> For the left singular vector matrix, three types of orthogonal
!> matrices are involved:
!>
!> (1L) Givens rotations: the number of such rotations is GIVPTR; the
!>      pairs of columns/rows they were applied to are stored in GIVCOL;
!>      and the C- and S-values of these rotations are stored in GIVNUM.
!>
!> (2L) Permutation. The (NL+1)-st row of B is to be moved to the first
!>      row, and for J=2:N, PERM(J)-th row of B is to be moved to the
!>      J-th row.
!>
!> (3L) The left singular vector matrix of the remaining matrix.
!>
!> For the right singular vector matrix, four types of orthogonal
!> matrices are involved:
!>
!> (1R) The right singular vector matrix of the remaining matrix.
!>
!> (2R) If SQRE = 1, one extra Givens rotation to generate the right
!>      null space.
!>
!> (3R) The inverse transformation of (2L).
!>
!> (4R) The inverse transformation of (1L).
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] ICOMPQ
!> \verbatim
!>          ICOMPQ is INTEGER
!>         Specifies whether singular vectors are to be computed in
!>         factored form:
!>         = 0: Left singular vector matrix.
!>         = 1: Right singular vector matrix.
!> \endverbatim
!>
!> \param[in] NL
!> \verbatim
!>          NL is INTEGER
!>         The row dimension of the upper block. NL >= 1.
!> \endverbatim
!>
!> \param[in] NR
!> \verbatim
!>          NR is INTEGER
!>         The row dimension of the lower block. NR >= 1.
!> \endverbatim
!>
!> \param[in] SQRE
!> \verbatim
!>          SQRE is INTEGER
!>         = 0: the lower block is an NR-by-NR square matrix.
!>         = 1: the lower block is an NR-by-(NR+1) rectangular matrix.
!>
!>         The bidiagonal matrix has row dimension N = NL + NR + 1,
!>         and column dimension M = N + SQRE.
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
!>          B is DOUBLE PRECISION array, dimension ( LDB, NRHS )
!>         On input, B contains the right hand sides of the least
!>         squares problem in rows 1 through M. On output, B contains
!>         the solution X in rows 1 through N.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>         The leading dimension of B. LDB must be at least
!>         max(1,MAX( M, N ) ).
!> \endverbatim
!>
!> \param[out] BX
!> \verbatim
!>          BX is DOUBLE PRECISION array, dimension ( LDBX, NRHS )
!> \endverbatim
!>
!> \param[in] LDBX
!> \verbatim
!>          LDBX is INTEGER
!>         The leading dimension of BX.
!> \endverbatim
!>
!> \param[in] PERM
!> \verbatim
!>          PERM is INTEGER array, dimension ( N )
!>         The permutations (from deflation and sorting) applied
!>         to the two blocks.
!> \endverbatim
!>
!> \param[in] GIVPTR
!> \verbatim
!>          GIVPTR is INTEGER
!>         The number of Givens rotations which took place in this
!>         subproblem.
!> \endverbatim
!>
!> \param[in] GIVCOL
!> \verbatim
!>          GIVCOL is INTEGER array, dimension ( LDGCOL, 2 )
!>         Each pair of numbers indicates a pair of rows/columns
!>         involved in a Givens rotation.
!> \endverbatim
!>
!> \param[in] LDGCOL
!> \verbatim
!>          LDGCOL is INTEGER
!>         The leading dimension of GIVCOL, must be at least N.
!> \endverbatim
!>
!> \param[in] GIVNUM
!> \verbatim
!>          GIVNUM is DOUBLE PRECISION array, dimension ( LDGNUM, 2 )
!>         Each number indicates the C or S value used in the
!>         corresponding Givens rotation.
!> \endverbatim
!>
!> \param[in] LDGNUM
!> \verbatim
!>          LDGNUM is INTEGER
!>         The leading dimension of arrays DIFR, POLES and
!>         GIVNUM, must be at least K.
!> \endverbatim
!>
!> \param[in] POLES
!> \verbatim
!>          POLES is DOUBLE PRECISION array, dimension ( LDGNUM, 2 )
!>         On entry, POLES(1:K, 1) contains the new singular
!>         values obtained from solving the secular equation, and
!>         POLES(1:K, 2) is an array containing the poles in the secular
!>         equation.
!> \endverbatim
!>
!> \param[in] DIFL
!> \verbatim
!>          DIFL is DOUBLE PRECISION array, dimension ( K ).
!>         On entry, DIFL(I) is the distance between I-th updated
!>         (undeflated) singular value and the I-th (undeflated) old
!>         singular value.
!> \endverbatim
!>
!> \param[in] DIFR
!> \verbatim
!>          DIFR is DOUBLE PRECISION array, dimension ( LDGNUM, 2 ).
!>         On entry, DIFR(I, 1) contains the distances between I-th
!>         updated (undeflated) singular value and the I+1-th
!>         (undeflated) old singular value. And DIFR(I, 2) is the
!>         normalizing factor for the I-th right singular vector.
!> \endverbatim
!>
!> \param[in] Z
!> \verbatim
!>          Z is DOUBLE PRECISION array, dimension ( K )
!>         Contain the components of the deflation-adjusted updating row
!>         vector.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>         Contains the dimension of the non-deflated matrix,
!>         This is the order of the related secular equation. 1 <= K <=N.
!> \endverbatim
!>
!> \param[in] C
!> \verbatim
!>          C is DOUBLE PRECISION
!>         C contains garbage if SQRE =0 and the C-value of a Givens
!>         rotation related to the right null space if SQRE = 1.
!> \endverbatim
!>
!> \param[in] S
!> \verbatim
!>          S is DOUBLE PRECISION
!>         S contains garbage if SQRE =0 and the S-value of a Givens
!>         rotation related to the right null space if SQRE = 1.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension ( K )
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
!> \date December 2016
!
!> \ingroup doubleOTHERcomputational
!
!> \par Contributors:
!  ==================
!>
!>     Ming Gu and Ren-Cang Li, Computer Science Division, University of
!>       California at Berkeley, USA \n
!>     Osni Marques, LBNL/NERSC, USA \n
!
!  =====================================================================
      SUBROUTINE DLALS0(Icompq,Nl,Nr,Sqre,Nrhs,B,Ldb,Bx,Ldbx,Perm,      &
     &                  Givptr,Givcol,Ldgcol,Givnum,Ldgnum,Poles,Difl,  &
     &                  Difr,Z,K,C,S,Work,Info)
      USE F77KINDS                        
      USE S_DCOPY
      USE S_DGEMV
      USE S_DLACPY
      USE S_DLAMC3
      USE S_DLASCL
      USE S_DNRM2
      USE S_DROT
      USE S_DSCAL
      USE S_XERBLA
      IMPLICIT NONE
!*--DLALS0282
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D0 , ZERO = 0.0D0 ,        &
     &                              NEGONE = -1.0D0
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: Icompq
      INTEGER , INTENT(IN) :: Nl
      INTEGER , INTENT(IN) :: Nr
      INTEGER , INTENT(IN) :: Sqre
      INTEGER :: Nrhs
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , DIMENSION(Ldbx,*) :: Bx
      INTEGER :: Ldbx
      INTEGER , INTENT(IN) , DIMENSION(*) :: Perm
      INTEGER , INTENT(IN) :: Givptr
      INTEGER , INTENT(IN) , DIMENSION(Ldgcol,*) :: Givcol
      INTEGER , INTENT(IN) :: Ldgcol
      REAL(R8KIND) , DIMENSION(Ldgnum,*) :: Givnum
      INTEGER , INTENT(IN) :: Ldgnum
      REAL(R8KIND) , DIMENSION(Ldgnum,*) :: Poles
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: Difl
      REAL(R8KIND) , INTENT(IN) , DIMENSION(Ldgnum,*) :: Difr
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: Z
      INTEGER :: K
      REAL(R8KIND) :: C
      REAL(R8KIND) :: S
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: diflj , difrj , dj , dsigj , dsigjp , temp
      INTEGER :: i , j , m , n , nlp1
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
!     .. External Functions ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      Info = 0
      n = Nl + Nr + 1
!
      IF ( (Icompq<0) .OR. (Icompq>1) ) THEN
         Info = -1
      ELSEIF ( Nl<1 ) THEN
         Info = -2
      ELSEIF ( Nr<1 ) THEN
         Info = -3
      ELSEIF ( (Sqre<0) .OR. (Sqre>1) ) THEN
         Info = -4
      ELSEIF ( Nrhs<1 ) THEN
         Info = -5
      ELSEIF ( Ldb<n ) THEN
         Info = -7
      ELSEIF ( Ldbx<n ) THEN
         Info = -9
      ELSEIF ( Givptr<0 ) THEN
         Info = -11
      ELSEIF ( Ldgcol<n ) THEN
         Info = -13
      ELSEIF ( Ldgnum<n ) THEN
         Info = -15
      ELSEIF ( K<1 ) THEN
         Info = -20
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('DLALS0',-Info)
         RETURN
      ENDIF
!
      m = n + Sqre
      nlp1 = Nl + 1
!
      IF ( Icompq==0 ) THEN
!
!        Apply back orthogonal transformations from the left.
!
!        Step (1L): apply back the Givens rotations performed.
!
         DO i = 1 , Givptr
            CALL DROT(Nrhs,B(Givcol(i,2),1),Ldb,B(Givcol(i,1),1),Ldb,   &
     &                Givnum(i,2),Givnum(i,1))
         ENDDO
!
!        Step (2L): permute rows of B.
!
         CALL DCOPY(Nrhs,B(nlp1,1),Ldb,Bx(1,1),Ldbx)
         DO i = 2 , n
            CALL DCOPY(Nrhs,B(Perm(i),1),Ldb,Bx(i,1),Ldbx)
         ENDDO
!
!        Step (3L): apply the inverse of the left singular vector
!        matrix to BX.
!
         IF ( K==1 ) THEN
            CALL DCOPY(Nrhs,Bx,Ldbx,B,Ldb)
            IF ( Z(1)<ZERO ) CALL DSCAL(Nrhs,NEGONE,B,Ldb)
         ELSE
            DO j = 1 , K
               diflj = Difl(j)
               dj = Poles(j,1)
               dsigj = -Poles(j,2)
               IF ( j<K ) THEN
                  difrj = -Difr(j,1)
                  dsigjp = -Poles(j+1,2)
               ENDIF
               IF ( (Z(j)==ZERO) .OR. (Poles(j,2)==ZERO) ) THEN
                  Work(j) = ZERO
               ELSE
                  Work(j) = -Poles(j,2)*Z(j)/diflj/(Poles(j,2)+dj)
               ENDIF
               DO i = 1 , j - 1
                  IF ( (Z(i)==ZERO) .OR. (Poles(i,2)==ZERO) ) THEN
                     Work(i) = ZERO
                  ELSE
                     Work(i) = Poles(i,2)*Z(i)                          &
     &                         /(DLAMC3(Poles(i,2),dsigj)-diflj)        &
     &                         /(Poles(i,2)+dj)
                  ENDIF
               ENDDO
               DO i = j + 1 , K
                  IF ( (Z(i)==ZERO) .OR. (Poles(i,2)==ZERO) ) THEN
                     Work(i) = ZERO
                  ELSE
                     Work(i) = Poles(i,2)*Z(i)                          &
     &                         /(DLAMC3(Poles(i,2),dsigjp)+difrj)       &
     &                         /(Poles(i,2)+dj)
                  ENDIF
               ENDDO
               Work(1) = NEGONE
               temp = DNRM2(K,Work,1)
               CALL DGEMV('T',K,Nrhs,ONE,Bx,Ldbx,Work,1,ZERO,B(j,1),Ldb)
               CALL DLASCL('G',0,0,temp,ONE,1,Nrhs,B(j,1),Ldb,Info)
            ENDDO
         ENDIF
!
!        Move the deflated rows of BX to B also.
!
         IF ( K<MAX(m,n) ) CALL DLACPY('A',n-K,Nrhs,Bx(K+1,1),Ldbx,     &
     &                                 B(K+1,1),Ldb)
      ELSE
!
!        Apply back the right orthogonal transformations.
!
!        Step (1R): apply back the new right singular vector matrix
!        to B.
!
         IF ( K==1 ) THEN
            CALL DCOPY(Nrhs,B,Ldb,Bx,Ldbx)
         ELSE
            DO j = 1 , K
               dsigj = Poles(j,2)
               IF ( Z(j)==ZERO ) THEN
                  Work(j) = ZERO
               ELSE
                  Work(j) = -Z(j)/Difl(j)/(dsigj+Poles(j,1))/Difr(j,2)
               ENDIF
               DO i = 1 , j - 1
                  IF ( Z(j)==ZERO ) THEN
                     Work(i) = ZERO
                  ELSE
                     Work(i) = Z(j)                                     &
     &                         /(DLAMC3(dsigj,-Poles(i+1,2))-Difr(i,1)) &
     &                         /(dsigj+Poles(i,1))/Difr(i,2)
                  ENDIF
               ENDDO
               DO i = j + 1 , K
                  IF ( Z(j)==ZERO ) THEN
                     Work(i) = ZERO
                  ELSE
                     Work(i) = Z(j)/(DLAMC3(dsigj,-Poles(i,2))-Difl(i)) &
     &                         /(dsigj+Poles(i,1))/Difr(i,2)
                  ENDIF
               ENDDO
               CALL DGEMV('T',K,Nrhs,ONE,B,Ldb,Work,1,ZERO,Bx(j,1),Ldbx)
            ENDDO
         ENDIF
!
!        Step (2R): if SQRE = 1, apply back the rotation that is
!        related to the right null space of the subproblem.
!
         IF ( Sqre==1 ) THEN
            CALL DCOPY(Nrhs,B(m,1),Ldb,Bx(m,1),Ldbx)
            CALL DROT(Nrhs,Bx(1,1),Ldbx,Bx(m,1),Ldbx,C,S)
         ENDIF
         IF ( K<MAX(m,n) ) CALL DLACPY('A',n-K,Nrhs,B(K+1,1),Ldb,       &
     &                                 Bx(K+1,1),Ldbx)
!
!        Step (3R): permute rows of B.
!
         CALL DCOPY(Nrhs,Bx(1,1),Ldbx,B(nlp1,1),Ldb)
         IF ( Sqre==1 ) CALL DCOPY(Nrhs,Bx(m,1),Ldbx,B(m,1),Ldb)
         DO i = 2 , n
            CALL DCOPY(Nrhs,Bx(i,1),Ldbx,B(Perm(i),1),Ldb)
         ENDDO
!
!        Step (4R): apply back the Givens rotations performed.
!
         DO i = Givptr , 1 , -1
            CALL DROT(Nrhs,B(Givcol(i,2),1),Ldb,B(Givcol(i,1),1),Ldb,   &
     &                Givnum(i,2),-Givnum(i,1))
         ENDDO
      ENDIF
!
!
!     End of DLALS0
!
      END SUBROUTINE DLALS0
