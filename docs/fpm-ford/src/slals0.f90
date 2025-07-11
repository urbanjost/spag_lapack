!*==slals0.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b SLALS0 applies back multiplying factors in solving the least squares problem using divide and conquer SVD approach. Used by sgelsd.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SLALS0 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slals0.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slals0.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slals0.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLALS0( ICOMPQ, NL, NR, SQRE, NRHS, B, LDB, BX, LDBX,
!                          PERM, GIVPTR, GIVCOL, LDGCOL, GIVNUM, LDGNUM,
!                          POLES, DIFL, DIFR, Z, K, C, S, WORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            GIVPTR, ICOMPQ, INFO, K, LDB, LDBX, LDGCOL,
!      $                   LDGNUM, NL, NR, NRHS, SQRE
!       REAL               C, S
!       ..
!       .. Array Arguments ..
!       INTEGER            GIVCOL( LDGCOL, * ), PERM( * )
!       REAL               B( LDB, * ), BX( LDBX, * ), DIFL( * ),
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
!> SLALS0 applies back the multiplying factors of either the left or the
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
!>          B is REAL array, dimension ( LDB, NRHS )
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
!>          BX is REAL array, dimension ( LDBX, NRHS )
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
!>          GIVNUM is REAL array, dimension ( LDGNUM, 2 )
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
!>          POLES is REAL array, dimension ( LDGNUM, 2 )
!>         On entry, POLES(1:K, 1) contains the new singular
!>         values obtained from solving the secular equation, and
!>         POLES(1:K, 2) is an array containing the poles in the secular
!>         equation.
!> \endverbatim
!>
!> \param[in] DIFL
!> \verbatim
!>          DIFL is REAL array, dimension ( K ).
!>         On entry, DIFL(I) is the distance between I-th updated
!>         (undeflated) singular value and the I-th (undeflated) old
!>         singular value.
!> \endverbatim
!>
!> \param[in] DIFR
!> \verbatim
!>          DIFR is REAL array, dimension ( LDGNUM, 2 ).
!>         On entry, DIFR(I, 1) contains the distances between I-th
!>         updated (undeflated) singular value and the I+1-th
!>         (undeflated) old singular value. And DIFR(I, 2) is the
!>         normalizing factor for the I-th right singular vector.
!> \endverbatim
!>
!> \param[in] Z
!> \verbatim
!>          Z is REAL array, dimension ( K )
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
!>          C is REAL
!>         C contains garbage if SQRE =0 and the C-value of a Givens
!>         rotation related to the right null space if SQRE = 1.
!> \endverbatim
!>
!> \param[in] S
!> \verbatim
!>          S is REAL
!>         S contains garbage if SQRE =0 and the S-value of a Givens
!>         rotation related to the right null space if SQRE = 1.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension ( K )
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
      SUBROUTINE SLALS0(Icompq,Nl,Nr,Sqre,Nrhs,B,Ldb,Bx,Ldbx,Perm,      &
     &                  Givptr,Givcol,Ldgcol,Givnum,Ldgnum,Poles,Difl,  &
     &                  Difr,Z,K,C,S,Work,Info)
      IMPLICIT NONE
!*--SLALS0272
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Givptr , Icompq , Info , K , Ldb , Ldbx , Ldgcol ,        &
     &        Ldgnum , Nl , Nr , Nrhs , Sqre
      REAL C , S
!     ..
!     .. Array Arguments ..
      INTEGER Givcol(Ldgcol,*) , Perm(*)
      REAL B(Ldb,*) , Bx(Ldbx,*) , Difl(*) , Difr(Ldgnum,*) ,           &
     &     Givnum(Ldgnum,*) , Poles(Ldgnum,*) , Work(*) , Z(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ONE , ZERO , NEGONE
      PARAMETER (ONE=1.0E0,ZERO=0.0E0,NEGONE=-1.0E0)
!     ..
!     .. Local Scalars ..
      INTEGER i , j , m , n , nlp1
      REAL diflj , difrj , dj , dsigj , dsigjp , temp
!     ..
!     .. External Subroutines ..
      EXTERNAL SCOPY , SGEMV , SLACPY , SLASCL , SROT , SSCAL , XERBLA
!     ..
!     .. External Functions ..
      REAL SLAMC3 , SNRM2
      EXTERNAL SLAMC3 , SNRM2
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX
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
         CALL XERBLA('SLALS0',-Info)
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
            CALL SROT(Nrhs,B(Givcol(i,2),1),Ldb,B(Givcol(i,1),1),Ldb,   &
     &                Givnum(i,2),Givnum(i,1))
         ENDDO
!
!        Step (2L): permute rows of B.
!
         CALL SCOPY(Nrhs,B(nlp1,1),Ldb,Bx(1,1),Ldbx)
         DO i = 2 , n
            CALL SCOPY(Nrhs,B(Perm(i),1),Ldb,Bx(i,1),Ldbx)
         ENDDO
!
!        Step (3L): apply the inverse of the left singular vector
!        matrix to BX.
!
         IF ( K==1 ) THEN
            CALL SCOPY(Nrhs,Bx,Ldbx,B,Ldb)
            IF ( Z(1)<ZERO ) CALL SSCAL(Nrhs,NEGONE,B,Ldb)
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
     &                         /(SLAMC3(Poles(i,2),dsigj)-diflj)        &
     &                         /(Poles(i,2)+dj)
                  ENDIF
               ENDDO
               DO i = j + 1 , K
                  IF ( (Z(i)==ZERO) .OR. (Poles(i,2)==ZERO) ) THEN
                     Work(i) = ZERO
                  ELSE
                     Work(i) = Poles(i,2)*Z(i)                          &
     &                         /(SLAMC3(Poles(i,2),dsigjp)+difrj)       &
     &                         /(Poles(i,2)+dj)
                  ENDIF
               ENDDO
               Work(1) = NEGONE
               temp = SNRM2(K,Work,1)
               CALL SGEMV('T',K,Nrhs,ONE,Bx,Ldbx,Work,1,ZERO,B(j,1),Ldb)
               CALL SLASCL('G',0,0,temp,ONE,1,Nrhs,B(j,1),Ldb,Info)
            ENDDO
         ENDIF
!
!        Move the deflated rows of BX to B also.
!
         IF ( K<MAX(m,n) ) CALL SLACPY('A',n-K,Nrhs,Bx(K+1,1),Ldbx,     &
     &                                 B(K+1,1),Ldb)
      ELSE
!
!        Apply back the right orthogonal transformations.
!
!        Step (1R): apply back the new right singular vector matrix
!        to B.
!
         IF ( K==1 ) THEN
            CALL SCOPY(Nrhs,B,Ldb,Bx,Ldbx)
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
     &                         /(SLAMC3(dsigj,-Poles(i+1,2))-Difr(i,1)) &
     &                         /(dsigj+Poles(i,1))/Difr(i,2)
                  ENDIF
               ENDDO
               DO i = j + 1 , K
                  IF ( Z(j)==ZERO ) THEN
                     Work(i) = ZERO
                  ELSE
                     Work(i) = Z(j)/(SLAMC3(dsigj,-Poles(i,2))-Difl(i)) &
     &                         /(dsigj+Poles(i,1))/Difr(i,2)
                  ENDIF
               ENDDO
               CALL SGEMV('T',K,Nrhs,ONE,B,Ldb,Work,1,ZERO,Bx(j,1),Ldbx)
            ENDDO
         ENDIF
!
!        Step (2R): if SQRE = 1, apply back the rotation that is
!        related to the right null space of the subproblem.
!
         IF ( Sqre==1 ) THEN
            CALL SCOPY(Nrhs,B(m,1),Ldb,Bx(m,1),Ldbx)
            CALL SROT(Nrhs,Bx(1,1),Ldbx,Bx(m,1),Ldbx,C,S)
         ENDIF
         IF ( K<MAX(m,n) ) CALL SLACPY('A',n-K,Nrhs,B(K+1,1),Ldb,       &
     &                                 Bx(K+1,1),Ldbx)
!
!        Step (3R): permute rows of B.
!
         CALL SCOPY(Nrhs,Bx(1,1),Ldbx,B(nlp1,1),Ldb)
         IF ( Sqre==1 ) CALL SCOPY(Nrhs,Bx(m,1),Ldbx,B(m,1),Ldb)
         DO i = 2 , n
            CALL SCOPY(Nrhs,Bx(i,1),Ldbx,B(Perm(i),1),Ldb)
         ENDDO
!
!        Step (4R): apply back the Givens rotations performed.
!
         DO i = Givptr , 1 , -1
            CALL SROT(Nrhs,B(Givcol(i,2),1),Ldb,B(Givcol(i,1),1),Ldb,   &
     &                Givnum(i,2),-Givnum(i,1))
         ENDDO
      ENDIF
!
!
!     End of SLALS0
!
      END SUBROUTINE SLALS0
