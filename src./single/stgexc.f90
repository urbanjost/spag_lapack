!*==stgexc.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b STGEXC
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download STGEXC + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/stgexc.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/stgexc.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/stgexc.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE STGEXC( WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z,
!                          LDZ, IFST, ILST, WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       LOGICAL            WANTQ, WANTZ
!       INTEGER            IFST, ILST, INFO, LDA, LDB, LDQ, LDZ, LWORK, N
!       ..
!       .. Array Arguments ..
!       REAL               A( LDA, * ), B( LDB, * ), Q( LDQ, * ),
!      $                   WORK( * ), Z( LDZ, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> STGEXC reorders the generalized real Schur decomposition of a real
!> matrix pair (A,B) using an orthogonal equivalence transformation
!>
!>                (A, B) = Q * (A, B) * Z**T,
!>
!> so that the diagonal block of (A, B) with row index IFST is moved
!> to row ILST.
!>
!> (A, B) must be in generalized real Schur canonical form (as returned
!> by SGGES), i.e. A is block upper triangular with 1-by-1 and 2-by-2
!> diagonal blocks. B is upper triangular.
!>
!> Optionally, the matrices Q and Z of generalized Schur vectors are
!> updated.
!>
!>        Q(in) * A(in) * Z(in)**T = Q(out) * A(out) * Z(out)**T
!>        Q(in) * B(in) * Z(in)**T = Q(out) * B(out) * Z(out)**T
!>
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] WANTQ
!> \verbatim
!>          WANTQ is LOGICAL
!>          .TRUE. : update the left transformation matrix Q;
!>          .FALSE.: do not update Q.
!> \endverbatim
!>
!> \param[in] WANTZ
!> \verbatim
!>          WANTZ is LOGICAL
!>          .TRUE. : update the right transformation matrix Z;
!>          .FALSE.: do not update Z.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrices A and B. N >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is REAL array, dimension (LDA,N)
!>          On entry, the matrix A in generalized real Schur canonical
!>          form.
!>          On exit, the updated matrix A, again in generalized
!>          real Schur canonical form.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A. LDA >= max(1,N).
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is REAL array, dimension (LDB,N)
!>          On entry, the matrix B in generalized real Schur canonical
!>          form (A,B).
!>          On exit, the updated matrix B, again in generalized
!>          real Schur canonical form (A,B).
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B. LDB >= max(1,N).
!> \endverbatim
!>
!> \param[in,out] Q
!> \verbatim
!>          Q is REAL array, dimension (LDQ,N)
!>          On entry, if WANTQ = .TRUE., the orthogonal matrix Q.
!>          On exit, the updated matrix Q.
!>          If WANTQ = .FALSE., Q is not referenced.
!> \endverbatim
!>
!> \param[in] LDQ
!> \verbatim
!>          LDQ is INTEGER
!>          The leading dimension of the array Q. LDQ >= 1.
!>          If WANTQ = .TRUE., LDQ >= N.
!> \endverbatim
!>
!> \param[in,out] Z
!> \verbatim
!>          Z is REAL array, dimension (LDZ,N)
!>          On entry, if WANTZ = .TRUE., the orthogonal matrix Z.
!>          On exit, the updated matrix Z.
!>          If WANTZ = .FALSE., Z is not referenced.
!> \endverbatim
!>
!> \param[in] LDZ
!> \verbatim
!>          LDZ is INTEGER
!>          The leading dimension of the array Z. LDZ >= 1.
!>          If WANTZ = .TRUE., LDZ >= N.
!> \endverbatim
!>
!> \param[in,out] IFST
!> \verbatim
!>          IFST is INTEGER
!> \endverbatim
!>
!> \param[in,out] ILST
!> \verbatim
!>          ILST is INTEGER
!>          Specify the reordering of the diagonal blocks of (A, B).
!>          The block with row index IFST is moved to row ILST, by a
!>          sequence of swapping between adjacent blocks.
!>          On exit, if IFST pointed on entry to the second row of
!>          a 2-by-2 block, it is changed to point to the first row;
!>          ILST always points to the first row of the block in its
!>          final position (which may differ from its input value by
!>          +1 or -1). 1 <= IFST, ILST <= N.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK.
!>          LWORK >= 1 when N <= 1, otherwise LWORK >= 4*N + 16.
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the WORK array, returns
!>          this value as the first entry of the WORK array, and no error
!>          message related to LWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>           =0:  successful exit.
!>           <0:  if INFO = -i, the i-th argument had an illegal value.
!>           =1:  The transformed matrix pair (A, B) would be too far
!>                from generalized Schur form; the problem is ill-
!>                conditioned. (A, B) may have been partially reordered,
!>                and ILST points to the first row of the current
!>                position of the block being moved.
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
!> \ingroup realGEcomputational
!
!> \par Contributors:
!  ==================
!>
!>     Bo Kagstrom and Peter Poromaa, Department of Computing Science,
!>     Umea University, S-901 87 Umea, Sweden.
!
!> \par References:
!  ================
!>
!> \verbatim
!>
!>  [1] B. Kagstrom; A Direct Method for Reordering Eigenvalues in the
!>      Generalized Real Schur Form of a Regular Matrix Pair (A, B), in
!>      M.S. Moonen et al (eds), Linear Algebra for Large Scale and
!>      Real-Time Applications, Kluwer Academic Publ. 1993, pp 195-218.
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE STGEXC(Wantq,Wantz,N,A,Lda,B,Ldb,Q,Ldq,Z,Ldz,Ifst,Ilst,&
     &                  Work,Lwork,Info)
      USE S_STGEX2
      USE S_XERBLA
      IMPLICIT NONE
!*--STGEXC226
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0
!
! Dummy argument declarations rewritten by SPAG
!
      LOGICAL :: Wantq
      LOGICAL :: Wantz
      INTEGER :: N
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      REAL , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      INTEGER , INTENT(INOUT) :: Ifst
      INTEGER , INTENT(INOUT) :: Ilst
      REAL , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: here , lwmin , nbf , nbl , nbnext
      LOGICAL :: lquery
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
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
!     Decode and test input arguments.
!
      Info = 0
      lquery = (Lwork==-1)
      IF ( N<0 ) THEN
         Info = -3
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -5
      ELSEIF ( Ldb<MAX(1,N) ) THEN
         Info = -7
      ELSEIF ( Ldq<1 .OR. Wantq .AND. (Ldq<MAX(1,N)) ) THEN
         Info = -9
      ELSEIF ( Ldz<1 .OR. Wantz .AND. (Ldz<MAX(1,N)) ) THEN
         Info = -11
      ELSEIF ( Ifst<1 .OR. Ifst>N ) THEN
         Info = -12
      ELSEIF ( Ilst<1 .OR. Ilst>N ) THEN
         Info = -13
      ENDIF
!
      IF ( Info==0 ) THEN
         IF ( N<=1 ) THEN
            lwmin = 1
         ELSE
            lwmin = 4*N + 16
         ENDIF
         Work(1) = lwmin
!
         IF ( Lwork<lwmin .AND. .NOT.lquery ) Info = -15
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('STGEXC',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N<=1 ) RETURN
!
!     Determine the first row of the specified block and find out
!     if it is 1-by-1 or 2-by-2.
!
      IF ( Ifst>1 ) THEN
         IF ( A(Ifst,Ifst-1)/=ZERO ) Ifst = Ifst - 1
      ENDIF
      nbf = 1
      IF ( Ifst<N ) THEN
         IF ( A(Ifst+1,Ifst)/=ZERO ) nbf = 2
      ENDIF
!
!     Determine the first row of the final block
!     and find out if it is 1-by-1 or 2-by-2.
!
      IF ( Ilst>1 ) THEN
         IF ( A(Ilst,Ilst-1)/=ZERO ) Ilst = Ilst - 1
      ENDIF
      nbl = 1
      IF ( Ilst<N ) THEN
         IF ( A(Ilst+1,Ilst)/=ZERO ) nbl = 2
      ENDIF
      IF ( Ifst==Ilst ) RETURN
!
      IF ( Ifst<Ilst ) THEN
!
!        Update ILST.
!
         IF ( nbf==2 .AND. nbl==1 ) Ilst = Ilst - 1
         IF ( nbf==1 .AND. nbl==2 ) Ilst = Ilst + 1
!
         here = Ifst
         DO
!
!
!        Swap with next one below.
!
            IF ( nbf==1 .OR. nbf==2 ) THEN
!
!           Current block either 1-by-1 or 2-by-2.
!
               nbnext = 1
               IF ( here+nbf+1<=N ) THEN
                  IF ( A(here+nbf+1,here+nbf)/=ZERO ) nbnext = 2
               ENDIF
               CALL STGEX2(Wantq,Wantz,N,A,Lda,B,Ldb,Q,Ldq,Z,Ldz,here,  &
     &                     nbf,nbnext,Work,Lwork,Info)
               IF ( Info/=0 ) THEN
                  Ilst = here
                  RETURN
               ENDIF
               here = here + nbnext
!
!           Test if 2-by-2 block breaks into two 1-by-1 blocks.
!
               IF ( nbf==2 ) THEN
                  IF ( A(here+1,here)==ZERO ) nbf = 3
               ENDIF
!
            ELSE
!
!           Current block consists of two 1-by-1 blocks, each of which
!           must be swapped individually.
!
               nbnext = 1
               IF ( here+3<=N ) THEN
                  IF ( A(here+3,here+2)/=ZERO ) nbnext = 2
               ENDIF
               CALL STGEX2(Wantq,Wantz,N,A,Lda,B,Ldb,Q,Ldq,Z,Ldz,here+1,&
     &                     1,nbnext,Work,Lwork,Info)
               IF ( Info/=0 ) THEN
                  Ilst = here
                  RETURN
               ENDIF
               IF ( nbnext==1 ) THEN
!
!              Swap two 1-by-1 blocks.
!
                  CALL STGEX2(Wantq,Wantz,N,A,Lda,B,Ldb,Q,Ldq,Z,Ldz,    &
     &                        here,1,1,Work,Lwork,Info)
                  IF ( Info/=0 ) THEN
                     Ilst = here
                     RETURN
                  ENDIF
                  here = here + 1
!
               ELSE
!
!              Recompute NBNEXT in case of 2-by-2 split.
!
                  IF ( A(here+2,here+1)==ZERO ) nbnext = 1
                  IF ( nbnext==2 ) THEN
!
!                 2-by-2 block did not split.
!
                     CALL STGEX2(Wantq,Wantz,N,A,Lda,B,Ldb,Q,Ldq,Z,Ldz, &
     &                           here,1,nbnext,Work,Lwork,Info)
                     IF ( Info/=0 ) THEN
                        Ilst = here
                        RETURN
                     ENDIF
                     here = here + 2
                  ELSE
!
!                 2-by-2 block did split.
!
                     CALL STGEX2(Wantq,Wantz,N,A,Lda,B,Ldb,Q,Ldq,Z,Ldz, &
     &                           here,1,1,Work,Lwork,Info)
                     IF ( Info/=0 ) THEN
                        Ilst = here
                        RETURN
                     ENDIF
                     here = here + 1
                     CALL STGEX2(Wantq,Wantz,N,A,Lda,B,Ldb,Q,Ldq,Z,Ldz, &
     &                           here,1,1,Work,Lwork,Info)
                     IF ( Info/=0 ) THEN
                        Ilst = here
                        RETURN
                     ENDIF
                     here = here + 1
                  ENDIF
!
               ENDIF
            ENDIF
            IF ( here>=Ilst ) EXIT
         ENDDO
      ELSE
         here = Ifst
         DO
!
!
!        Swap with next one below.
!
            IF ( nbf==1 .OR. nbf==2 ) THEN
!
!           Current block either 1-by-1 or 2-by-2.
!
               nbnext = 1
               IF ( here>=3 ) THEN
                  IF ( A(here-1,here-2)/=ZERO ) nbnext = 2
               ENDIF
               CALL STGEX2(Wantq,Wantz,N,A,Lda,B,Ldb,Q,Ldq,Z,Ldz,       &
     &                     here-nbnext,nbnext,nbf,Work,Lwork,Info)
               IF ( Info/=0 ) THEN
                  Ilst = here
                  RETURN
               ENDIF
               here = here - nbnext
!
!           Test if 2-by-2 block breaks into two 1-by-1 blocks.
!
               IF ( nbf==2 ) THEN
                  IF ( A(here+1,here)==ZERO ) nbf = 3
               ENDIF
!
            ELSE
!
!           Current block consists of two 1-by-1 blocks, each of which
!           must be swapped individually.
!
               nbnext = 1
               IF ( here>=3 ) THEN
                  IF ( A(here-1,here-2)/=ZERO ) nbnext = 2
               ENDIF
               CALL STGEX2(Wantq,Wantz,N,A,Lda,B,Ldb,Q,Ldq,Z,Ldz,       &
     &                     here-nbnext,nbnext,1,Work,Lwork,Info)
               IF ( Info/=0 ) THEN
                  Ilst = here
                  RETURN
               ENDIF
               IF ( nbnext==1 ) THEN
!
!              Swap two 1-by-1 blocks.
!
                  CALL STGEX2(Wantq,Wantz,N,A,Lda,B,Ldb,Q,Ldq,Z,Ldz,    &
     &                        here,nbnext,1,Work,Lwork,Info)
                  IF ( Info/=0 ) THEN
                     Ilst = here
                     RETURN
                  ENDIF
                  here = here - 1
               ELSE
!
!             Recompute NBNEXT in case of 2-by-2 split.
!
                  IF ( A(here,here-1)==ZERO ) nbnext = 1
                  IF ( nbnext==2 ) THEN
!
!                 2-by-2 block did not split.
!
                     CALL STGEX2(Wantq,Wantz,N,A,Lda,B,Ldb,Q,Ldq,Z,Ldz, &
     &                           here-1,2,1,Work,Lwork,Info)
                     IF ( Info/=0 ) THEN
                        Ilst = here
                        RETURN
                     ENDIF
                     here = here - 2
                  ELSE
!
!                 2-by-2 block did split.
!
                     CALL STGEX2(Wantq,Wantz,N,A,Lda,B,Ldb,Q,Ldq,Z,Ldz, &
     &                           here,1,1,Work,Lwork,Info)
                     IF ( Info/=0 ) THEN
                        Ilst = here
                        RETURN
                     ENDIF
                     here = here - 1
                     CALL STGEX2(Wantq,Wantz,N,A,Lda,B,Ldb,Q,Ldq,Z,Ldz, &
     &                           here,1,1,Work,Lwork,Info)
                     IF ( Info/=0 ) THEN
                        Ilst = here
                        RETURN
                     ENDIF
                     here = here - 1
                  ENDIF
               ENDIF
            ENDIF
            IF ( here<=Ilst ) EXIT
         ENDDO
      ENDIF
      Ilst = here
      Work(1) = lwmin
!
!     End of STGEXC
!
      END SUBROUTINE STGEXC
