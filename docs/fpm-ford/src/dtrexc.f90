!*==dtrexc.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b DTREXC
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DTREXC + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dtrexc.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dtrexc.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dtrexc.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DTREXC( COMPQ, N, T, LDT, Q, LDQ, IFST, ILST, WORK,
!                          INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          COMPQ
!       INTEGER            IFST, ILST, INFO, LDQ, LDT, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   Q( LDQ, * ), T( LDT, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DTREXC reorders the real Schur factorization of a real matrix
!> A = Q*T*Q**T, so that the diagonal block of T with row index IFST is
!> moved to row ILST.
!>
!> The real Schur form T is reordered by an orthogonal similarity
!> transformation Z**T*T*Z, and optionally the matrix Q of Schur vectors
!> is updated by postmultiplying it with Z.
!>
!> T must be in Schur canonical form (as returned by DHSEQR), that is,
!> block upper triangular with 1-by-1 and 2-by-2 diagonal blocks; each
!> 2-by-2 diagonal block has its diagonal elements equal and its
!> off-diagonal elements of opposite sign.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] COMPQ
!> \verbatim
!>          COMPQ is CHARACTER*1
!>          = 'V':  update the matrix Q of Schur vectors;
!>          = 'N':  do not update Q.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix T. N >= 0.
!>          If N == 0 arguments ILST and IFST may be any value.
!> \endverbatim
!>
!> \param[in,out] T
!> \verbatim
!>          T is DOUBLE PRECISION array, dimension (LDT,N)
!>          On entry, the upper quasi-triangular matrix T, in Schur
!>          Schur canonical form.
!>          On exit, the reordered upper quasi-triangular matrix, again
!>          in Schur canonical form.
!> \endverbatim
!>
!> \param[in] LDT
!> \verbatim
!>          LDT is INTEGER
!>          The leading dimension of the array T. LDT >= max(1,N).
!> \endverbatim
!>
!> \param[in,out] Q
!> \verbatim
!>          Q is DOUBLE PRECISION array, dimension (LDQ,N)
!>          On entry, if COMPQ = 'V', the matrix Q of Schur vectors.
!>          On exit, if COMPQ = 'V', Q has been postmultiplied by the
!>          orthogonal transformation matrix Z which reorders T.
!>          If COMPQ = 'N', Q is not referenced.
!> \endverbatim
!>
!> \param[in] LDQ
!> \verbatim
!>          LDQ is INTEGER
!>          The leading dimension of the array Q.  LDQ >= 1, and if
!>          COMPQ = 'V', LDQ >= max(1,N).
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
!>
!>          Specify the reordering of the diagonal blocks of T.
!>          The block with row index IFST is moved to row ILST, by a
!>          sequence of transpositions between adjacent blocks.
!>          On exit, if IFST pointed on entry to the second row of a
!>          2-by-2 block, it is changed to point to the first row; ILST
!>          always points to the first row of the block in its final
!>          position (which may differ from its input value by +1 or -1).
!>          1 <= IFST <= N; 1 <= ILST <= N.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (N)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!>          = 1:  two adjacent blocks were too close to swap (the problem
!>                is very ill-conditioned); T may have been partially
!>                reordered, and ILST points to the first row of the
!>                current position of the block being moved.
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
!  =====================================================================
      SUBROUTINE DTREXC(Compq,N,T,Ldt,Q,Ldq,Ifst,Ilst,Work,Info)
      IMPLICIT NONE
!*--DTREXC151
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Compq
      INTEGER Ifst , Ilst , Info , Ldq , Ldt , N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION Q(Ldq,*) , T(Ldt,*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D+0)
!     ..
!     .. Local Scalars ..
      LOGICAL wantq
      INTEGER here , nbf , nbl , nbnext
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL DLAEXC , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX
!     ..
!     .. Executable Statements ..
!
!     Decode and test the input arguments.
!
      Info = 0
      wantq = LSAME(Compq,'V')
      IF ( .NOT.wantq .AND. .NOT.LSAME(Compq,'N') ) THEN
         Info = -1
      ELSEIF ( N<0 ) THEN
         Info = -2
      ELSEIF ( Ldt<MAX(1,N) ) THEN
         Info = -4
      ELSEIF ( Ldq<1 .OR. (wantq .AND. Ldq<MAX(1,N)) ) THEN
         Info = -6
      ELSEIF ( (Ifst<1 .OR. Ifst>N) .AND. (N>0) ) THEN
         Info = -7
      ELSEIF ( (Ilst<1 .OR. Ilst>N) .AND. (N>0) ) THEN
         Info = -8
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('DTREXC',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N<=1 ) RETURN
!
!     Determine the first row of specified block
!     and find out it is 1 by 1 or 2 by 2.
!
      IF ( Ifst>1 ) THEN
         IF ( T(Ifst,Ifst-1)/=ZERO ) Ifst = Ifst - 1
      ENDIF
      nbf = 1
      IF ( Ifst<N ) THEN
         IF ( T(Ifst+1,Ifst)/=ZERO ) nbf = 2
      ENDIF
!
!     Determine the first row of the final block
!     and find out it is 1 by 1 or 2 by 2.
!
      IF ( Ilst>1 ) THEN
         IF ( T(Ilst,Ilst-1)/=ZERO ) Ilst = Ilst - 1
      ENDIF
      nbl = 1
      IF ( Ilst<N ) THEN
         IF ( T(Ilst+1,Ilst)/=ZERO ) nbl = 2
      ENDIF
!
      IF ( Ifst==Ilst ) RETURN
!
      IF ( Ifst<Ilst ) THEN
!
!        Update ILST
!
         IF ( nbf==2 .AND. nbl==1 ) Ilst = Ilst - 1
         IF ( nbf==1 .AND. nbl==2 ) Ilst = Ilst + 1
!
         here = Ifst
         DO
!
!
!        Swap block with next one below
!
            IF ( nbf==1 .OR. nbf==2 ) THEN
!
!           Current block either 1 by 1 or 2 by 2
!
               nbnext = 1
               IF ( here+nbf+1<=N ) THEN
                  IF ( T(here+nbf+1,here+nbf)/=ZERO ) nbnext = 2
               ENDIF
               CALL DLAEXC(wantq,N,T,Ldt,Q,Ldq,here,nbf,nbnext,Work,    &
     &                     Info)
               IF ( Info/=0 ) THEN
                  Ilst = here
                  RETURN
               ENDIF
               here = here + nbnext
!
!           Test if 2 by 2 block breaks into two 1 by 1 blocks
!
               IF ( nbf==2 ) THEN
                  IF ( T(here+1,here)==ZERO ) nbf = 3
               ENDIF
!
            ELSE
!
!           Current block consists of two 1 by 1 blocks each of which
!           must be swapped individually
!
               nbnext = 1
               IF ( here+3<=N ) THEN
                  IF ( T(here+3,here+2)/=ZERO ) nbnext = 2
               ENDIF
               CALL DLAEXC(wantq,N,T,Ldt,Q,Ldq,here+1,1,nbnext,Work,    &
     &                     Info)
               IF ( Info/=0 ) THEN
                  Ilst = here
                  RETURN
               ENDIF
               IF ( nbnext==1 ) THEN
!
!              Swap two 1 by 1 blocks, no problems possible
!
                  CALL DLAEXC(wantq,N,T,Ldt,Q,Ldq,here,1,nbnext,Work,   &
     &                        Info)
                  here = here + 1
               ELSE
!
!              Recompute NBNEXT in case 2 by 2 split
!
                  IF ( T(here+2,here+1)==ZERO ) nbnext = 1
                  IF ( nbnext==2 ) THEN
!
!                 2 by 2 Block did not split
!
                     CALL DLAEXC(wantq,N,T,Ldt,Q,Ldq,here,1,nbnext,Work,&
     &                           Info)
                     IF ( Info/=0 ) THEN
                        Ilst = here
                        RETURN
                     ENDIF
                     here = here + 2
                  ELSE
!
!                 2 by 2 Block did split
!
                     CALL DLAEXC(wantq,N,T,Ldt,Q,Ldq,here,1,1,Work,Info)
                     CALL DLAEXC(wantq,N,T,Ldt,Q,Ldq,here+1,1,1,Work,   &
     &                           Info)
                     here = here + 2
                  ENDIF
               ENDIF
            ENDIF
            IF ( here>=Ilst ) THEN
               Ilst = here
               EXIT
            ENDIF
         ENDDO
!
      ELSE
!
         here = Ifst
         DO
!
!        Swap block with next one above
!
            IF ( nbf==1 .OR. nbf==2 ) THEN
!
!           Current block either 1 by 1 or 2 by 2
!
               nbnext = 1
               IF ( here>=3 ) THEN
                  IF ( T(here-1,here-2)/=ZERO ) nbnext = 2
               ENDIF
               CALL DLAEXC(wantq,N,T,Ldt,Q,Ldq,here-nbnext,nbnext,nbf,  &
     &                     Work,Info)
               IF ( Info/=0 ) THEN
                  Ilst = here
                  RETURN
               ENDIF
               here = here - nbnext
!
!           Test if 2 by 2 block breaks into two 1 by 1 blocks
!
               IF ( nbf==2 ) THEN
                  IF ( T(here+1,here)==ZERO ) nbf = 3
               ENDIF
!
            ELSE
!
!           Current block consists of two 1 by 1 blocks each of which
!           must be swapped individually
!
               nbnext = 1
               IF ( here>=3 ) THEN
                  IF ( T(here-1,here-2)/=ZERO ) nbnext = 2
               ENDIF
               CALL DLAEXC(wantq,N,T,Ldt,Q,Ldq,here-nbnext,nbnext,1,    &
     &                     Work,Info)
               IF ( Info/=0 ) THEN
                  Ilst = here
                  RETURN
               ENDIF
               IF ( nbnext==1 ) THEN
!
!              Swap two 1 by 1 blocks, no problems possible
!
                  CALL DLAEXC(wantq,N,T,Ldt,Q,Ldq,here,nbnext,1,Work,   &
     &                        Info)
                  here = here - 1
               ELSE
!
!              Recompute NBNEXT in case 2 by 2 split
!
                  IF ( T(here,here-1)==ZERO ) nbnext = 1
                  IF ( nbnext==2 ) THEN
!
!                 2 by 2 Block did not split
!
                     CALL DLAEXC(wantq,N,T,Ldt,Q,Ldq,here-1,2,1,Work,   &
     &                           Info)
                     IF ( Info/=0 ) THEN
                        Ilst = here
                        RETURN
                     ENDIF
                     here = here - 2
                  ELSE
!
!                 2 by 2 Block did split
!
                     CALL DLAEXC(wantq,N,T,Ldt,Q,Ldq,here,1,1,Work,Info)
                     CALL DLAEXC(wantq,N,T,Ldt,Q,Ldq,here-1,1,1,Work,   &
     &                           Info)
                     here = here - 2
                  ENDIF
               ENDIF
            ENDIF
            IF ( here<=Ilst ) THEN
               Ilst = here
               EXIT
            ENDIF
         ENDDO
      ENDIF
!
!
!     End of DTREXC
!
      END SUBROUTINE DTREXC
