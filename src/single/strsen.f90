!*==strsen.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b STRSEN
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download STRSEN + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/strsen.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/strsen.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/strsen.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE STRSEN( JOB, COMPQ, SELECT, N, T, LDT, Q, LDQ, WR, WI,
!                          M, S, SEP, WORK, LWORK, IWORK, LIWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          COMPQ, JOB
!       INTEGER            INFO, LDQ, LDT, LIWORK, LWORK, M, N
!       REAL               S, SEP
!       ..
!       .. Array Arguments ..
!       LOGICAL            SELECT( * )
!       INTEGER            IWORK( * )
!       REAL               Q( LDQ, * ), T( LDT, * ), WI( * ), WORK( * ),
!      $                   WR( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> STRSEN reorders the real Schur factorization of a real matrix
!> A = Q*T*Q**T, so that a selected cluster of eigenvalues appears in
!> the leading diagonal blocks of the upper quasi-triangular matrix T,
!> and the leading columns of Q form an orthonormal basis of the
!> corresponding right invariant subspace.
!>
!> Optionally the routine computes the reciprocal condition numbers of
!> the cluster of eigenvalues and/or the invariant subspace.
!>
!> T must be in Schur canonical form (as returned by SHSEQR), that is,
!> block upper triangular with 1-by-1 and 2-by-2 diagonal blocks; each
!> 2-by-2 diagonal block has its diagonal elements equal and its
!> off-diagonal elements of opposite sign.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] JOB
!> \verbatim
!>          JOB is CHARACTER*1
!>          Specifies whether condition numbers are required for the
!>          cluster of eigenvalues (S) or the invariant subspace (SEP):
!>          = 'N': none;
!>          = 'E': for eigenvalues only (S);
!>          = 'V': for invariant subspace only (SEP);
!>          = 'B': for both eigenvalues and invariant subspace (S and
!>                 SEP).
!> \endverbatim
!>
!> \param[in] COMPQ
!> \verbatim
!>          COMPQ is CHARACTER*1
!>          = 'V': update the matrix Q of Schur vectors;
!>          = 'N': do not update Q.
!> \endverbatim
!>
!> \param[in] SELECT
!> \verbatim
!>          SELECT is LOGICAL array, dimension (N)
!>          SELECT specifies the eigenvalues in the selected cluster. To
!>          select a real eigenvalue w(j), SELECT(j) must be set to
!>          .TRUE.. To select a complex conjugate pair of eigenvalues
!>          w(j) and w(j+1), corresponding to a 2-by-2 diagonal block,
!>          either SELECT(j) or SELECT(j+1) or both must be set to
!>          .TRUE.; a complex conjugate pair of eigenvalues must be
!>          either both included in the cluster or both excluded.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix T. N >= 0.
!> \endverbatim
!>
!> \param[in,out] T
!> \verbatim
!>          T is REAL array, dimension (LDT,N)
!>          On entry, the upper quasi-triangular matrix T, in Schur
!>          canonical form.
!>          On exit, T is overwritten by the reordered matrix T, again in
!>          Schur canonical form, with the selected eigenvalues in the
!>          leading diagonal blocks.
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
!>          Q is REAL array, dimension (LDQ,N)
!>          On entry, if COMPQ = 'V', the matrix Q of Schur vectors.
!>          On exit, if COMPQ = 'V', Q has been postmultiplied by the
!>          orthogonal transformation matrix which reorders T; the
!>          leading M columns of Q form an orthonormal basis for the
!>          specified invariant subspace.
!>          If COMPQ = 'N', Q is not referenced.
!> \endverbatim
!>
!> \param[in] LDQ
!> \verbatim
!>          LDQ is INTEGER
!>          The leading dimension of the array Q.
!>          LDQ >= 1; and if COMPQ = 'V', LDQ >= N.
!> \endverbatim
!>
!> \param[out] WR
!> \verbatim
!>          WR is REAL array, dimension (N)
!> \endverbatim
!>
!> \param[out] WI
!> \verbatim
!>          WI is REAL array, dimension (N)
!>
!>          The real and imaginary parts, respectively, of the reordered
!>          eigenvalues of T. The eigenvalues are stored in the same
!>          order as on the diagonal of T, with WR(i) = T(i,i) and, if
!>          T(i:i+1,i:i+1) is a 2-by-2 diagonal block, WI(i) > 0 and
!>          WI(i+1) = -WI(i). Note that if a complex eigenvalue is
!>          sufficiently ill-conditioned, then its value may differ
!>          significantly from its value before reordering.
!> \endverbatim
!>
!> \param[out] M
!> \verbatim
!>          M is INTEGER
!>          The dimension of the specified invariant subspace.
!>          0 < = M <= N.
!> \endverbatim
!>
!> \param[out] S
!> \verbatim
!>          S is REAL
!>          If JOB = 'E' or 'B', S is a lower bound on the reciprocal
!>          condition number for the selected cluster of eigenvalues.
!>          S cannot underestimate the true reciprocal condition number
!>          by more than a factor of sqrt(N). If M = 0 or N, S = 1.
!>          If JOB = 'N' or 'V', S is not referenced.
!> \endverbatim
!>
!> \param[out] SEP
!> \verbatim
!>          SEP is REAL
!>          If JOB = 'V' or 'B', SEP is the estimated reciprocal
!>          condition number of the specified invariant subspace. If
!>          M = 0 or N, SEP = norm(T).
!>          If JOB = 'N' or 'E', SEP is not referenced.
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
!>          If JOB = 'N', LWORK >= max(1,N);
!>          if JOB = 'E', LWORK >= max(1,M*(N-M));
!>          if JOB = 'V' or 'B', LWORK >= max(1,2*M*(N-M)).
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the WORK array, returns
!>          this value as the first entry of the WORK array, and no error
!>          message related to LWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (MAX(1,LIWORK))
!>          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.
!> \endverbatim
!>
!> \param[in] LIWORK
!> \verbatim
!>          LIWORK is INTEGER
!>          The dimension of the array IWORK.
!>          If JOB = 'N' or 'E', LIWORK >= 1;
!>          if JOB = 'V' or 'B', LIWORK >= max(1,M*(N-M)).
!>
!>          If LIWORK = -1, then a workspace query is assumed; the
!>          routine only calculates the optimal size of the IWORK array,
!>          returns this value as the first entry of the IWORK array, and
!>          no error message related to LIWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: successful exit
!>          < 0: if INFO = -i, the i-th argument had an illegal value
!>          = 1: reordering of T failed because some eigenvalues are too
!>               close to separate (the problem is very ill-conditioned);
!>               T may have been partially reordered, and WR and WI
!>               contain the eigenvalues in the same order as in T; S and
!>               SEP (if requested) are set to zero.
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
!> \date April 2012
!
!> \ingroup realOTHERcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  STRSEN first collects the selected eigenvalues by computing an
!>  orthogonal transformation Z to move them to the top left corner of T.
!>  In other words, the selected eigenvalues are the eigenvalues of T11
!>  in:
!>
!>          Z**T * T * Z = ( T11 T12 ) n1
!>                         (  0  T22 ) n2
!>                            n1  n2
!>
!>  where N = n1+n2 and Z**T means the transpose of Z. The first n1 columns
!>  of Z span the specified invariant subspace of T.
!>
!>  If T has been obtained from the real Schur factorization of a matrix
!>  A = Q*T*Q**T, then the reordered real Schur factorization of A is given
!>  by A = (Q*Z)*(Z**T*T*Z)*(Q*Z)**T, and the first n1 columns of Q*Z span
!>  the corresponding invariant subspace of A.
!>
!>  The reciprocal condition number of the average of the eigenvalues of
!>  T11 may be returned in S. S lies between 0 (very badly conditioned)
!>  and 1 (very well conditioned). It is computed as follows. First we
!>  compute R so that
!>
!>                         P = ( I  R ) n1
!>                             ( 0  0 ) n2
!>                               n1 n2
!>
!>  is the projector on the invariant subspace associated with T11.
!>  R is the solution of the Sylvester equation:
!>
!>                        T11*R - R*T22 = T12.
!>
!>  Let F-norm(M) denote the Frobenius-norm of M and 2-norm(M) denote
!>  the two-norm of M. Then S is computed as the lower bound
!>
!>                      (1 + F-norm(R)**2)**(-1/2)
!>
!>  on the reciprocal of 2-norm(P), the true reciprocal condition number.
!>  S cannot underestimate 1 / 2-norm(P) by more than a factor of
!>  sqrt(N).
!>
!>  An approximate error bound for the computed average of the
!>  eigenvalues of T11 is
!>
!>                         EPS * norm(T) / S
!>
!>  where EPS is the machine precision.
!>
!>  The reciprocal condition number of the right invariant subspace
!>  spanned by the first n1 columns of Z (or of Q*Z) is returned in SEP.
!>  SEP is defined as the separation of T11 and T22:
!>
!>                     sep( T11, T22 ) = sigma-min( C )
!>
!>  where sigma-min(C) is the smallest singular value of the
!>  n1*n2-by-n1*n2 matrix
!>
!>     C  = kprod( I(n2), T11 ) - kprod( transpose(T22), I(n1) )
!>
!>  I(m) is an m by m identity matrix, and kprod denotes the Kronecker
!>  product. We estimate sigma-min(C) by the reciprocal of an estimate of
!>  the 1-norm of inverse(C). The true reciprocal 1-norm of inverse(C)
!>  cannot differ from sigma-min(C) by more than a factor of sqrt(n1*n2).
!>
!>  When SEP is small, small changes in T can cause large changes in
!>  the invariant subspace. An approximate bound on the maximum angular
!>  error in the computed right invariant subspace is
!>
!>                      EPS * norm(T) / SEP
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE STRSEN(Job,Compq,Select,N,T,Ldt,Q,Ldq,Wr,Wi,M,S,Sep,   &
     &                  Work,Lwork,Iwork,Liwork,Info)
      IMPLICIT NONE
!*--STRSEN318
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     April 2012
!
!     .. Scalar Arguments ..
      CHARACTER Compq , Job
      INTEGER Info , Ldq , Ldt , Liwork , Lwork , M , N
      REAL S , Sep
!     ..
!     .. Array Arguments ..
      LOGICAL Select(*)
      INTEGER Iwork(*)
      REAL Q(Ldq,*) , T(Ldt,*) , Wi(*) , Work(*) , Wr(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE
      PARAMETER (ZERO=0.0E+0,ONE=1.0E+0)
!     ..
!     .. Local Scalars ..
      LOGICAL lquery , pair , swap , wantbh , wantq , wants , wantsp
      INTEGER ierr , k , kase , kk , ks , liwmin , lwmin , n1 , n2 , nn
      REAL est , rnorm , scale
!     ..
!     .. Local Arrays ..
      INTEGER isave(3)
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      REAL SLANGE
      EXTERNAL LSAME , SLANGE
!     ..
!     .. External Subroutines ..
      EXTERNAL SLACN2 , SLACPY , STREXC , STRSYL , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX , SQRT
!     ..
!     .. Executable Statements ..
!
!     Decode and test the input parameters
!
      wantbh = LSAME(Job,'B')
      wants = LSAME(Job,'E') .OR. wantbh
      wantsp = LSAME(Job,'V') .OR. wantbh
      wantq = LSAME(Compq,'V')
!
      Info = 0
      lquery = (Lwork==-1)
      IF ( .NOT.LSAME(Job,'N') .AND. .NOT.wants .AND. .NOT.wantsp ) THEN
         Info = -1
      ELSEIF ( .NOT.LSAME(Compq,'N') .AND. .NOT.wantq ) THEN
         Info = -2
      ELSEIF ( N<0 ) THEN
         Info = -4
      ELSEIF ( Ldt<MAX(1,N) ) THEN
         Info = -6
      ELSEIF ( Ldq<1 .OR. (wantq .AND. Ldq<N) ) THEN
         Info = -8
      ELSE
!
!        Set M to the dimension of the specified invariant subspace,
!        and test LWORK and LIWORK.
!
         M = 0
         pair = .FALSE.
         DO k = 1 , N
            IF ( pair ) THEN
               pair = .FALSE.
            ELSEIF ( k<N ) THEN
               IF ( T(k+1,k)==ZERO ) THEN
                  IF ( Select(k) ) M = M + 1
               ELSE
                  pair = .TRUE.
                  IF ( Select(k) .OR. Select(k+1) ) M = M + 2
               ENDIF
            ELSE
               IF ( Select(N) ) M = M + 1
            ENDIF
         ENDDO
!
         n1 = M
         n2 = N - M
         nn = n1*n2
!
         IF ( wantsp ) THEN
            lwmin = MAX(1,2*nn)
            liwmin = MAX(1,nn)
         ELSEIF ( LSAME(Job,'N') ) THEN
            lwmin = MAX(1,N)
            liwmin = 1
         ELSEIF ( LSAME(Job,'E') ) THEN
            lwmin = MAX(1,nn)
            liwmin = 1
         ENDIF
!
         IF ( Lwork<lwmin .AND. .NOT.lquery ) THEN
            Info = -15
         ELSEIF ( Liwork<liwmin .AND. .NOT.lquery ) THEN
            Info = -17
         ENDIF
      ENDIF
!
      IF ( Info==0 ) THEN
         Work(1) = lwmin
         Iwork(1) = liwmin
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('STRSEN',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
         RETURN
      ENDIF
!
!     Quick return if possible.
!
      IF ( M==N .OR. M==0 ) THEN
         IF ( wants ) S = ONE
         IF ( wantsp ) Sep = SLANGE('1',N,N,T,Ldt,Work)
         GOTO 100
      ENDIF
!
!     Collect the selected blocks at the top-left corner of T.
!
      ks = 0
      pair = .FALSE.
      DO k = 1 , N
         IF ( pair ) THEN
            pair = .FALSE.
         ELSE
            swap = Select(k)
            IF ( k<N ) THEN
               IF ( T(k+1,k)/=ZERO ) THEN
                  pair = .TRUE.
                  swap = swap .OR. Select(k+1)
               ENDIF
            ENDIF
            IF ( swap ) THEN
               ks = ks + 1
!
!              Swap the K-th block to position KS.
!
               ierr = 0
               kk = k
               IF ( k/=ks ) CALL STREXC(Compq,N,T,Ldt,Q,Ldq,kk,ks,Work, &
     &                                  ierr)
               IF ( ierr==1 .OR. ierr==2 ) THEN
!
!                 Blocks too close to swap: exit.
!
                  Info = 1
                  IF ( wants ) S = ZERO
                  IF ( wantsp ) Sep = ZERO
                  GOTO 100
               ENDIF
               IF ( pair ) ks = ks + 1
            ENDIF
         ENDIF
      ENDDO
!
      IF ( wants ) THEN
!
!        Solve Sylvester equation for R:
!
!           T11*R - R*T22 = scale*T12
!
         CALL SLACPY('F',n1,n2,T(1,n1+1),Ldt,Work,n1)
         CALL STRSYL('N','N',-1,n1,n2,T,Ldt,T(n1+1,n1+1),Ldt,Work,n1,   &
     &               scale,ierr)
!
!        Estimate the reciprocal of the condition number of the cluster
!        of eigenvalues.
!
         rnorm = SLANGE('F',n1,n2,Work,n1,Work)
         IF ( rnorm==ZERO ) THEN
            S = ONE
         ELSE
            S = scale/(SQRT(scale*scale/rnorm+rnorm)*SQRT(rnorm))
         ENDIF
      ENDIF
!
      IF ( wantsp ) THEN
!
!        Estimate sep(T11,T22).
!
         est = ZERO
         kase = 0
         DO
            CALL SLACN2(nn,Work(nn+1),Work,Iwork,est,kase,isave)
            IF ( kase/=0 ) THEN
               IF ( kase==1 ) THEN
!
!              Solve  T11*R - R*T22 = scale*X.
!
                  CALL STRSYL('N','N',-1,n1,n2,T,Ldt,T(n1+1,n1+1),Ldt,  &
     &                        Work,n1,scale,ierr)
               ELSE
!
!              Solve T11**T*R - R*T22**T = scale*X.
!
                  CALL STRSYL('T','T',-1,n1,n2,T,Ldt,T(n1+1,n1+1),Ldt,  &
     &                        Work,n1,scale,ierr)
               ENDIF
               CYCLE
            ENDIF
!
            Sep = scale/est
            EXIT
         ENDDO
      ENDIF
!
!
!     Store the output eigenvalues in WR and WI.
!
 100  DO k = 1 , N
         Wr(k) = T(k,k)
         Wi(k) = ZERO
      ENDDO
      DO k = 1 , N - 1
         IF ( T(k+1,k)/=ZERO ) THEN
            Wi(k) = SQRT(ABS(T(k,k+1)))*SQRT(ABS(T(k+1,k)))
            Wi(k+1) = -Wi(k)
         ENDIF
      ENDDO
!
      Work(1) = lwmin
      Iwork(1) = liwmin
!
!
!     End of STRSEN
!
      END SUBROUTINE STRSEN
