!*==ztrsen.f90  processed by SPAG 7.51RB at 20:09 on  3 Mar 2022
!> \brief \b ZTRSEN
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZTRSEN + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztrsen.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztrsen.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztrsen.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZTRSEN( JOB, COMPQ, SELECT, N, T, LDT, Q, LDQ, W, M, S,
!                          SEP, WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          COMPQ, JOB
!       INTEGER            INFO, LDQ, LDT, LWORK, M, N
!       DOUBLE PRECISION   S, SEP
!       ..
!       .. Array Arguments ..
!       LOGICAL            SELECT( * )
!       COMPLEX*16         Q( LDQ, * ), T( LDT, * ), W( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZTRSEN reorders the Schur factorization of a complex matrix
!> A = Q*T*Q**H, so that a selected cluster of eigenvalues appears in
!> the leading positions on the diagonal of the upper triangular matrix
!> T, and the leading columns of Q form an orthonormal basis of the
!> corresponding right invariant subspace.
!>
!> Optionally the routine computes the reciprocal condition numbers of
!> the cluster of eigenvalues and/or the invariant subspace.
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
!>          select the j-th eigenvalue, SELECT(j) must be set to .TRUE..
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
!>          T is COMPLEX*16 array, dimension (LDT,N)
!>          On entry, the upper triangular matrix T.
!>          On exit, T is overwritten by the reordered matrix T, with the
!>          selected eigenvalues as the leading diagonal elements.
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
!>          Q is COMPLEX*16 array, dimension (LDQ,N)
!>          On entry, if COMPQ = 'V', the matrix Q of Schur vectors.
!>          On exit, if COMPQ = 'V', Q has been postmultiplied by the
!>          unitary transformation matrix which reorders T; the leading M
!>          columns of Q form an orthonormal basis for the specified
!>          invariant subspace.
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
!> \param[out] W
!> \verbatim
!>          W is COMPLEX*16 array, dimension (N)
!>          The reordered eigenvalues of T, in the same order as they
!>          appear on the diagonal of T.
!> \endverbatim
!>
!> \param[out] M
!> \verbatim
!>          M is INTEGER
!>          The dimension of the specified invariant subspace.
!>          0 <= M <= N.
!> \endverbatim
!>
!> \param[out] S
!> \verbatim
!>          S is DOUBLE PRECISION
!>          If JOB = 'E' or 'B', S is a lower bound on the reciprocal
!>          condition number for the selected cluster of eigenvalues.
!>          S cannot underestimate the true reciprocal condition number
!>          by more than a factor of sqrt(N). If M = 0 or N, S = 1.
!>          If JOB = 'N' or 'V', S is not referenced.
!> \endverbatim
!>
!> \param[out] SEP
!> \verbatim
!>          SEP is DOUBLE PRECISION
!>          If JOB = 'V' or 'B', SEP is the estimated reciprocal
!>          condition number of the specified invariant subspace. If
!>          M = 0 or N, SEP = norm(T).
!>          If JOB = 'N' or 'E', SEP is not referenced.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK.
!>          If JOB = 'N', LWORK >= 1;
!>          if JOB = 'E', LWORK = max(1,M*(N-M));
!>          if JOB = 'V' or 'B', LWORK >= max(1,2*M*(N-M)).
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
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
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
!> \ingroup complex16OTHERcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  ZTRSEN first collects the selected eigenvalues by computing a unitary
!>  transformation Z to move them to the top left corner of T. In other
!>  words, the selected eigenvalues are the eigenvalues of T11 in:
!>
!>          Z**H * T * Z = ( T11 T12 ) n1
!>                         (  0  T22 ) n2
!>                            n1  n2
!>
!>  where N = n1+n2. The first
!>  n1 columns of Z span the specified invariant subspace of T.
!>
!>  If T has been obtained from the Schur factorization of a matrix
!>  A = Q*T*Q**H, then the reordered Schur factorization of A is given by
!>  A = (Q*Z)*(Z**H*T*Z)*(Q*Z)**H, and the first n1 columns of Q*Z span the
!>  corresponding invariant subspace of A.
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
      SUBROUTINE ZTRSEN(Job,Compq,Select,N,T,Ldt,Q,Ldq,W,M,S,Sep,Work,  &
     &                  Lwork,Info)
      IMPLICIT NONE
!*--ZTRSEN268
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Compq , Job
      INTEGER Info , Ldq , Ldt , Lwork , M , N
      DOUBLE PRECISION S , Sep
!     ..
!     .. Array Arguments ..
      LOGICAL Select(*)
      COMPLEX*16 Q(Ldq,*) , T(Ldt,*) , W(*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D+0,ONE=1.0D+0)
!     ..
!     .. Local Scalars ..
      LOGICAL lquery , wantbh , wantq , wants , wantsp
      INTEGER ierr , k , kase , ks , lwmin , n1 , n2 , nn
      DOUBLE PRECISION est , rnorm , scale
!     ..
!     .. Local Arrays ..
      INTEGER isave(3)
      DOUBLE PRECISION rwork(1)
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      DOUBLE PRECISION ZLANGE
      EXTERNAL LSAME , ZLANGE
!     ..
!     .. External Subroutines ..
      EXTERNAL XERBLA , ZLACN2 , ZLACPY , ZTREXC , ZTRSYL
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX , SQRT
!     ..
!     .. Executable Statements ..
!
!     Decode and test the input parameters.
!
      wantbh = LSAME(Job,'B')
      wants = LSAME(Job,'E') .OR. wantbh
      wantsp = LSAME(Job,'V') .OR. wantbh
      wantq = LSAME(Compq,'V')
!
!     Set M to the number of selected eigenvalues.
!
      M = 0
      DO k = 1 , N
         IF ( Select(k) ) M = M + 1
      ENDDO
!
      n1 = M
      n2 = N - M
      nn = n1*n2
!
      Info = 0
      lquery = (Lwork==-1)
!
      IF ( wantsp ) THEN
         lwmin = MAX(1,2*nn)
      ELSEIF ( LSAME(Job,'N') ) THEN
         lwmin = 1
      ELSEIF ( LSAME(Job,'E') ) THEN
         lwmin = MAX(1,nn)
      ENDIF
!
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
      ELSEIF ( Lwork<lwmin .AND. .NOT.lquery ) THEN
         Info = -14
      ENDIF
!
      IF ( Info==0 ) Work(1) = lwmin
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('ZTRSEN',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( M==N .OR. M==0 ) THEN
         IF ( wants ) S = ONE
         IF ( wantsp ) Sep = ZLANGE('1',N,N,T,Ldt,rwork)
         GOTO 100
      ENDIF
!
!     Collect the selected eigenvalues at the top left corner of T.
!
      ks = 0
      DO k = 1 , N
         IF ( Select(k) ) THEN
            ks = ks + 1
!
!           Swap the K-th eigenvalue to position KS.
!
            IF ( k/=ks ) CALL ZTREXC(Compq,N,T,Ldt,Q,Ldq,k,ks,ierr)
         ENDIF
      ENDDO
!
      IF ( wants ) THEN
!
!        Solve the Sylvester equation for R:
!
!           T11*R - R*T22 = scale*T12
!
         CALL ZLACPY('F',n1,n2,T(1,n1+1),Ldt,Work,n1)
         CALL ZTRSYL('N','N',-1,n1,n2,T,Ldt,T(n1+1,n1+1),Ldt,Work,n1,   &
     &               scale,ierr)
!
!        Estimate the reciprocal of the condition number of the cluster
!        of eigenvalues.
!
         rnorm = ZLANGE('F',n1,n2,Work,n1,rwork)
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
            CALL ZLACN2(nn,Work(nn+1),Work,est,kase,isave)
            IF ( kase/=0 ) THEN
               IF ( kase==1 ) THEN
!
!              Solve T11*R - R*T22 = scale*X.
!
                  CALL ZTRSYL('N','N',-1,n1,n2,T,Ldt,T(n1+1,n1+1),Ldt,  &
     &                        Work,n1,scale,ierr)
               ELSE
!
!              Solve T11**H*R - R*T22**H = scale*X.
!
                  CALL ZTRSYL('C','C',-1,n1,n2,T,Ldt,T(n1+1,n1+1),Ldt,  &
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
!     Copy reordered eigenvalues to W.
!
 100  DO k = 1 , N
         W(k) = T(k,k)
      ENDDO
!
      Work(1) = lwmin
!
!
!     End of ZTRSEN
!
      END SUBROUTINE ZTRSEN
