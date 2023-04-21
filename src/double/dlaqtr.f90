!*==dlaqtr.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b DLAQTR solves a real quasi-triangular system of equations, or a complex quasi-triangular system of special form, in real arithmetic.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLAQTR + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaqtr.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaqtr.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaqtr.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLAQTR( LTRAN, LREAL, N, T, LDT, B, W, SCALE, X, WORK,
!                          INFO )
!
!       .. Scalar Arguments ..
!       LOGICAL            LREAL, LTRAN
!       INTEGER            INFO, LDT, N
!       DOUBLE PRECISION   SCALE, W
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   B( * ), T( LDT, * ), WORK( * ), X( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLAQTR solves the real quasi-triangular system
!>
!>              op(T)*p = scale*c,               if LREAL = .TRUE.
!>
!> or the complex quasi-triangular systems
!>
!>            op(T + iB)*(p+iq) = scale*(c+id),  if LREAL = .FALSE.
!>
!> in real arithmetic, where T is upper quasi-triangular.
!> If LREAL = .FALSE., then the first diagonal block of T must be
!> 1 by 1, B is the specially structured matrix
!>
!>                B = [ b(1) b(2) ... b(n) ]
!>                    [       w            ]
!>                    [           w        ]
!>                    [              .     ]
!>                    [                 w  ]
!>
!> op(A) = A or A**T, A**T denotes the transpose of
!> matrix A.
!>
!> On input, X = [ c ].  On output, X = [ p ].
!>               [ d ]                  [ q ]
!>
!> This subroutine is designed for the condition number estimation
!> in routine DTRSNA.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] LTRAN
!> \verbatim
!>          LTRAN is LOGICAL
!>          On entry, LTRAN specifies the option of conjugate transpose:
!>             = .FALSE.,    op(T+i*B) = T+i*B,
!>             = .TRUE.,     op(T+i*B) = (T+i*B)**T.
!> \endverbatim
!>
!> \param[in] LREAL
!> \verbatim
!>          LREAL is LOGICAL
!>          On entry, LREAL specifies the input matrix structure:
!>             = .FALSE.,    the input is complex
!>             = .TRUE.,     the input is real
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          On entry, N specifies the order of T+i*B. N >= 0.
!> \endverbatim
!>
!> \param[in] T
!> \verbatim
!>          T is DOUBLE PRECISION array, dimension (LDT,N)
!>          On entry, T contains a matrix in Schur canonical form.
!>          If LREAL = .FALSE., then the first diagonal block of T mu
!>          be 1 by 1.
!> \endverbatim
!>
!> \param[in] LDT
!> \verbatim
!>          LDT is INTEGER
!>          The leading dimension of the matrix T. LDT >= max(1,N).
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is DOUBLE PRECISION array, dimension (N)
!>          On entry, B contains the elements to form the matrix
!>          B as described above.
!>          If LREAL = .TRUE., B is not referenced.
!> \endverbatim
!>
!> \param[in] W
!> \verbatim
!>          W is DOUBLE PRECISION
!>          On entry, W is the diagonal element of the matrix B.
!>          If LREAL = .TRUE., W is not referenced.
!> \endverbatim
!>
!> \param[out] SCALE
!> \verbatim
!>          SCALE is DOUBLE PRECISION
!>          On exit, SCALE is the scale factor.
!> \endverbatim
!>
!> \param[in,out] X
!> \verbatim
!>          X is DOUBLE PRECISION array, dimension (2*N)
!>          On entry, X contains the right hand side of the system.
!>          On exit, X is overwritten by the solution.
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
!>          On exit, INFO is set to
!>             0: successful exit.
!>               1: the some diagonal 1 by 1 block has been perturbed by
!>                  a small number SMIN to keep nonsingularity.
!>               2: the some diagonal 2 by 2 block has been perturbed by
!>                  a small number in DLALN2 to keep nonsingularity.
!>          NOTE: In the interests of speed, this routine does not
!>                check the inputs for errors.
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
!> \ingroup doubleOTHERauxiliary
!
!  =====================================================================
      SUBROUTINE DLAQTR(Ltran,Lreal,N,T,Ldt,B,W,Scale,X,Work,Info)
      IMPLICIT NONE
!*--DLAQTR168
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      LOGICAL Lreal , Ltran
      INTEGER Info , Ldt , N
      DOUBLE PRECISION Scale , W
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION B(*) , T(Ldt,*) , Work(*) , X(*)
!     ..
!
! =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D+0,ONE=1.0D+0)
!     ..
!     .. Local Scalars ..
      LOGICAL notran
      INTEGER i , ierr , j , j1 , j2 , jnext , k , n1 , n2
      DOUBLE PRECISION bignum , eps , rec , scaloc , si , smin , sminw ,&
     &                 smlnum , sr , tjj , tmp , xj , xmax , xnorm , z
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION d(2,2) , v(2,2)
!     ..
!     .. External Functions ..
      INTEGER IDAMAX
      DOUBLE PRECISION DASUM , DDOT , DLAMCH , DLANGE
      EXTERNAL IDAMAX , DASUM , DDOT , DLAMCH , DLANGE
!     ..
!     .. External Subroutines ..
      EXTERNAL DAXPY , DLADIV , DLALN2 , DSCAL
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX
!     ..
!     .. Executable Statements ..
!
!     Do not test the input parameters for errors
!
      notran = .NOT.Ltran
      Info = 0
!
!     Quick return if possible
!
      IF ( N==0 ) RETURN
!
!     Set constants to control overflow
!
      eps = DLAMCH('P')
      smlnum = DLAMCH('S')/eps
      bignum = ONE/smlnum
!
      xnorm = DLANGE('M',N,N,T,Ldt,d)
      IF ( .NOT.Lreal ) xnorm = MAX(xnorm,ABS(W),DLANGE('M',N,1,B,N,d))
      smin = MAX(smlnum,eps*xnorm)
!
!     Compute 1-norm of each column of strictly upper triangular
!     part of T to control overflow in triangular solver.
!
      Work(1) = ZERO
      DO j = 2 , N
         Work(j) = DASUM(j-1,T(1,j),1)
      ENDDO
!
      IF ( .NOT.Lreal ) THEN
         DO i = 2 , N
            Work(i) = Work(i) + ABS(B(i))
         ENDDO
      ENDIF
!
      n2 = 2*N
      n1 = N
      IF ( .NOT.Lreal ) n1 = n2
      k = IDAMAX(n1,X,1)
      xmax = ABS(X(k))
      Scale = ONE
!
      IF ( xmax>bignum ) THEN
         Scale = bignum/xmax
         CALL DSCAL(n1,Scale,X,1)
         xmax = bignum
      ENDIF
!
      IF ( .NOT.(Lreal) ) THEN
!
         sminw = MAX(eps*ABS(W),smin)
         IF ( notran ) THEN
!
!           Solve (T + iB)*(p+iq) = c+id
!
            jnext = N
            DO j = N , 1 , -1
               IF ( j<=jnext ) THEN
                  j1 = j
                  j2 = j
                  jnext = j - 1
                  IF ( j>1 ) THEN
                     IF ( T(j,j-1)/=ZERO ) THEN
                        j1 = j - 1
                        jnext = j - 2
                     ENDIF
                  ENDIF
!
                  IF ( j1==j2 ) THEN
!
!                 1 by 1 diagonal block
!
!                 Scale if necessary to avoid overflow in division
!
                     z = W
                     IF ( j1==1 ) z = B(1)
                     xj = ABS(X(j1)) + ABS(X(N+j1))
                     tjj = ABS(T(j1,j1)) + ABS(z)
                     tmp = T(j1,j1)
                     IF ( tjj<sminw ) THEN
                        tmp = sminw
                        tjj = sminw
                        Info = 1
                     ENDIF
!
                     IF ( xj/=ZERO ) THEN
!
                        IF ( tjj<ONE ) THEN
                           IF ( xj>bignum*tjj ) THEN
                              rec = ONE/xj
                              CALL DSCAL(n2,rec,X,1)
                              Scale = Scale*rec
                              xmax = xmax*rec
                           ENDIF
                        ENDIF
                        CALL DLADIV(X(j1),X(N+j1),tmp,z,sr,si)
                        X(j1) = sr
                        X(N+j1) = si
                        xj = ABS(X(j1)) + ABS(X(N+j1))
!
!                 Scale x if necessary to avoid overflow when adding a
!                 multiple of column j1 of T.
!
                        IF ( xj>ONE ) THEN
                           rec = ONE/xj
                           IF ( Work(j1)>(bignum-xmax)*rec ) THEN
                              CALL DSCAL(n2,rec,X,1)
                              Scale = Scale*rec
                           ENDIF
                        ENDIF
!
                        IF ( j1>1 ) THEN
                           CALL DAXPY(j1-1,-X(j1),T(1,j1),1,X,1)
                           CALL DAXPY(j1-1,-X(N+j1),T(1,j1),1,X(N+1),1)
!
                           X(1) = X(1) + B(j1)*X(N+j1)
                           X(N+1) = X(N+1) - B(j1)*X(j1)
!
                           xmax = ZERO
                           DO k = 1 , j1 - 1
                              xmax = MAX(xmax,ABS(X(k))+ABS(X(k+N)))
                           ENDDO
                        ENDIF
                     ENDIF
!
                  ELSE
!
!                 Meet 2 by 2 diagonal block
!
                     d(1,1) = X(j1)
                     d(2,1) = X(j2)
                     d(1,2) = X(N+j1)
                     d(2,2) = X(N+j2)
                     CALL DLALN2(.FALSE.,2,2,sminw,ONE,T(j1,j1),Ldt,ONE,&
     &                           ONE,d,2,ZERO,-W,v,2,scaloc,xnorm,ierr)
                     IF ( ierr/=0 ) Info = 2
!
                     IF ( scaloc/=ONE ) THEN
                        CALL DSCAL(2*N,scaloc,X,1)
                        Scale = scaloc*Scale
                     ENDIF
                     X(j1) = v(1,1)
                     X(j2) = v(2,1)
                     X(N+j1) = v(1,2)
                     X(N+j2) = v(2,2)
!
!                 Scale X(J1), .... to avoid overflow in
!                 updating right hand side.
!
                     xj = MAX(ABS(v(1,1))+ABS(v(1,2)),ABS(v(2,1))       &
     &                    +ABS(v(2,2)))
                     IF ( xj>ONE ) THEN
                        rec = ONE/xj
                        IF ( MAX(Work(j1),Work(j2))>(bignum-xmax)*rec ) &
     &                       THEN
                           CALL DSCAL(n2,rec,X,1)
                           Scale = Scale*rec
                        ENDIF
                     ENDIF
!
!                 Update the right-hand side.
!
                     IF ( j1>1 ) THEN
                        CALL DAXPY(j1-1,-X(j1),T(1,j1),1,X,1)
                        CALL DAXPY(j1-1,-X(j2),T(1,j2),1,X,1)
!
                        CALL DAXPY(j1-1,-X(N+j1),T(1,j1),1,X(N+1),1)
                        CALL DAXPY(j1-1,-X(N+j2),T(1,j2),1,X(N+1),1)
!
                        X(1) = X(1) + B(j1)*X(N+j1) + B(j2)*X(N+j2)
                        X(N+1) = X(N+1) - B(j1)*X(j1) - B(j2)*X(j2)
!
                        xmax = ZERO
                        DO k = 1 , j1 - 1
                           xmax = MAX(ABS(X(k))+ABS(X(k+N)),xmax)
                        ENDDO
                     ENDIF
!
                  ENDIF
               ENDIF
            ENDDO
!
         ELSE
!
!           Solve (T + iB)**T*(p+iq) = c+id
!
            jnext = 1
            DO j = 1 , N
               IF ( j>=jnext ) THEN
                  j1 = j
                  j2 = j
                  jnext = j + 1
                  IF ( j<N ) THEN
                     IF ( T(j+1,j)/=ZERO ) THEN
                        j2 = j + 1
                        jnext = j + 2
                     ENDIF
                  ENDIF
!
                  IF ( j1==j2 ) THEN
!
!                 1 by 1 diagonal block
!
!                 Scale if necessary to avoid overflow in forming the
!                 right-hand side element by inner product.
!
                     xj = ABS(X(j1)) + ABS(X(j1+N))
                     IF ( xmax>ONE ) THEN
                        rec = ONE/xmax
                        IF ( Work(j1)>(bignum-xj)*rec ) THEN
                           CALL DSCAL(n2,rec,X,1)
                           Scale = Scale*rec
                           xmax = xmax*rec
                        ENDIF
                     ENDIF
!
                     X(j1) = X(j1) - DDOT(j1-1,T(1,j1),1,X,1)
                     X(N+j1) = X(N+j1) - DDOT(j1-1,T(1,j1),1,X(N+1),1)
                     IF ( j1>1 ) THEN
                        X(j1) = X(j1) - B(j1)*X(N+1)
                        X(N+j1) = X(N+j1) + B(j1)*X(1)
                     ENDIF
                     xj = ABS(X(j1)) + ABS(X(j1+N))
!
                     z = W
                     IF ( j1==1 ) z = B(1)
!
!                 Scale if necessary to avoid overflow in
!                 complex division
!
                     tjj = ABS(T(j1,j1)) + ABS(z)
                     tmp = T(j1,j1)
                     IF ( tjj<sminw ) THEN
                        tmp = sminw
                        tjj = sminw
                        Info = 1
                     ENDIF
!
                     IF ( tjj<ONE ) THEN
                        IF ( xj>bignum*tjj ) THEN
                           rec = ONE/xj
                           CALL DSCAL(n2,rec,X,1)
                           Scale = Scale*rec
                           xmax = xmax*rec
                        ENDIF
                     ENDIF
                     CALL DLADIV(X(j1),X(N+j1),tmp,-z,sr,si)
                     X(j1) = sr
                     X(j1+N) = si
                     xmax = MAX(ABS(X(j1))+ABS(X(j1+N)),xmax)
!
                  ELSE
!
!                 2 by 2 diagonal block
!
!                 Scale if necessary to avoid overflow in forming the
!                 right-hand side element by inner product.
!
                     xj = MAX(ABS(X(j1))+ABS(X(N+j1)),ABS(X(j2))        &
     &                    +ABS(X(N+j2)))
                     IF ( xmax>ONE ) THEN
                        rec = ONE/xmax
                        IF ( MAX(Work(j1),Work(j2))>(bignum-xj)/xmax )  &
     &                       THEN
                           CALL DSCAL(n2,rec,X,1)
                           Scale = Scale*rec
                           xmax = xmax*rec
                        ENDIF
                     ENDIF
!
                     d(1,1) = X(j1) - DDOT(j1-1,T(1,j1),1,X,1)
                     d(2,1) = X(j2) - DDOT(j1-1,T(1,j2),1,X,1)
                     d(1,2) = X(N+j1) - DDOT(j1-1,T(1,j1),1,X(N+1),1)
                     d(2,2) = X(N+j2) - DDOT(j1-1,T(1,j2),1,X(N+1),1)
                     d(1,1) = d(1,1) - B(j1)*X(N+1)
                     d(2,1) = d(2,1) - B(j2)*X(N+1)
                     d(1,2) = d(1,2) + B(j1)*X(1)
                     d(2,2) = d(2,2) + B(j2)*X(1)
!
                     CALL DLALN2(.TRUE.,2,2,sminw,ONE,T(j1,j1),Ldt,ONE, &
     &                           ONE,d,2,ZERO,W,v,2,scaloc,xnorm,ierr)
                     IF ( ierr/=0 ) Info = 2
!
                     IF ( scaloc/=ONE ) THEN
                        CALL DSCAL(n2,scaloc,X,1)
                        Scale = scaloc*Scale
                     ENDIF
                     X(j1) = v(1,1)
                     X(j2) = v(2,1)
                     X(N+j1) = v(1,2)
                     X(N+j2) = v(2,2)
                     xmax = MAX(ABS(X(j1))+ABS(X(N+j1)),ABS(X(j2))      &
     &                      +ABS(X(N+j2)),xmax)
!
                  ENDIF
               ENDIF
!
            ENDDO
!
         ENDIF
!
      ELSEIF ( notran ) THEN
!
!           Solve T*p = scale*c
!
         jnext = N
         DO j = N , 1 , -1
            IF ( j<=jnext ) THEN
               j1 = j
               j2 = j
               jnext = j - 1
               IF ( j>1 ) THEN
                  IF ( T(j,j-1)/=ZERO ) THEN
                     j1 = j - 1
                     jnext = j - 2
                  ENDIF
               ENDIF
!
               IF ( j1==j2 ) THEN
!
!                 Meet 1 by 1 diagonal block
!
!                 Scale to avoid overflow when computing
!                     x(j) = b(j)/T(j,j)
!
                  xj = ABS(X(j1))
                  tjj = ABS(T(j1,j1))
                  tmp = T(j1,j1)
                  IF ( tjj<smin ) THEN
                     tmp = smin
                     tjj = smin
                     Info = 1
                  ENDIF
!
                  IF ( xj/=ZERO ) THEN
!
                     IF ( tjj<ONE ) THEN
                        IF ( xj>bignum*tjj ) THEN
                           rec = ONE/xj
                           CALL DSCAL(N,rec,X,1)
                           Scale = Scale*rec
                           xmax = xmax*rec
                        ENDIF
                     ENDIF
                     X(j1) = X(j1)/tmp
                     xj = ABS(X(j1))
!
!                 Scale x if necessary to avoid overflow when adding a
!                 multiple of column j1 of T.
!
                     IF ( xj>ONE ) THEN
                        rec = ONE/xj
                        IF ( Work(j1)>(bignum-xmax)*rec ) THEN
                           CALL DSCAL(N,rec,X,1)
                           Scale = Scale*rec
                        ENDIF
                     ENDIF
                     IF ( j1>1 ) THEN
                        CALL DAXPY(j1-1,-X(j1),T(1,j1),1,X,1)
                        k = IDAMAX(j1-1,X,1)
                        xmax = ABS(X(k))
                     ENDIF
                  ENDIF
!
               ELSE
!
!                 Meet 2 by 2 diagonal block
!
!                 Call 2 by 2 linear system solve, to take
!                 care of possible overflow by scaling factor.
!
                  d(1,1) = X(j1)
                  d(2,1) = X(j2)
                  CALL DLALN2(.FALSE.,2,1,smin,ONE,T(j1,j1),Ldt,ONE,ONE,&
     &                        d,2,ZERO,ZERO,v,2,scaloc,xnorm,ierr)
                  IF ( ierr/=0 ) Info = 2
!
                  IF ( scaloc/=ONE ) THEN
                     CALL DSCAL(N,scaloc,X,1)
                     Scale = Scale*scaloc
                  ENDIF
                  X(j1) = v(1,1)
                  X(j2) = v(2,1)
!
!                 Scale V(1,1) (= X(J1)) and/or V(2,1) (=X(J2))
!                 to avoid overflow in updating right-hand side.
!
                  xj = MAX(ABS(v(1,1)),ABS(v(2,1)))
                  IF ( xj>ONE ) THEN
                     rec = ONE/xj
                     IF ( MAX(Work(j1),Work(j2))>(bignum-xmax)*rec )    &
     &                    THEN
                        CALL DSCAL(N,rec,X,1)
                        Scale = Scale*rec
                     ENDIF
                  ENDIF
!
!                 Update right-hand side
!
                  IF ( j1>1 ) THEN
                     CALL DAXPY(j1-1,-X(j1),T(1,j1),1,X,1)
                     CALL DAXPY(j1-1,-X(j2),T(1,j2),1,X,1)
                     k = IDAMAX(j1-1,X,1)
                     xmax = ABS(X(k))
                  ENDIF
!
               ENDIF
            ENDIF
!
         ENDDO
!
      ELSE
!
!           Solve T**T*p = scale*c
!
         jnext = 1
         DO j = 1 , N
            IF ( j>=jnext ) THEN
               j1 = j
               j2 = j
               jnext = j + 1
               IF ( j<N ) THEN
                  IF ( T(j+1,j)/=ZERO ) THEN
                     j2 = j + 1
                     jnext = j + 2
                  ENDIF
               ENDIF
!
               IF ( j1==j2 ) THEN
!
!                 1 by 1 diagonal block
!
!                 Scale if necessary to avoid overflow in forming the
!                 right-hand side element by inner product.
!
                  xj = ABS(X(j1))
                  IF ( xmax>ONE ) THEN
                     rec = ONE/xmax
                     IF ( Work(j1)>(bignum-xj)*rec ) THEN
                        CALL DSCAL(N,rec,X,1)
                        Scale = Scale*rec
                        xmax = xmax*rec
                     ENDIF
                  ENDIF
!
                  X(j1) = X(j1) - DDOT(j1-1,T(1,j1),1,X,1)
!
                  xj = ABS(X(j1))
                  tjj = ABS(T(j1,j1))
                  tmp = T(j1,j1)
                  IF ( tjj<smin ) THEN
                     tmp = smin
                     tjj = smin
                     Info = 1
                  ENDIF
!
                  IF ( tjj<ONE ) THEN
                     IF ( xj>bignum*tjj ) THEN
                        rec = ONE/xj
                        CALL DSCAL(N,rec,X,1)
                        Scale = Scale*rec
                        xmax = xmax*rec
                     ENDIF
                  ENDIF
                  X(j1) = X(j1)/tmp
                  xmax = MAX(xmax,ABS(X(j1)))
!
               ELSE
!
!                 2 by 2 diagonal block
!
!                 Scale if necessary to avoid overflow in forming the
!                 right-hand side elements by inner product.
!
                  xj = MAX(ABS(X(j1)),ABS(X(j2)))
                  IF ( xmax>ONE ) THEN
                     rec = ONE/xmax
                     IF ( MAX(Work(j2),Work(j1))>(bignum-xj)*rec ) THEN
                        CALL DSCAL(N,rec,X,1)
                        Scale = Scale*rec
                        xmax = xmax*rec
                     ENDIF
                  ENDIF
!
                  d(1,1) = X(j1) - DDOT(j1-1,T(1,j1),1,X,1)
                  d(2,1) = X(j2) - DDOT(j1-1,T(1,j2),1,X,1)
!
                  CALL DLALN2(.TRUE.,2,1,smin,ONE,T(j1,j1),Ldt,ONE,ONE, &
     &                        d,2,ZERO,ZERO,v,2,scaloc,xnorm,ierr)
                  IF ( ierr/=0 ) Info = 2
!
                  IF ( scaloc/=ONE ) THEN
                     CALL DSCAL(N,scaloc,X,1)
                     Scale = Scale*scaloc
                  ENDIF
                  X(j1) = v(1,1)
                  X(j2) = v(2,1)
                  xmax = MAX(ABS(X(j1)),ABS(X(j2)),xmax)
!
               ENDIF
            ENDIF
         ENDDO
!
      ENDIF
!
!
!     End of DLAQTR
!
      END SUBROUTINE DLAQTR
