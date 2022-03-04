!*==zunm22.f90  processed by SPAG 7.51RB at 20:09 on  3 Mar 2022
!> \brief \b ZUNM22 multiplies a general matrix by a banded unitary matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZUNM22 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zunm22.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zunm22.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zunm22.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!     SUBROUTINE ZUNM22( SIDE, TRANS, M, N, N1, N2, Q, LDQ, C, LDC,
!    $                   WORK, LWORK, INFO )
!
!     .. Scalar Arguments ..
!     CHARACTER          SIDE, TRANS
!     INTEGER            M, N, N1, N2, LDQ, LDC, LWORK, INFO
!     ..
!     .. Array Arguments ..
!     COMPLEX*16            Q( LDQ, * ), C( LDC, * ), WORK( * )
!     ..
!
!> \par Purpose
!  ============
!>
!> \verbatim
!>
!>  ZUNM22 overwrites the general complex M-by-N matrix C with
!>
!>                  SIDE = 'L'     SIDE = 'R'
!>  TRANS = 'N':      Q * C          C * Q
!>  TRANS = 'C':      Q**H * C       C * Q**H
!>
!>  where Q is a complex unitary matrix of order NQ, with NQ = M if
!>  SIDE = 'L' and NQ = N if SIDE = 'R'.
!>  The unitary matrix Q processes a 2-by-2 block structure
!>
!>         [  Q11  Q12  ]
!>     Q = [            ]
!>         [  Q21  Q22  ],
!>
!>  where Q12 is an N1-by-N1 lower triangular matrix and Q21 is an
!>  N2-by-N2 upper triangular matrix.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SIDE
!> \verbatim
!>          SIDE is CHARACTER*1
!>          = 'L': apply Q or Q**H from the Left;
!>          = 'R': apply Q or Q**H from the Right.
!> \endverbatim
!>
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>          = 'N':  apply Q (No transpose);
!>          = 'C':  apply Q**H (Conjugate transpose).
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix C. M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix C. N >= 0.
!> \endverbatim
!>
!> \param[in] N1
!> \param[in] N2
!> \verbatim
!>          N1 is INTEGER
!>          N2 is INTEGER
!>          The dimension of Q12 and Q21, respectively. N1, N2 >= 0.
!>          The following requirement must be satisfied:
!>          N1 + N2 = M if SIDE = 'L' and N1 + N2 = N if SIDE = 'R'.
!> \endverbatim
!>
!> \param[in] Q
!> \verbatim
!>          Q is COMPLEX*16 array, dimension
!>                              (LDQ,M) if SIDE = 'L'
!>                              (LDQ,N) if SIDE = 'R'
!> \endverbatim
!>
!> \param[in] LDQ
!> \verbatim
!>          LDQ is INTEGER
!>          The leading dimension of the array Q.
!>          LDQ >= max(1,M) if SIDE = 'L'; LDQ >= max(1,N) if SIDE = 'R'.
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is COMPLEX*16 array, dimension (LDC,N)
!>          On entry, the M-by-N matrix C.
!>          On exit, C is overwritten by Q*C or Q**H*C or C*Q**H or C*Q.
!> \endverbatim
!>
!> \param[in] LDC
!> \verbatim
!>          LDC is INTEGER
!>          The leading dimension of the array C. LDC >= max(1,M).
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
!>          If SIDE = 'L', LWORK >= max(1,N);
!>          if SIDE = 'R', LWORK >= max(1,M).
!>          For optimum performance LWORK >= M*N.
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
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \date January 2015
!
!> \ingroup complexOTHERcomputational
!
!  =====================================================================
      SUBROUTINE ZUNM22(Side,Trans,M,N,N1,N2,Q,Ldq,C,Ldc,Work,Lwork,    &
     &                  Info)
!
!  -- LAPACK computational routine (version 3.7.1) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     January 2015
!
      IMPLICIT NONE
!*--ZUNM22172
!
!     .. Scalar Arguments ..
      CHARACTER Side , Trans
      INTEGER M , N , N1 , N2 , Ldq , Ldc , Lwork , Info
!     ..
!     .. Array Arguments ..
      COMPLEX*16 Q(Ldq,*) , C(Ldc,*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX*16 ONE
      PARAMETER (ONE=(1.0D+0,0.0D+0))
!
!     .. Local Scalars ..
      LOGICAL left , lquery , notran
      INTEGER i , ldwork , len , lwkopt , nb , nq , nw
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL ZGEMM , ZLACPY , ZTRMM , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DCMPLX , MAX , MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      Info = 0
      left = LSAME(Side,'L')
      notran = LSAME(Trans,'N')
      lquery = (Lwork==-1)
!
!     NQ is the order of Q;
!     NW is the minimum dimension of WORK.
!
      IF ( left ) THEN
         nq = M
      ELSE
         nq = N
      ENDIF
      nw = nq
      IF ( N1==0 .OR. N2==0 ) nw = 1
      IF ( .NOT.left .AND. .NOT.LSAME(Side,'R') ) THEN
         Info = -1
      ELSEIF ( .NOT.LSAME(Trans,'N') .AND. .NOT.LSAME(Trans,'C') ) THEN
         Info = -2
      ELSEIF ( M<0 ) THEN
         Info = -3
      ELSEIF ( N<0 ) THEN
         Info = -4
      ELSEIF ( N1<0 .OR. N1+N2/=nq ) THEN
         Info = -5
      ELSEIF ( N2<0 ) THEN
         Info = -6
      ELSEIF ( Ldq<MAX(1,nq) ) THEN
         Info = -8
      ELSEIF ( Ldc<MAX(1,M) ) THEN
         Info = -10
      ELSEIF ( Lwork<nw .AND. .NOT.lquery ) THEN
         Info = -12
      ENDIF
!
      IF ( Info==0 ) THEN
         lwkopt = M*N
         Work(1) = DCMPLX(lwkopt)
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('ZUNM22',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( M==0 .OR. N==0 ) THEN
         Work(1) = 1
         RETURN
      ENDIF
!
!     Degenerate cases (N1 = 0 or N2 = 0) are handled using ZTRMM.
!
      IF ( N1==0 ) THEN
         CALL ZTRMM(Side,'Upper',Trans,'Non-Unit',M,N,ONE,Q,Ldq,C,Ldc)
         Work(1) = ONE
         RETURN
      ELSEIF ( N2==0 ) THEN
         CALL ZTRMM(Side,'Lower',Trans,'Non-Unit',M,N,ONE,Q,Ldq,C,Ldc)
         Work(1) = ONE
         RETURN
      ENDIF
!
!     Compute the largest chunk size available from the workspace.
!
      nb = MAX(1,MIN(Lwork,lwkopt)/nq)
!
      IF ( left ) THEN
         IF ( notran ) THEN
            DO i = 1 , N , nb
               len = MIN(nb,N-i+1)
               ldwork = M
!
!              Multiply bottom part of C by Q12.
!
               CALL ZLACPY('All',N1,len,C(N2+1,i),Ldc,Work,ldwork)
               CALL ZTRMM('Left','Lower','No Transpose','Non-Unit',N1,  &
     &                    len,ONE,Q(1,N2+1),Ldq,Work,ldwork)
!
!              Multiply top part of C by Q11.
!
               CALL ZGEMM('No Transpose','No Transpose',N1,len,N2,ONE,Q,&
     &                    Ldq,C(1,i),Ldc,ONE,Work,ldwork)
!
!              Multiply top part of C by Q21.
!
               CALL ZLACPY('All',N2,len,C(1,i),Ldc,Work(N1+1),ldwork)
               CALL ZTRMM('Left','Upper','No Transpose','Non-Unit',N2,  &
     &                    len,ONE,Q(N1+1,1),Ldq,Work(N1+1),ldwork)
!
!              Multiply bottom part of C by Q22.
!
               CALL ZGEMM('No Transpose','No Transpose',N2,len,N1,ONE,  &
     &                    Q(N1+1,N2+1),Ldq,C(N2+1,i),Ldc,ONE,Work(N1+1),&
     &                    ldwork)
!
!              Copy everything back.
!
               CALL ZLACPY('All',M,len,Work,ldwork,C(1,i),Ldc)
            ENDDO
         ELSE
            DO i = 1 , N , nb
               len = MIN(nb,N-i+1)
               ldwork = M
!
!              Multiply bottom part of C by Q21**H.
!
               CALL ZLACPY('All',N2,len,C(N1+1,i),Ldc,Work,ldwork)
               CALL ZTRMM('Left','Upper','Conjugate','Non-Unit',N2,len, &
     &                    ONE,Q(N1+1,1),Ldq,Work,ldwork)
!
!              Multiply top part of C by Q11**H.
!
               CALL ZGEMM('Conjugate','No Transpose',N2,len,N1,ONE,Q,   &
     &                    Ldq,C(1,i),Ldc,ONE,Work,ldwork)
!
!              Multiply top part of C by Q12**H.
!
               CALL ZLACPY('All',N1,len,C(1,i),Ldc,Work(N2+1),ldwork)
               CALL ZTRMM('Left','Lower','Conjugate','Non-Unit',N1,len, &
     &                    ONE,Q(1,N2+1),Ldq,Work(N2+1),ldwork)
!
!              Multiply bottom part of C by Q22**H.
!
               CALL ZGEMM('Conjugate','No Transpose',N1,len,N2,ONE,     &
     &                    Q(N1+1,N2+1),Ldq,C(N1+1,i),Ldc,ONE,Work(N2+1),&
     &                    ldwork)
!
!              Copy everything back.
!
               CALL ZLACPY('All',M,len,Work,ldwork,C(1,i),Ldc)
            ENDDO
         ENDIF
      ELSEIF ( notran ) THEN
         DO i = 1 , M , nb
            len = MIN(nb,M-i+1)
            ldwork = len
!
!              Multiply right part of C by Q21.
!
            CALL ZLACPY('All',len,N2,C(i,N1+1),Ldc,Work,ldwork)
            CALL ZTRMM('Right','Upper','No Transpose','Non-Unit',len,N2,&
     &                 ONE,Q(N1+1,1),Ldq,Work,ldwork)
!
!              Multiply left part of C by Q11.
!
            CALL ZGEMM('No Transpose','No Transpose',len,N2,N1,ONE,     &
     &                 C(i,1),Ldc,Q,Ldq,ONE,Work,ldwork)
!
!              Multiply left part of C by Q12.
!
            CALL ZLACPY('All',len,N1,C(i,1),Ldc,Work(1+N2*ldwork),      &
     &                  ldwork)
            CALL ZTRMM('Right','Lower','No Transpose','Non-Unit',len,N1,&
     &                 ONE,Q(1,N2+1),Ldq,Work(1+N2*ldwork),ldwork)
!
!              Multiply right part of C by Q22.
!
            CALL ZGEMM('No Transpose','No Transpose',len,N1,N2,ONE,     &
     &                 C(i,N1+1),Ldc,Q(N1+1,N2+1),Ldq,ONE,              &
     &                 Work(1+N2*ldwork),ldwork)
!
!              Copy everything back.
!
            CALL ZLACPY('All',len,N,Work,ldwork,C(i,1),Ldc)
         ENDDO
      ELSE
         DO i = 1 , M , nb
            len = MIN(nb,M-i+1)
            ldwork = len
!
!              Multiply right part of C by Q12**H.
!
            CALL ZLACPY('All',len,N1,C(i,N2+1),Ldc,Work,ldwork)
            CALL ZTRMM('Right','Lower','Conjugate','Non-Unit',len,N1,   &
     &                 ONE,Q(1,N2+1),Ldq,Work,ldwork)
!
!              Multiply left part of C by Q11**H.
!
            CALL ZGEMM('No Transpose','Conjugate',len,N1,N2,ONE,C(i,1), &
     &                 Ldc,Q,Ldq,ONE,Work,ldwork)
!
!              Multiply left part of C by Q21**H.
!
            CALL ZLACPY('All',len,N2,C(i,1),Ldc,Work(1+N1*ldwork),      &
     &                  ldwork)
            CALL ZTRMM('Right','Upper','Conjugate','Non-Unit',len,N2,   &
     &                 ONE,Q(N1+1,1),Ldq,Work(1+N1*ldwork),ldwork)
!
!              Multiply right part of C by Q22**H.
!
            CALL ZGEMM('No Transpose','Conjugate',len,N2,N1,ONE,        &
     &                 C(i,N2+1),Ldc,Q(N1+1,N2+1),Ldq,ONE,              &
     &                 Work(1+N1*ldwork),ldwork)
!
!              Copy everything back.
!
            CALL ZLACPY('All',len,N,Work,ldwork,C(i,1),Ldc)
         ENDDO
      ENDIF
!
      Work(1) = DCMPLX(lwkopt)
!
!     End of ZUNM22
!
      END SUBROUTINE ZUNM22
