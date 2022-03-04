!*==sgehrd.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b SGEHRD
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SGEHRD + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgehrd.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgehrd.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgehrd.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SGEHRD( N, ILO, IHI, A, LDA, TAU, WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            IHI, ILO, INFO, LDA, LWORK, N
!       ..
!       .. Array Arguments ..
!       REAL              A( LDA, * ), TAU( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SGEHRD reduces a real general matrix A to upper Hessenberg form H by
!> an orthogonal similarity transformation:  Q**T * A * Q = H .
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] ILO
!> \verbatim
!>          ILO is INTEGER
!> \endverbatim
!>
!> \param[in] IHI
!> \verbatim
!>          IHI is INTEGER
!>
!>          It is assumed that A is already upper triangular in rows
!>          and columns 1:ILO-1 and IHI+1:N. ILO and IHI are normally
!>          set by a previous call to SGEBAL; otherwise they should be
!>          set to 1 and N respectively. See Further Details.
!>          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is REAL array, dimension (LDA,N)
!>          On entry, the N-by-N general matrix to be reduced.
!>          On exit, the upper triangle and the first subdiagonal of A
!>          are overwritten with the upper Hessenberg matrix H, and the
!>          elements below the first subdiagonal, with the array TAU,
!>          represent the orthogonal matrix Q as a product of elementary
!>          reflectors. See Further Details.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[out] TAU
!> \verbatim
!>          TAU is REAL array, dimension (N-1)
!>          The scalar factors of the elementary reflectors (see Further
!>          Details). Elements 1:ILO-1 and IHI:N-1 of TAU are set to
!>          zero.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (LWORK)
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The length of the array WORK.  LWORK >= max(1,N).
!>          For good performance, LWORK should generally be larger.
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
!> \ingroup realGEcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The matrix Q is represented as a product of (ihi-ilo) elementary
!>  reflectors
!>
!>     Q = H(ilo) H(ilo+1) . . . H(ihi-1).
!>
!>  Each H(i) has the form
!>
!>     H(i) = I - tau * v * v**T
!>
!>  where tau is a real scalar, and v is a real vector with
!>  v(1:i) = 0, v(i+1) = 1 and v(ihi+1:n) = 0; v(i+2:ihi) is stored on
!>  exit in A(i+2:ihi,i), and tau in TAU(i).
!>
!>  The contents of A are illustrated by the following example, with
!>  n = 7, ilo = 2 and ihi = 6:
!>
!>  on entry,                        on exit,
!>
!>  ( a   a   a   a   a   a   a )    (  a   a   h   h   h   h   a )
!>  (     a   a   a   a   a   a )    (      a   h   h   h   h   a )
!>  (     a   a   a   a   a   a )    (      h   h   h   h   h   h )
!>  (     a   a   a   a   a   a )    (      v2  h   h   h   h   h )
!>  (     a   a   a   a   a   a )    (      v2  v3  h   h   h   h )
!>  (     a   a   a   a   a   a )    (      v2  v3  v4  h   h   h )
!>  (                         a )    (                          a )
!>
!>  where a denotes an element of the original matrix A, h denotes a
!>  modified element of the upper Hessenberg matrix H, and vi denotes an
!>  element of the vector defining H(i).
!>
!>  This file is a slight modification of LAPACK-3.0's DGEHRD
!>  subroutine incorporating improvements proposed by Quintana-Orti and
!>  Van de Geijn (2006). (See DLAHR2.)
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE SGEHRD(N,Ilo,Ihi,A,Lda,Tau,Work,Lwork,Info)
      IMPLICIT NONE
!*--SGEHRD171
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Ihi , Ilo , Info , Lda , Lwork , N
!     ..
!     .. Array Arguments ..
      REAL A(Lda,*) , Tau(*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER NBMAX , LDT , TSIZE
      PARAMETER (NBMAX=64,LDT=NBMAX+1,TSIZE=LDT*NBMAX)
      REAL ZERO , ONE
      PARAMETER (ZERO=0.0E+0,ONE=1.0E+0)
!     ..
!     .. Local Scalars ..
      LOGICAL lquery
      INTEGER i , ib , iinfo , iwt , j , ldwork , lwkopt , nb , nbmin , &
     &        nh , nx
      REAL ei
!     ..
!     .. External Subroutines ..
      EXTERNAL SAXPY , SGEHD2 , SGEMM , SLAHR2 , SLARFB , STRMM , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX , MIN
!     ..
!     .. External Functions ..
      INTEGER ILAENV
      EXTERNAL ILAENV
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters
!
      Info = 0
      lquery = (Lwork==-1)
      IF ( N<0 ) THEN
         Info = -1
      ELSEIF ( Ilo<1 .OR. Ilo>MAX(1,N) ) THEN
         Info = -2
      ELSEIF ( Ihi<MIN(Ilo,N) .OR. Ihi>N ) THEN
         Info = -3
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -5
      ELSEIF ( Lwork<MAX(1,N) .AND. .NOT.lquery ) THEN
         Info = -8
      ENDIF
!
      IF ( Info==0 ) THEN
!
!       Compute the workspace requirements
!
         nb = MIN(NBMAX,ILAENV(1,'SGEHRD',' ',N,Ilo,Ihi,-1))
         lwkopt = N*nb + TSIZE
         Work(1) = lwkopt
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('SGEHRD',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
         RETURN
      ENDIF
!
!     Set elements 1:ILO-1 and IHI:N-1 of TAU to zero
!
      DO i = 1 , Ilo - 1
         Tau(i) = ZERO
      ENDDO
      DO i = MAX(1,Ihi) , N - 1
         Tau(i) = ZERO
      ENDDO
!
!     Quick return if possible
!
      nh = Ihi - Ilo + 1
      IF ( nh<=1 ) THEN
         Work(1) = 1
         RETURN
      ENDIF
!
!     Determine the block size
!
      nb = MIN(NBMAX,ILAENV(1,'SGEHRD',' ',N,Ilo,Ihi,-1))
      nbmin = 2
      IF ( nb>1 .AND. nb<nh ) THEN
!
!        Determine when to cross over from blocked to unblocked code
!        (last block is always handled by unblocked code)
!
         nx = MAX(nb,ILAENV(3,'SGEHRD',' ',N,Ilo,Ihi,-1))
         IF ( nx<nh ) THEN
!
!           Determine if workspace is large enough for blocked code
!
            IF ( Lwork<N*nb+TSIZE ) THEN
!
!              Not enough workspace to use optimal NB:  determine the
!              minimum value of NB, and reduce NB or force use of
!              unblocked code
!
               nbmin = MAX(2,ILAENV(2,'SGEHRD',' ',N,Ilo,Ihi,-1))
               IF ( Lwork>=(N*nbmin+TSIZE) ) THEN
                  nb = (Lwork-TSIZE)/N
               ELSE
                  nb = 1
               ENDIF
            ENDIF
         ENDIF
      ENDIF
      ldwork = N
!
      IF ( nb<nbmin .OR. nb>=nh ) THEN
!
!        Use unblocked code below
!
         i = Ilo
!
      ELSE
!
!        Use blocked code
!
         iwt = 1 + N*nb
         DO i = Ilo , Ihi - 1 - nx , nb
            ib = MIN(nb,Ihi-i)
!
!           Reduce columns i:i+ib-1 to Hessenberg form, returning the
!           matrices V and T of the block reflector H = I - V*T*V**T
!           which performs the reduction, and also the matrix Y = A*V*T
!
            CALL SLAHR2(Ihi,i,ib,A(1,i),Lda,Tau(i),Work(iwt),LDT,Work,  &
     &                  ldwork)
!
!           Apply the block reflector H to A(1:ihi,i+ib:ihi) from the
!           right, computing  A := A - Y * V**T. V(i+ib,ib-1) must be set
!           to 1
!
            ei = A(i+ib,i+ib-1)
            A(i+ib,i+ib-1) = ONE
            CALL SGEMM('No transpose','Transpose',Ihi,Ihi-i-ib+1,ib,    &
     &                 -ONE,Work,ldwork,A(i+ib,i),Lda,ONE,A(1,i+ib),Lda)
            A(i+ib,i+ib-1) = ei
!
!           Apply the block reflector H to A(1:i,i+1:i+ib-1) from the
!           right
!
            CALL STRMM('Right','Lower','Transpose','Unit',i,ib-1,ONE,   &
     &                 A(i+1,i),Lda,Work,ldwork)
            DO j = 0 , ib - 2
               CALL SAXPY(i,-ONE,Work(ldwork*j+1),1,A(1,i+j+1),1)
            ENDDO
!
!           Apply the block reflector H to A(i+1:ihi,i+ib:n) from the
!           left
!
            CALL SLARFB('Left','Transpose','Forward','Columnwise',Ihi-i,&
     &                  N-i-ib+1,ib,A(i+1,i),Lda,Work(iwt),LDT,         &
     &                  A(i+1,i+ib),Lda,Work,ldwork)
         ENDDO
      ENDIF
!
!     Use unblocked code to reduce the rest of the matrix
!
      CALL SGEHD2(N,i,Ihi,A,Lda,Tau,Work,iinfo)
      Work(1) = lwkopt
!
!
!     End of SGEHRD
!
      END SUBROUTINE SGEHRD
