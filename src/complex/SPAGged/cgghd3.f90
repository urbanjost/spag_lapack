!*==cgghd3.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b CGGHD3
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CGGHD3 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgghd3.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgghd3.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgghd3.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!        SUBROUTINE CGGHD3( COMPQ, COMPZ, N, ILO, IHI, A, LDA, B, LDB, Q,
!       $                   LDQ, Z, LDZ, WORK, LWORK, INFO )
!
!        .. Scalar Arguments ..
!        CHARACTER          COMPQ, COMPZ
!        INTEGER            IHI, ILO, INFO, LDA, LDB, LDQ, LDZ, N, LWORK
!        ..
!        .. Array Arguments ..
!        COMPLEX            A( LDA, * ), B( LDB, * ), Q( LDQ, * ),
!       $                   Z( LDZ, * ), WORK( * )
!        ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>
!> CGGHD3 reduces a pair of complex matrices (A,B) to generalized upper
!> Hessenberg form using unitary transformations, where A is a
!> general matrix and B is upper triangular.  The form of the
!> generalized eigenvalue problem is
!>    A*x = lambda*B*x,
!> and B is typically made upper triangular by computing its QR
!> factorization and moving the unitary matrix Q to the left side
!> of the equation.
!>
!> This subroutine simultaneously reduces A to a Hessenberg matrix H:
!>    Q**H*A*Z = H
!> and transforms B to another upper triangular matrix T:
!>    Q**H*B*Z = T
!> in order to reduce the problem to its standard form
!>    H*y = lambda*T*y
!> where y = Z**H*x.
!>
!> The unitary matrices Q and Z are determined as products of Givens
!> rotations.  They may either be formed explicitly, or they may be
!> postmultiplied into input matrices Q1 and Z1, so that
!>
!>      Q1 * A * Z1**H = (Q1*Q) * H * (Z1*Z)**H
!>
!>      Q1 * B * Z1**H = (Q1*Q) * T * (Z1*Z)**H
!>
!> If Q1 is the unitary matrix from the QR factorization of B in the
!> original equation A*x = lambda*B*x, then CGGHD3 reduces the original
!> problem to generalized Hessenberg form.
!>
!> This is a blocked variant of CGGHRD, using matrix-matrix
!> multiplications for parts of the computation to enhance performance.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] COMPQ
!> \verbatim
!>          COMPQ is CHARACTER*1
!>          = 'N': do not compute Q;
!>          = 'I': Q is initialized to the unit matrix, and the
!>                 unitary matrix Q is returned;
!>          = 'V': Q must contain a unitary matrix Q1 on entry,
!>                 and the product Q1*Q is returned.
!> \endverbatim
!>
!> \param[in] COMPZ
!> \verbatim
!>          COMPZ is CHARACTER*1
!>          = 'N': do not compute Z;
!>          = 'I': Z is initialized to the unit matrix, and the
!>                 unitary matrix Z is returned;
!>          = 'V': Z must contain a unitary matrix Z1 on entry,
!>                 and the product Z1*Z is returned.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrices A and B.  N >= 0.
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
!>          ILO and IHI mark the rows and columns of A which are to be
!>          reduced.  It is assumed that A is already upper triangular
!>          in rows and columns 1:ILO-1 and IHI+1:N.  ILO and IHI are
!>          normally set by a previous call to CGGBAL; otherwise they
!>          should be set to 1 and N respectively.
!>          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA, N)
!>          On entry, the N-by-N general matrix to be reduced.
!>          On exit, the upper triangle and the first subdiagonal of A
!>          are overwritten with the upper Hessenberg matrix H, and the
!>          rest is set to zero.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is COMPLEX array, dimension (LDB, N)
!>          On entry, the N-by-N upper triangular matrix B.
!>          On exit, the upper triangular matrix T = Q**H B Z.  The
!>          elements below the diagonal are set to zero.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B.  LDB >= max(1,N).
!> \endverbatim
!>
!> \param[in,out] Q
!> \verbatim
!>          Q is COMPLEX array, dimension (LDQ, N)
!>          On entry, if COMPQ = 'V', the unitary matrix Q1, typically
!>          from the QR factorization of B.
!>          On exit, if COMPQ='I', the unitary matrix Q, and if
!>          COMPQ = 'V', the product Q1*Q.
!>          Not referenced if COMPQ='N'.
!> \endverbatim
!>
!> \param[in] LDQ
!> \verbatim
!>          LDQ is INTEGER
!>          The leading dimension of the array Q.
!>          LDQ >= N if COMPQ='V' or 'I'; LDQ >= 1 otherwise.
!> \endverbatim
!>
!> \param[in,out] Z
!> \verbatim
!>          Z is COMPLEX array, dimension (LDZ, N)
!>          On entry, if COMPZ = 'V', the unitary matrix Z1.
!>          On exit, if COMPZ='I', the unitary matrix Z, and if
!>          COMPZ = 'V', the product Z1*Z.
!>          Not referenced if COMPZ='N'.
!> \endverbatim
!>
!> \param[in] LDZ
!> \verbatim
!>          LDZ is INTEGER
!>          The leading dimension of the array Z.
!>          LDZ >= N if COMPZ='V' or 'I'; LDZ >= 1 otherwise.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension (LWORK)
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in]  LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The length of the array WORK.  LWORK >= 1.
!>          For optimum performance LWORK >= 6*N*NB, where NB is the
!>          optimal blocksize.
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
!> \date January 2015
!
!> \ingroup complexOTHERcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  This routine reduces A to Hessenberg form and maintains B in triangular form
!>  using a blocked variant of Moler and Stewart's original algorithm,
!>  as described by Kagstrom, Kressner, Quintana-Orti, and Quintana-Orti
!>  (BIT 2008).
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE CGGHD3(Compq,Compz,N,Ilo,Ihi,A,Lda,B,Ldb,Q,Ldq,Z,Ldz,  &
     &                  Work,Lwork,Info)
!
!  -- LAPACK computational routine (version 3.8.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     January 2015
!
!
      USE F77KINDS                        
      USE S_CGEMM
      USE S_CGEMV
      USE S_CGGHRD
      USE S_CLACPY
      USE S_CLARTG
      USE S_CLASET
      USE S_CROT
      USE S_CTRMV
      USE S_CUNM22
      USE S_ILAENV
      USE S_LSAME
      USE S_XERBLA
      IMPLICIT NONE
!*--CGGHD3255
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX , PARAMETER  ::  CONE = (1.0E+0,0.0E+0) ,                 &
     &                         CZERO = (0.0E+0,0.0E+0)
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Compq
      CHARACTER :: Compz
      INTEGER :: N
      INTEGER :: Ilo
      INTEGER :: Ihi
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      COMPLEX , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      LOGICAL :: blk22 , initq , initz , lquery , wantq , wantz
      REAL :: c
      COMPLEX :: c1 , c2 , ctemp , s , s1 , s2 , temp , temp1 , temp2 , &
     &           temp3
      INTEGER :: cola , i , ierr , j , j0 , jcol , jj , jrow , k ,      &
     &           kacc22 , len , lwkopt , n2nb , nb , nblst , nbmin ,    &
     &           nh , nnb , nx , ppw , ppwo , pw , top , topq
      CHARACTER(1) :: compq2 , compz2
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
!     .. External Functions ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
!     Decode and test the input parameters.
!
      Info = 0
      nb = ILAENV(1,'CGGHD3',' ',N,Ilo,Ihi,-1)
      lwkopt = MAX(6*N*nb,1)
      Work(1) = CMPLX(lwkopt)
      initq = LSAME(Compq,'I')
      wantq = initq .OR. LSAME(Compq,'V')
      initz = LSAME(Compz,'I')
      wantz = initz .OR. LSAME(Compz,'V')
      lquery = (Lwork==-1)
!
      IF ( .NOT.LSAME(Compq,'N') .AND. .NOT.wantq ) THEN
         Info = -1
      ELSEIF ( .NOT.LSAME(Compz,'N') .AND. .NOT.wantz ) THEN
         Info = -2
      ELSEIF ( N<0 ) THEN
         Info = -3
      ELSEIF ( Ilo<1 ) THEN
         Info = -4
      ELSEIF ( Ihi>N .OR. Ihi<Ilo-1 ) THEN
         Info = -5
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -7
      ELSEIF ( Ldb<MAX(1,N) ) THEN
         Info = -9
      ELSEIF ( (wantq .AND. Ldq<N) .OR. Ldq<1 ) THEN
         Info = -11
      ELSEIF ( (wantz .AND. Ldz<N) .OR. Ldz<1 ) THEN
         Info = -13
      ELSEIF ( Lwork<1 .AND. .NOT.lquery ) THEN
         Info = -15
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('CGGHD3',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
         RETURN
      ENDIF
!
!     Initialize Q and Z if desired.
!
      IF ( initq ) CALL CLASET('All',N,N,CZERO,CONE,Q,Ldq)
      IF ( initz ) CALL CLASET('All',N,N,CZERO,CONE,Z,Ldz)
!
!     Zero out lower triangle of B.
!
      IF ( N>1 ) CALL CLASET('Lower',N-1,N-1,CZERO,CZERO,B(2,1),Ldb)
!
!     Quick return if possible
!
      nh = Ihi - Ilo + 1
      IF ( nh<=1 ) THEN
         Work(1) = CONE
         RETURN
      ENDIF
!
!     Determine the blocksize.
!
      nbmin = ILAENV(2,'CGGHD3',' ',N,Ilo,Ihi,-1)
      IF ( nb>1 .AND. nb<nh ) THEN
!
!        Determine when to use unblocked instead of blocked code.
!
         nx = MAX(nb,ILAENV(3,'CGGHD3',' ',N,Ilo,Ihi,-1))
         IF ( nx<nh ) THEN
!
!           Determine if workspace is large enough for blocked code.
!
            IF ( Lwork<lwkopt ) THEN
!
!              Not enough workspace to use optimal NB:  determine the
!              minimum value of NB, and reduce NB or force use of
!              unblocked code.
!
               nbmin = MAX(2,ILAENV(2,'CGGHD3',' ',N,Ilo,Ihi,-1))
               IF ( Lwork>=6*N*nbmin ) THEN
                  nb = Lwork/(6*N)
               ELSE
                  nb = 1
               ENDIF
            ENDIF
         ENDIF
      ENDIF
!
      IF ( nb<nbmin .OR. nb>=nh ) THEN
!
!        Use unblocked code below
!
         jcol = Ilo
!
      ELSE
!
!        Use blocked code
!
         kacc22 = ILAENV(16,'CGGHD3',' ',N,Ilo,Ihi,-1)
         blk22 = kacc22==2
         DO jcol = Ilo , Ihi - 2 , nb
            nnb = MIN(nb,Ihi-jcol-1)
!
!           Initialize small unitary factors that will hold the
!           accumulated Givens rotations in workspace.
!           N2NB   denotes the number of 2*NNB-by-2*NNB factors
!           NBLST  denotes the (possibly smaller) order of the last
!                  factor.
!
            n2nb = (Ihi-jcol-1)/nnb - 1
            nblst = Ihi - jcol - n2nb*nnb
            CALL CLASET('All',nblst,nblst,CZERO,CONE,Work,nblst)
            pw = nblst*nblst + 1
            DO i = 1 , n2nb
               CALL CLASET('All',2*nnb,2*nnb,CZERO,CONE,Work(pw),2*nnb)
               pw = pw + 4*nnb*nnb
            ENDDO
!
!           Reduce columns JCOL:JCOL+NNB-1 of A to Hessenberg form.
!
            DO j = jcol , jcol + nnb - 1
!
!              Reduce Jth column of A. Store cosines and sines in Jth
!              column of A and B, respectively.
!
               DO i = Ihi , j + 2 , -1
                  temp = A(i-1,j)
                  CALL CLARTG(temp,A(i,j),c,s,A(i-1,j))
                  A(i,j) = CMPLX(c)
                  B(i,j) = s
               ENDDO
!
!              Accumulate Givens rotations into workspace array.
!
               ppw = (nblst+1)*(nblst-2) - j + jcol + 1
               len = 2 + j - jcol
               jrow = j + n2nb*nnb + 2
               DO i = Ihi , jrow , -1
                  ctemp = A(i,j)
                  s = B(i,j)
                  DO jj = ppw , ppw + len - 1
                     temp = Work(jj+nblst)
                     Work(jj+nblst) = ctemp*temp - s*Work(jj)
                     Work(jj) = CONJG(s)*temp + ctemp*Work(jj)
                  ENDDO
                  len = len + 1
                  ppw = ppw - nblst - 1
               ENDDO
!
               ppwo = nblst*nblst + (nnb+j-jcol-1)*2*nnb + nnb
               j0 = jrow - nnb
               DO jrow = j0 , j + 2 , -nnb
                  ppw = ppwo
                  len = 2 + j - jcol
                  DO i = jrow + nnb - 1 , jrow , -1
                     ctemp = A(i,j)
                     s = B(i,j)
                     DO jj = ppw , ppw + len - 1
                        temp = Work(jj+2*nnb)
                        Work(jj+2*nnb) = ctemp*temp - s*Work(jj)
                        Work(jj) = CONJG(s)*temp + ctemp*Work(jj)
                     ENDDO
                     len = len + 1
                     ppw = ppw - 2*nnb - 1
                  ENDDO
                  ppwo = ppwo + 4*nnb*nnb
               ENDDO
!
!              TOP denotes the number of top rows in A and B that will
!              not be updated during the next steps.
!
               IF ( jcol<=2 ) THEN
                  top = 0
               ELSE
                  top = jcol
               ENDIF
!
!              Propagate transformations through B and replace stored
!              left sines/cosines by right sines/cosines.
!
               DO jj = N , j + 1 , -1
!
!                 Update JJth column of B.
!
                  DO i = MIN(jj+1,Ihi) , j + 2 , -1
                     ctemp = A(i,j)
                     s = B(i,j)
                     temp = B(i,jj)
                     B(i,jj) = ctemp*temp - CONJG(s)*B(i-1,jj)
                     B(i-1,jj) = s*temp + ctemp*B(i-1,jj)
                  ENDDO
!
!                 Annihilate B( JJ+1, JJ ).
!
                  IF ( jj<Ihi ) THEN
                     temp = B(jj+1,jj+1)
                     CALL CLARTG(temp,B(jj+1,jj),c,s,B(jj+1,jj+1))
                     B(jj+1,jj) = CZERO
                     CALL CROT(jj-top,B(top+1,jj+1),1,B(top+1,jj),1,c,s)
                     A(jj+1,j) = CMPLX(c)
                     B(jj+1,j) = -CONJG(s)
                  ENDIF
               ENDDO
!
!              Update A by transformations from right.
!
               jj = MOD(Ihi-j-1,3)
               DO i = Ihi - j - 3 , jj + 1 , -3
                  ctemp = A(j+1+i,j)
                  s = -B(j+1+i,j)
                  c1 = A(j+2+i,j)
                  s1 = -B(j+2+i,j)
                  c2 = A(j+3+i,j)
                  s2 = -B(j+3+i,j)
!
                  DO k = top + 1 , Ihi
                     temp = A(k,j+i)
                     temp1 = A(k,j+i+1)
                     temp2 = A(k,j+i+2)
                     temp3 = A(k,j+i+3)
                     A(k,j+i+3) = c2*temp3 + CONJG(s2)*temp2
                     temp2 = -s2*temp3 + c2*temp2
                     A(k,j+i+2) = c1*temp2 + CONJG(s1)*temp1
                     temp1 = -s1*temp2 + c1*temp1
                     A(k,j+i+1) = ctemp*temp1 + CONJG(s)*temp
                     A(k,j+i) = -s*temp1 + ctemp*temp
                  ENDDO
               ENDDO
!
               IF ( jj>0 ) THEN
                  DO i = jj , 1 , -1
                     c = DBLE(A(j+1+i,j))
                     CALL CROT(Ihi-top,A(top+1,j+i+1),1,A(top+1,j+i),1, &
     &                         c,-CONJG(B(j+1+i,j)))
                  ENDDO
               ENDIF
!
!              Update (J+1)th column of A by transformations from left.
!
               IF ( j<jcol+nnb-1 ) THEN
                  len = 1 + j - jcol
!
!                 Multiply with the trailing accumulated unitary
!                 matrix, which takes the form
!
!                        [  U11  U12  ]
!                    U = [            ],
!                        [  U21  U22  ]
!
!                 where U21 is a LEN-by-LEN matrix and U12 is lower
!                 triangular.
!
                  jrow = Ihi - nblst + 1
                  CALL CGEMV('Conjugate',nblst,len,CONE,Work,nblst,     &
     &                       A(jrow,j+1),1,CZERO,Work(pw),1)
                  ppw = pw + len
                  DO i = jrow , jrow + nblst - len - 1
                     Work(ppw) = A(i,j+1)
                     ppw = ppw + 1
                  ENDDO
                  CALL CTRMV('Lower','Conjugate','Non-unit',nblst-len,  &
     &                       Work(len*nblst+1),nblst,Work(pw+len),1)
                  CALL CGEMV('Conjugate',len,nblst-len,CONE,            &
     &                       Work((len+1)*nblst-len+1),nblst,           &
     &                       A(jrow+nblst-len,j+1),1,CONE,Work(pw+len), &
     &                       1)
                  ppw = pw
                  DO i = jrow , jrow + nblst - 1
                     A(i,j+1) = Work(ppw)
                     ppw = ppw + 1
                  ENDDO
!
!                 Multiply with the other accumulated unitary
!                 matrices, which take the form
!
!                        [  U11  U12   0  ]
!                        [                ]
!                    U = [  U21  U22   0  ],
!                        [                ]
!                        [   0    0    I  ]
!
!                 where I denotes the (NNB-LEN)-by-(NNB-LEN) identity
!                 matrix, U21 is a LEN-by-LEN upper triangular matrix
!                 and U12 is an NNB-by-NNB lower triangular matrix.
!
                  ppwo = 1 + nblst*nblst
                  j0 = jrow - nnb
                  DO jrow = j0 , jcol + 1 , -nnb
                     ppw = pw + len
                     DO i = jrow , jrow + nnb - 1
                        Work(ppw) = A(i,j+1)
                        ppw = ppw + 1
                     ENDDO
                     ppw = pw
                     DO i = jrow + nnb , jrow + nnb + len - 1
                        Work(ppw) = A(i,j+1)
                        ppw = ppw + 1
                     ENDDO
                     CALL CTRMV('Upper','Conjugate','Non-unit',len,     &
     &                          Work(ppwo+nnb),2*nnb,Work(pw),1)
                     CALL CTRMV('Lower','Conjugate','Non-unit',nnb,     &
     &                          Work(ppwo+2*len*nnb),2*nnb,Work(pw+len),&
     &                          1)
                     CALL CGEMV('Conjugate',nnb,len,CONE,Work(ppwo),    &
     &                          2*nnb,A(jrow,j+1),1,CONE,Work(pw),1)
                     CALL CGEMV('Conjugate',len,nnb,CONE,               &
     &                          Work(ppwo+2*len*nnb+nnb),2*nnb,         &
     &                          A(jrow+nnb,j+1),1,CONE,Work(pw+len),1)
                     ppw = pw
                     DO i = jrow , jrow + len + nnb - 1
                        A(i,j+1) = Work(ppw)
                        ppw = ppw + 1
                     ENDDO
                     ppwo = ppwo + 4*nnb*nnb
                  ENDDO
               ENDIF
            ENDDO
!
!           Apply accumulated unitary matrices to A.
!
            cola = N - jcol - nnb + 1
            j = Ihi - nblst + 1
            CALL CGEMM('Conjugate','No Transpose',nblst,cola,nblst,CONE,&
     &                 Work,nblst,A(j,jcol+nnb),Lda,CZERO,Work(pw),     &
     &                 nblst)
            CALL CLACPY('All',nblst,cola,Work(pw),nblst,A(j,jcol+nnb),  &
     &                  Lda)
            ppwo = nblst*nblst + 1
            j0 = j - nnb
            DO j = j0 , jcol + 1 , -nnb
               IF ( blk22 ) THEN
!
!                 Exploit the structure of
!
!                        [  U11  U12  ]
!                    U = [            ]
!                        [  U21  U22  ],
!
!                 where all blocks are NNB-by-NNB, U21 is upper
!                 triangular and U12 is lower triangular.
!
                  CALL CUNM22('Left','Conjugate',2*nnb,cola,nnb,nnb,    &
     &                        Work(ppwo),2*nnb,A(j,jcol+nnb),Lda,       &
     &                        Work(pw),Lwork-pw+1,ierr)
               ELSE
!
!                 Ignore the structure of U.
!
                  CALL CGEMM('Conjugate','No Transpose',2*nnb,cola,     &
     &                       2*nnb,CONE,Work(ppwo),2*nnb,A(j,jcol+nnb), &
     &                       Lda,CZERO,Work(pw),2*nnb)
                  CALL CLACPY('All',2*nnb,cola,Work(pw),2*nnb,          &
     &                        A(j,jcol+nnb),Lda)
               ENDIF
               ppwo = ppwo + 4*nnb*nnb
            ENDDO
!
!           Apply accumulated unitary matrices to Q.
!
            IF ( wantq ) THEN
               j = Ihi - nblst + 1
               IF ( initq ) THEN
                  topq = MAX(2,j-jcol+1)
                  nh = Ihi - topq + 1
               ELSE
                  topq = 1
                  nh = N
               ENDIF
               CALL CGEMM('No Transpose','No Transpose',nh,nblst,nblst, &
     &                    CONE,Q(topq,j),Ldq,Work,nblst,CZERO,Work(pw), &
     &                    nh)
               CALL CLACPY('All',nh,nblst,Work(pw),nh,Q(topq,j),Ldq)
               ppwo = nblst*nblst + 1
               j0 = j - nnb
               DO j = j0 , jcol + 1 , -nnb
                  IF ( initq ) THEN
                     topq = MAX(2,j-jcol+1)
                     nh = Ihi - topq + 1
                  ENDIF
                  IF ( blk22 ) THEN
!
!                    Exploit the structure of U.
!
                     CALL CUNM22('Right','No Transpose',nh,2*nnb,nnb,   &
     &                           nnb,Work(ppwo),2*nnb,Q(topq,j),Ldq,    &
     &                           Work(pw),Lwork-pw+1,ierr)
                  ELSE
!
!                    Ignore the structure of U.
!
                     CALL CGEMM('No Transpose','No Transpose',nh,2*nnb, &
     &                          2*nnb,CONE,Q(topq,j),Ldq,Work(ppwo),    &
     &                          2*nnb,CZERO,Work(pw),nh)
                     CALL CLACPY('All',nh,2*nnb,Work(pw),nh,Q(topq,j),  &
     &                           Ldq)
                  ENDIF
                  ppwo = ppwo + 4*nnb*nnb
               ENDDO
            ENDIF
!
!           Accumulate right Givens rotations if required.
!
            IF ( wantz .OR. top>0 ) THEN
!
!              Initialize small unitary factors that will hold the
!              accumulated Givens rotations in workspace.
!
               CALL CLASET('All',nblst,nblst,CZERO,CONE,Work,nblst)
               pw = nblst*nblst + 1
               DO i = 1 , n2nb
                  CALL CLASET('All',2*nnb,2*nnb,CZERO,CONE,Work(pw),    &
     &                        2*nnb)
                  pw = pw + 4*nnb*nnb
               ENDDO
!
!              Accumulate Givens rotations into workspace array.
!
               DO j = jcol , jcol + nnb - 1
                  ppw = (nblst+1)*(nblst-2) - j + jcol + 1
                  len = 2 + j - jcol
                  jrow = j + n2nb*nnb + 2
                  DO i = Ihi , jrow , -1
                     ctemp = A(i,j)
                     A(i,j) = CZERO
                     s = B(i,j)
                     B(i,j) = CZERO
                     DO jj = ppw , ppw + len - 1
                        temp = Work(jj+nblst)
                        Work(jj+nblst) = ctemp*temp - CONJG(s)*Work(jj)
                        Work(jj) = s*temp + ctemp*Work(jj)
                     ENDDO
                     len = len + 1
                     ppw = ppw - nblst - 1
                  ENDDO
!
                  ppwo = nblst*nblst + (nnb+j-jcol-1)*2*nnb + nnb
                  j0 = jrow - nnb
                  DO jrow = j0 , j + 2 , -nnb
                     ppw = ppwo
                     len = 2 + j - jcol
                     DO i = jrow + nnb - 1 , jrow , -1
                        ctemp = A(i,j)
                        A(i,j) = CZERO
                        s = B(i,j)
                        B(i,j) = CZERO
                        DO jj = ppw , ppw + len - 1
                           temp = Work(jj+2*nnb)
                           Work(jj+2*nnb) = ctemp*temp - CONJG(s)       &
     &                        *Work(jj)
                           Work(jj) = s*temp + ctemp*Work(jj)
                        ENDDO
                        len = len + 1
                        ppw = ppw - 2*nnb - 1
                     ENDDO
                     ppwo = ppwo + 4*nnb*nnb
                  ENDDO
               ENDDO
            ELSE
!
               CALL CLASET('Lower',Ihi-jcol-1,nnb,CZERO,CZERO,          &
     &                     A(jcol+2,jcol),Lda)
               CALL CLASET('Lower',Ihi-jcol-1,nnb,CZERO,CZERO,          &
     &                     B(jcol+2,jcol),Ldb)
            ENDIF
!
!           Apply accumulated unitary matrices to A and B.
!
            IF ( top>0 ) THEN
               j = Ihi - nblst + 1
               CALL CGEMM('No Transpose','No Transpose',top,nblst,nblst,&
     &                    CONE,A(1,j),Lda,Work,nblst,CZERO,Work(pw),top)
               CALL CLACPY('All',top,nblst,Work(pw),top,A(1,j),Lda)
               ppwo = nblst*nblst + 1
               j0 = j - nnb
               DO j = j0 , jcol + 1 , -nnb
                  IF ( blk22 ) THEN
!
!                    Exploit the structure of U.
!
                     CALL CUNM22('Right','No Transpose',top,2*nnb,nnb,  &
     &                           nnb,Work(ppwo),2*nnb,A(1,j),Lda,       &
     &                           Work(pw),Lwork-pw+1,ierr)
                  ELSE
!
!                    Ignore the structure of U.
!
                     CALL CGEMM('No Transpose','No Transpose',top,2*nnb,&
     &                          2*nnb,CONE,A(1,j),Lda,Work(ppwo),2*nnb, &
     &                          CZERO,Work(pw),top)
                     CALL CLACPY('All',top,2*nnb,Work(pw),top,A(1,j),   &
     &                           Lda)
                  ENDIF
                  ppwo = ppwo + 4*nnb*nnb
               ENDDO
!
               j = Ihi - nblst + 1
               CALL CGEMM('No Transpose','No Transpose',top,nblst,nblst,&
     &                    CONE,B(1,j),Ldb,Work,nblst,CZERO,Work(pw),top)
               CALL CLACPY('All',top,nblst,Work(pw),top,B(1,j),Ldb)
               ppwo = nblst*nblst + 1
               j0 = j - nnb
               DO j = j0 , jcol + 1 , -nnb
                  IF ( blk22 ) THEN
!
!                    Exploit the structure of U.
!
                     CALL CUNM22('Right','No Transpose',top,2*nnb,nnb,  &
     &                           nnb,Work(ppwo),2*nnb,B(1,j),Ldb,       &
     &                           Work(pw),Lwork-pw+1,ierr)
                  ELSE
!
!                    Ignore the structure of U.
!
                     CALL CGEMM('No Transpose','No Transpose',top,2*nnb,&
     &                          2*nnb,CONE,B(1,j),Ldb,Work(ppwo),2*nnb, &
     &                          CZERO,Work(pw),top)
                     CALL CLACPY('All',top,2*nnb,Work(pw),top,B(1,j),   &
     &                           Ldb)
                  ENDIF
                  ppwo = ppwo + 4*nnb*nnb
               ENDDO
            ENDIF
!
!           Apply accumulated unitary matrices to Z.
!
            IF ( wantz ) THEN
               j = Ihi - nblst + 1
               IF ( initq ) THEN
                  topq = MAX(2,j-jcol+1)
                  nh = Ihi - topq + 1
               ELSE
                  topq = 1
                  nh = N
               ENDIF
               CALL CGEMM('No Transpose','No Transpose',nh,nblst,nblst, &
     &                    CONE,Z(topq,j),Ldz,Work,nblst,CZERO,Work(pw), &
     &                    nh)
               CALL CLACPY('All',nh,nblst,Work(pw),nh,Z(topq,j),Ldz)
               ppwo = nblst*nblst + 1
               j0 = j - nnb
               DO j = j0 , jcol + 1 , -nnb
                  IF ( initq ) THEN
                     topq = MAX(2,j-jcol+1)
                     nh = Ihi - topq + 1
                  ENDIF
                  IF ( blk22 ) THEN
!
!                    Exploit the structure of U.
!
                     CALL CUNM22('Right','No Transpose',nh,2*nnb,nnb,   &
     &                           nnb,Work(ppwo),2*nnb,Z(topq,j),Ldz,    &
     &                           Work(pw),Lwork-pw+1,ierr)
                  ELSE
!
!                    Ignore the structure of U.
!
                     CALL CGEMM('No Transpose','No Transpose',nh,2*nnb, &
     &                          2*nnb,CONE,Z(topq,j),Ldz,Work(ppwo),    &
     &                          2*nnb,CZERO,Work(pw),nh)
                     CALL CLACPY('All',nh,2*nnb,Work(pw),nh,Z(topq,j),  &
     &                           Ldz)
                  ENDIF
                  ppwo = ppwo + 4*nnb*nnb
               ENDDO
            ENDIF
         ENDDO
      ENDIF
!
!     Use unblocked code to reduce the rest of the matrix
!     Avoid re-initialization of modified Q and Z.
!
      compq2 = Compq
      compz2 = Compz
      IF ( jcol/=Ilo ) THEN
         IF ( wantq ) compq2 = 'V'
         IF ( wantz ) compz2 = 'V'
      ENDIF
!
      IF ( jcol<Ihi ) CALL CGGHRD(compq2,compz2,N,jcol,Ihi,A,Lda,B,Ldb, &
     &                            Q,Ldq,Z,Ldz,ierr)
      Work(1) = CMPLX(lwkopt)
!
!
!     End of CGGHD3
!
      END SUBROUTINE CGGHD3
