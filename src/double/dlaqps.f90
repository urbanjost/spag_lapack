!*==dlaqps.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b DLAQPS computes a step of QR factorization with column pivoting of a real m-by-n matrix A by using BLAS level 3.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLAQPS + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaqps.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaqps.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaqps.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLAQPS( M, N, OFFSET, NB, KB, A, LDA, JPVT, TAU, VN1,
!                          VN2, AUXV, F, LDF )
!
!       .. Scalar Arguments ..
!       INTEGER            KB, LDA, LDF, M, N, NB, OFFSET
!       ..
!       .. Array Arguments ..
!       INTEGER            JPVT( * )
!       DOUBLE PRECISION   A( LDA, * ), AUXV( * ), F( LDF, * ), TAU( * ),
!      $                   VN1( * ), VN2( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLAQPS computes a step of QR factorization with column pivoting
!> of a real M-by-N matrix A by using Blas-3.  It tries to factorize
!> NB columns from A starting from the row OFFSET+1, and updates all
!> of the matrix with Blas-3 xGEMM.
!>
!> In some cases, due to catastrophic cancellations, it cannot
!> factorize NB columns.  Hence, the actual number of factorized
!> columns is returned in KB.
!>
!> Block A(1:OFFSET,1:N) is accordingly pivoted, but not factorized.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A. M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A. N >= 0
!> \endverbatim
!>
!> \param[in] OFFSET
!> \verbatim
!>          OFFSET is INTEGER
!>          The number of rows of A that have been factorized in
!>          previous steps.
!> \endverbatim
!>
!> \param[in] NB
!> \verbatim
!>          NB is INTEGER
!>          The number of columns to factorize.
!> \endverbatim
!>
!> \param[out] KB
!> \verbatim
!>          KB is INTEGER
!>          The number of columns actually factorized.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          On entry, the M-by-N matrix A.
!>          On exit, block A(OFFSET+1:M,1:KB) is the triangular
!>          factor obtained and block A(1:OFFSET,1:N) has been
!>          accordingly pivoted, but no factorized.
!>          The rest of the matrix, block A(OFFSET+1:M,KB+1:N) has
!>          been updated.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A. LDA >= max(1,M).
!> \endverbatim
!>
!> \param[in,out] JPVT
!> \verbatim
!>          JPVT is INTEGER array, dimension (N)
!>          JPVT(I) = K <==> Column K of the full matrix A has been
!>          permuted into position I in AP.
!> \endverbatim
!>
!> \param[out] TAU
!> \verbatim
!>          TAU is DOUBLE PRECISION array, dimension (KB)
!>          The scalar factors of the elementary reflectors.
!> \endverbatim
!>
!> \param[in,out] VN1
!> \verbatim
!>          VN1 is DOUBLE PRECISION array, dimension (N)
!>          The vector with the partial column norms.
!> \endverbatim
!>
!> \param[in,out] VN2
!> \verbatim
!>          VN2 is DOUBLE PRECISION array, dimension (N)
!>          The vector with the exact column norms.
!> \endverbatim
!>
!> \param[in,out] AUXV
!> \verbatim
!>          AUXV is DOUBLE PRECISION array, dimension (NB)
!>          Auxiliary vector.
!> \endverbatim
!>
!> \param[in,out] F
!> \verbatim
!>          F is DOUBLE PRECISION array, dimension (LDF,NB)
!>          Matrix F**T = L*Y**T*A.
!> \endverbatim
!>
!> \param[in] LDF
!> \verbatim
!>          LDF is INTEGER
!>          The leading dimension of the array F. LDF >= max(1,N).
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
!> \par Contributors:
!  ==================
!>
!>    G. Quintana-Orti, Depto. de Informatica, Universidad Jaime I, Spain
!>    X. Sun, Computer Science Dept., Duke University, USA
!> \n
!>  Partial column norm updating strategy modified on April 2011
!>    Z. Drmac and Z. Bujanovic, Dept. of Mathematics,
!>    University of Zagreb, Croatia.
!
!> \par References:
!  ================
!>
!> LAPACK Working Note 176
!
!> \htmlonly
!> <a href="http://www.netlib.org/lapack/lawnspdf/lawn176.pdf">[PDF]</a>
!> \endhtmlonly
!
!  =====================================================================
      SUBROUTINE DLAQPS(M,N,Offset,Nb,Kb,A,Lda,Jpvt,Tau,Vn1,Vn2,Auxv,F, &
     &                  Ldf)
      IMPLICIT NONE
!*--DLAQPS181
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Kb , Lda , Ldf , M , N , Nb , Offset
!     ..
!     .. Array Arguments ..
      INTEGER Jpvt(*)
      DOUBLE PRECISION A(Lda,*) , Auxv(*) , F(Ldf,*) , Tau(*) , Vn1(*) ,&
     &                 Vn2(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D+0,ONE=1.0D+0)
!     ..
!     .. Local Scalars ..
      INTEGER itemp , j , k , lastrk , lsticc , pvt , rk
      DOUBLE PRECISION akk , temp , temp2 , tol3z
!     ..
!     .. External Subroutines ..
      EXTERNAL DGEMM , DGEMV , DLARFG , DSWAP
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , DBLE , MAX , MIN , NINT , SQRT
!     ..
!     .. External Functions ..
      INTEGER IDAMAX
      DOUBLE PRECISION DLAMCH , DNRM2
      EXTERNAL IDAMAX , DLAMCH , DNRM2
!     ..
!     .. Executable Statements ..
!
      lastrk = MIN(M,N+Offset)
      lsticc = 0
      k = 0
      tol3z = SQRT(DLAMCH('Epsilon'))
      DO
!
!     Beginning of while loop.
!
         IF ( (k<Nb) .AND. (lsticc==0) ) THEN
            k = k + 1
            rk = Offset + k
!
!        Determine ith pivot column and swap if necessary
!
            pvt = (k-1) + IDAMAX(N-k+1,Vn1(k),1)
            IF ( pvt/=k ) THEN
               CALL DSWAP(M,A(1,pvt),1,A(1,k),1)
               CALL DSWAP(k-1,F(pvt,1),Ldf,F(k,1),Ldf)
               itemp = Jpvt(pvt)
               Jpvt(pvt) = Jpvt(k)
               Jpvt(k) = itemp
               Vn1(pvt) = Vn1(k)
               Vn2(pvt) = Vn2(k)
            ENDIF
!
!        Apply previous Householder reflectors to column K:
!        A(RK:M,K) := A(RK:M,K) - A(RK:M,1:K-1)*F(K,1:K-1)**T.
!
            IF ( k>1 ) CALL DGEMV('No transpose',M-rk+1,k-1,-ONE,A(rk,1)&
     &                            ,Lda,F(k,1),Ldf,ONE,A(rk,k),1)
!
!        Generate elementary reflector H(k).
!
            IF ( rk<M ) THEN
               CALL DLARFG(M-rk+1,A(rk,k),A(rk+1,k),1,Tau(k))
            ELSE
               CALL DLARFG(1,A(rk,k),A(rk,k),1,Tau(k))
            ENDIF
!
            akk = A(rk,k)
            A(rk,k) = ONE
!
!        Compute Kth column of F:
!
!        Compute  F(K+1:N,K) := tau(K)*A(RK:M,K+1:N)**T*A(RK:M,K).
!
            IF ( k<N ) CALL DGEMV('Transpose',M-rk+1,N-k,Tau(k),        &
     &                            A(rk,k+1),Lda,A(rk,k),1,ZERO,F(k+1,k),&
     &                            1)
!
!        Padding F(1:K,K) with zeros.
!
            DO j = 1 , k
               F(j,k) = ZERO
            ENDDO
!
!        Incremental updating of F:
!        F(1:N,K) := F(1:N,K) - tau(K)*F(1:N,1:K-1)*A(RK:M,1:K-1)**T
!                    *A(RK:M,K).
!
            IF ( k>1 ) THEN
               CALL DGEMV('Transpose',M-rk+1,k-1,-Tau(k),A(rk,1),Lda,   &
     &                    A(rk,k),1,ZERO,Auxv(1),1)
!
               CALL DGEMV('No transpose',N,k-1,ONE,F(1,1),Ldf,Auxv(1),1,&
     &                    ONE,F(1,k),1)
            ENDIF
!
!        Update the current row of A:
!        A(RK,K+1:N) := A(RK,K+1:N) - A(RK,1:K)*F(K+1:N,1:K)**T.
!
            IF ( k<N ) CALL DGEMV('No transpose',N-k,k,-ONE,F(k+1,1),   &
     &                            Ldf,A(rk,1),Lda,ONE,A(rk,k+1),Lda)
!
!        Update partial column norms.
!
            IF ( rk<lastrk ) THEN
               DO j = k + 1 , N
                  IF ( Vn1(j)/=ZERO ) THEN
!
!                 NOTE: The following 4 lines follow from the analysis in
!                 Lapack Working Note 176.
!
                     temp = ABS(A(rk,j))/Vn1(j)
                     temp = MAX(ZERO,(ONE+temp)*(ONE-temp))
                     temp2 = temp*(Vn1(j)/Vn2(j))**2
                     IF ( temp2<=tol3z ) THEN
                        Vn2(j) = DBLE(lsticc)
                        lsticc = j
                     ELSE
                        Vn1(j) = Vn1(j)*SQRT(temp)
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF
!
            A(rk,k) = akk
!
!        End of while loop.
!
            CYCLE
         ENDIF
         Kb = k
         rk = Offset + Kb
!
!     Apply the block reflector to the rest of the matrix:
!     A(OFFSET+KB+1:M,KB+1:N) := A(OFFSET+KB+1:M,KB+1:N) -
!                         A(OFFSET+KB+1:M,1:KB)*F(KB+1:N,1:KB)**T.
!
         IF ( Kb<MIN(N,M-Offset) )                                      &
     &        CALL DGEMM('No transpose','Transpose',M-rk,N-Kb,Kb,-ONE,  &
     &        A(rk+1,1),Lda,F(Kb+1,1),Ldf,ONE,A(rk+1,Kb+1),Lda)
         EXIT
      ENDDO
!
!     Recomputation of difficult columns.
!
      DO WHILE ( lsticc>0 )
         itemp = NINT(Vn2(lsticc))
         Vn1(lsticc) = DNRM2(M-rk,A(rk+1,lsticc),1)
!
!        NOTE: The computation of VN1( LSTICC ) relies on the fact that
!        SNRM2 does not fail on vectors with norm below the value of
!        SQRT(DLAMCH('S'))
!
         Vn2(lsticc) = Vn1(lsticc)
         lsticc = itemp
      ENDDO
!
!
!     End of DLAQPS
!
      END SUBROUTINE DLAQPS
