!*==claqps.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b CLAQPS computes a step of QR factorization with column pivoting of a real m-by-n matrix A by using BLAS level 3.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CLAQPS + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/claqps.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/claqps.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/claqps.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CLAQPS( M, N, OFFSET, NB, KB, A, LDA, JPVT, TAU, VN1,
!                          VN2, AUXV, F, LDF )
!
!       .. Scalar Arguments ..
!       INTEGER            KB, LDA, LDF, M, N, NB, OFFSET
!       ..
!       .. Array Arguments ..
!       INTEGER            JPVT( * )
!       REAL               VN1( * ), VN2( * )
!       COMPLEX            A( LDA, * ), AUXV( * ), F( LDF, * ), TAU( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CLAQPS computes a step of QR factorization with column pivoting
!> of a complex M-by-N matrix A by using Blas-3.  It tries to factorize
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
!>          A is COMPLEX array, dimension (LDA,N)
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
!>          TAU is COMPLEX array, dimension (KB)
!>          The scalar factors of the elementary reflectors.
!> \endverbatim
!>
!> \param[in,out] VN1
!> \verbatim
!>          VN1 is REAL array, dimension (N)
!>          The vector with the partial column norms.
!> \endverbatim
!>
!> \param[in,out] VN2
!> \verbatim
!>          VN2 is REAL array, dimension (N)
!>          The vector with the exact column norms.
!> \endverbatim
!>
!> \param[in,out] AUXV
!> \verbatim
!>          AUXV is COMPLEX array, dimension (NB)
!>          Auxiliary vector.
!> \endverbatim
!>
!> \param[in,out] F
!> \verbatim
!>          F is COMPLEX array, dimension (LDF,NB)
!>          Matrix  F**H = L * Y**H * A.
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
!> \ingroup complexOTHERauxiliary
!
!> \par Contributors:
!  ==================
!>
!>    G. Quintana-Orti, Depto. de Informatica, Universidad Jaime I, Spain
!>    X. Sun, Computer Science Dept., Duke University, USA
!>
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
      SUBROUTINE CLAQPS(M,N,Offset,Nb,Kb,A,Lda,Jpvt,Tau,Vn1,Vn2,Auxv,F, &
     &                  Ldf)
      USE S_CGEMM
      USE S_CGEMV
      USE S_CLARFG
      USE S_CSWAP
      USE S_ISAMAX
      USE S_SCNRM2
      USE S_SLAMCH
      IMPLICIT NONE
!*--CLAQPS189
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      COMPLEX , PARAMETER  ::  CZERO = (0.0E+0,0.0E+0) ,                &
     &                         CONE = (1.0E+0,0.0E+0)
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER :: M
      INTEGER :: N
      INTEGER , INTENT(IN) :: Offset
      INTEGER , INTENT(IN) :: Nb
      INTEGER , INTENT(INOUT) :: Kb
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Jpvt
      COMPLEX , DIMENSION(*) :: Tau
      REAL , INTENT(INOUT) , DIMENSION(*) :: Vn1
      REAL , INTENT(INOUT) , DIMENSION(*) :: Vn2
      COMPLEX , DIMENSION(*) :: Auxv
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldf,*) :: F
      INTEGER :: Ldf
!
! Local variable declarations rewritten by SPAG
!
      COMPLEX :: akk
      INTEGER :: itemp , j , k , lastrk , lsticc , pvt , rk
      REAL :: temp , temp2 , tol3z
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
!     .. External Functions ..
!     ..
!     .. Executable Statements ..
!
      lastrk = MIN(M,N+Offset)
      lsticc = 0
      k = 0
      tol3z = SQRT(SLAMCH('Epsilon'))
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
            pvt = (k-1) + ISAMAX(N-k+1,Vn1(k),1)
            IF ( pvt/=k ) THEN
               CALL CSWAP(M,A(1,pvt),1,A(1,k),1)
               CALL CSWAP(k-1,F(pvt,1),Ldf,F(k,1),Ldf)
               itemp = Jpvt(pvt)
               Jpvt(pvt) = Jpvt(k)
               Jpvt(k) = itemp
               Vn1(pvt) = Vn1(k)
               Vn2(pvt) = Vn2(k)
            ENDIF
!
!        Apply previous Householder reflectors to column K:
!        A(RK:M,K) := A(RK:M,K) - A(RK:M,1:K-1)*F(K,1:K-1)**H.
!
            IF ( k>1 ) THEN
               DO j = 1 , k - 1
                  F(k,j) = CONJG(F(k,j))
               ENDDO
               CALL CGEMV('No transpose',M-rk+1,k-1,-CONE,A(rk,1),Lda,  &
     &                    F(k,1),Ldf,CONE,A(rk,k),1)
               DO j = 1 , k - 1
                  F(k,j) = CONJG(F(k,j))
               ENDDO
            ENDIF
!
!        Generate elementary reflector H(k).
!
            IF ( rk<M ) THEN
               CALL CLARFG(M-rk+1,A(rk,k),A(rk+1,k),1,Tau(k))
            ELSE
               CALL CLARFG(1,A(rk,k),A(rk,k),1,Tau(k))
            ENDIF
!
            akk = A(rk,k)
            A(rk,k) = CONE
!
!        Compute Kth column of F:
!
!        Compute  F(K+1:N,K) := tau(K)*A(RK:M,K+1:N)**H*A(RK:M,K).
!
            IF ( k<N ) CALL CGEMV('Conjugate transpose',M-rk+1,N-k,     &
     &                            Tau(k),A(rk,k+1),Lda,A(rk,k),1,CZERO, &
     &                            F(k+1,k),1)
!
!        Padding F(1:K,K) with zeros.
!
            DO j = 1 , k
               F(j,k) = CZERO
            ENDDO
!
!        Incremental updating of F:
!        F(1:N,K) := F(1:N,K) - tau(K)*F(1:N,1:K-1)*A(RK:M,1:K-1)**H
!                    *A(RK:M,K).
!
            IF ( k>1 ) THEN
               CALL CGEMV('Conjugate transpose',M-rk+1,k-1,-Tau(k),     &
     &                    A(rk,1),Lda,A(rk,k),1,CZERO,Auxv(1),1)
!
               CALL CGEMV('No transpose',N,k-1,CONE,F(1,1),Ldf,Auxv(1), &
     &                    1,CONE,F(1,k),1)
            ENDIF
!
!        Update the current row of A:
!        A(RK,K+1:N) := A(RK,K+1:N) - A(RK,1:K)*F(K+1:N,1:K)**H.
!
            IF ( k<N ) CALL CGEMM('No transpose','Conjugate transpose', &
     &                            1,N-k,k,-CONE,A(rk,1),Lda,F(k+1,1),   &
     &                            Ldf,CONE,A(rk,k+1),Lda)
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
                        Vn2(j) = REAL(lsticc)
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
!                         A(OFFSET+KB+1:M,1:KB)*F(KB+1:N,1:KB)**H.
!
         IF ( Kb<MIN(N,M-Offset) )                                      &
     &         CALL CGEMM('No transpose','Conjugate transpose',M-rk,    &
     &        N-Kb,Kb,-CONE,A(rk+1,1),Lda,F(Kb+1,1),Ldf,CONE,           &
     &        A(rk+1,Kb+1),Lda)
         EXIT
      ENDDO
!
!     Recomputation of difficult columns.
!
      DO WHILE ( lsticc>0 )
         itemp = NINT(Vn2(lsticc))
         Vn1(lsticc) = SCNRM2(M-rk,A(rk+1,lsticc),1)
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
!     End of CLAQPS
!
      END SUBROUTINE CLAQPS
