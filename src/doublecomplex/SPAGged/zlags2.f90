!*==zlags2.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b ZLAGS2
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZLAGS2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlags2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlags2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlags2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZLAGS2( UPPER, A1, A2, A3, B1, B2, B3, CSU, SNU, CSV,
!                          SNV, CSQ, SNQ )
!
!       .. Scalar Arguments ..
!       LOGICAL            UPPER
!       DOUBLE PRECISION   A1, A3, B1, B3, CSQ, CSU, CSV
!       COMPLEX*16         A2, B2, SNQ, SNU, SNV
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZLAGS2 computes 2-by-2 unitary matrices U, V and Q, such
!> that if ( UPPER ) then
!>
!>           U**H *A*Q = U**H *( A1 A2 )*Q = ( x  0  )
!>                             ( 0  A3 )     ( x  x  )
!> and
!>           V**H*B*Q = V**H *( B1 B2 )*Q = ( x  0  )
!>                            ( 0  B3 )     ( x  x  )
!>
!> or if ( .NOT.UPPER ) then
!>
!>           U**H *A*Q = U**H *( A1 0  )*Q = ( x  x  )
!>                             ( A2 A3 )     ( 0  x  )
!> and
!>           V**H *B*Q = V**H *( B1 0  )*Q = ( x  x  )
!>                             ( B2 B3 )     ( 0  x  )
!> where
!>
!>   U = (   CSU    SNU ), V = (  CSV    SNV ),
!>       ( -SNU**H  CSU )      ( -SNV**H CSV )
!>
!>   Q = (   CSQ    SNQ )
!>       ( -SNQ**H  CSQ )
!>
!> The rows of the transformed A and B are parallel. Moreover, if the
!> input 2-by-2 matrix A is not zero, then the transformed (1,1) entry
!> of A is not zero. If the input matrices A and B are both not zero,
!> then the transformed (2,2) element of B is not zero, except when the
!> first rows of input A and B are parallel and the second rows are
!> zero.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPPER
!> \verbatim
!>          UPPER is LOGICAL
!>          = .TRUE.: the input matrices A and B are upper triangular.
!>          = .FALSE.: the input matrices A and B are lower triangular.
!> \endverbatim
!>
!> \param[in] A1
!> \verbatim
!>          A1 is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[in] A2
!> \verbatim
!>          A2 is COMPLEX*16
!> \endverbatim
!>
!> \param[in] A3
!> \verbatim
!>          A3 is DOUBLE PRECISION
!>          On entry, A1, A2 and A3 are elements of the input 2-by-2
!>          upper (lower) triangular matrix A.
!> \endverbatim
!>
!> \param[in] B1
!> \verbatim
!>          B1 is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[in] B2
!> \verbatim
!>          B2 is COMPLEX*16
!> \endverbatim
!>
!> \param[in] B3
!> \verbatim
!>          B3 is DOUBLE PRECISION
!>          On entry, B1, B2 and B3 are elements of the input 2-by-2
!>          upper (lower) triangular matrix B.
!> \endverbatim
!>
!> \param[out] CSU
!> \verbatim
!>          CSU is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[out] SNU
!> \verbatim
!>          SNU is COMPLEX*16
!>          The desired unitary matrix U.
!> \endverbatim
!>
!> \param[out] CSV
!> \verbatim
!>          CSV is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[out] SNV
!> \verbatim
!>          SNV is COMPLEX*16
!>          The desired unitary matrix V.
!> \endverbatim
!>
!> \param[out] CSQ
!> \verbatim
!>          CSQ is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[out] SNQ
!> \verbatim
!>          SNQ is COMPLEX*16
!>          The desired unitary matrix Q.
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
!> \ingroup complex16OTHERauxiliary
!
!  =====================================================================
      SUBROUTINE ZLAGS2(Upper,A1,A2,A3,B1,B2,B3,Csu,Snu,Csv,Snv,Csq,Snq)
      USE F77KINDS                        
      USE S_DLASV2
      USE S_ZLARTG
      IMPLICIT NONE
!*--ZLAGS2164
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
!
! Dummy argument declarations rewritten by SPAG
!
      LOGICAL , INTENT(IN) :: Upper
      REAL(R8KIND) , INTENT(IN) :: A1
      COMPLEX(CX16KIND) , INTENT(IN) :: A2
      REAL(R8KIND) , INTENT(IN) :: A3
      REAL(R8KIND) , INTENT(IN) :: B1
      COMPLEX(CX16KIND) , INTENT(IN) :: B2
      REAL(R8KIND) , INTENT(IN) :: B3
      REAL(R8KIND) , INTENT(OUT) :: Csu
      COMPLEX(CX16KIND) , INTENT(OUT) :: Snu
      REAL(R8KIND) , INTENT(OUT) :: Csv
      COMPLEX(CX16KIND) , INTENT(OUT) :: Snv
      REAL(R8KIND) :: Csq
      COMPLEX(CX16KIND) :: Snq
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: a , aua11 , aua12 , aua21 , aua22 , avb11 ,       &
     &                avb12 , avb21 , avb22 , csl , csr , d , fb , fc , &
     &                s1 , s2 , snl , snr , ua11r , ua22r , vb11r ,     &
     &                vb22r
      REAL(R8KIND) :: ABS1
      COMPLEX(CX16KIND) :: b , c , d1 , r , t , ua11 , ua12 , ua21 ,    &
     &                     ua22 , vb11 , vb12 , vb21 , vb22
!
! End of declarations rewritten by SPAG
!
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
!     .. Statement Functions ..
!     ..
!     .. Statement Function definitions ..
      ABS1(t) = ABS(DBLE(t)) + ABS(DIMAG(t))
!     ..
!     .. Executable Statements ..
!
      IF ( Upper ) THEN
!
!        Input matrices A and B are upper triangular matrices
!
!        Form matrix C = A*adj(B) = ( a b )
!                                   ( 0 d )
!
         a = A1*B3
         d = A3*B1
         b = A2*B1 - A1*B2
         fb = ABS(b)
!
!        Transform complex 2-by-2 matrix C to real matrix by unitary
!        diagonal matrix diag(1,D1).
!
         d1 = ONE
         IF ( fb/=ZERO ) d1 = b/fb
!
!        The SVD of real 2 by 2 triangular C
!
!         ( CSL -SNL )*( A B )*(  CSR  SNR ) = ( R 0 )
!         ( SNL  CSL ) ( 0 D ) ( -SNR  CSR )   ( 0 T )
!
         CALL DLASV2(a,fb,d,s1,s2,snr,csr,snl,csl)
!
         IF ( ABS(csl)>=ABS(snl) .OR. ABS(csr)>=ABS(snr) ) THEN
!
!           Compute the (1,1) and (1,2) elements of U**H *A and V**H *B,
!           and (1,2) element of |U|**H *|A| and |V|**H *|B|.
!
            ua11r = csl*A1
            ua12 = csl*A2 + d1*snl*A3
!
            vb11r = csr*B1
            vb12 = csr*B2 + d1*snr*B3
!
            aua12 = ABS(csl)*ABS1(A2) + ABS(snl)*ABS(A3)
            avb12 = ABS(csr)*ABS1(B2) + ABS(snr)*ABS(B3)
!
!           zero (1,2) elements of U**H *A and V**H *B
!
            IF ( (ABS(ua11r)+ABS1(ua12))==ZERO ) THEN
               CALL ZLARTG(-DCMPLX(vb11r),DCONJG(vb12),Csq,Snq,r)
            ELSEIF ( (ABS(vb11r)+ABS1(vb12))==ZERO ) THEN
               CALL ZLARTG(-DCMPLX(ua11r),DCONJG(ua12),Csq,Snq,r)
            ELSEIF ( aua12/(ABS(ua11r)+ABS1(ua12))                      &
     &               <=avb12/(ABS(vb11r)+ABS1(vb12)) ) THEN
               CALL ZLARTG(-DCMPLX(ua11r),DCONJG(ua12),Csq,Snq,r)
            ELSE
               CALL ZLARTG(-DCMPLX(vb11r),DCONJG(vb12),Csq,Snq,r)
            ENDIF
!
            Csu = csl
            Snu = -d1*snl
            Csv = csr
            Snv = -d1*snr
!
         ELSE
!
!           Compute the (2,1) and (2,2) elements of U**H *A and V**H *B,
!           and (2,2) element of |U|**H *|A| and |V|**H *|B|.
!
            ua21 = -DCONJG(d1)*snl*A1
            ua22 = -DCONJG(d1)*snl*A2 + csl*A3
!
            vb21 = -DCONJG(d1)*snr*B1
            vb22 = -DCONJG(d1)*snr*B2 + csr*B3
!
            aua22 = ABS(snl)*ABS1(A2) + ABS(csl)*ABS(A3)
            avb22 = ABS(snr)*ABS1(B2) + ABS(csr)*ABS(B3)
!
!           zero (2,2) elements of U**H *A and V**H *B, and then swap.
!
            IF ( (ABS1(ua21)+ABS1(ua22))==ZERO ) THEN
               CALL ZLARTG(-DCONJG(vb21),DCONJG(vb22),Csq,Snq,r)
            ELSEIF ( (ABS1(vb21)+ABS(vb22))==ZERO ) THEN
               CALL ZLARTG(-DCONJG(ua21),DCONJG(ua22),Csq,Snq,r)
            ELSEIF ( aua22/(ABS1(ua21)+ABS1(ua22))                      &
     &               <=avb22/(ABS1(vb21)+ABS1(vb22)) ) THEN
               CALL ZLARTG(-DCONJG(ua21),DCONJG(ua22),Csq,Snq,r)
            ELSE
               CALL ZLARTG(-DCONJG(vb21),DCONJG(vb22),Csq,Snq,r)
            ENDIF
!
            Csu = snl
            Snu = d1*csl
            Csv = snr
            Snv = d1*csr
!
         ENDIF
!
      ELSE
!
!        Input matrices A and B are lower triangular matrices
!
!        Form matrix C = A*adj(B) = ( a 0 )
!                                   ( c d )
!
         a = A1*B3
         d = A3*B1
         c = A2*B3 - A3*B2
         fc = ABS(c)
!
!        Transform complex 2-by-2 matrix C to real matrix by unitary
!        diagonal matrix diag(d1,1).
!
         d1 = ONE
         IF ( fc/=ZERO ) d1 = c/fc
!
!        The SVD of real 2 by 2 triangular C
!
!         ( CSL -SNL )*( A 0 )*(  CSR  SNR ) = ( R 0 )
!         ( SNL  CSL ) ( C D ) ( -SNR  CSR )   ( 0 T )
!
         CALL DLASV2(a,fc,d,s1,s2,snr,csr,snl,csl)
!
         IF ( ABS(csr)>=ABS(snr) .OR. ABS(csl)>=ABS(snl) ) THEN
!
!           Compute the (2,1) and (2,2) elements of U**H *A and V**H *B,
!           and (2,1) element of |U|**H *|A| and |V|**H *|B|.
!
            ua21 = -d1*snr*A1 + csr*A2
            ua22r = csr*A3
!
            vb21 = -d1*snl*B1 + csl*B2
            vb22r = csl*B3
!
            aua21 = ABS(snr)*ABS(A1) + ABS(csr)*ABS1(A2)
            avb21 = ABS(snl)*ABS(B1) + ABS(csl)*ABS1(B2)
!
!           zero (2,1) elements of U**H *A and V**H *B.
!
            IF ( (ABS1(ua21)+ABS(ua22r))==ZERO ) THEN
               CALL ZLARTG(DCMPLX(vb22r),vb21,Csq,Snq,r)
            ELSEIF ( (ABS1(vb21)+ABS(vb22r))==ZERO ) THEN
               CALL ZLARTG(DCMPLX(ua22r),ua21,Csq,Snq,r)
            ELSEIF ( aua21/(ABS1(ua21)+ABS(ua22r))                      &
     &               <=avb21/(ABS1(vb21)+ABS(vb22r)) ) THEN
               CALL ZLARTG(DCMPLX(ua22r),ua21,Csq,Snq,r)
            ELSE
               CALL ZLARTG(DCMPLX(vb22r),vb21,Csq,Snq,r)
            ENDIF
!
            Csu = csr
            Snu = -DCONJG(d1)*snr
            Csv = csl
            Snv = -DCONJG(d1)*snl
!
         ELSE
!
!           Compute the (1,1) and (1,2) elements of U**H *A and V**H *B,
!           and (1,1) element of |U|**H *|A| and |V|**H *|B|.
!
            ua11 = csr*A1 + DCONJG(d1)*snr*A2
            ua12 = DCONJG(d1)*snr*A3
!
            vb11 = csl*B1 + DCONJG(d1)*snl*B2
            vb12 = DCONJG(d1)*snl*B3
!
            aua11 = ABS(csr)*ABS(A1) + ABS(snr)*ABS1(A2)
            avb11 = ABS(csl)*ABS(B1) + ABS(snl)*ABS1(B2)
!
!           zero (1,1) elements of U**H *A and V**H *B, and then swap.
!
            IF ( (ABS1(ua11)+ABS1(ua12))==ZERO ) THEN
               CALL ZLARTG(vb12,vb11,Csq,Snq,r)
            ELSEIF ( (ABS1(vb11)+ABS1(vb12))==ZERO ) THEN
               CALL ZLARTG(ua12,ua11,Csq,Snq,r)
            ELSEIF ( aua11/(ABS1(ua11)+ABS1(ua12))                      &
     &               <=avb11/(ABS1(vb11)+ABS1(vb12)) ) THEN
               CALL ZLARTG(ua12,ua11,Csq,Snq,r)
            ELSE
               CALL ZLARTG(vb12,vb11,Csq,Snq,r)
            ENDIF
!
            Csu = snr
            Snu = DCONJG(d1)*csr
            Csv = snl
            Snv = DCONJG(d1)*csl
!
         ENDIF
!
      ENDIF
!
!
!     End of ZLAGS2
!
      END SUBROUTINE ZLAGS2
