!*==slags2.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b SLAGS2 computes 2-by-2 orthogonal matrices U, V, and Q, and applies them to matrices A and B such that the rows of the transformed A and B are parallel.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SLAGS2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slags2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slags2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slags2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLAGS2( UPPER, A1, A2, A3, B1, B2, B3, CSU, SNU, CSV,
!                          SNV, CSQ, SNQ )
!
!       .. Scalar Arguments ..
!       LOGICAL            UPPER
!       REAL               A1, A2, A3, B1, B2, B3, CSQ, CSU, CSV, SNQ,
!      $                   SNU, SNV
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLAGS2 computes 2-by-2 orthogonal matrices U, V and Q, such
!> that if ( UPPER ) then
!>
!>           U**T *A*Q = U**T *( A1 A2 )*Q = ( x  0  )
!>                             ( 0  A3 )     ( x  x  )
!> and
!>           V**T*B*Q = V**T *( B1 B2 )*Q = ( x  0  )
!>                            ( 0  B3 )     ( x  x  )
!>
!> or if ( .NOT.UPPER ) then
!>
!>           U**T *A*Q = U**T *( A1 0  )*Q = ( x  x  )
!>                             ( A2 A3 )     ( 0  x  )
!> and
!>           V**T*B*Q = V**T*( B1 0  )*Q = ( x  x  )
!>                           ( B2 B3 )     ( 0  x  )
!>
!> The rows of the transformed A and B are parallel, where
!>
!>   U = (  CSU  SNU ), V = (  CSV SNV ), Q = (  CSQ   SNQ )
!>       ( -SNU  CSU )      ( -SNV CSV )      ( -SNQ   CSQ )
!>
!> Z**T denotes the transpose of Z.
!>
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
!>          A1 is REAL
!> \endverbatim
!>
!> \param[in] A2
!> \verbatim
!>          A2 is REAL
!> \endverbatim
!>
!> \param[in] A3
!> \verbatim
!>          A3 is REAL
!>          On entry, A1, A2 and A3 are elements of the input 2-by-2
!>          upper (lower) triangular matrix A.
!> \endverbatim
!>
!> \param[in] B1
!> \verbatim
!>          B1 is REAL
!> \endverbatim
!>
!> \param[in] B2
!> \verbatim
!>          B2 is REAL
!> \endverbatim
!>
!> \param[in] B3
!> \verbatim
!>          B3 is REAL
!>          On entry, B1, B2 and B3 are elements of the input 2-by-2
!>          upper (lower) triangular matrix B.
!> \endverbatim
!>
!> \param[out] CSU
!> \verbatim
!>          CSU is REAL
!> \endverbatim
!>
!> \param[out] SNU
!> \verbatim
!>          SNU is REAL
!>          The desired orthogonal matrix U.
!> \endverbatim
!>
!> \param[out] CSV
!> \verbatim
!>          CSV is REAL
!> \endverbatim
!>
!> \param[out] SNV
!> \verbatim
!>          SNV is REAL
!>          The desired orthogonal matrix V.
!> \endverbatim
!>
!> \param[out] CSQ
!> \verbatim
!>          CSQ is REAL
!> \endverbatim
!>
!> \param[out] SNQ
!> \verbatim
!>          SNQ is REAL
!>          The desired orthogonal matrix Q.
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
!> \ingroup realOTHERauxiliary
!
!  =====================================================================
      SUBROUTINE SLAGS2(Upper,A1,A2,A3,B1,B2,B3,Csu,Snu,Csv,Snv,Csq,Snq)
      IMPLICIT NONE
!*--SLAGS2155
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      LOGICAL Upper
      REAL A1 , A2 , A3 , B1 , B2 , B3 , Csq , Csu , Csv , Snq , Snu ,  &
     &     Snv
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO
      PARAMETER (ZERO=0.0E+0)
!     ..
!     .. Local Scalars ..
      REAL a , aua11 , aua12 , aua21 , aua22 , avb11 , avb12 , avb21 ,  &
     &     avb22 , csl , csr , d , s1 , s2 , snl , snr , ua11r , ua22r ,&
     &     vb11r , vb22r , b , c , r , ua11 , ua12 , ua21 , ua22 ,      &
     &     vb11 , vb12 , vb21 , vb22
!     ..
!     .. External Subroutines ..
      EXTERNAL SLARTG , SLASV2
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS
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
!
!        The SVD of real 2-by-2 triangular C
!
!         ( CSL -SNL )*( A B )*(  CSR  SNR ) = ( R 0 )
!         ( SNL  CSL ) ( 0 D ) ( -SNR  CSR )   ( 0 T )
!
         CALL SLASV2(a,b,d,s1,s2,snr,csr,snl,csl)
!
         IF ( ABS(csl)>=ABS(snl) .OR. ABS(csr)>=ABS(snr) ) THEN
!
!           Compute the (1,1) and (1,2) elements of U**T *A and V**T *B,
!           and (1,2) element of |U|**T *|A| and |V|**T *|B|.
!
            ua11r = csl*A1
            ua12 = csl*A2 + snl*A3
!
            vb11r = csr*B1
            vb12 = csr*B2 + snr*B3
!
            aua12 = ABS(csl)*ABS(A2) + ABS(snl)*ABS(A3)
            avb12 = ABS(csr)*ABS(B2) + ABS(snr)*ABS(B3)
!
!           zero (1,2) elements of U**T *A and V**T *B
!
            IF ( (ABS(ua11r)+ABS(ua12))==ZERO ) THEN
               CALL SLARTG(-vb11r,vb12,Csq,Snq,r)
            ELSEIF ( aua12/(ABS(ua11r)+ABS(ua12))                       &
     &               <=avb12/(ABS(vb11r)+ABS(vb12)) ) THEN
               CALL SLARTG(-ua11r,ua12,Csq,Snq,r)
            ELSE
               CALL SLARTG(-vb11r,vb12,Csq,Snq,r)
            ENDIF
!
            Csu = csl
            Snu = -snl
            Csv = csr
            Snv = -snr
!
         ELSE
!
!           Compute the (2,1) and (2,2) elements of U**T *A and V**T *B,
!           and (2,2) element of |U|**T *|A| and |V|**T *|B|.
!
            ua21 = -snl*A1
            ua22 = -snl*A2 + csl*A3
!
            vb21 = -snr*B1
            vb22 = -snr*B2 + csr*B3
!
            aua22 = ABS(snl)*ABS(A2) + ABS(csl)*ABS(A3)
            avb22 = ABS(snr)*ABS(B2) + ABS(csr)*ABS(B3)
!
!           zero (2,2) elements of U**T*A and V**T*B, and then swap.
!
            IF ( (ABS(ua21)+ABS(ua22))==ZERO ) THEN
               CALL SLARTG(-vb21,vb22,Csq,Snq,r)
            ELSEIF ( aua22/(ABS(ua21)+ABS(ua22))                        &
     &               <=avb22/(ABS(vb21)+ABS(vb22)) ) THEN
               CALL SLARTG(-ua21,ua22,Csq,Snq,r)
            ELSE
               CALL SLARTG(-vb21,vb22,Csq,Snq,r)
            ENDIF
!
            Csu = snl
            Snu = csl
            Csv = snr
            Snv = csr
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
!
!        The SVD of real 2-by-2 triangular C
!
!         ( CSL -SNL )*( A 0 )*(  CSR  SNR ) = ( R 0 )
!         ( SNL  CSL ) ( C D ) ( -SNR  CSR )   ( 0 T )
!
         CALL SLASV2(a,c,d,s1,s2,snr,csr,snl,csl)
!
         IF ( ABS(csr)>=ABS(snr) .OR. ABS(csl)>=ABS(snl) ) THEN
!
!           Compute the (2,1) and (2,2) elements of U**T *A and V**T *B,
!           and (2,1) element of |U|**T *|A| and |V|**T *|B|.
!
            ua21 = -snr*A1 + csr*A2
            ua22r = csr*A3
!
            vb21 = -snl*B1 + csl*B2
            vb22r = csl*B3
!
            aua21 = ABS(snr)*ABS(A1) + ABS(csr)*ABS(A2)
            avb21 = ABS(snl)*ABS(B1) + ABS(csl)*ABS(B2)
!
!           zero (2,1) elements of U**T *A and V**T *B.
!
            IF ( (ABS(ua21)+ABS(ua22r))==ZERO ) THEN
               CALL SLARTG(vb22r,vb21,Csq,Snq,r)
            ELSEIF ( aua21/(ABS(ua21)+ABS(ua22r))                       &
     &               <=avb21/(ABS(vb21)+ABS(vb22r)) ) THEN
               CALL SLARTG(ua22r,ua21,Csq,Snq,r)
            ELSE
               CALL SLARTG(vb22r,vb21,Csq,Snq,r)
            ENDIF
!
            Csu = csr
            Snu = -snr
            Csv = csl
            Snv = -snl
!
         ELSE
!
!           Compute the (1,1) and (1,2) elements of U**T *A and V**T *B,
!           and (1,1) element of |U|**T *|A| and |V|**T *|B|.
!
            ua11 = csr*A1 + snr*A2
            ua12 = snr*A3
!
            vb11 = csl*B1 + snl*B2
            vb12 = snl*B3
!
            aua11 = ABS(csr)*ABS(A1) + ABS(snr)*ABS(A2)
            avb11 = ABS(csl)*ABS(B1) + ABS(snl)*ABS(B2)
!
!           zero (1,1) elements of U**T*A and V**T*B, and then swap.
!
            IF ( (ABS(ua11)+ABS(ua12))==ZERO ) THEN
               CALL SLARTG(vb12,vb11,Csq,Snq,r)
            ELSEIF ( aua11/(ABS(ua11)+ABS(ua12))                        &
     &               <=avb11/(ABS(vb11)+ABS(vb12)) ) THEN
               CALL SLARTG(ua12,ua11,Csq,Snq,r)
            ELSE
               CALL SLARTG(vb12,vb11,Csq,Snq,r)
            ENDIF
!
            Csu = snr
            Snu = csr
            Csv = snl
            Snv = csl
!
         ENDIF
!
      ENDIF
!
!
!     End of SLAGS2
!
      END SUBROUTINE SLAGS2
