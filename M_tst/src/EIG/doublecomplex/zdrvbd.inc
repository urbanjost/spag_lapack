!*==zdrvbd.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b zdrvbd
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE zdrvbd( NSIZES, MM, NN, NTYPES, DOTYPE, ISEED, THRESH,
!                          A, LDA, U, LDU, VT, LDVT, ASAV, USAV, VTSAV, S,
!                          SSAV, E, WORK, LWORK, RWORK, IWORK, NOUNIT,
!                          INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, LDU, LDVT, LWORK, NOUNIT, NSIZES,
!      $                   NTYPES
!       DOUBLE PRECISION   THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            DOTYPE( * )
!       INTEGER            ISEED( 4 ), IWORK( * ), MM( * ), NN( * )
!       DOUBLE PRECISION   E( * ), RWORK( * ), S( * ), SSAV( * )
!       COMPLEX*16         A( LDA, * ), ASAV( LDA, * ), U( LDU, * ),
!      $                   USAV( LDU, * ), VT( LDVT, * ),
!      $                   VTSAV( LDVT, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> zdrvbd checks the singular value decomposition (SVD) driver ZGESVD,
!> ZGESDD, ZGESVJ, ZGEJSV, ZGESVDX, and ZGESVDQ.
!>
!> ZGESVD and ZGESDD factors A = U diag(S) VT, where U and VT are
!> unitary and diag(S) is diagonal with the entries of the array S on
!> its diagonal. The entries of S are the singular values, nonnegative
!> and stored in decreasing order.  U and VT can be optionally not
!> computed, overwritten on A, or computed partially.
!>
!> A is M by N. Let MNMIN = min( M, N ). S has dimension MNMIN.
!> U can be M by M or M by MNMIN. VT can be N by N or MNMIN by N.
!>
!> When zdrvbd is called, a number of matrix "sizes" (M's and N's)
!> and a number of matrix "types" are specified.  For each size (M,N)
!> and each type of matrix, and for the minimal workspace as well as
!> workspace adequate to permit blocking, an  M x N  matrix "A" will be
!> generated and used to test the SVD routines.  For each matrix, A will
!> be factored as A = U diag(S) VT and the following 12 tests computed:
!>
!> Test for ZGESVD:
!>
!> (1)   | A - U diag(S) VT | / ( |A| max(M,N) ulp )
!>
!> (2)   | I - U'U | / ( M ulp )
!>
!> (3)   | I - VT VT' | / ( N ulp )
!>
!> (4)   S contains MNMIN nonnegative values in decreasing order.
!>       (Return 0 if true, 1/ULP if false.)
!>
!> (5)   | U - Upartial | / ( M ulp ) where Upartial is a partially
!>       computed U.
!>
!> (6)   | VT - VTpartial | / ( N ulp ) where VTpartial is a partially
!>       computed VT.
!>
!> (7)   | S - Spartial | / ( MNMIN ulp |S| ) where Spartial is the
!>       vector of singular values from the partial SVD
!>
!> Test for ZGESDD:
!>
!> (8)   | A - U diag(S) VT | / ( |A| max(M,N) ulp )
!>
!> (9)   | I - U'U | / ( M ulp )
!>
!> (10)  | I - VT VT' | / ( N ulp )
!>
!> (11)  S contains MNMIN nonnegative values in decreasing order.
!>       (Return 0 if true, 1/ULP if false.)
!>
!> (12)  | U - Upartial | / ( M ulp ) where Upartial is a partially
!>       computed U.
!>
!> (13)  | VT - VTpartial | / ( N ulp ) where VTpartial is a partially
!>       computed VT.
!>
!> (14)  | S - Spartial | / ( MNMIN ulp |S| ) where Spartial is the
!>       vector of singular values from the partial SVD
!>
!> Test for ZGESVDQ:
!>
!> (36)  | A - U diag(S) VT | / ( |A| max(M,N) ulp )
!>
!> (37)  | I - U'U | / ( M ulp )
!>
!> (38)  | I - VT VT' | / ( N ulp )
!>
!> (39)  S contains MNMIN nonnegative values in decreasing order.
!>       (Return 0 if true, 1/ULP if false.)
!>
!> Test for ZGESVJ:
!>
!> (15)  | A - U diag(S) VT | / ( |A| max(M,N) ulp )
!>
!> (16)  | I - U'U | / ( M ulp )
!>
!> (17)  | I - VT VT' | / ( N ulp )
!>
!> (18)  S contains MNMIN nonnegative values in decreasing order.
!>       (Return 0 if true, 1/ULP if false.)
!>
!> Test for ZGEJSV:
!>
!> (19)  | A - U diag(S) VT | / ( |A| max(M,N) ulp )
!>
!> (20)  | I - U'U | / ( M ulp )
!>
!> (21)  | I - VT VT' | / ( N ulp )
!>
!> (22)  S contains MNMIN nonnegative values in decreasing order.
!>        (Return 0 if true, 1/ULP if false.)
!>
!> Test for ZGESVDX( 'V', 'V', 'A' )/ZGESVDX( 'N', 'N', 'A' )
!>
!> (23)  | A - U diag(S) VT | / ( |A| max(M,N) ulp )
!>
!> (24)  | I - U'U | / ( M ulp )
!>
!> (25)  | I - VT VT' | / ( N ulp )
!>
!> (26)  S contains MNMIN nonnegative values in decreasing order.
!>       (Return 0 if true, 1/ULP if false.)
!>
!> (27)  | U - Upartial | / ( M ulp ) where Upartial is a partially
!>       computed U.
!>
!> (28)  | VT - VTpartial | / ( N ulp ) where VTpartial is a partially
!>       computed VT.
!>
!> (29)  | S - Spartial | / ( MNMIN ulp |S| ) where Spartial is the
!>       vector of singular values from the partial SVD
!>
!> Test for ZGESVDX( 'V', 'V', 'I' )
!>
!> (30)  | U' A VT''' - diag(S) | / ( |A| max(M,N) ulp )
!>
!> (31)  | I - U'U | / ( M ulp )
!>
!> (32)  | I - VT VT' | / ( N ulp )
!>
!> Test for ZGESVDX( 'V', 'V', 'V' )
!>
!> (33)   | U' A VT''' - diag(S) | / ( |A| max(M,N) ulp )
!>
!> (34)   | I - U'U | / ( M ulp )
!>
!> (35)   | I - VT VT' | / ( N ulp )
!>
!> The "sizes" are specified by the arrays MM(1:NSIZES) and
!> NN(1:NSIZES); the value of each element pair (MM(j),NN(j))
!> specifies one size.  The "types" are specified by a logical array
!> DOTYPE( 1:NTYPES ); if DOTYPE(j) is .TRUE., then matrix type "j"
!> will be generated.
!> Currently, the list of possible types is:
!>
!> (1)  The zero matrix.
!> (2)  The identity matrix.
!> (3)  A matrix of the form  U D V, where U and V are unitary and
!>      D has evenly spaced entries 1, ..., ULP with random signs
!>      on the diagonal.
!> (4)  Same as (3), but multiplied by the underflow-threshold / ULP.
!> (5)  Same as (3), but multiplied by the overflow-threshold * ULP.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] NSIZES
!> \verbatim
!>          NSIZES is INTEGER
!>          The number of sizes of matrices to use.  If it is zero,
!>          zdrvbd does nothing.  It must be at least zero.
!> \endverbatim
!>
!> \param[in] MM
!> \verbatim
!>          MM is INTEGER array, dimension (NSIZES)
!>          An array containing the matrix "heights" to be used.  For
!>          each j=1,...,NSIZES, if MM(j) is zero, then MM(j) and NN(j)
!>          will be ignored.  The MM(j) values must be at least zero.
!> \endverbatim
!>
!> \param[in] NN
!> \verbatim
!>          NN is INTEGER array, dimension (NSIZES)
!>          An array containing the matrix "widths" to be used.  For
!>          each j=1,...,NSIZES, if NN(j) is zero, then MM(j) and NN(j)
!>          will be ignored.  The NN(j) values must be at least zero.
!> \endverbatim
!>
!> \param[in] NTYPES
!> \verbatim
!>          NTYPES is INTEGER
!>          The number of elements in DOTYPE.   If it is zero, zdrvbd
!>          does nothing.  It must be at least zero.  If it is MAXTYP+1
!>          and NSIZES is 1, then an additional type, MAXTYP+1 is
!>          defined, which is to use whatever matrices are in A and B.
!>          This is only useful if DOTYPE(1:MAXTYP) is .FALSE. and
!>          DOTYPE(MAXTYP+1) is .TRUE. .
!> \endverbatim
!>
!> \param[in] DOTYPE
!> \verbatim
!>          DOTYPE is LOGICAL array, dimension (NTYPES)
!>          If DOTYPE(j) is .TRUE., then for each size (m,n), a matrix
!>          of type j will be generated.  If NTYPES is smaller than the
!>          maximum number of types defined (PARAMETER MAXTYP), then
!>          types NTYPES+1 through MAXTYP will not be generated.  If
!>          NTYPES is larger than MAXTYP, DOTYPE(MAXTYP+1) through
!>          DOTYPE(NTYPES) will be ignored.
!> \endverbatim
!>
!> \param[in,out] ISEED
!> \verbatim
!>          ISEED is INTEGER array, dimension (4)
!>          On entry ISEED specifies the seed of the random number
!>          generator. The array elements should be between 0 and 4095;
!>          if not they will be reduced mod 4096.  Also, ISEED(4) must
!>          be odd.  The random number generator uses a linear
!>          congruential sequence limited to small integers, and so
!>          should produce machine independent random numbers. The
!>          values of ISEED are changed on exit, and can be used in the
!>          next call to zdrvbd to continue the same random number
!>          sequence.
!> \endverbatim
!>
!> \param[in] THRESH
!> \verbatim
!>          THRESH is DOUBLE PRECISION
!>          A test will count as "failed" if the "error", computed as
!>          described above, exceeds THRESH.  Note that the error
!>          is scaled to be O(1), so THRESH should be a reasonably
!>          small multiple of 1, e.g., 10 or 100.  In particular,
!>          it should not depend on the precision (single vs. double)
!>          or the size of the matrix.  It must be at least zero.
!> \endverbatim
!>
!> \param[out] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,max(NN))
!>          Used to hold the matrix whose singular values are to be
!>          computed.  On exit, A contains the last matrix actually
!>          used.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of A.  It must be at
!>          least 1 and at least max( MM ).
!> \endverbatim
!>
!> \param[out] U
!> \verbatim
!>          U is COMPLEX*16 array, dimension (LDU,max(MM))
!>          Used to hold the computed matrix of right singular vectors.
!>          On exit, U contains the last such vectors actually computed.
!> \endverbatim
!>
!> \param[in] LDU
!> \verbatim
!>          LDU is INTEGER
!>          The leading dimension of U.  It must be at
!>          least 1 and at least max( MM ).
!> \endverbatim
!>
!> \param[out] VT
!> \verbatim
!>          VT is COMPLEX*16 array, dimension (LDVT,max(NN))
!>          Used to hold the computed matrix of left singular vectors.
!>          On exit, VT contains the last such vectors actually computed.
!> \endverbatim
!>
!> \param[in] LDVT
!> \verbatim
!>          LDVT is INTEGER
!>          The leading dimension of VT.  It must be at
!>          least 1 and at least max( NN ).
!> \endverbatim
!>
!> \param[out] ASAV
!> \verbatim
!>          ASAV is COMPLEX*16 array, dimension (LDA,max(NN))
!>          Used to hold a different copy of the matrix whose singular
!>          values are to be computed.  On exit, A contains the last
!>          matrix actually used.
!> \endverbatim
!>
!> \param[out] USAV
!> \verbatim
!>          USAV is COMPLEX*16 array, dimension (LDU,max(MM))
!>          Used to hold a different copy of the computed matrix of
!>          right singular vectors. On exit, USAV contains the last such
!>          vectors actually computed.
!> \endverbatim
!>
!> \param[out] VTSAV
!> \verbatim
!>          VTSAV is COMPLEX*16 array, dimension (LDVT,max(NN))
!>          Used to hold a different copy of the computed matrix of
!>          left singular vectors. On exit, VTSAV contains the last such
!>          vectors actually computed.
!> \endverbatim
!>
!> \param[out] S
!> \verbatim
!>          S is DOUBLE PRECISION array, dimension (max(min(MM,NN)))
!>          Contains the computed singular values.
!> \endverbatim
!>
!> \param[out] SSAV
!> \verbatim
!>          SSAV is DOUBLE PRECISION array, dimension (max(min(MM,NN)))
!>          Contains another copy of the computed singular values.
!> \endverbatim
!>
!> \param[out] E
!> \verbatim
!>          E is DOUBLE PRECISION array, dimension (max(min(MM,NN)))
!>          Workspace for ZGESVD.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension (LWORK)
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The number of entries in WORK.  This must be at least
!>          MAX(3*MIN(M,N)+MAX(M,N)**2,5*MIN(M,N),3*MAX(M,N)) for all
!>          pairs  (M,N)=(MM(j),NN(j))
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array,
!>                      dimension ( 5*max(max(MM,NN)) )
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension at least 8*min(M,N)
!> \endverbatim
!>
!> \param[in] NOUNIT
!> \verbatim
!>          NOUNIT is INTEGER
!>          The FORTRAN unit number for printing out error messages
!>          (e.g., if a routine returns IINFO not equal to 0.)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          If 0, then everything ran OK.
!>           -1: NSIZES < 0
!>           -2: Some MM(j) < 0
!>           -3: Some NN(j) < 0
!>           -4: NTYPES < 0
!>           -7: THRESH < 0
!>          -10: LDA < 1 or LDA < MMAX, where MMAX is max( MM(j) ).
!>          -12: LDU < 1 or LDU < MMAX.
!>          -14: LDVT < 1 or LDVT < NMAX, where NMAX is max( NN(j) ).
!>          -21: LWORK too small.
!>          If  zlatms, or ZGESVD returns an error code, the
!>              absolute value of it is returned.
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
!> \date June 2016
!
!> \ingroup complex16_eig
!
!  =====================================================================
      SUBROUTINE zdrvbd(Nsizes,Mm,Nn,Ntypes,Dotype,Iseed,Thresh,A,Lda,U,&
     &                  Ldu,Vt,Ldvt,Asav,Usav,Vtsav,S,Ssav,E,Work,Lwork,&
     &                  Rwork,Iwork,Nounit,Info)
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2016
!
      IMPLICIT NONE
!*--ZDRVBD410
!
!     .. Scalar Arguments ..
      INTEGER Info , Lda , Ldu , Ldvt , Lwork , Nounit , Nsizes , Ntypes
      DOUBLE PRECISION Thresh
!     ..
!     .. Array Arguments ..
      LOGICAL Dotype(*)
      INTEGER Iseed(4) , Iwork(*) , Mm(*) , Nn(*)
      DOUBLE PRECISION E(*) , Rwork(*) , S(*) , Ssav(*)
      COMPLEX*16 A(Lda,*) , Asav(Lda,*) , U(Ldu,*) , Usav(Ldu,*) ,      &
     &           Vt(Ldvt,*) , Vtsav(Ldvt,*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE , TWO , HALF
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,HALF=0.5D0)
      COMPLEX*16 CZERO , CONE
      PARAMETER (CZERO=(0.0D+0,0.0D+0),CONE=(1.0D+0,0.0D+0))
      INTEGER MAXTYP
      PARAMETER (MAXTYP=5)
!     ..
!     .. Local Scalars ..
      LOGICAL badmm , badnn
      CHARACTER jobq , jobu , jobvt , range
      INTEGER i , iinfo , ijq , iju , ijvt , il , iu , itemp , iwspc ,  &
     &        iwtmp , j , jsize , jtype , lswork , m , minwrk , mmax ,  &
     &        mnmax , mnmin , mtypes , n , nerrs , nfail , nmax , ns ,  &
     &        nsi , nsv , ntest , ntestf , ntestt , lrwork
      DOUBLE PRECISION anorm , dif , div , ovfl , rtunfl , ulp ,        &
     &                 ulpinv , unfl , vl , vu
!     ..
!     .. Local Scalars for ZGESVDQ ..
      INTEGER liwork , numrank
!     ..
!     .. Local Arrays ..
      CHARACTER cjob(4) , cjobr(3) , cjobv(2)
      INTEGER ioldsd(4) , iseed2(4)
      DOUBLE PRECISION result(39)
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMCH 
      EXTERNAL DLAMCH 
!     ..
!     .. External Subroutines ..
      EXTERNAL ZGESDD , ZGESVD 
     EXTERNAL         ZGESVDQ , ZGESVJ , ZGEJSV , ZGESVDX , ZLACPY , ZLASET
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , DBLE , MAX , MIN
!     ..
!     .. Scalars in Common ..
      CHARACTER*32 SRNamt
!     ..
!     .. Common blocks ..
      COMMON /SRNAMC/ SRNamt
!     ..
!     .. Data statements ..
      DATA cjob/'N' , 'O' , 'S' , 'A'/
      DATA cjobr/'A' , 'V' , 'I'/
      DATA cjobv/'N' , 'V'/
!     ..
!     .. Executable Statements ..
!
!     Check for errors
!
      Info = 0
!
!     Important constants
!
      nerrs = 0
      ntestt = 0
      ntestf = 0
      badmm = .FALSE.
      badnn = .FALSE.
      mmax = 1
      nmax = 1
      mnmax = 1
      minwrk = 1
      DO j = 1 , Nsizes
         mmax = MAX(mmax,Mm(j))
         IF ( Mm(j)<0 ) badmm = .TRUE.
         nmax = MAX(nmax,Nn(j))
         IF ( Nn(j)<0 ) badnn = .TRUE.
         mnmax = MAX(mnmax,MIN(Mm(j),Nn(j)))
         minwrk = MAX(minwrk,MAX(3*MIN(Mm(j),Nn(j))+MAX(Mm(j),Nn(j))**2,&
     &            5*MIN(Mm(j),Nn(j)),3*MAX(Mm(j),Nn(j))))
      ENDDO
!
!     Check for errors
!
      IF ( Nsizes<0 ) THEN
         Info = -1
      ELSEIF ( badmm ) THEN
         Info = -2
      ELSEIF ( badnn ) THEN
         Info = -3
      ELSEIF ( Ntypes<0 ) THEN
         Info = -4
      ELSEIF ( Lda<MAX(1,mmax) ) THEN
         Info = -10
      ELSEIF ( Ldu<MAX(1,mmax) ) THEN
         Info = -12
      ELSEIF ( Ldvt<MAX(1,nmax) ) THEN
         Info = -14
      ELSEIF ( minwrk>Lwork ) THEN
         Info = -21
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL xerbla('zdrvbd',-Info)
         RETURN
      ENDIF
!
!     Quick return if nothing to do
!
      IF ( Nsizes==0 .OR. Ntypes==0 ) RETURN
!
!     More Important constants
!
      unfl = DLAMCH('S')
      ovfl = ONE/unfl
      ulp = DLAMCH('E')
      ulpinv = ONE/ulp
      rtunfl = SQRT(unfl)
!
!     Loop over sizes, types
!
      nerrs = 0
!
      DO jsize = 1 , Nsizes
         m = Mm(jsize)
         n = Nn(jsize)
         mnmin = MIN(m,n)
!
         IF ( Nsizes/=1 ) THEN
            mtypes = MIN(MAXTYP,Ntypes)
         ELSE
            mtypes = MIN(MAXTYP+1,Ntypes)
         ENDIF
!
         DO jtype = 1 , mtypes
            IF ( Dotype(jtype) ) THEN
               ntest = 0
!
               DO j = 1 , 4
                  ioldsd(j) = Iseed(j)
               ENDDO
!
!           Compute "A"
!
               IF ( mtypes<=MAXTYP ) THEN
!
                  IF ( jtype==1 ) THEN
!
!              Zero matrix
!
                     CALL ZLASET('Full',m,n,CZERO,CZERO,A,Lda)
                     DO i = 1 , MIN(m,n)
                        S(i) = ZERO
                     ENDDO
!
                  ELSEIF ( jtype==2 ) THEN
!
!              Identity matrix
!
                     CALL ZLASET('Full',m,n,CZERO,CONE,A,Lda)
                     DO i = 1 , MIN(m,n)
                        S(i) = ONE
                     ENDDO
!
                  ELSE
!
!              (Scaled) random matrix
!
                     IF ( jtype==3 ) anorm = ONE
                     IF ( jtype==4 ) anorm = unfl/ulp
                     IF ( jtype==5 ) anorm = ovfl*ulp
                     CALL zlatms(m,n,'U',Iseed,'N',S,4,DBLE(mnmin),     &
     &                           anorm,m-1,n-1,'N',A,Lda,Work,iinfo)
                     IF ( iinfo/=0 ) THEN
                        WRITE (Nounit,FMT=99004) 'Generator' , iinfo ,  &
     &                         m , n , jtype , ioldsd
                        Info = ABS(iinfo)
                        RETURN
                     ENDIF
                  ENDIF
               ENDIF
!
               CALL ZLACPY('F',m,n,A,Lda,Asav,Lda)
!
!           Do for minimal and adequate (for blocking) workspace
!
               DO iwspc = 1 , 4
!
!              Test for ZGESVD
!
                  iwtmp = 2*MIN(m,n) + MAX(m,n)
                  lswork = iwtmp + (iwspc-1)*(Lwork-iwtmp)/3
                  lswork = MIN(lswork,Lwork)
                  lswork = MAX(lswork,1)
                  IF ( iwspc==4 ) lswork = Lwork
!
                  DO j = 1 , 35
                     result(j) = -ONE
                  ENDDO
!
!              Factorize A
!
                  IF ( iwspc>1 ) CALL ZLACPY('F',m,n,Asav,Lda,A,Lda)
                  SRNamt = 'ZGESVD'
                  CALL ZGESVD('A','A',m,n,A,Lda,Ssav,Usav,Ldu,Vtsav,    &
     &                        Ldvt,Work,lswork,Rwork,iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99005) 'GESVD' , iinfo , m , n , &
     &                      jtype , lswork , ioldsd
                     Info = ABS(iinfo)
                     RETURN
                  ENDIF
!
!              Do tests 1--4
!
                  CALL zbdt01(m,n,0,Asav,Lda,Usav,Ldu,Ssav,E,Vtsav,Ldvt,&
     &                        Work,Rwork,result(1))
                  IF ( m/=0 .AND. n/=0 ) THEN
                     CALL zunt01('Columns',mnmin,m,Usav,Ldu,Work,Lwork, &
     &                           Rwork,result(2))
                     CALL zunt01('Rows',mnmin,n,Vtsav,Ldvt,Work,Lwork,  &
     &                           Rwork,result(3))
                  ENDIF
                  result(4) = 0
                  DO i = 1 , mnmin - 1
                     IF ( Ssav(i)<Ssav(i+1) ) result(4) = ulpinv
                     IF ( Ssav(i)<ZERO ) result(4) = ulpinv
                  ENDDO
                  IF ( mnmin>=1 ) THEN
                     IF ( Ssav(mnmin)<ZERO ) result(4) = ulpinv
                  ENDIF
!
!              Do partial SVDs, comparing to SSAV, USAV, and VTSAV
!
                  result(5) = ZERO
                  result(6) = ZERO
                  result(7) = ZERO
                  DO iju = 0 , 3
                     DO ijvt = 0 , 3
                        IF ( .NOT.((iju==3 .AND. ijvt==3) .OR. (iju==1  &
     &                       .AND. ijvt==1)) ) THEN
                           jobu = cjob(iju+1)
                           jobvt = cjob(ijvt+1)
                           CALL ZLACPY('F',m,n,Asav,Lda,A,Lda)
                           SRNamt = 'ZGESVD'
                           CALL ZGESVD(jobu,jobvt,m,n,A,Lda,S,U,Ldu,Vt, &
     &                                 Ldvt,Work,lswork,Rwork,iinfo)
!
!                    Compare U
!
                           dif = ZERO
                           IF ( m>0 .AND. n>0 ) THEN
                              IF ( iju==1 ) THEN
                                 CALL zunt03('C',m,mnmin,m,mnmin,Usav,  &
     &                              Ldu,A,Lda,Work,Lwork,Rwork,dif,     &
     &                              iinfo)
                              ELSEIF ( iju==2 ) THEN
                                 CALL zunt03('C',m,mnmin,m,mnmin,Usav,  &
     &                              Ldu,U,Ldu,Work,Lwork,Rwork,dif,     &
     &                              iinfo)
                              ELSEIF ( iju==3 ) THEN
                                 CALL zunt03('C',m,m,m,mnmin,Usav,Ldu,U,&
     &                              Ldu,Work,Lwork,Rwork,dif,iinfo)
                              ENDIF
                           ENDIF
                           result(5) = MAX(result(5),dif)
!
!                    Compare VT
!
                           dif = ZERO
                           IF ( m>0 .AND. n>0 ) THEN
                              IF ( ijvt==1 ) THEN
                                 CALL zunt03('R',n,mnmin,n,mnmin,Vtsav, &
     &                              Ldvt,A,Lda,Work,Lwork,Rwork,dif,    &
     &                              iinfo)
                              ELSEIF ( ijvt==2 ) THEN
                                 CALL zunt03('R',n,mnmin,n,mnmin,Vtsav, &
     &                              Ldvt,Vt,Ldvt,Work,Lwork,Rwork,dif,  &
     &                              iinfo)
                              ELSEIF ( ijvt==3 ) THEN
                                 CALL zunt03('R',n,n,n,mnmin,Vtsav,Ldvt,&
     &                              Vt,Ldvt,Work,Lwork,Rwork,dif,iinfo)
                              ENDIF
                           ENDIF
                           result(6) = MAX(result(6),dif)
!
!                    Compare S
!
                           dif = ZERO
                           div = MAX(DBLE(mnmin)*ulp*S(1),              &
     &                           DLAMCH('Safe minimum'))
                           DO i = 1 , mnmin - 1
                              IF ( Ssav(i)<Ssav(i+1) ) dif = ulpinv
                              IF ( Ssav(i)<ZERO ) dif = ulpinv
                              dif = MAX(dif,ABS(Ssav(i)-S(i))/div)
                           ENDDO
                           result(7) = MAX(result(7),dif)
                        ENDIF
                     ENDDO
                  ENDDO
!
!              Test for ZGESDD
!
                  iwtmp = 2*mnmin*mnmin + 2*mnmin + MAX(m,n)
                  lswork = iwtmp + (iwspc-1)*(Lwork-iwtmp)/3
                  lswork = MIN(lswork,Lwork)
                  lswork = MAX(lswork,1)
                  IF ( iwspc==4 ) lswork = Lwork
!
!              Factorize A
!
                  CALL ZLACPY('F',m,n,Asav,Lda,A,Lda)
                  SRNamt = 'ZGESDD'
                  CALL ZGESDD('A',m,n,A,Lda,Ssav,Usav,Ldu,Vtsav,Ldvt,   &
     &                        Work,lswork,Rwork,Iwork,iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99005) 'GESDD' , iinfo , m , n , &
     &                      jtype , lswork , ioldsd
                     Info = ABS(iinfo)
                     RETURN
                  ENDIF
!
!              Do tests 1--4
!
                  CALL zbdt01(m,n,0,Asav,Lda,Usav,Ldu,Ssav,E,Vtsav,Ldvt,&
     &                        Work,Rwork,result(8))
                  IF ( m/=0 .AND. n/=0 ) THEN
                     CALL zunt01('Columns',mnmin,m,Usav,Ldu,Work,Lwork, &
     &                           Rwork,result(9))
                     CALL zunt01('Rows',mnmin,n,Vtsav,Ldvt,Work,Lwork,  &
     &                           Rwork,result(10))
                  ENDIF
                  result(11) = 0
                  DO i = 1 , mnmin - 1
                     IF ( Ssav(i)<Ssav(i+1) ) result(11) = ulpinv
                     IF ( Ssav(i)<ZERO ) result(11) = ulpinv
                  ENDDO
                  IF ( mnmin>=1 ) THEN
                     IF ( Ssav(mnmin)<ZERO ) result(11) = ulpinv
                  ENDIF
!
!              Do partial SVDs, comparing to SSAV, USAV, and VTSAV
!
                  result(12) = ZERO
                  result(13) = ZERO
                  result(14) = ZERO
                  DO ijq = 0 , 2
                     jobq = cjob(ijq+1)
                     CALL ZLACPY('F',m,n,Asav,Lda,A,Lda)
                     SRNamt = 'ZGESDD'
                     CALL ZGESDD(jobq,m,n,A,Lda,S,U,Ldu,Vt,Ldvt,Work,   &
     &                           lswork,Rwork,Iwork,iinfo)
!
!                 Compare U
!
                     dif = ZERO
                     IF ( m>0 .AND. n>0 ) THEN
                        IF ( ijq==1 ) THEN
                           IF ( m>=n ) THEN
                              CALL zunt03('C',m,mnmin,m,mnmin,Usav,Ldu, &
     &                           A,Lda,Work,Lwork,Rwork,dif,iinfo)
                           ELSE
                              CALL zunt03('C',m,mnmin,m,mnmin,Usav,Ldu, &
     &                           U,Ldu,Work,Lwork,Rwork,dif,iinfo)
                           ENDIF
                        ELSEIF ( ijq==2 ) THEN
                           CALL zunt03('C',m,mnmin,m,mnmin,Usav,Ldu,U,  &
     &                                 Ldu,Work,Lwork,Rwork,dif,iinfo)
                        ENDIF
                     ENDIF
                     result(12) = MAX(result(12),dif)
!
!                 Compare VT
!
                     dif = ZERO
                     IF ( m>0 .AND. n>0 ) THEN
                        IF ( ijq==1 ) THEN
                           IF ( m>=n ) THEN
                              CALL zunt03('R',n,mnmin,n,mnmin,Vtsav,    &
     &                           Ldvt,Vt,Ldvt,Work,Lwork,Rwork,dif,     &
     &                           iinfo)
                           ELSE
                              CALL zunt03('R',n,mnmin,n,mnmin,Vtsav,    &
     &                           Ldvt,A,Lda,Work,Lwork,Rwork,dif,iinfo)
                           ENDIF
                        ELSEIF ( ijq==2 ) THEN
                           CALL zunt03('R',n,mnmin,n,mnmin,Vtsav,Ldvt,  &
     &                                 Vt,Ldvt,Work,Lwork,Rwork,dif,    &
     &                                 iinfo)
                        ENDIF
                     ENDIF
                     result(13) = MAX(result(13),dif)
!
!                 Compare S
!
                     dif = ZERO
                     div = MAX(DBLE(mnmin)*ulp*S(1),                    &
     &                     DLAMCH('Safe minimum'))
                     DO i = 1 , mnmin - 1
                        IF ( Ssav(i)<Ssav(i+1) ) dif = ulpinv
                        IF ( Ssav(i)<ZERO ) dif = ulpinv
                        dif = MAX(dif,ABS(Ssav(i)-S(i))/div)
                     ENDDO
                     result(14) = MAX(result(14),dif)
                  ENDDO
!
!              Test ZGESVDQ
!              Note: ZGESVDQ only works for M >= N
!
                  result(36) = ZERO
                  result(37) = ZERO
                  result(38) = ZERO
                  result(39) = ZERO
!
                  IF ( m>=n ) THEN
                     iwtmp = 2*mnmin*mnmin + 2*mnmin + MAX(m,n)
                     lswork = iwtmp + (iwspc-1)*(Lwork-iwtmp)/3
                     lswork = MIN(lswork,Lwork)
                     lswork = MAX(lswork,1)
                     IF ( iwspc==4 ) lswork = Lwork
!
                     CALL ZLACPY('F',m,n,Asav,Lda,A,Lda)
                     SRNamt = 'ZGESVDQ'
!
                     lrwork = MAX(2,m,5*n)
                     liwork = MAX(n,1)
                     CALL ZGESVDQ('H','N','N','A','A',m,n,A,Lda,Ssav,   &
     &                            Usav,Ldu,Vtsav,Ldvt,numrank,Iwork,    &
     &                            liwork,Work,Lwork,Rwork,lrwork,iinfo)
!
                     IF ( iinfo/=0 ) THEN
                        WRITE (Nounit,FMT=99005) 'ZGESVDQ' , iinfo , m ,&
     &                         n , jtype , lswork , ioldsd
                        Info = ABS(iinfo)
                        RETURN
                     ENDIF
!
!                 Do tests 36--39
!
                     CALL zbdt01(m,n,0,Asav,Lda,Usav,Ldu,Ssav,E,Vtsav,  &
     &                           Ldvt,Work,Rwork,result(36))
                     IF ( m/=0 .AND. n/=0 ) THEN
                        CALL zunt01('Columns',m,m,Usav,Ldu,Work,Lwork,  &
     &                              Rwork,result(37))
                        CALL zunt01('Rows',n,n,Vtsav,Ldvt,Work,Lwork,   &
     &                              Rwork,result(38))
                     ENDIF
                     result(39) = ZERO
                     DO i = 1 , mnmin - 1
                        IF ( Ssav(i)<Ssav(i+1) ) result(39) = ulpinv
                        IF ( Ssav(i)<ZERO ) result(39) = ulpinv
                     ENDDO
                     IF ( mnmin>=1 ) THEN
                        IF ( Ssav(mnmin)<ZERO ) result(39) = ulpinv
                     ENDIF
                  ENDIF
!
!              Test ZGESVJ
!              Note: ZGESVJ only works for M >= N
!
                  result(15) = ZERO
                  result(16) = ZERO
                  result(17) = ZERO
                  result(18) = ZERO
!
                  IF ( m>=n ) THEN
                     iwtmp = 2*mnmin*mnmin + 2*mnmin + MAX(m,n)
                     lswork = iwtmp + (iwspc-1)*(Lwork-iwtmp)/3
                     lswork = MIN(lswork,Lwork)
                     lswork = MAX(lswork,1)
                     lrwork = MAX(6,n)
                     IF ( iwspc==4 ) lswork = Lwork
!
                     CALL ZLACPY('F',m,n,Asav,Lda,Usav,Lda)
                     SRNamt = 'ZGESVJ'
                     CALL ZGESVJ('G','U','V',m,n,Usav,Lda,Ssav,0,A,Ldvt,&
     &                           Work,Lwork,Rwork,lrwork,iinfo)
!
!                 ZGESVJ returns V not VH
!
                     DO j = 1 , n
                        DO i = 1 , n
                           Vtsav(j,i) = CONJG(A(i,j))
                        ENDDO
                     ENDDO
!
                     IF ( iinfo/=0 ) THEN
                        WRITE (Nounit,FMT=99005) 'GESVJ' , iinfo , m ,  &
     &                         n , jtype , lswork , ioldsd
                        Info = ABS(iinfo)
                        RETURN
                     ENDIF
!
!                 Do tests 15--18
!
                     CALL zbdt01(m,n,0,Asav,Lda,Usav,Ldu,Ssav,E,Vtsav,  &
     &                           Ldvt,Work,Rwork,result(15))
                     IF ( m/=0 .AND. n/=0 ) THEN
                        CALL zunt01('Columns',m,m,Usav,Ldu,Work,Lwork,  &
     &                              Rwork,result(16))
                        CALL zunt01('Rows',n,n,Vtsav,Ldvt,Work,Lwork,   &
     &                              Rwork,result(17))
                     ENDIF
                     result(18) = ZERO
                     DO i = 1 , mnmin - 1
                        IF ( Ssav(i)<Ssav(i+1) ) result(18) = ulpinv
                        IF ( Ssav(i)<ZERO ) result(18) = ulpinv
                     ENDDO
                     IF ( mnmin>=1 ) THEN
                        IF ( Ssav(mnmin)<ZERO ) result(18) = ulpinv
                     ENDIF
                  ENDIF
!
!              Test ZGEJSV
!              Note: ZGEJSV only works for M >= N
!
                  result(19) = ZERO
                  result(20) = ZERO
                  result(21) = ZERO
                  result(22) = ZERO
                  IF ( m>=n ) THEN
                     iwtmp = 2*mnmin*mnmin + 2*mnmin + MAX(m,n)
                     lswork = iwtmp + (iwspc-1)*(Lwork-iwtmp)/3
                     lswork = MIN(lswork,Lwork)
                     lswork = MAX(lswork,1)
                     IF ( iwspc==4 ) lswork = Lwork
                     lrwork = MAX(7,n+2*m)
!
                     CALL ZLACPY('F',m,n,Asav,Lda,Vtsav,Lda)
                     SRNamt = 'ZGEJSV'
                     CALL ZGEJSV('G','U','V','R','N','N',m,n,Vtsav,Lda, &
     &                           Ssav,Usav,Ldu,A,Ldvt,Work,Lwork,Rwork, &
     &                           lrwork,Iwork,iinfo)
!
!                 ZGEJSV returns V not VH
!
                     DO j = 1 , n
                        DO i = 1 , n
                           Vtsav(j,i) = CONJG(A(i,j))
                        ENDDO
                     ENDDO
!
                     IF ( iinfo/=0 ) THEN
                        WRITE (Nounit,FMT=99005) 'GEJSV' , iinfo , m ,  &
     &                         n , jtype , lswork , ioldsd
                        Info = ABS(iinfo)
                        RETURN
                     ENDIF
!
!                 Do tests 19--22
!
                     CALL zbdt01(m,n,0,Asav,Lda,Usav,Ldu,Ssav,E,Vtsav,  &
     &                           Ldvt,Work,Rwork,result(19))
                     IF ( m/=0 .AND. n/=0 ) THEN
                        CALL zunt01('Columns',m,m,Usav,Ldu,Work,Lwork,  &
     &                              Rwork,result(20))
                        CALL zunt01('Rows',n,n,Vtsav,Ldvt,Work,Lwork,   &
     &                              Rwork,result(21))
                     ENDIF
                     result(22) = ZERO
                     DO i = 1 , mnmin - 1
                        IF ( Ssav(i)<Ssav(i+1) ) result(22) = ulpinv
                        IF ( Ssav(i)<ZERO ) result(22) = ulpinv
                     ENDDO
                     IF ( mnmin>=1 ) THEN
                        IF ( Ssav(mnmin)<ZERO ) result(22) = ulpinv
                     ENDIF
                  ENDIF
!
!              Test ZGESVDX
!
!              Factorize A
!
                  CALL ZLACPY('F',m,n,Asav,Lda,A,Lda)
                  SRNamt = 'ZGESVDX'
                  CALL ZGESVDX('V','V','A',m,n,A,Lda,vl,vu,il,iu,ns,    &
     &                         Ssav,Usav,Ldu,Vtsav,Ldvt,Work,Lwork,     &
     &                         Rwork,Iwork,iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99005) 'GESVDX' , iinfo , m , n ,&
     &                      jtype , lswork , ioldsd
                     Info = ABS(iinfo)
                     RETURN
                  ENDIF
!
!              Do tests 1--4
!
                  result(23) = ZERO
                  result(24) = ZERO
                  result(25) = ZERO
                  CALL zbdt01(m,n,0,Asav,Lda,Usav,Ldu,Ssav,E,Vtsav,Ldvt,&
     &                        Work,Rwork,result(23))
                  IF ( m/=0 .AND. n/=0 ) THEN
                     CALL zunt01('Columns',mnmin,m,Usav,Ldu,Work,Lwork, &
     &                           Rwork,result(24))
                     CALL zunt01('Rows',mnmin,n,Vtsav,Ldvt,Work,Lwork,  &
     &                           Rwork,result(25))
                  ENDIF
                  result(26) = ZERO
                  DO i = 1 , mnmin - 1
                     IF ( Ssav(i)<Ssav(i+1) ) result(26) = ulpinv
                     IF ( Ssav(i)<ZERO ) result(26) = ulpinv
                  ENDDO
                  IF ( mnmin>=1 ) THEN
                     IF ( Ssav(mnmin)<ZERO ) result(26) = ulpinv
                  ENDIF
!
!              Do partial SVDs, comparing to SSAV, USAV, and VTSAV
!
                  result(27) = ZERO
                  result(28) = ZERO
                  result(29) = ZERO
                  DO iju = 0 , 1
                     DO ijvt = 0 , 1
                        IF ( .NOT.((iju==0 .AND. ijvt==0) .OR. (iju==1  &
     &                       .AND. ijvt==1)) ) THEN
                           jobu = cjobv(iju+1)
                           jobvt = cjobv(ijvt+1)
                           range = cjobr(1)
                           CALL ZLACPY('F',m,n,Asav,Lda,A,Lda)
                           SRNamt = 'ZGESVDX'
                           CALL ZGESVDX(jobu,jobvt,'A',m,n,A,Lda,vl,vu, &
     &                                  il,iu,ns,Ssav,U,Ldu,Vt,Ldvt,    &
     &                                  Work,Lwork,Rwork,Iwork,iinfo)
!
!                    Compare U
!
                           dif = ZERO
                           IF ( m>0 .AND. n>0 ) THEN
                              IF ( iju==1 )                             &
     &                             CALL zunt03('C',m,mnmin,m,mnmin,Usav,&
     &                             Ldu,U,Ldu,Work,Lwork,Rwork,dif,iinfo)
                           ENDIF
                           result(27) = MAX(result(27),dif)
!
!                    Compare VT
!
                           dif = ZERO
                           IF ( m>0 .AND. n>0 ) THEN
                              IF ( ijvt==1 )                            &
     &                             CALL zunt03('R',n,mnmin,n,mnmin,     &
     &                             Vtsav,Ldvt,Vt,Ldvt,Work,Lwork,Rwork, &
     &                             dif,iinfo)
                           ENDIF
                           result(28) = MAX(result(28),dif)
!
!                    Compare S
!
                           dif = ZERO
                           div = MAX(DBLE(mnmin)*ulp*S(1),              &
     &                           DLAMCH('Safe minimum'))
                           DO i = 1 , mnmin - 1
                              IF ( Ssav(i)<Ssav(i+1) ) dif = ulpinv
                              IF ( Ssav(i)<ZERO ) dif = ulpinv
                              dif = MAX(dif,ABS(Ssav(i)-S(i))/div)
                           ENDDO
                           result(29) = MAX(result(29),dif)
                        ENDIF
                     ENDDO
                  ENDDO
!
!              Do tests 8--10
!
                  DO i = 1 , 4
                     iseed2(i) = Iseed(i)
                  ENDDO
                  IF ( mnmin<=1 ) THEN
                     il = 1
                     iu = MAX(1,mnmin)
                  ELSE
                     il = 1 + INT((mnmin-1)*dlarnd(1,iseed2))
                     iu = 1 + INT((mnmin-1)*dlarnd(1,iseed2))
                     IF ( iu<il ) THEN
                        itemp = iu
                        iu = il
                        il = itemp
                     ENDIF
                  ENDIF
                  CALL ZLACPY('F',m,n,Asav,Lda,A,Lda)
                  SRNamt = 'ZGESVDX'
                  CALL ZGESVDX('V','V','I',m,n,A,Lda,vl,vu,il,iu,nsi,S, &
     &                         U,Ldu,Vt,Ldvt,Work,Lwork,Rwork,Iwork,    &
     &                         iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99005) 'GESVDX' , iinfo , m , n ,&
     &                      jtype , lswork , ioldsd
                     Info = ABS(iinfo)
                     RETURN
                  ENDIF
!
                  result(30) = ZERO
                  result(31) = ZERO
                  result(32) = ZERO
                  CALL zbdt05(m,n,Asav,Lda,S,nsi,U,Ldu,Vt,Ldvt,Work,    &
     &                        result(30))
                  IF ( m/=0 .AND. n/=0 ) THEN
                     CALL zunt01('Columns',m,nsi,U,Ldu,Work,Lwork,Rwork,&
     &                           result(31))
                     CALL zunt01('Rows',nsi,n,Vt,Ldvt,Work,Lwork,Rwork, &
     &                           result(32))
                  ENDIF
!
!              Do tests 11--13
!
                  IF ( mnmin>0 .AND. nsi>1 ) THEN
                     IF ( il/=1 ) THEN
                        vu = Ssav(il)                                   &
     &                       + MAX(HALF*ABS(Ssav(il)-Ssav(il-1)),       &
     &                       ulp*anorm,TWO*rtunfl)
                     ELSE
                        vu = Ssav(1)                                    &
     &                       + MAX(HALF*ABS(Ssav(ns)-Ssav(1)),ulp*anorm,&
     &                       TWO*rtunfl)
                     ENDIF
                     IF ( iu/=ns ) THEN
                        vl = Ssav(iu)                                   &
     &                       - MAX(ulp*anorm,TWO*rtunfl,HALF*ABS        &
     &                       (Ssav(iu+1)-Ssav(iu)))
                     ELSE
                        vl = Ssav(ns)                                   &
     &                       - MAX(ulp*anorm,TWO*rtunfl,HALF*ABS        &
     &                       (Ssav(ns)-Ssav(1)))
                     ENDIF
                     vl = MAX(vl,ZERO)
                     vu = MAX(vu,ZERO)
                     IF ( vl>=vu ) vu = MAX(vu*2,vu+vl+HALF)
                  ELSE
                     vl = ZERO
                     vu = ONE
                  ENDIF
                  CALL ZLACPY('F',m,n,Asav,Lda,A,Lda)
                  SRNamt = 'ZGESVDX'
                  CALL ZGESVDX('V','V','V',m,n,A,Lda,vl,vu,il,iu,nsv,S, &
     &                         U,Ldu,Vt,Ldvt,Work,Lwork,Rwork,Iwork,    &
     &                         iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99005) 'GESVDX' , iinfo , m , n ,&
     &                      jtype , lswork , ioldsd
                     Info = ABS(iinfo)
                     RETURN
                  ENDIF
!
                  result(33) = ZERO
                  result(34) = ZERO
                  result(35) = ZERO
                  CALL zbdt05(m,n,Asav,Lda,S,nsv,U,Ldu,Vt,Ldvt,Work,    &
     &                        result(33))
                  IF ( m/=0 .AND. n/=0 ) THEN
                     CALL zunt01('Columns',m,nsv,U,Ldu,Work,Lwork,Rwork,&
     &                           result(34))
                     CALL zunt01('Rows',nsv,n,Vt,Ldvt,Work,Lwork,Rwork, &
     &                           result(35))
                  ENDIF
!
!              End of Loop -- Check for RESULT(j) > THRESH
!
                  ntest = 0
                  nfail = 0
                  DO j = 1 , 39
                     IF ( result(j)>=ZERO ) ntest = ntest + 1
                     IF ( result(j)>=Thresh ) nfail = nfail + 1
                  ENDDO
!
                  IF ( nfail>0 ) ntestf = ntestf + 1
                  IF ( ntestf==1 ) THEN
                     WRITE (Nounit,FMT=99001)
                     WRITE (Nounit,FMT=99002) Thresh
                     ntestf = 2
                  ENDIF
!
                  DO j = 1 , 39
                     IF ( result(j)>=Thresh ) WRITE (Nounit,FMT=99003)  &
     &                    m , n , jtype , iwspc , ioldsd , j , result(j)
                  ENDDO
!
                  nerrs = nerrs + nfail
                  ntestt = ntestt + ntest
!
               ENDDO
            ENDIF
!
         ENDDO
      ENDDO
!
!     Summary
!
      CALL alasvm('ZBD',Nounit,nerrs,ntestt,0)
!
99001 FORMAT (' SVD -- Complex Singular Value Decomposition Driver ',   &
     &        /' Matrix types (see zdrvbd for details):',               &
     &        //' 1 = Zero matrix',/' 2 = Identity matrix',             &
     &        /' 3 = Evenly spaced singular values near 1',             &
     &        /' 4 = Evenly spaced singular values near underflow',     &
     &        /' 5 = Evenly spaced singular values near overflow',      &
     &        //' Tests performed: ( A is dense, U and V are unitary,', &
     &        /19X,' S is an array, and Upartial, VTpartial, and',/19X, &
     &        ' Spartial are partially computed U, VT and S),',/)
99002 FORMAT (' Tests performed with Test Threshold = ',F8.2,           &
     &        /' ZGESVD: ',                                             &
     &        /' 1 = | A - U diag(S) VT | / ( |A| max(M,N) ulp ) ',     &
     &        /' 2 = | I - U**T U | / ( M ulp ) ',                      &
     &        /' 3 = | I - VT VT**T | / ( N ulp ) ',                    &
     &        /' 4 = 0 if S contains min(M,N) nonnegative values in',   &
     &        ' decreasing order, else 1/ulp',                          &
     &        /' 5 = | U - Upartial | / ( M ulp )',                     &
     &        /' 6 = | VT - VTpartial | / ( N ulp )',                   &
     &        /' 7 = | S - Spartial | / ( min(M,N) ulp |S| )',          &
     &        /' ZGESDD: ',                                             &
     &        /' 8 = | A - U diag(S) VT | / ( |A| max(M,N) ulp ) ',     &
     &        /' 9 = | I - U**T U | / ( M ulp ) ',                      &
     &        /'10 = | I - VT VT**T | / ( N ulp ) ',                    &
     &        /'11 = 0 if S contains min(M,N) nonnegative values in',   &
     &        ' decreasing order, else 1/ulp',                          &
     &        /'12 = | U - Upartial | / ( M ulp )',                     &
     &        /'13 = | VT - VTpartial | / ( N ulp )',                   &
     &        /'14 = | S - Spartial | / ( min(M,N) ulp |S| )',          &
     &        /' ZGESVJ: ',                                             &
     &        //'15 = | A - U diag(S) VT | / ( |A| max(M,N) ulp ) ',    &
     &        /'16 = | I - U**T U | / ( M ulp ) ',                      &
     &        /'17 = | I - VT VT**T | / ( N ulp ) ',                    &
     &        /'18 = 0 if S contains min(M,N) nonnegative values in',   &
     &        ' decreasing order, else 1/ulp',/' ZGESJV: ',             &
     &        //'19 = | A - U diag(S) VT | / ( |A| max(M,N) ulp )',     &
     &        /'20 = | I - U**T U | / ( M ulp ) ',                      &
     &        /'21 = | I - VT VT**T | / ( N ulp ) ',                    &
     &        /'22 = 0 if S contains min(M,N) nonnegative values in',   &
     &        ' decreasing order, else 1/ulp',/' ZGESVDX(V,V,A): ',     &
     &        /'23 = | A - U diag(S) VT | / ( |A| max(M,N) ulp ) ',     &
     &        /'24 = | I - U**T U | / ( M ulp ) ',                      &
     &        /'25 = | I - VT VT**T | / ( N ulp ) ',                    &
     &        /'26 = 0 if S contains min(M,N) nonnegative values in',   &
     &        ' decreasing order, else 1/ulp',                          &
     &        /'27 = | U - Upartial | / ( M ulp )',                     &
     &        /'28 = | VT - VTpartial | / ( N ulp )',                   &
     &        /'29 = | S - Spartial | / ( min(M,N) ulp |S| )',          &
     &        /' ZGESVDX(V,V,I): ',                                     &
     &        /'30 = | U**T A VT**T - diag(S) | / ( |A| max(M,N) ulp )',&
     &        /'31 = | I - U**T U | / ( M ulp ) ',                      &
     &        /'32 = | I - VT VT**T | / ( N ulp ) ',/' ZGESVDX(V,V,V) ',&
     &        /'33 = | U**T A VT**T - diag(S) | / ( |A| max(M,N) ulp )',&
     &        /'34 = | I - U**T U | / ( M ulp ) ',                      &
     &        /'35 = | I - VT VT**T | / ( N ulp ) ',                    &
     &        ' ZGESVDQ(H,N,N,A,A',                                     &
     &        /'36 = | A - U diag(S) VT | / ( |A| max(M,N) ulp ) ',     &
     &        /'37 = | I - U**T U | / ( M ulp ) ',                      &
     &        /'38 = | I - VT VT**T | / ( N ulp ) ',                    &
     &        /'39 = 0 if S contains min(M,N) nonnegative values in',   &
     &        ' decreasing order, else 1/ulp',//)
99003 FORMAT (' M=',I5,', N=',I5,', type ',I1,', IWS=',I1,', seed=',    &
     &        4(I4,','),' test(',I2,')=',G11.4)
99004 FORMAT (' zdrvbd: ',A,' returned INFO=',I6,'.',/9X,'M=',I6,', N=',&
     &        I6,', JTYPE=',I6,', ISEED=(',3(I5,','),I5,')')
99005 FORMAT (' zdrvbd: ',A,' returned INFO=',I6,'.',/9X,'M=',I6,', N=',&
     &        I6,', JTYPE=',I6,', LSWORK=',I6,/9X,'ISEED=(',3(I5,','),  &
     &        I5,')')
!
!
!     End of zdrvbd
!
      END SUBROUTINE zdrvbd
