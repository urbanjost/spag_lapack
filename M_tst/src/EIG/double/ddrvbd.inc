!*==ddrvbd.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b ddrvbd
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ddrvbd( NSIZES, MM, NN, NTYPES, DOTYPE, ISEED, THRESH,
!                          A, LDA, U, LDU, VT, LDVT, ASAV, USAV, VTSAV, S,
!                          SSAV, E, WORK, LWORK, IWORK, NOUT, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, LDU, LDVT, LWORK, NOUT, NSIZES,
!      $                   NTYPES
!       DOUBLE PRECISION   THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            DOTYPE( * )
!       INTEGER            ISEED( 4 ), IWORK( * ), MM( * ), NN( * )
!       DOUBLE PRECISION   A( LDA, * ), ASAV( LDA, * ), E( * ), S( * ),
!      $                   SSAV( * ), U( LDU, * ), USAV( LDU, * ),
!      $                   VT( LDVT, * ), VTSAV( LDVT, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ddrvbd checks the singular value decomposition (SVD) drivers
!> DGESVD, DGESDD, DGESVDQ, DGESVJ, DGEJSV, and DGESVDX.
!>
!> Both DGESVD and DGESDD factor A = U diag(S) VT, where U and VT are
!> orthogonal and diag(S) is diagonal with the entries of the array S
!> on its diagonal. The entries of S are the singular values,
!> nonnegative and stored in decreasing order.  U and VT can be
!> optionally not computed, overwritten on A, or computed partially.
!>
!> A is M by N. Let MNMIN = min( M, N ). S has dimension MNMIN.
!> U can be M by M or M by MNMIN. VT can be N by N or MNMIN by N.
!>
!> When ddrvbd is called, a number of matrix "sizes" (M's and N's)
!> and a number of matrix "types" are specified.  For each size (M,N)
!> and each type of matrix, and for the minimal workspace as well as
!> workspace adequate to permit blocking, an  M x N  matrix "A" will be
!> generated and used to test the SVD routines.  For each matrix, A will
!> be factored as A = U diag(S) VT and the following 12 tests computed:
!>
!> Test for DGESVD:
!>
!> (1)    | A - U diag(S) VT | / ( |A| max(M,N) ulp )
!>
!> (2)    | I - U'U | / ( M ulp )
!>
!> (3)    | I - VT VT' | / ( N ulp )
!>
!> (4)    S contains MNMIN nonnegative values in decreasing order.
!>        (Return 0 if true, 1/ULP if false.)
!>
!> (5)    | U - Upartial | / ( M ulp ) where Upartial is a partially
!>        computed U.
!>
!> (6)    | VT - VTpartial | / ( N ulp ) where VTpartial is a partially
!>        computed VT.
!>
!> (7)    | S - Spartial | / ( MNMIN ulp |S| ) where Spartial is the
!>        vector of singular values from the partial SVD
!>
!> Test for DGESDD:
!>
!> (8)    | A - U diag(S) VT | / ( |A| max(M,N) ulp )
!>
!> (9)    | I - U'U | / ( M ulp )
!>
!> (10)   | I - VT VT' | / ( N ulp )
!>
!> (11)   S contains MNMIN nonnegative values in decreasing order.
!>        (Return 0 if true, 1/ULP if false.)
!>
!> (12)   | U - Upartial | / ( M ulp ) where Upartial is a partially
!>        computed U.
!>
!> (13)   | VT - VTpartial | / ( N ulp ) where VTpartial is a partially
!>        computed VT.
!>
!> (14)   | S - Spartial | / ( MNMIN ulp |S| ) where Spartial is the
!>        vector of singular values from the partial SVD
!>
!> Test for DGESVDQ:
!>
!> (36)   | A - U diag(S) VT | / ( |A| max(M,N) ulp )
!>
!> (37)   | I - U'U | / ( M ulp )
!>
!> (38)   | I - VT VT' | / ( N ulp )
!>
!> (39)   S contains MNMIN nonnegative values in decreasing order.
!>        (Return 0 if true, 1/ULP if false.)
!>
!> Test for DGESVJ:
!>
!> (15)   | A - U diag(S) VT | / ( |A| max(M,N) ulp )
!>
!> (16)   | I - U'U | / ( M ulp )
!>
!> (17)   | I - VT VT' | / ( N ulp )
!>
!> (18)   S contains MNMIN nonnegative values in decreasing order.
!>        (Return 0 if true, 1/ULP if false.)
!>
!> Test for DGEJSV:
!>
!> (19)   | A - U diag(S) VT | / ( |A| max(M,N) ulp )
!>
!> (20)   | I - U'U | / ( M ulp )
!>
!> (21)   | I - VT VT' | / ( N ulp )
!>
!> (22)   S contains MNMIN nonnegative values in decreasing order.
!>        (Return 0 if true, 1/ULP if false.)
!>
!> Test for DGESVDX( 'V', 'V', 'A' )/DGESVDX( 'N', 'N', 'A' )
!>
!> (23)   | A - U diag(S) VT | / ( |A| max(M,N) ulp )
!>
!> (24)   | I - U'U | / ( M ulp )
!>
!> (25)   | I - VT VT' | / ( N ulp )
!>
!> (26)   S contains MNMIN nonnegative values in decreasing order.
!>        (Return 0 if true, 1/ULP if false.)
!>
!> (27)   | U - Upartial | / ( M ulp ) where Upartial is a partially
!>        computed U.
!>
!> (28)   | VT - VTpartial | / ( N ulp ) where VTpartial is a partially
!>        computed VT.
!>
!> (29)   | S - Spartial | / ( MNMIN ulp |S| ) where Spartial is the
!>        vector of singular values from the partial SVD
!>
!> Test for DGESVDX( 'V', 'V', 'I' )
!>
!> (30)   | U' A VT''' - diag(S) | / ( |A| max(M,N) ulp )
!>
!> (31)   | I - U'U | / ( M ulp )
!>
!> (32)   | I - VT VT' | / ( N ulp )
!>
!> Test for DGESVDX( 'V', 'V', 'V' )
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
!> (3)  A matrix of the form  U D V, where U and V are orthogonal and
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
!>          The number of matrix sizes (M,N) contained in the vectors
!>          MM and NN.
!> \endverbatim
!>
!> \param[in] MM
!> \verbatim
!>          MM is INTEGER array, dimension (NSIZES)
!>          The values of the matrix row dimension M.
!> \endverbatim
!>
!> \param[in] NN
!> \verbatim
!>          NN is INTEGER array, dimension (NSIZES)
!>          The values of the matrix column dimension N.
!> \endverbatim
!>
!> \param[in] NTYPES
!> \verbatim
!>          NTYPES is INTEGER
!>          The number of elements in DOTYPE.   If it is zero, ddrvbd
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
!>          On entry, the seed of the random number generator.  The array
!>          elements should be between 0 and 4095; if not they will be
!>          reduced mod 4096.  Also, ISEED(4) must be odd.
!>          On exit, ISEED is changed and can be used in the next call to
!>          ddrvbd to continue the same random number sequence.
!> \endverbatim
!>
!> \param[in] THRESH
!> \verbatim
!>          THRESH is DOUBLE PRECISION
!>          The threshold value for the test ratios.  A result is
!>          included in the output file if RESULT >= THRESH.  The test
!>          ratios are scaled to be O(1), so THRESH should be a small
!>          multiple of 1, e.g., 10 or 100.  To have every test ratio
!>          printed, use THRESH = 0.
!> \endverbatim
!>
!> \param[out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,NMAX)
!>          where NMAX is the maximum value of N in NN.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,MMAX),
!>          where MMAX is the maximum value of M in MM.
!> \endverbatim
!>
!> \param[out] U
!> \verbatim
!>          U is DOUBLE PRECISION array, dimension (LDU,MMAX)
!> \endverbatim
!>
!> \param[in] LDU
!> \verbatim
!>          LDU is INTEGER
!>          The leading dimension of the array U.  LDU >= max(1,MMAX).
!> \endverbatim
!>
!> \param[out] VT
!> \verbatim
!>          VT is DOUBLE PRECISION array, dimension (LDVT,NMAX)
!> \endverbatim
!>
!> \param[in] LDVT
!> \verbatim
!>          LDVT is INTEGER
!>          The leading dimension of the array VT.  LDVT >= max(1,NMAX).
!> \endverbatim
!>
!> \param[out] ASAV
!> \verbatim
!>          ASAV is DOUBLE PRECISION array, dimension (LDA,NMAX)
!> \endverbatim
!>
!> \param[out] USAV
!> \verbatim
!>          USAV is DOUBLE PRECISION array, dimension (LDU,MMAX)
!> \endverbatim
!>
!> \param[out] VTSAV
!> \verbatim
!>          VTSAV is DOUBLE PRECISION array, dimension (LDVT,NMAX)
!> \endverbatim
!>
!> \param[out] S
!> \verbatim
!>          S is DOUBLE PRECISION array, dimension
!>                      (max(min(MM,NN)))
!> \endverbatim
!>
!> \param[out] SSAV
!> \verbatim
!>          SSAV is DOUBLE PRECISION array, dimension
!>                      (max(min(MM,NN)))
!> \endverbatim
!>
!> \param[out] E
!> \verbatim
!>          E is DOUBLE PRECISION array, dimension
!>                      (max(min(MM,NN)))
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (LWORK)
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The number of entries in WORK.  This must be at least
!>          max(3*MN+MX,5*MN-4)+2*MN**2 for all pairs
!>          pairs  (MN,MX)=( min(MM(j),NN(j), max(MM(j),NN(j)) )
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension at least 8*min(M,N)
!> \endverbatim
!>
!> \param[in] NOUT
!> \verbatim
!>          NOUT is INTEGER
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
!>          If  dlatms, or DGESVD returns an error code, the
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
!> \ingroup double_eig
!
!  =====================================================================
      SUBROUTINE ddrvbd(Nsizes,Mm,Nn,Ntypes,Dotype,Iseed,Thresh,A,Lda,U,&
     &                  Ldu,Vt,Ldvt,Asav,Usav,Vtsav,S,Ssav,E,Work,Lwork,&
     &                  Iwork,Nout,Info)
!
      IMPLICIT NONE
!*--DDRVBD371
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2016
!
!     .. Scalar Arguments ..
      INTEGER Info , Lda , Ldu , Ldvt , Lwork , Nout , Nsizes , Ntypes
      DOUBLE PRECISION Thresh
!     ..
!     .. Array Arguments ..
      LOGICAL Dotype(*)
      INTEGER Iseed(4) , Iwork(*) , Mm(*) , Nn(*)
      DOUBLE PRECISION A(Lda,*) , Asav(Lda,*) , E(*) , S(*) , Ssav(*) , &
     &                 U(Ldu,*) , Usav(Ldu,*) , Vt(Ldvt,*) ,            &
     &                 Vtsav(Ldvt,*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE , TWO , HALF
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,HALF=0.5D0)
      INTEGER MAXTYP
      PARAMETER (MAXTYP=5)
!     ..
!     .. Local Scalars ..
      LOGICAL badmm , badnn
      CHARACTER jobq , jobu , jobvt , range
      CHARACTER*3 path
      INTEGER i , iinfo , ijq , iju , ijvt , il , iu , iws , iwtmp ,    &
     &        itemp , j , jsize , jtype , lswork , m , minwrk , mmax ,  &
     &        mnmax , mnmin , mtypes , n , nfail , nmax , ns , nsi ,    &
     &        nsv , ntest
      DOUBLE PRECISION anorm , dif , div , ovfl , rtunfl , ulp ,        &
     &                 ulpinv , unfl , vl , vu
!     ..
!     .. Local Scalars for DGESVDQ ..
      INTEGER liwork , lrwork , numrank
!     ..
!     .. Local Arrays for DGESVDQ ..
      DOUBLE PRECISION rwork(2)
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
      EXTERNAL DGEJSV , DGESDD , DGESVD , DGESVDQ 
     EXTERNAL         DGESVDX , DGESVJ , DLABAD , DLACPY , DLASET 
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , DBLE , INT , MAX , MIN
!     ..
!     .. Scalars in Common ..
      LOGICAL LERr , OK
      CHARACTER*32 SRNamt
      INTEGER INFot , NUNit
!     ..
!     .. Common blocks ..
      COMMON /INFOC / INFot , NUNit , OK , LERr
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
         minwrk = MAX(minwrk,                                           &
     &            MAX(3*MIN(Mm(j),Nn(j))+MAX(Mm(j),Nn(j)),5*MIN(Mm(j),  &
     &            Nn(j)-4))+2*MIN(Mm(j),Nn(j))**2)
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
         CALL xerbla('ddrvbd',-Info)
         RETURN
      ENDIF
!
!     Initialize constants
!
      path(1:1) = 'Double precision'
      path(2:3) = 'BD'
      nfail = 0
      ntest = 0
      unfl = DLAMCH('Safe minimum')
      ovfl = ONE/unfl
      CALL DLABAD(unfl,ovfl)
      ulp = DLAMCH('Precision')
      rtunfl = SQRT(unfl)
      ulpinv = ONE/ulp
      INFot = 0
!
!     Loop over sizes, types
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
                     CALL DLASET('Full',m,n,ZERO,ZERO,A,Lda)
!
                  ELSEIF ( jtype==2 ) THEN
!
!              Identity matrix
!
                     CALL DLASET('Full',m,n,ZERO,ONE,A,Lda)
!
                  ELSE
!
!              (Scaled) random matrix
!
                     IF ( jtype==3 ) anorm = ONE
                     IF ( jtype==4 ) anorm = unfl/ulp
                     IF ( jtype==5 ) anorm = ovfl*ulp
                     CALL dlatms(m,n,'U',Iseed,'N',S,4,DBLE(mnmin),     &
     &                           anorm,m-1,n-1,'N',A,Lda,Work,iinfo)
                     IF ( iinfo/=0 ) THEN
                        WRITE (Nout,FMT=99004) 'Generator' , iinfo , m ,&
     &                         n , jtype , ioldsd
                        Info = ABS(iinfo)
                        RETURN
                     ENDIF
                  ENDIF
               ENDIF
!
               CALL DLACPY('F',m,n,A,Lda,Asav,Lda)
!
!           Do for minimal and adequate (for blocking) workspace
!
               DO iws = 1 , 4
!
                  DO j = 1 , 32
                     result(j) = -ONE
                  ENDDO
!
!              Test DGESVD: Factorize A
!
                  iwtmp = MAX(3*MIN(m,n)+MAX(m,n),5*MIN(m,n))
                  lswork = iwtmp + (iws-1)*(Lwork-iwtmp)/3
                  lswork = MIN(lswork,Lwork)
                  lswork = MAX(lswork,1)
                  IF ( iws==4 ) lswork = Lwork
!
                  IF ( iws>1 ) CALL DLACPY('F',m,n,Asav,Lda,A,Lda)
                  SRNamt = 'DGESVD'
                  CALL DGESVD('A','A',m,n,A,Lda,Ssav,Usav,Ldu,Vtsav,    &
     &                        Ldvt,Work,lswork,iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nout,FMT=99005) 'GESVD' , iinfo , m , n ,   &
     &                      jtype , lswork , ioldsd
                     Info = ABS(iinfo)
                     RETURN
                  ENDIF
!
!              Do tests 1--4
!
                  CALL dbdt01(m,n,0,Asav,Lda,Usav,Ldu,Ssav,E,Vtsav,Ldvt,&
     &                        Work,result(1))
                  IF ( m/=0 .AND. n/=0 ) THEN
                     CALL dort01('Columns',m,m,Usav,Ldu,Work,Lwork,     &
     &                           result(2))
                     CALL dort01('Rows',n,n,Vtsav,Ldvt,Work,Lwork,      &
     &                           result(3))
                  ENDIF
                  result(4) = ZERO
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
                           CALL DLACPY('F',m,n,Asav,Lda,A,Lda)
                           SRNamt = 'DGESVD'
                           CALL DGESVD(jobu,jobvt,m,n,A,Lda,S,U,Ldu,Vt, &
     &                                 Ldvt,Work,lswork,iinfo)
!
!                    Compare U
!
                           dif = ZERO
                           IF ( m>0 .AND. n>0 ) THEN
                              IF ( iju==1 ) THEN
                                 CALL dort03('C',m,mnmin,m,mnmin,Usav,  &
     &                              Ldu,A,Lda,Work,Lwork,dif,iinfo)
                              ELSEIF ( iju==2 ) THEN
                                 CALL dort03('C',m,mnmin,m,mnmin,Usav,  &
     &                              Ldu,U,Ldu,Work,Lwork,dif,iinfo)
                              ELSEIF ( iju==3 ) THEN
                                 CALL dort03('C',m,m,m,mnmin,Usav,Ldu,U,&
     &                              Ldu,Work,Lwork,dif,iinfo)
                              ENDIF
                           ENDIF
                           result(5) = MAX(result(5),dif)
!
!                    Compare VT
!
                           dif = ZERO
                           IF ( m>0 .AND. n>0 ) THEN
                              IF ( ijvt==1 ) THEN
                                 CALL dort03('R',n,mnmin,n,mnmin,Vtsav, &
     &                              Ldvt,A,Lda,Work,Lwork,dif,iinfo)
                              ELSEIF ( ijvt==2 ) THEN
                                 CALL dort03('R',n,mnmin,n,mnmin,Vtsav, &
     &                              Ldvt,Vt,Ldvt,Work,Lwork,dif,iinfo)
                              ELSEIF ( ijvt==3 ) THEN
                                 CALL dort03('R',n,n,n,mnmin,Vtsav,Ldvt,&
     &                              Vt,Ldvt,Work,Lwork,dif,iinfo)
                              ENDIF
                           ENDIF
                           result(6) = MAX(result(6),dif)
!
!                    Compare S
!
                           dif = ZERO
                           div = MAX(mnmin*ulp*S(1),unfl)
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
!              Test DGESDD: Factorize A
!
                  iwtmp = 5*mnmin*mnmin + 9*mnmin + MAX(m,n)
                  lswork = iwtmp + (iws-1)*(Lwork-iwtmp)/3
                  lswork = MIN(lswork,Lwork)
                  lswork = MAX(lswork,1)
                  IF ( iws==4 ) lswork = Lwork
!
                  CALL DLACPY('F',m,n,Asav,Lda,A,Lda)
                  SRNamt = 'DGESDD'
                  CALL DGESDD('A',m,n,A,Lda,Ssav,Usav,Ldu,Vtsav,Ldvt,   &
     &                        Work,lswork,Iwork,iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nout,FMT=99005) 'GESDD' , iinfo , m , n ,   &
     &                      jtype , lswork , ioldsd
                     Info = ABS(iinfo)
                     RETURN
                  ENDIF
!
!              Do tests 8--11
!
                  CALL dbdt01(m,n,0,Asav,Lda,Usav,Ldu,Ssav,E,Vtsav,Ldvt,&
     &                        Work,result(8))
                  IF ( m/=0 .AND. n/=0 ) THEN
                     CALL dort01('Columns',m,m,Usav,Ldu,Work,Lwork,     &
     &                           result(9))
                     CALL dort01('Rows',n,n,Vtsav,Ldvt,Work,Lwork,      &
     &                           result(10))
                  ENDIF
                  result(11) = ZERO
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
                     CALL DLACPY('F',m,n,Asav,Lda,A,Lda)
                     SRNamt = 'DGESDD'
                     CALL DGESDD(jobq,m,n,A,Lda,S,U,Ldu,Vt,Ldvt,Work,   &
     &                           lswork,Iwork,iinfo)
!
!                 Compare U
!
                     dif = ZERO
                     IF ( m>0 .AND. n>0 ) THEN
                        IF ( ijq==1 ) THEN
                           IF ( m>=n ) THEN
                              CALL dort03('C',m,mnmin,m,mnmin,Usav,Ldu, &
     &                           A,Lda,Work,Lwork,dif,Info)
                           ELSE
                              CALL dort03('C',m,mnmin,m,mnmin,Usav,Ldu, &
     &                           U,Ldu,Work,Lwork,dif,Info)
                           ENDIF
                        ELSEIF ( ijq==2 ) THEN
                           CALL dort03('C',m,mnmin,m,mnmin,Usav,Ldu,U,  &
     &                                 Ldu,Work,Lwork,dif,Info)
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
                              CALL dort03('R',n,mnmin,n,mnmin,Vtsav,    &
     &                           Ldvt,Vt,Ldvt,Work,Lwork,dif,Info)
                           ELSE
                              CALL dort03('R',n,mnmin,n,mnmin,Vtsav,    &
     &                           Ldvt,A,Lda,Work,Lwork,dif,Info)
                           ENDIF
                        ELSEIF ( ijq==2 ) THEN
                           CALL dort03('R',n,mnmin,n,mnmin,Vtsav,Ldvt,  &
     &                                 Vt,Ldvt,Work,Lwork,dif,Info)
                        ENDIF
                     ENDIF
                     result(13) = MAX(result(13),dif)
!
!                 Compare S
!
                     dif = ZERO
                     div = MAX(mnmin*ulp*S(1),unfl)
                     DO i = 1 , mnmin - 1
                        IF ( Ssav(i)<Ssav(i+1) ) dif = ulpinv
                        IF ( Ssav(i)<ZERO ) dif = ulpinv
                        dif = MAX(dif,ABS(Ssav(i)-S(i))/div)
                     ENDDO
                     result(14) = MAX(result(14),dif)
                  ENDDO
!
!              Test DGESVDQ
!              Note: DGESVDQ only works for M >= N
!
                  result(36) = ZERO
                  result(37) = ZERO
                  result(38) = ZERO
                  result(39) = ZERO
!
                  IF ( m>=n ) THEN
                     iwtmp = 5*mnmin*mnmin + 9*mnmin + MAX(m,n)
                     lswork = iwtmp + (iws-1)*(Lwork-iwtmp)/3
                     lswork = MIN(lswork,Lwork)
                     lswork = MAX(lswork,1)
                     IF ( iws==4 ) lswork = Lwork
!
                     CALL DLACPY('F',m,n,Asav,Lda,A,Lda)
                     SRNamt = 'DGESVDQ'
!
                     lrwork = 2
                     liwork = MAX(n,1)
                     CALL DGESVDQ('H','N','N','A','A',m,n,A,Lda,Ssav,   &
     &                            Usav,Ldu,Vtsav,Ldvt,numrank,Iwork,    &
     &                            liwork,Work,Lwork,rwork,lrwork,iinfo)
!
                     IF ( iinfo/=0 ) THEN
                        WRITE (Nout,FMT=99005) 'DGESVDQ' , iinfo , m ,  &
     &                         n , jtype , lswork , ioldsd
                        Info = ABS(iinfo)
                        RETURN
                     ENDIF
!
!                 Do tests 36--39
!
                     CALL dbdt01(m,n,0,Asav,Lda,Usav,Ldu,Ssav,E,Vtsav,  &
     &                           Ldvt,Work,result(36))
                     IF ( m/=0 .AND. n/=0 ) THEN
                        CALL dort01('Columns',m,m,Usav,Ldu,Work,Lwork,  &
     &                              result(37))
                        CALL dort01('Rows',n,n,Vtsav,Ldvt,Work,Lwork,   &
     &                              result(38))
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
!              Test DGESVJ
!              Note: DGESVJ only works for M >= N
!
                  result(15) = ZERO
                  result(16) = ZERO
                  result(17) = ZERO
                  result(18) = ZERO
!
                  IF ( m>=n ) THEN
                     iwtmp = 5*mnmin*mnmin + 9*mnmin + MAX(m,n)
                     lswork = iwtmp + (iws-1)*(Lwork-iwtmp)/3
                     lswork = MIN(lswork,Lwork)
                     lswork = MAX(lswork,1)
                     IF ( iws==4 ) lswork = Lwork
!
                     CALL DLACPY('F',m,n,Asav,Lda,Usav,Lda)
                     SRNamt = 'DGESVJ'
                     CALL DGESVJ('G','U','V',m,n,Usav,Lda,Ssav,0,A,Ldvt,&
     &                           Work,Lwork,Info)
!
!                 DGESVJ returns V not VT
!
                     DO j = 1 , n
                        DO i = 1 , n
                           Vtsav(j,i) = A(i,j)
                        ENDDO
                     ENDDO
!
                     IF ( iinfo/=0 ) THEN
                        WRITE (Nout,FMT=99005) 'GESVJ' , iinfo , m , n ,&
     &                         jtype , lswork , ioldsd
                        Info = ABS(iinfo)
                        RETURN
                     ENDIF
!
!                 Do tests 15--18
!
                     CALL dbdt01(m,n,0,Asav,Lda,Usav,Ldu,Ssav,E,Vtsav,  &
     &                           Ldvt,Work,result(15))
                     IF ( m/=0 .AND. n/=0 ) THEN
                        CALL dort01('Columns',m,m,Usav,Ldu,Work,Lwork,  &
     &                              result(16))
                        CALL dort01('Rows',n,n,Vtsav,Ldvt,Work,Lwork,   &
     &                              result(17))
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
!              Test DGEJSV
!              Note: DGEJSV only works for M >= N
!
                  result(19) = ZERO
                  result(20) = ZERO
                  result(21) = ZERO
                  result(22) = ZERO
                  IF ( m>=n ) THEN
                     iwtmp = 5*mnmin*mnmin + 9*mnmin + MAX(m,n)
                     lswork = iwtmp + (iws-1)*(Lwork-iwtmp)/3
                     lswork = MIN(lswork,Lwork)
                     lswork = MAX(lswork,1)
                     IF ( iws==4 ) lswork = Lwork
!
                     CALL DLACPY('F',m,n,Asav,Lda,Vtsav,Lda)
                     SRNamt = 'DGEJSV'
                     CALL DGEJSV('G','U','V','R','N','N',m,n,Vtsav,Lda, &
     &                           Ssav,Usav,Ldu,A,Ldvt,Work,Lwork,Iwork, &
     &                           Info)
!
!                 DGEJSV returns V not VT
!
                     DO j = 1 , n
                        DO i = 1 , n
                           Vtsav(j,i) = A(i,j)
                        ENDDO
                     ENDDO
!
                     IF ( iinfo/=0 ) THEN
                        WRITE (Nout,FMT=99005) 'GEJSV' , iinfo , m , n ,&
     &                         jtype , lswork , ioldsd
                        Info = ABS(iinfo)
                        RETURN
                     ENDIF
!
!                 Do tests 19--22
!
                     CALL dbdt01(m,n,0,Asav,Lda,Usav,Ldu,Ssav,E,Vtsav,  &
     &                           Ldvt,Work,result(19))
                     IF ( m/=0 .AND. n/=0 ) THEN
                        CALL dort01('Columns',m,m,Usav,Ldu,Work,Lwork,  &
     &                              result(20))
                        CALL dort01('Rows',n,n,Vtsav,Ldvt,Work,Lwork,   &
     &                              result(21))
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
!              Test DGESVDX
!
                  CALL DLACPY('F',m,n,Asav,Lda,A,Lda)
                  CALL DGESVDX('V','V','A',m,n,A,Lda,vl,vu,il,iu,ns,    &
     &                         Ssav,Usav,Ldu,Vtsav,Ldvt,Work,Lwork,     &
     &                         Iwork,iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nout,FMT=99005) 'GESVDX' , iinfo , m , n ,  &
     &                      jtype , lswork , ioldsd
                     Info = ABS(iinfo)
                     RETURN
                  ENDIF
!
!              Do tests 23--29
!
                  result(23) = ZERO
                  result(24) = ZERO
                  result(25) = ZERO
                  CALL dbdt01(m,n,0,Asav,Lda,Usav,Ldu,Ssav,E,Vtsav,Ldvt,&
     &                        Work,result(23))
                  IF ( m/=0 .AND. n/=0 ) THEN
                     CALL dort01('Columns',m,m,Usav,Ldu,Work,Lwork,     &
     &                           result(24))
                     CALL dort01('Rows',n,n,Vtsav,Ldvt,Work,Lwork,      &
     &                           result(25))
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
                           CALL DLACPY('F',m,n,Asav,Lda,A,Lda)
                           CALL DGESVDX(jobu,jobvt,range,m,n,A,Lda,vl,  &
     &                                  vu,il,iu,ns,S,U,Ldu,Vt,Ldvt,    &
     &                                  Work,Lwork,Iwork,iinfo)
!
!                    Compare U
!
                           dif = ZERO
                           IF ( m>0 .AND. n>0 ) THEN
                              IF ( iju==1 )                             &
     &                             CALL dort03('C',m,mnmin,m,mnmin,Usav,&
     &                             Ldu,U,Ldu,Work,Lwork,dif,iinfo)
                           ENDIF
                           result(27) = MAX(result(27),dif)
!
!                    Compare VT
!
                           dif = ZERO
                           IF ( m>0 .AND. n>0 ) THEN
                              IF ( ijvt==1 )                            &
     &                             CALL dort03('R',n,mnmin,n,mnmin,     &
     &                             Vtsav,Ldvt,Vt,Ldvt,Work,Lwork,dif,   &
     &                             iinfo)
                           ENDIF
                           result(28) = MAX(result(28),dif)
!
!                    Compare S
!
                           dif = ZERO
                           div = MAX(mnmin*ulp*S(1),unfl)
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
!              Do tests 30--32: DGESVDX( 'V', 'V', 'I' )
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
                  CALL DLACPY('F',m,n,Asav,Lda,A,Lda)
                  CALL DGESVDX('V','V','I',m,n,A,Lda,vl,vu,il,iu,nsi,S, &
     &                         U,Ldu,Vt,Ldvt,Work,Lwork,Iwork,iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nout,FMT=99005) 'GESVDX' , iinfo , m , n ,  &
     &                      jtype , lswork , ioldsd
                     Info = ABS(iinfo)
                     RETURN
                  ENDIF
!
                  result(30) = ZERO
                  result(31) = ZERO
                  result(32) = ZERO
                  CALL dbdt05(m,n,Asav,Lda,S,nsi,U,Ldu,Vt,Ldvt,Work,    &
     &                        result(30))
                  CALL dort01('Columns',m,nsi,U,Ldu,Work,Lwork,         &
     &                        result(31))
                  CALL dort01('Rows',nsi,n,Vt,Ldvt,Work,Lwork,result(32)&
     &                        )
!
!              Do tests 33--35: DGESVDX( 'V', 'V', 'V' )
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
                  CALL DLACPY('F',m,n,Asav,Lda,A,Lda)
                  CALL DGESVDX('V','V','V',m,n,A,Lda,vl,vu,il,iu,nsv,S, &
     &                         U,Ldu,Vt,Ldvt,Work,Lwork,Iwork,iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nout,FMT=99005) 'GESVDX' , iinfo , m , n ,  &
     &                      jtype , lswork , ioldsd
                     Info = ABS(iinfo)
                     RETURN
                  ENDIF
!
                  result(33) = ZERO
                  result(34) = ZERO
                  result(35) = ZERO
                  CALL dbdt05(m,n,Asav,Lda,S,nsv,U,Ldu,Vt,Ldvt,Work,    &
     &                        result(33))
                  CALL dort01('Columns',m,nsv,U,Ldu,Work,Lwork,         &
     &                        result(34))
                  CALL dort01('Rows',nsv,n,Vt,Ldvt,Work,Lwork,result(35)&
     &                        )
!
!              End of Loop -- Check for RESULT(j) > THRESH
!
                  DO j = 1 , 39
                     IF ( result(j)>=Thresh ) THEN
                        IF ( nfail==0 ) THEN
                           WRITE (Nout,FMT=99001)
                           WRITE (Nout,FMT=99002)
                        ENDIF
                        WRITE (Nout,FMT=99003) m , n , jtype , iws ,    &
     &                         ioldsd , j , result(j)
                        nfail = nfail + 1
                     ENDIF
                  ENDDO
                  ntest = ntest + 39
               ENDDO
            ENDIF
         ENDDO
      ENDDO
!
!     Summary
!
      CALL alasvm(path,Nout,nfail,ntest,0)
!
99001 FORMAT (' SVD -- Real Singular Value Decomposition Driver ',      &
     &        /' Matrix types (see ddrvbd for details):',               &
     &        //' 1 = Zero matrix',/' 2 = Identity matrix',             &
     &        /' 3 = Evenly spaced singular values near 1',             &
     &        /' 4 = Evenly spaced singular values near underflow',     &
     &        /' 5 = Evenly spaced singular values near overflow',//    &
     &        ' Tests performed: ( A is dense, U and V are orthogonal,',&
     &        /19X,' S is an array, and Upartial, VTpartial, and',/19X, &
     &        ' Spartial are partially computed U, VT and S),',/)
99002 FORMAT (' 1 = | A - U diag(S) VT | / ( |A| max(M,N) ulp ) ',      &
     &        /' 2 = | I - U**T U | / ( M ulp ) ',                      &
     &        /' 3 = | I - VT VT**T | / ( N ulp ) ',                    &
     &        /' 4 = 0 if S contains min(M,N) nonnegative values in',   &
     &        ' decreasing order, else 1/ulp',                          &
     &        /' 5 = | U - Upartial | / ( M ulp )',                     &
     &        /' 6 = | VT - VTpartial | / ( N ulp )',                   &
     &        /' 7 = | S - Spartial | / ( min(M,N) ulp |S| )',          &
     &        /' 8 = | A - U diag(S) VT | / ( |A| max(M,N) ulp ) ',     &
     &        /' 9 = | I - U**T U | / ( M ulp ) ',                      &
     &        /'10 = | I - VT VT**T | / ( N ulp ) ',                    &
     &        /'11 = 0 if S contains min(M,N) nonnegative values in',   &
     &        ' decreasing order, else 1/ulp',                          &
     &        /'12 = | U - Upartial | / ( M ulp )',                     &
     &        /'13 = | VT - VTpartial | / ( N ulp )',                   &
     &        /'14 = | S - Spartial | / ( min(M,N) ulp |S| )',          &
     &        /'15 = | A - U diag(S) VT | / ( |A| max(M,N) ulp ) ',     &
     &        /'16 = | I - U**T U | / ( M ulp ) ',                      &
     &        /'17 = | I - VT VT**T | / ( N ulp ) ',                    &
     &        /'18 = 0 if S contains min(M,N) nonnegative values in',   &
     &        ' decreasing order, else 1/ulp',                          &
     &        /'19 = | U - Upartial | / ( M ulp )',                     &
     &        /'20 = | VT - VTpartial | / ( N ulp )',                   &
     &        /'21 = | S - Spartial | / ( min(M,N) ulp |S| )',          &
     &        /'22 = 0 if S contains min(M,N) nonnegative values in',   &
     &        ' decreasing order, else 1/ulp',                          &
     &        /'23 = | A - U diag(S) VT | / ( |A| max(M,N) ulp ),',     &
     &        ' DGESVDX(V,V,A) ',/'24 = | I - U**T U | / ( M ulp ) ',   &
     &        /'25 = | I - VT VT**T | / ( N ulp ) ',                    &
     &        /'26 = 0 if S contains min(M,N) nonnegative values in',   &
     &        ' decreasing order, else 1/ulp',                          &
     &        /'27 = | U - Upartial | / ( M ulp )',                     &
     &        /'28 = | VT - VTpartial | / ( N ulp )',                   &
     &        /'29 = | S - Spartial | / ( min(M,N) ulp |S| )',          &
     &        /'30 = | U**T A VT**T - diag(S) | / ( |A| max(M,N) ulp ),'&
     &        ,' DGESVDX(V,V,I) ',/'31 = | I - U**T U | / ( M ulp ) ',  &
     &        /'32 = | I - VT VT**T | / ( N ulp ) ',                    &
     &        /'33 = | U**T A VT**T - diag(S) | / ( |A| max(M,N) ulp ),'&
     &        ,' DGESVDX(V,V,V) ',/'34 = | I - U**T U | / ( M ulp ) ',  &
     &        /'35 = | I - VT VT**T | / ( N ulp ) ',                    &
     &        ' DGESVDQ(H,N,N,A,A',                                     &
     &        /'36 = | A - U diag(S) VT | / ( |A| max(M,N) ulp ) ',     &
     &        /'37 = | I - U**T U | / ( M ulp ) ',                      &
     &        /'38 = | I - VT VT**T | / ( N ulp ) ',                    &
     &        /'39 = 0 if S contains min(M,N) nonnegative values in',   &
     &        ' decreasing order, else 1/ulp',//)
99003 FORMAT (' M=',I5,', N=',I5,', type ',I1,', IWS=',I1,', seed=',    &
     &        4(I4,','),' test(',I2,')=',G11.4)
99004 FORMAT (' ddrvbd: ',A,' returned INFO=',I6,'.',/9X,'M=',I6,', N=',&
     &        I6,', JTYPE=',I6,', ISEED=(',3(I5,','),I5,')')
99005 FORMAT (' ddrvbd: ',A,' returned INFO=',I6,'.',/9X,'M=',I6,', N=',&
     &        I6,', JTYPE=',I6,', LSWORK=',I6,/9X,'ISEED=(',3(I5,','),  &
     &        I5,')')
!
!
!     End of ddrvbd
!
      END SUBROUTINE ddrvbd
