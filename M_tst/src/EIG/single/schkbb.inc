!*==schkbb.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b schkbb
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE schkbb( NSIZES, MVAL, NVAL, NWDTHS, KK, NTYPES, DOTYPE,
!                          NRHS, ISEED, THRESH, NOUNIT, A, LDA, AB, LDAB,
!                          BD, BE, Q, LDQ, P, LDP, C, LDC, CC, WORK,
!                          LWORK, RESULT, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, LDAB, LDC, LDP, LDQ, LWORK, NOUNIT,
!      $                   NRHS, NSIZES, NTYPES, NWDTHS
!       REAL               THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            DOTYPE( * )
!       INTEGER            ISEED( 4 ), KK( * ), MVAL( * ), NVAL( * )
!       REAL               A( LDA, * ), AB( LDAB, * ), BD( * ), BE( * ),
!      $                   C( LDC, * ), CC( LDC, * ), P( LDP, * ),
!      $                   Q( LDQ, * ), RESULT( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> schkbb tests the reduction of a general real rectangular band
!> matrix to bidiagonal form.
!>
!> SGBBRD factors a general band matrix A as  Q B P* , where * means
!> transpose, B is upper bidiagonal, and Q and P are orthogonal;
!> SGBBRD can also overwrite a given matrix C with Q* C .
!>
!> For each pair of matrix dimensions (M,N) and each selected matrix
!> type, an M by N matrix A and an M by NRHS matrix C are generated.
!> The problem dimensions are as follows
!>    A:          M x N
!>    Q:          M x M
!>    P:          N x N
!>    B:          min(M,N) x min(M,N)
!>    C:          M x NRHS
!>
!> For each generated matrix, 4 tests are performed:
!>
!> (1)   | A - Q B PT | / ( |A| max(M,N) ulp ), PT = P'
!>
!> (2)   | I - Q' Q | / ( M ulp )
!>
!> (3)   | I - PT PT' | / ( N ulp )
!>
!> (4)   | Y - Q' C | / ( |Y| max(M,NRHS) ulp ), where Y = Q' C.
!>
!> The "types" are specified by a logical array DOTYPE( 1:NTYPES );
!> if DOTYPE(j) is .TRUE., then matrix type "j" will be generated.
!> Currently, the list of possible types is:
!>
!> The possible matrix types are
!>
!> (1)  The zero matrix.
!> (2)  The identity matrix.
!>
!> (3)  A diagonal matrix with evenly spaced entries
!>      1, ..., ULP  and random signs.
!>      (ULP = (first number larger than 1) - 1 )
!> (4)  A diagonal matrix with geometrically spaced entries
!>      1, ..., ULP  and random signs.
!> (5)  A diagonal matrix with "clustered" entries 1, ULP, ..., ULP
!>      and random signs.
!>
!> (6)  Same as (3), but multiplied by SQRT( overflow threshold )
!> (7)  Same as (3), but multiplied by SQRT( underflow threshold )
!>
!> (8)  A matrix of the form  U D V, where U and V are orthogonal and
!>      D has evenly spaced entries 1, ..., ULP with random signs
!>      on the diagonal.
!>
!> (9)  A matrix of the form  U D V, where U and V are orthogonal and
!>      D has geometrically spaced entries 1, ..., ULP with random
!>      signs on the diagonal.
!>
!> (10) A matrix of the form  U D V, where U and V are orthogonal and
!>      D has "clustered" entries 1, ULP,..., ULP with random
!>      signs on the diagonal.
!>
!> (11) Same as (8), but multiplied by SQRT( overflow threshold )
!> (12) Same as (8), but multiplied by SQRT( underflow threshold )
!>
!> (13) Rectangular matrix with random entries chosen from (-1,1).
!> (14) Same as (13), but multiplied by SQRT( overflow threshold )
!> (15) Same as (13), but multiplied by SQRT( underflow threshold )
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] NSIZES
!> \verbatim
!>          NSIZES is INTEGER
!>          The number of values of M and N contained in the vectors
!>          MVAL and NVAL.  The matrix sizes are used in pairs (M,N).
!>          If NSIZES is zero, schkbb does nothing.  NSIZES must be at
!>          least zero.
!> \endverbatim
!>
!> \param[in] MVAL
!> \verbatim
!>          MVAL is INTEGER array, dimension (NSIZES)
!>          The values of the matrix row dimension M.
!> \endverbatim
!>
!> \param[in] NVAL
!> \verbatim
!>          NVAL is INTEGER array, dimension (NSIZES)
!>          The values of the matrix column dimension N.
!> \endverbatim
!>
!> \param[in] NWDTHS
!> \verbatim
!>          NWDTHS is INTEGER
!>          The number of bandwidths to use.  If it is zero,
!>          schkbb does nothing.  It must be at least zero.
!> \endverbatim
!>
!> \param[in] KK
!> \verbatim
!>          KK is INTEGER array, dimension (NWDTHS)
!>          An array containing the bandwidths to be used for the band
!>          matrices.  The values must be at least zero.
!> \endverbatim
!>
!> \param[in] NTYPES
!> \verbatim
!>          NTYPES is INTEGER
!>          The number of elements in DOTYPE.   If it is zero, schkbb
!>          does nothing.  It must be at least zero.  If it is MAXTYP+1
!>          and NSIZES is 1, then an additional type, MAXTYP+1 is
!>          defined, which is to use whatever matrix is in A.  This
!>          is only useful if DOTYPE(1:MAXTYP) is .FALSE. and
!>          DOTYPE(MAXTYP+1) is .TRUE. .
!> \endverbatim
!>
!> \param[in] DOTYPE
!> \verbatim
!>          DOTYPE is LOGICAL array, dimension (NTYPES)
!>          If DOTYPE(j) is .TRUE., then for each size in NN a
!>          matrix of that size and of type j will be generated.
!>          If NTYPES is smaller than the maximum number of types
!>          defined (PARAMETER MAXTYP), then types NTYPES+1 through
!>          MAXTYP will not be generated.  If NTYPES is larger
!>          than MAXTYP, DOTYPE(MAXTYP+1) through DOTYPE(NTYPES)
!>          will be ignored.
!> \endverbatim
!>
!> \param[in] NRHS
!> \verbatim
!>          NRHS is INTEGER
!>          The number of columns in the "right-hand side" matrix C.
!>          If NRHS = 0, then the operations on the right-hand side will
!>          not be tested. NRHS must be at least 0.
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
!>          next call to schkbb to continue the same random number
!>          sequence.
!> \endverbatim
!>
!> \param[in] THRESH
!> \verbatim
!>          THRESH is REAL
!>          A test will count as "failed" if the "error", computed as
!>          described above, exceeds THRESH.  Note that the error
!>          is scaled to be O(1), so THRESH should be a reasonably
!>          small multiple of 1, e.g., 10 or 100.  In particular,
!>          it should not depend on the precision (single vs. double)
!>          or the size of the matrix.  It must be at least zero.
!> \endverbatim
!>
!> \param[in] NOUNIT
!> \verbatim
!>          NOUNIT is INTEGER
!>          The FORTRAN unit number for printing out error messages
!>          (e.g., if a routine returns IINFO not equal to 0.)
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is REAL array, dimension
!>                            (LDA, max(NN))
!>          Used to hold the matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of A.  It must be at least 1
!>          and at least max( NN ).
!> \endverbatim
!>
!> \param[out] AB
!> \verbatim
!>          AB is REAL array, dimension (LDAB, max(NN))
!>          Used to hold A in band storage format.
!> \endverbatim
!>
!> \param[in] LDAB
!> \verbatim
!>          LDAB is INTEGER
!>          The leading dimension of AB.  It must be at least 2 (not 1!)
!>          and at least max( KK )+1.
!> \endverbatim
!>
!> \param[out] BD
!> \verbatim
!>          BD is REAL array, dimension (max(NN))
!>          Used to hold the diagonal of the bidiagonal matrix computed
!>          by SGBBRD.
!> \endverbatim
!>
!> \param[out] BE
!> \verbatim
!>          BE is REAL array, dimension (max(NN))
!>          Used to hold the off-diagonal of the bidiagonal matrix
!>          computed by SGBBRD.
!> \endverbatim
!>
!> \param[out] Q
!> \verbatim
!>          Q is REAL array, dimension (LDQ, max(NN))
!>          Used to hold the orthogonal matrix Q computed by SGBBRD.
!> \endverbatim
!>
!> \param[in] LDQ
!> \verbatim
!>          LDQ is INTEGER
!>          The leading dimension of Q.  It must be at least 1
!>          and at least max( NN ).
!> \endverbatim
!>
!> \param[out] P
!> \verbatim
!>          P is REAL array, dimension (LDP, max(NN))
!>          Used to hold the orthogonal matrix P computed by SGBBRD.
!> \endverbatim
!>
!> \param[in] LDP
!> \verbatim
!>          LDP is INTEGER
!>          The leading dimension of P.  It must be at least 1
!>          and at least max( NN ).
!> \endverbatim
!>
!> \param[out] C
!> \verbatim
!>          C is REAL array, dimension (LDC, max(NN))
!>          Used to hold the matrix C updated by SGBBRD.
!> \endverbatim
!>
!> \param[in] LDC
!> \verbatim
!>          LDC is INTEGER
!>          The leading dimension of U.  It must be at least 1
!>          and at least max( NN ).
!> \endverbatim
!>
!> \param[out] CC
!> \verbatim
!>          CC is REAL array, dimension (LDC, max(NN))
!>          Used to hold a copy of the matrix C.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (LWORK)
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The number of entries in WORK.  This must be at least
!>          max( LDA+1, max(NN)+1 )*max(NN).
!> \endverbatim
!>
!> \param[out] RESULT
!> \verbatim
!>          RESULT is REAL array, dimension (4)
!>          The values computed by the tests described above.
!>          The values are currently limited to 1/ulp, to avoid
!>          overflow.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          If 0, then everything ran OK.
!>
!>-----------------------------------------------------------------------
!>
!>       Some Local Variables and Parameters:
!>       ---- ----- --------- --- ----------
!>       ZERO, ONE       Real 0 and 1.
!>       MAXTYP          The number of types defined.
!>       NTEST           The number of tests performed, or which can
!>                       be performed so far, for the current matrix.
!>       NTESTT          The total number of tests performed so far.
!>       NMAX            Largest value in NN.
!>       NMATS           The number of matrices generated so far.
!>       NERRS           The number of tests which have exceeded THRESH
!>                       so far.
!>       COND, IMODE     Values to be passed to the matrix generators.
!>       ANORM           Norm of A; passed to matrix generators.
!>
!>       OVFL, UNFL      Overflow and underflow thresholds.
!>       ULP, ULPINV     Finest relative precision and its inverse.
!>       RTOVFL, RTUNFL  Square roots of the previous 2 values.
!>               The following four arrays decode JTYPE:
!>       KTYPE(j)        The general type (1-10) for type "j".
!>       KMODE(j)        The MODE value to be passed to the matrix
!>                       generator for type "j".
!>       KMAGN(j)        The order of magnitude ( O(1),
!>                       O(overflow^(1/2) ), O(underflow^(1/2) )
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
!> \ingroup single_eig
!
!  =====================================================================
      SUBROUTINE schkbb(Nsizes,Mval,Nval,Nwdths,Kk,Ntypes,Dotype,Nrhs,  &
     &                  Iseed,Thresh,Nounit,A,Lda,Ab,Ldab,Bd,Be,Q,Ldq,P,&
     &                  Ldp,C,Ldc,Cc,Work,Lwork,Result,Info)
      IMPLICIT NONE
!*--SCHKBB358
!
!  -- LAPACK test routine (input) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Info , Lda , Ldab , Ldc , Ldp , Ldq , Lwork , Nounit ,    &
     &        Nrhs , Nsizes , Ntypes , Nwdths
      REAL Thresh
!     ..
!     .. Array Arguments ..
      LOGICAL Dotype(*)
      INTEGER Iseed(4) , Kk(*) , Mval(*) , Nval(*)
      REAL A(Lda,*) , Ab(Ldab,*) , Bd(*) , Be(*) , C(Ldc,*) , Cc(Ldc,*) &
     &     , P(Ldp,*) , Q(Ldq,*) , Result(*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE
      PARAMETER (ZERO=0.0E0,ONE=1.0E0)
      INTEGER MAXTYP
      PARAMETER (MAXTYP=15)
!     ..
!     .. Local Scalars ..
      LOGICAL badmm , badnn , badnnb
      INTEGER i , iinfo , imode , itype , j , jcol , jr , jsize ,       &
     &        jtype , jwidth , k , kl , kmax , ku , m , mmax , mnmax ,  &
     &        mnmin , mtypes , n , nerrs , nmats , nmax , ntest , ntestt
      REAL amninv , anorm , cond , ovfl , rtovfl , rtunfl , ulp ,       &
     &     ulpinv , unfl
!     ..
!     .. Local Arrays ..
      INTEGER idumma(1) , ioldsd(4) , kmagn(MAXTYP) , kmode(MAXTYP) ,   &
     &        ktype(MAXTYP)
!     ..
!     .. External Functions ..
      REAL SLAMCH
      EXTERNAL SLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL SGBBRD , SLACPY , SLASET
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX , MIN , REAL , SQRT
!     ..
!     .. Data statements ..
      DATA ktype/1 , 2 , 5*4 , 5*6 , 3*9/
      DATA kmagn/2*1 , 3*1 , 2 , 3 , 3*1 , 2 , 3 , 1 , 2 , 3/
      DATA kmode/2*0 , 4 , 3 , 1 , 4 , 4 , 4 , 3 , 1 , 4 , 4 , 0 , 0 ,  &
     &     0/
!     ..
!     .. Executable Statements ..
!
!     Check for errors
!
      ntestt = 0
      Info = 0
!
!     Important constants
!
      badmm = .FALSE.
      badnn = .FALSE.
      mmax = 1
      nmax = 1
      mnmax = 1
      DO j = 1 , Nsizes
         mmax = MAX(mmax,Mval(j))
         IF ( Mval(j)<0 ) badmm = .TRUE.
         nmax = MAX(nmax,Nval(j))
         IF ( Nval(j)<0 ) badnn = .TRUE.
         mnmax = MAX(mnmax,MIN(Mval(j),Nval(j)))
      ENDDO
!
      badnnb = .FALSE.
      kmax = 0
      DO j = 1 , Nwdths
         kmax = MAX(kmax,Kk(j))
         IF ( Kk(j)<0 ) badnnb = .TRUE.
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
      ELSEIF ( Nwdths<0 ) THEN
         Info = -4
      ELSEIF ( badnnb ) THEN
         Info = -5
      ELSEIF ( Ntypes<0 ) THEN
         Info = -6
      ELSEIF ( Nrhs<0 ) THEN
         Info = -8
      ELSEIF ( Lda<nmax ) THEN
         Info = -13
      ELSEIF ( Ldab<2*kmax+1 ) THEN
         Info = -15
      ELSEIF ( Ldq<nmax ) THEN
         Info = -19
      ELSEIF ( Ldp<nmax ) THEN
         Info = -21
      ELSEIF ( Ldc<nmax ) THEN
         Info = -23
      ELSEIF ( (MAX(Lda,nmax)+1)*nmax>Lwork ) THEN
         Info = -26
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL xerbla('schkbb',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( Nsizes==0 .OR. Ntypes==0 .OR. Nwdths==0 ) RETURN
!
!     More Important constants
!
      unfl = SLAMCH('Safe minimum')
      ovfl = ONE/unfl
      ulp = SLAMCH('Epsilon')*SLAMCH('Base')
      ulpinv = ONE/ulp
      rtunfl = SQRT(unfl)
      rtovfl = SQRT(ovfl)
!
!     Loop over sizes, widths, types
!
      nerrs = 0
      nmats = 0
!
      DO jsize = 1 , Nsizes
         m = Mval(jsize)
         n = Nval(jsize)
         mnmin = MIN(m,n)
         amninv = ONE/REAL(MAX(1,m,n))
!
         DO jwidth = 1 , Nwdths
            k = Kk(jwidth)
            IF ( k<m .OR. k<n ) THEN
               kl = MAX(0,MIN(m-1,k))
               ku = MAX(0,MIN(n-1,k))
!
               IF ( Nsizes/=1 ) THEN
                  mtypes = MIN(MAXTYP,Ntypes)
               ELSE
                  mtypes = MIN(MAXTYP+1,Ntypes)
               ENDIF
!
               DO jtype = 1 , mtypes
                  IF ( Dotype(jtype) ) THEN
                     nmats = nmats + 1
                     ntest = 0
!
                     DO j = 1 , 4
                        ioldsd(j) = Iseed(j)
                     ENDDO
!
!              Compute "A".
!
!              Control parameters:
!
!                  KMAGN  KMODE        KTYPE
!              =1  O(1)   clustered 1  zero
!              =2  large  clustered 2  identity
!              =3  small  exponential  (none)
!              =4         arithmetic   diagonal, (w/ singular values)
!              =5         random log   (none)
!              =6         random       nonhermitian, w/ singular values
!              =7                      (none)
!              =8                      (none)
!              =9                      random nonhermitian
!
                     IF ( mtypes<=MAXTYP ) THEN
!
                        itype = ktype(jtype)
                        imode = kmode(jtype)
!
!              Compute norm
!
                        IF ( kmagn(jtype)==2 ) THEN
!
                           anorm = (rtovfl*ulp)*amninv
                        ELSEIF ( kmagn(jtype)==3 ) THEN
!
                           anorm = rtunfl*MAX(m,n)*ulpinv
                        ELSE
!
                           anorm = ONE
                        ENDIF
!
!
                        CALL SLASET('Full',Lda,n,ZERO,ZERO,A,Lda)
                        CALL SLASET('Full',Ldab,n,ZERO,ZERO,Ab,Ldab)
                        iinfo = 0
                        cond = ulpinv
!
!              Special Matrices -- Identity & Jordan block
!
!                 Zero
!
                        IF ( itype==1 ) THEN
                           iinfo = 0
!
                        ELSEIF ( itype==2 ) THEN
!
!                 Identity
!
                           DO jcol = 1 , n
                              A(jcol,jcol) = anorm
                           ENDDO
!
                        ELSEIF ( itype==4 ) THEN
!
!                 Diagonal Matrix, singular values specified
!
                           CALL slatms(m,n,'S',Iseed,'N',Work,imode,    &
     &                                 cond,anorm,0,0,'N',A,Lda,        &
     &                                 Work(m+1),iinfo)
!
                        ELSEIF ( itype==6 ) THEN
!
!                 Nonhermitian, singular values specified
!
                           CALL slatms(m,n,'S',Iseed,'N',Work,imode,    &
     &                                 cond,anorm,kl,ku,'N',A,Lda,      &
     &                                 Work(m+1),iinfo)
!
                        ELSEIF ( itype==9 ) THEN
!
!                 Nonhermitian, random entries
!
                           CALL slatmr(m,n,'S',Iseed,'N',Work,6,ONE,ONE,&
     &                                 'T','N',Work(n+1),1,ONE,         &
     &                                 Work(2*n+1),1,ONE,'N',idumma,kl, &
     &                                 ku,ZERO,anorm,'N',A,Lda,idumma,  &
     &                                 iinfo)
!
                        ELSE
!
                           iinfo = 1
                        ENDIF
!
!              Generate Right-Hand Side
!
                        CALL slatmr(m,Nrhs,'S',Iseed,'N',Work,6,ONE,ONE,&
     &                              'T','N',Work(m+1),1,ONE,Work(2*m+1),&
     &                              1,ONE,'N',idumma,m,Nrhs,ZERO,ONE,   &
     &                              'NO',C,Ldc,idumma,iinfo)
!
                        IF ( iinfo/=0 ) THEN
                           WRITE (Nounit,FMT=99001) 'Generator' ,       &
     &                            iinfo , n , jtype , ioldsd
                           Info = ABS(iinfo)
                           RETURN
                        ENDIF
                     ENDIF
!
!
!              Copy A to band storage.
!
                     DO j = 1 , n
                        DO i = MAX(1,j-ku) , MIN(m,j+kl)
                           Ab(ku+1+i-j,j) = A(i,j)
                        ENDDO
                     ENDDO
!
!              Copy C
!
                     CALL SLACPY('Full',m,Nrhs,C,Ldc,Cc,Ldc)
!
!              Call SGBBRD to compute B, Q and P, and to update C.
!
                     CALL SGBBRD('B',m,n,Nrhs,kl,ku,Ab,Ldab,Bd,Be,Q,Ldq,&
     &                           P,Ldp,Cc,Ldc,Work,iinfo)
!
                     IF ( iinfo/=0 ) THEN
                        WRITE (Nounit,FMT=99001) 'SGBBRD' , iinfo , n , &
     &                         jtype , ioldsd
                        Info = ABS(iinfo)
                        IF ( iinfo<0 ) THEN
                           RETURN
                        ELSE
                           Result(1) = ulpinv
                           GOTO 2
                        ENDIF
                     ENDIF
!
!              Test 1:  Check the decomposition A := Q * B * P'
!                   2:  Check the orthogonality of Q
!                   3:  Check the orthogonality of P
!                   4:  Check the computation of Q' * C
!
                     CALL sbdt01(m,n,-1,A,Lda,Q,Ldq,Bd,Be,P,Ldp,Work,   &
     &                           Result(1))
                     CALL sort01('Columns',m,m,Q,Ldq,Work,Lwork,        &
     &                           Result(2))
                     CALL sort01('Rows',n,n,P,Ldp,Work,Lwork,Result(3))
                     CALL sbdt02(m,Nrhs,C,Ldc,Cc,Ldc,Q,Ldq,Work,        &
     &                           Result(4))
!
!              End of Loop -- Check for RESULT(j) > THRESH
!
                     ntest = 4
 2                   ntestt = ntestt + ntest
!
!              Print out tests which fail.
!
                     DO jr = 1 , ntest
                        IF ( Result(jr)>=Thresh ) THEN
                           IF ( nerrs==0 ) CALL slahd2(Nounit,'SBB')
                           nerrs = nerrs + 1
                           WRITE (Nounit,FMT=99002) m , n , k , ioldsd ,&
     &                            jtype , jr , Result(jr)
                        ENDIF
                     ENDDO
                  ENDIF
!
               ENDDO
            ENDIF
         ENDDO
      ENDDO
!
!     Summary
!
      CALL slasum('SBB',Nounit,nerrs,ntestt)
      RETURN
!
99001 FORMAT (' schkbb: ',A,' returned INFO=',I5,'.',/9X,'M=',I5,' N=', &
     &        I5,' K=',I5,', JTYPE=',I5,', ISEED=(',3(I5,','),I5,')')
99002 FORMAT (' M =',I4,' N=',I4,', K=',I3,', seed=',4(I4,','),' type ',&
     &        I2,', test(',I2,')=',G10.3)
!
!     End of schkbb
!
      END SUBROUTINE schkbb
