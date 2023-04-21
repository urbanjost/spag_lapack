!*==dsytrd_sb2st.f90  processed by SPAG 7.51RB at 20:06 on  3 Mar 2022
!> \brief \b DSYTRD_SB2ST reduces a real symmetric band matrix A to real symmetric tridiagonal form T
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DSYTRD_SB2ST + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsytrd_sb2st.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsytrd_sb2st.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsytrd_sb2st.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DSYTRD_SB2ST( STAGE1, VECT, UPLO, N, KD, AB, LDAB,
!                               D, E, HOUS, LHOUS, WORK, LWORK, INFO )
!
!       #if defined(_OPENMP)
!       use omp_lib
!       #endif
!
!       IMPLICIT NONE
!
!       .. Scalar Arguments ..
!       CHARACTER          STAGE1, UPLO, VECT
!       INTEGER            N, KD, IB, LDAB, LHOUS, LWORK, INFO
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   D( * ), E( * )
!       DOUBLE PRECISION   AB( LDAB, * ), HOUS( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DSYTRD_SB2ST reduces a real symmetric band matrix A to real symmetric
!> tridiagonal form T by a orthogonal similarity transformation:
!> Q**T * A * Q = T.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] STAGE1
!> \verbatim
!>          STAGE1 is CHARACTER*1
!>          = 'N':  "No": to mention that the stage 1 of the reduction
!>                  from dense to band using the dsytrd_sy2sb routine
!>                  was not called before this routine to reproduce AB.
!>                  In other term this routine is called as standalone.
!>          = 'Y':  "Yes": to mention that the stage 1 of the
!>                  reduction from dense to band using the dsytrd_sy2sb
!>                  routine has been called to produce AB (e.g., AB is
!>                  the output of dsytrd_sy2sb.
!> \endverbatim
!>
!> \param[in] VECT
!> \verbatim
!>          VECT is CHARACTER*1
!>          = 'N':  No need for the Housholder representation,
!>                  and thus LHOUS is of size max(1, 4*N);
!>          = 'V':  the Householder representation is needed to
!>                  either generate or to apply Q later on,
!>                  then LHOUS is to be queried and computed.
!>                  (NOT AVAILABLE IN THIS RELEASE).
!> \endverbatim
!>
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          = 'U':  Upper triangle of A is stored;
!>          = 'L':  Lower triangle of A is stored.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] KD
!> \verbatim
!>          KD is INTEGER
!>          The number of superdiagonals of the matrix A if UPLO = 'U',
!>          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.
!> \endverbatim
!>
!> \param[in,out] AB
!> \verbatim
!>          AB is DOUBLE PRECISION array, dimension (LDAB,N)
!>          On entry, the upper or lower triangle of the symmetric band
!>          matrix A, stored in the first KD+1 rows of the array.  The
!>          j-th column of A is stored in the j-th column of the array AB
!>          as follows:
!>          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
!>          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).
!>          On exit, the diagonal elements of AB are overwritten by the
!>          diagonal elements of the tridiagonal matrix T; if KD > 0, the
!>          elements on the first superdiagonal (if UPLO = 'U') or the
!>          first subdiagonal (if UPLO = 'L') are overwritten by the
!>          off-diagonal elements of T; the rest of AB is overwritten by
!>          values generated during the reduction.
!> \endverbatim
!>
!> \param[in] LDAB
!> \verbatim
!>          LDAB is INTEGER
!>          The leading dimension of the array AB.  LDAB >= KD+1.
!> \endverbatim
!>
!> \param[out] D
!> \verbatim
!>          D is DOUBLE PRECISION array, dimension (N)
!>          The diagonal elements of the tridiagonal matrix T.
!> \endverbatim
!>
!> \param[out] E
!> \verbatim
!>          E is DOUBLE PRECISION array, dimension (N-1)
!>          The off-diagonal elements of the tridiagonal matrix T:
!>          E(i) = T(i,i+1) if UPLO = 'U'; E(i) = T(i+1,i) if UPLO = 'L'.
!> \endverbatim
!>
!> \param[out] HOUS
!> \verbatim
!>          HOUS is DOUBLE PRECISION array, dimension LHOUS, that
!>          store the Householder representation.
!> \endverbatim
!>
!> \param[in] LHOUS
!> \verbatim
!>          LHOUS is INTEGER
!>          The dimension of the array HOUS. LHOUS = MAX(1, dimension)
!>          If LWORK = -1, or LHOUS=-1,
!>          then a query is assumed; the routine
!>          only calculates the optimal size of the HOUS array, returns
!>          this value as the first entry of the HOUS array, and no error
!>          message related to LHOUS is issued by XERBLA.
!>          LHOUS = MAX(1, dimension) where
!>          dimension = 4*N if VECT='N'
!>          not available now if VECT='H'
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK. LWORK = MAX(1, dimension)
!>          If LWORK = -1, or LHOUS=-1,
!>          then a workspace query is assumed; the routine
!>          only calculates the optimal size of the WORK array, returns
!>          this value as the first entry of the WORK array, and no error
!>          message related to LWORK is issued by XERBLA.
!>          LWORK = MAX(1, dimension) where
!>          dimension   = (2KD+1)*N + KD*NTHREADS
!>          where KD is the blocking size of the reduction,
!>          FACTOPTNB is the blocking used by the QR or LQ
!>          algorithm, usually FACTOPTNB=128 is a good choice
!>          NTHREADS is the number of threads used when
!>          openMP compilation is enabled, otherwise =1.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
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
!> \date November 2017
!
!> \ingroup real16OTHERcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  Implemented by Azzam Haidar.
!>
!>  All details are available on technical report, SC11, SC13 papers.
!>
!>  Azzam Haidar, Hatem Ltaief, and Jack Dongarra.
!>  Parallel reduction to condensed forms for symmetric eigenvalue problems
!>  using aggregated fine-grained and memory-aware kernels. In Proceedings
!>  of 2011 International Conference for High Performance Computing,
!>  Networking, Storage and Analysis (SC '11), New York, NY, USA,
!>  Article 8 , 11 pages.
!>  http://doi.acm.org/10.1145/2063384.2063394
!>
!>  A. Haidar, J. Kurzak, P. Luszczek, 2013.
!>  An improved parallel singular value algorithm and its implementation
!>  for multicore hardware, In Proceedings of 2013 International Conference
!>  for High Performance Computing, Networking, Storage and Analysis (SC '13).
!>  Denver, Colorado, USA, 2013.
!>  Article 90, 12 pages.
!>  http://doi.acm.org/10.1145/2503210.2503292
!>
!>  A. Haidar, R. Solca, S. Tomov, T. Schulthess and J. Dongarra.
!>  A novel hybrid CPU-GPU generalized eigensolver for electronic structure
!>  calculations based on fine-grained memory aware tasks.
!>  International Journal of High Performance Computing Applications.
!>  Volume 28 Issue 2, Pages 196-209, May 2014.
!>  http://hpc.sagepub.com/content/28/2/196
!>
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE DSYTRD_SB2ST(Stage1,Vect,Uplo,N,Kd,Ab,Ldab,D,E,Hous,   &
     &                        Lhous,Work,Lwork,Info)
!
#if defined(_OPENMP)
      USE OMP_LIB
#endif
!
      IMPLICIT NONE
!*--DSYTRD_SB2ST239
!
!  -- LAPACK computational routine (version 3.8.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2017
!
!     .. Scalar Arguments ..
      CHARACTER Stage1 , Uplo , Vect
      INTEGER N , Kd , Ldab , Lhous , Lwork , Info
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION D(*) , E(*)
      DOUBLE PRECISION Ab(Ldab,*) , Hous(*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION RZERO
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (RZERO=0.0D+0,ZERO=0.0D+0,ONE=1.0D+0)
!     ..
!     .. Local Scalars ..
      LOGICAL lquery , wantq , upper , afters1
      INTEGER i , m , k , ib , sweepid , myid , shift , stt , st , ed , &
     &        stind , edind , blklastind , colpt , thed , stepercol ,   &
     &        grsiz , thgrsiz , thgrnb , thgrid , nbtiles , ttype ,     &
     &        tid , nthreads , debug , abdpos , abofdpos , dpos ,       &
     &        ofdpos , awpos , inda , indw , apos , sizea , lda , indv ,&
     &        indtau , sidev , sizetau , ldv , lhmin , lwmin
!     ..
!     .. External Subroutines ..
      EXTERNAL DSB2ST_KERNELS , DLACPY , DLASET , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MIN , MAX , CEILING , REAL
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      INTEGER ILAENV2STAGE
      EXTERNAL LSAME , ILAENV2STAGE
!     ..
!     .. Executable Statements ..
!
!     Determine the minimal workspace size required.
!     Test the input parameters
!
      debug = 0
      Info = 0
      afters1 = LSAME(Stage1,'Y')
      wantq = LSAME(Vect,'V')
      upper = LSAME(Uplo,'U')
      lquery = (Lwork==-1) .OR. (Lhous==-1)
!
!     Determine the block size, the workspace size and the hous size.
!
      ib = ILAENV2STAGE(2,'DSYTRD_SB2ST',Vect,N,Kd,-1,-1)
      lhmin = ILAENV2STAGE(3,'DSYTRD_SB2ST',Vect,N,Kd,ib,-1)
      lwmin = ILAENV2STAGE(4,'DSYTRD_SB2ST',Vect,N,Kd,ib,-1)
!
      IF ( .NOT.afters1 .AND. .NOT.LSAME(Stage1,'N') ) THEN
         Info = -1
      ELSEIF ( .NOT.LSAME(Vect,'N') ) THEN
         Info = -2
      ELSEIF ( .NOT.upper .AND. .NOT.LSAME(Uplo,'L') ) THEN
         Info = -3
      ELSEIF ( N<0 ) THEN
         Info = -4
      ELSEIF ( Kd<0 ) THEN
         Info = -5
      ELSEIF ( Ldab<(Kd+1) ) THEN
         Info = -7
      ELSEIF ( Lhous<lhmin .AND. .NOT.lquery ) THEN
         Info = -11
      ELSEIF ( Lwork<lwmin .AND. .NOT.lquery ) THEN
         Info = -13
      ENDIF
!
      IF ( Info==0 ) THEN
         Hous(1) = lhmin
         Work(1) = lwmin
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('DSYTRD_SB2ST',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 ) THEN
         Hous(1) = 1
         Work(1) = 1
         RETURN
      ENDIF
!
!     Determine pointer position
!
      ldv = Kd + ib
      sizetau = 2*N
      sidev = 2*N
      indtau = 1
      indv = indtau + sizetau
      lda = 2*Kd + 1
      sizea = lda*N
      inda = 1
      indw = inda + sizea
      nthreads = 1
      tid = 0
!
      IF ( upper ) THEN
         apos = inda + Kd
         awpos = inda
         dpos = apos + Kd
         ofdpos = dpos - 1
         abdpos = Kd + 1
         abofdpos = Kd
      ELSE
         apos = inda
         awpos = inda + Kd + 1
         dpos = apos
         ofdpos = dpos + 1
         abdpos = 1
         abofdpos = 2
 
      ENDIF
!
!     Case KD=0:
!     The matrix is diagonal. We just copy it (convert to "real" for
!     real because D is double and the imaginary part should be 0)
!     and store it in D. A sequential code here is better or
!     in a parallel environment it might need two cores for D and E
!
      IF ( Kd==0 ) THEN
         DO i = 1 , N
            D(i) = (Ab(abdpos,i))
         ENDDO
         DO i = 1 , N - 1
            E(i) = RZERO
         ENDDO
!
         Hous(1) = 1
         Work(1) = 1
         RETURN
      ENDIF
!
!     Case KD=1:
!     The matrix is already Tridiagonal. We have to make diagonal
!     and offdiagonal elements real, and store them in D and E.
!     For that, for real precision just copy the diag and offdiag
!     to D and E while for the COMPLEX case the bulge chasing is
!     performed to convert the hermetian tridiagonal to symmetric
!     tridiagonal. A simpler conversion formula might be used, but then
!     updating the Q matrix will be required and based if Q is generated
!     or not this might complicate the story.
!
      IF ( Kd==1 ) THEN
         DO i = 1 , N
            D(i) = (Ab(abdpos,i))
         ENDDO
!
         IF ( upper ) THEN
            DO i = 1 , N - 1
               E(i) = (Ab(abofdpos,i+1))
            ENDDO
         ELSE
            DO i = 1 , N - 1
               E(i) = (Ab(abofdpos,i))
            ENDDO
         ENDIF
!
         Hous(1) = 1
         Work(1) = 1
         RETURN
      ENDIF
!
!     Main code start here.
!     Reduce the symmetric band of A to a tridiagonal matrix.
!
      thgrsiz = N
      grsiz = 1
      shift = 3
      nbtiles = CEILING(REAL(N)/REAL(Kd))
      stepercol = CEILING(REAL(shift)/REAL(grsiz))
      thgrnb = CEILING(REAL(N-1)/REAL(thgrsiz))
!
      CALL DLACPY("A",Kd+1,N,Ab,Ldab,Work(apos),lda)
      CALL DLASET("A",Kd,N,ZERO,ZERO,Work(awpos),lda)
!
!
!     openMP parallelisation start here
!
#if defined(_OPENMP)
!$OMP PARALLEL PRIVATE( TID, THGRID, BLKLASTIND )
!$OMP$         PRIVATE( THED, I, M, K, ST, ED, STT, SWEEPID )
!$OMP$         PRIVATE( MYID, TTYPE, COLPT, STIND, EDIND )
!$OMP$         SHARED ( UPLO, WANTQ, INDV, INDTAU, HOUS, WORK)
!$OMP$         SHARED ( N, KD, IB, NBTILES, LDA, LDV, INDA )
!$OMP$         SHARED ( STEPERCOL, THGRNB, THGRSIZ, GRSIZ, SHIFT )
!$OMP MASTER
#endif
!
!     main bulge chasing loop
!
      DO thgrid = 1 , thgrnb
         stt = (thgrid-1)*thgrsiz + 1
         thed = MIN((stt+thgrsiz-1),(N-1))
         DO i = stt , N - 1
            ed = MIN(i,thed)
            IF ( stt>ed ) EXIT
            DO m = 1 , stepercol
               st = stt
               DO sweepid = st , ed
                  DO k = 1 , grsiz
                     myid = (i-sweepid)*(stepercol*grsiz) + (m-1)       &
     &                      *grsiz + k
                     IF ( myid==1 ) THEN
                        ttype = 1
                     ELSE
                        ttype = MOD(myid,2) + 2
                     ENDIF
 
                     IF ( ttype==2 ) THEN
                        colpt = (myid/2)*Kd + sweepid
                        stind = colpt - Kd + 1
                        edind = MIN(colpt,N)
                        blklastind = colpt
                     ELSE
                        colpt = ((myid+1)/2)*Kd + sweepid
                        stind = colpt - Kd + 1
                        edind = MIN(colpt,N)
                        IF ( (stind>=edind-1) .AND. (edind==N) ) THEN
                           blklastind = N
                        ELSE
                           blklastind = 0
                        ENDIF
                     ENDIF
!
!                         Call the kernel
!
#if defined(_OPENMP)
                     IF ( ttype/=1 ) THEN
!$OMP TASK DEPEND(in:WORK(MYID+SHIFT-1))
!$OMP$     DEPEND(in:WORK(MYID-1))
!$OMP$     DEPEND(out:WORK(MYID))
                        tid = OMP_GET_THREAD_NUM()
                        CALL DSB2ST_KERNELS(Uplo,wantq,ttype,stind,     &
     &                     edind,sweepid,N,Kd,ib,Work(inda),lda,        &
     &                     Hous(indv),Hous(indtau),ldv,Work(indw+tid*Kd)&
     &                     )
!$OMP END TASK
                     ELSE
!$OMP TASK DEPEND(in:WORK(MYID+SHIFT-1))
!$OMP$     DEPEND(out:WORK(MYID))
                        tid = OMP_GET_THREAD_NUM()
                        CALL DSB2ST_KERNELS(Uplo,wantq,ttype,stind,     &
     &                     edind,sweepid,N,Kd,ib,Work(inda),lda,        &
     &                     Hous(indv),Hous(indtau),ldv,Work(indw+tid*Kd)&
     &                     )
!$OMP END TASK
                     ENDIF
#else
                     CALL DSB2ST_KERNELS(Uplo,wantq,ttype,stind,edind,  &
     &                  sweepid,N,Kd,ib,Work(inda),lda,Hous(indv),      &
     &                  Hous(indtau),ldv,Work(indw+tid*Kd))
#endif
                     IF ( blklastind>=(N-1) ) THEN
                        stt = stt + 1
                        EXIT
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
!
#if defined(_OPENMP)
!$OMP END MASTER
!$OMP END PARALLEL
#endif
!
!     Copy the diagonal from A to D. Note that D is REAL thus only
!     the Real part is needed, the imaginary part should be zero.
!
      DO i = 1 , N
         D(i) = (Work(dpos+(i-1)*lda))
      ENDDO
!
!     Copy the off diagonal from A to E. Note that E is REAL thus only
!     the Real part is needed, the imaginary part should be zero.
!
      IF ( upper ) THEN
         DO i = 1 , N - 1
            E(i) = (Work(ofdpos+i*lda))
         ENDDO
      ELSE
         DO i = 1 , N - 1
            E(i) = (Work(ofdpos+(i-1)*lda))
         ENDDO
      ENDIF
!
      Hous(1) = lhmin
      Work(1) = lwmin
!
!     End of DSYTRD_SB2ST
!
      END SUBROUTINE DSYTRD_SB2ST
