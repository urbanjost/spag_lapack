!*==dsb2st_kernels.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DSB2ST_KERNELS
!
!  @generated from zhb2st_kernels.f, fortran z -> d, Wed Dec  7 08:22:39 2016
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DSB2ST_KERNELS + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsb2st_kernels.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsb2st_kernels.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsb2st_kernels.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE  DSB2ST_KERNELS( UPLO, WANTZ, TTYPE,
!                                   ST, ED, SWEEP, N, NB, IB,
!                                   A, LDA, V, TAU, LDVT, WORK)
!
!       IMPLICIT NONE
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       LOGICAL            WANTZ
!       INTEGER            TTYPE, ST, ED, SWEEP, N, NB, IB, LDA, LDVT
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * ), V( * ),
!                          TAU( * ), WORK( * )
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DSB2ST_KERNELS is an internal routine used by the DSYTRD_SB2ST
!> subroutine.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!> \endverbatim
!>
!> \param[in] WANTZ
!> \verbatim
!>          WANTZ is LOGICAL which indicate if Eigenvalue are requested or both
!>          Eigenvalue/Eigenvectors.
!> \endverbatim
!>
!> \param[in] TTYPE
!> \verbatim
!>          TTYPE is INTEGER
!> \endverbatim
!>
!> \param[in] ST
!> \verbatim
!>          ST is INTEGER
!>          internal parameter for indices.
!> \endverbatim
!>
!> \param[in] ED
!> \verbatim
!>          ED is INTEGER
!>          internal parameter for indices.
!> \endverbatim
!>
!> \param[in] SWEEP
!> \verbatim
!>          SWEEP is INTEGER
!>          internal parameter for indices.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER. The order of the matrix A.
!> \endverbatim
!>
!> \param[in] NB
!> \verbatim
!>          NB is INTEGER. The size of the band.
!> \endverbatim
!>
!> \param[in] IB
!> \verbatim
!>          IB is INTEGER.
!> \endverbatim
!>
!> \param[in, out] A
!> \verbatim
!>          A is DOUBLE PRECISION array. A pointer to the matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER. The leading dimension of the matrix A.
!> \endverbatim
!>
!> \param[out] V
!> \verbatim
!>          V is DOUBLE PRECISION array, dimension 2*n if eigenvalues only are
!>          requested or to be queried for vectors.
!> \endverbatim
!>
!> \param[out] TAU
!> \verbatim
!>          TAU is DOUBLE PRECISION array, dimension (2*n).
!>          The scalar factors of the Householder reflectors are stored
!>          in this array.
!> \endverbatim
!>
!> \param[in] LDVT
!> \verbatim
!>          LDVT is INTEGER.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array. Workspace of size nb.
!> \endverbatim
!>
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
      SUBROUTINE DSB2ST_KERNELS(Uplo,Wantz,Ttype,St,Ed,Sweep,N,Nb,Ib,A, &
     &                          Lda,V,Tau,Ldvt,Work)
!
      USE F77KINDS                        
      USE S_DLARFG
      USE S_DLARFX
      USE S_DLARFY
      USE S_LSAME
      IMPLICIT NONE
!*--DSB2ST_KERNELS177
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Uplo
      LOGICAL , INTENT(IN) :: Wantz
      INTEGER , INTENT(IN) :: Ttype
      INTEGER , INTENT(IN) :: St
      INTEGER , INTENT(IN) :: Ed
      INTEGER , INTENT(IN) :: Sweep
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Nb
      INTEGER , INTENT(IN) :: Ib
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL(R8KIND) , DIMENSION(*) :: V
      REAL(R8KIND) , DIMENSION(*) :: Tau
      INTEGER , INTENT(IN) :: Ldvt
      REAL(R8KIND) , DIMENSION(*) :: Work
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: ajeter , dpos , i , j1 , j2 , lm , ln , ofdpos ,       &
     &           taupos , vpos
      REAL(R8KIND) :: ctmp
      LOGICAL :: upper
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
!     .. External Functions ..
!     ..
!     ..
!     .. Executable Statements ..
!
      ajeter = Ib + Ldvt
      upper = LSAME(Uplo,'U')
 
      IF ( upper ) THEN
         dpos = 2*Nb + 1
         ofdpos = 2*Nb
      ELSE
         dpos = 1
         ofdpos = 2
      ENDIF
 
!
!     Upper case
!
      IF ( upper ) THEN
!
         IF ( Wantz ) THEN
            vpos = MOD(Sweep-1,2)*N + St
            taupos = MOD(Sweep-1,2)*N + St
         ELSE
            vpos = MOD(Sweep-1,2)*N + St
            taupos = MOD(Sweep-1,2)*N + St
         ENDIF
!
         IF ( Ttype==1 ) THEN
            lm = Ed - St + 1
!
            V(vpos) = ONE
            DO i = 1 , lm - 1
               V(vpos+i) = (A(ofdpos-i,St+i))
               A(ofdpos-i,St+i) = ZERO
            ENDDO
            ctmp = (A(ofdpos,St))
            CALL DLARFG(lm,ctmp,V(vpos+1),1,Tau(taupos))
            A(ofdpos,St) = ctmp
!
            lm = Ed - St + 1
            CALL DLARFY(Uplo,lm,V(vpos),1,(Tau(taupos)),A(dpos,St),     &
     &                  Lda-1,Work)
         ENDIF
!
         IF ( Ttype==3 ) THEN
!
            lm = Ed - St + 1
            CALL DLARFY(Uplo,lm,V(vpos),1,(Tau(taupos)),A(dpos,St),     &
     &                  Lda-1,Work)
         ENDIF
!
         IF ( Ttype==2 ) THEN
            j1 = Ed + 1
            j2 = MIN(Ed+Nb,N)
            ln = Ed - St + 1
            lm = j2 - j1 + 1
            IF ( lm>0 ) THEN
               CALL DLARFX('Left',ln,lm,V(vpos),(Tau(taupos)),          &
     &                     A(dpos-Nb,j1),Lda-1,Work)
!
               IF ( Wantz ) THEN
                  vpos = MOD(Sweep-1,2)*N + j1
                  taupos = MOD(Sweep-1,2)*N + j1
               ELSE
                  vpos = MOD(Sweep-1,2)*N + j1
                  taupos = MOD(Sweep-1,2)*N + j1
               ENDIF
!
               V(vpos) = ONE
               DO i = 1 , lm - 1
                  V(vpos+i) = (A(dpos-Nb-i,j1+i))
                  A(dpos-Nb-i,j1+i) = ZERO
               ENDDO
               ctmp = (A(dpos-Nb,j1))
               CALL DLARFG(lm,ctmp,V(vpos+1),1,Tau(taupos))
               A(dpos-Nb,j1) = ctmp
!
               CALL DLARFX('Right',ln-1,lm,V(vpos),Tau(taupos),         &
     &                     A(dpos-Nb+1,j1),Lda-1,Work)
            ENDIF
         ENDIF
!
!     Lower case
!
      ELSE
!
         IF ( Wantz ) THEN
            vpos = MOD(Sweep-1,2)*N + St
            taupos = MOD(Sweep-1,2)*N + St
         ELSE
            vpos = MOD(Sweep-1,2)*N + St
            taupos = MOD(Sweep-1,2)*N + St
         ENDIF
!
         IF ( Ttype==1 ) THEN
            lm = Ed - St + 1
!
            V(vpos) = ONE
            DO i = 1 , lm - 1
               V(vpos+i) = A(ofdpos+i,St-1)
               A(ofdpos+i,St-1) = ZERO
            ENDDO
            CALL DLARFG(lm,A(ofdpos,St-1),V(vpos+1),1,Tau(taupos))
!
            lm = Ed - St + 1
!
            CALL DLARFY(Uplo,lm,V(vpos),1,(Tau(taupos)),A(dpos,St),     &
     &                  Lda-1,Work)
 
         ENDIF
!
         IF ( Ttype==3 ) THEN
            lm = Ed - St + 1
!
            CALL DLARFY(Uplo,lm,V(vpos),1,(Tau(taupos)),A(dpos,St),     &
     &                  Lda-1,Work)
 
         ENDIF
!
         IF ( Ttype==2 ) THEN
            j1 = Ed + 1
            j2 = MIN(Ed+Nb,N)
            ln = Ed - St + 1
            lm = j2 - j1 + 1
!
            IF ( lm>0 ) THEN
               CALL DLARFX('Right',lm,ln,V(vpos),Tau(taupos),           &
     &                     A(dpos+Nb,St),Lda-1,Work)
!
               IF ( Wantz ) THEN
                  vpos = MOD(Sweep-1,2)*N + j1
                  taupos = MOD(Sweep-1,2)*N + j1
               ELSE
                  vpos = MOD(Sweep-1,2)*N + j1
                  taupos = MOD(Sweep-1,2)*N + j1
               ENDIF
!
               V(vpos) = ONE
               DO i = 1 , lm - 1
                  V(vpos+i) = A(dpos+Nb+i,St)
                  A(dpos+Nb+i,St) = ZERO
               ENDDO
               CALL DLARFG(lm,A(dpos+Nb,St),V(vpos+1),1,Tau(taupos))
!
               CALL DLARFX('Left',lm,ln-1,V(vpos),(Tau(taupos)),        &
     &                     A(dpos+Nb-1,St+1),Lda-1,Work)
 
            ENDIF
         ENDIF
      ENDIF
!
!
!     END OF DSB2ST_KERNELS
!
      END SUBROUTINE DSB2ST_KERNELS
