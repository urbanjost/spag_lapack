!*==zstein.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b ZSTEIN
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZSTEIN + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zstein.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zstein.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zstein.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZSTEIN( N, D, E, M, W, IBLOCK, ISPLIT, Z, LDZ, WORK,
!                          IWORK, IFAIL, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDZ, M, N
!       ..
!       .. Array Arguments ..
!       INTEGER            IBLOCK( * ), IFAIL( * ), ISPLIT( * ),
!      $                   IWORK( * )
!       DOUBLE PRECISION   D( * ), E( * ), W( * ), WORK( * )
!       COMPLEX*16         Z( LDZ, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZSTEIN computes the eigenvectors of a real symmetric tridiagonal
!> matrix T corresponding to specified eigenvalues, using inverse
!> iteration.
!>
!> The maximum number of iterations allowed for each eigenvector is
!> specified by an internal parameter MAXITS (currently set to 5).
!>
!> Although the eigenvectors are real, they are stored in a complex
!> array, which may be passed to ZUNMTR or ZUPMTR for back
!> transformation to the eigenvectors of a complex Hermitian matrix
!> which was reduced to tridiagonal form.
!>
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix.  N >= 0.
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is DOUBLE PRECISION array, dimension (N)
!>          The n diagonal elements of the tridiagonal matrix T.
!> \endverbatim
!>
!> \param[in] E
!> \verbatim
!>          E is DOUBLE PRECISION array, dimension (N-1)
!>          The (n-1) subdiagonal elements of the tridiagonal matrix
!>          T, stored in elements 1 to N-1.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of eigenvectors to be found.  0 <= M <= N.
!> \endverbatim
!>
!> \param[in] W
!> \verbatim
!>          W is DOUBLE PRECISION array, dimension (N)
!>          The first M elements of W contain the eigenvalues for
!>          which eigenvectors are to be computed.  The eigenvalues
!>          should be grouped by split-off block and ordered from
!>          smallest to largest within the block.  ( The output array
!>          W from DSTEBZ with ORDER = 'B' is expected here. )
!> \endverbatim
!>
!> \param[in] IBLOCK
!> \verbatim
!>          IBLOCK is INTEGER array, dimension (N)
!>          The submatrix indices associated with the corresponding
!>          eigenvalues in W; IBLOCK(i)=1 if eigenvalue W(i) belongs to
!>          the first submatrix from the top, =2 if W(i) belongs to
!>          the second submatrix, etc.  ( The output array IBLOCK
!>          from DSTEBZ is expected here. )
!> \endverbatim
!>
!> \param[in] ISPLIT
!> \verbatim
!>          ISPLIT is INTEGER array, dimension (N)
!>          The splitting points, at which T breaks up into submatrices.
!>          The first submatrix consists of rows/columns 1 to
!>          ISPLIT( 1 ), the second of rows/columns ISPLIT( 1 )+1
!>          through ISPLIT( 2 ), etc.
!>          ( The output array ISPLIT from DSTEBZ is expected here. )
!> \endverbatim
!>
!> \param[out] Z
!> \verbatim
!>          Z is COMPLEX*16 array, dimension (LDZ, M)
!>          The computed eigenvectors.  The eigenvector associated
!>          with the eigenvalue W(i) is stored in the i-th column of
!>          Z.  Any vector which fails to converge is set to its current
!>          iterate after MAXITS iterations.
!>          The imaginary parts of the eigenvectors are set to zero.
!> \endverbatim
!>
!> \param[in] LDZ
!> \verbatim
!>          LDZ is INTEGER
!>          The leading dimension of the array Z.  LDZ >= max(1,N).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (5*N)
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (N)
!> \endverbatim
!>
!> \param[out] IFAIL
!> \verbatim
!>          IFAIL is INTEGER array, dimension (M)
!>          On normal exit, all elements of IFAIL are zero.
!>          If one or more eigenvectors fail to converge after
!>          MAXITS iterations, then their indices are stored in
!>          array IFAIL.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: successful exit
!>          < 0: if INFO = -i, the i-th argument had an illegal value
!>          > 0: if INFO = i, then i eigenvectors failed to converge
!>               in MAXITS iterations.  Their indices are stored in
!>               array IFAIL.
!> \endverbatim
!
!> \par Internal Parameters:
!  =========================
!>
!> \verbatim
!>  MAXITS  INTEGER, default = 5
!>          The maximum number of iterations performed.
!>
!>  EXTRA   INTEGER, default = 2
!>          The number of iterations performed after norm growth
!>          criterion is satisfied, should be at least 1.
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
!> \ingroup complex16OTHERcomputational
!
!  =====================================================================
      SUBROUTINE ZSTEIN(N,D,E,M,W,Iblock,Isplit,Z,Ldz,Work,Iwork,Ifail, &
     &                  Info)
      USE F77KINDS                        
      USE S_DCOPY
      USE S_DLAGTF
      USE S_DLAGTS
      USE S_DLAMCH
      USE S_DLARNV
      USE S_DNRM2
      USE S_DSCAL
      USE S_IDAMAX
      USE S_XERBLA
      IMPLICIT NONE
!*--ZSTEIN196
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX(CX16KIND) , PARAMETER  ::  CZERO = (0.0D+0,0.0D+0) ,      &
     &                 CONE = (1.0D+0,0.0D+0)
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0 ,      &
     &                              TEN = 1.0D+1 , ODM3 = 1.0D-3 ,      &
     &                              ODM1 = 1.0D-1
      INTEGER , PARAMETER  ::  MAXITS = 5 , EXTRA = 2
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , DIMENSION(*) :: D
      REAL(R8KIND) , DIMENSION(*) :: E
      INTEGER , INTENT(IN) :: M
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: W
      INTEGER , INTENT(IN) , DIMENSION(*) :: Iblock
      INTEGER , INTENT(IN) , DIMENSION(*) :: Isplit
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldz,*) :: Z
      INTEGER , INTENT(IN) :: Ldz
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Ifail
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: b1 , blksiz , bn , gpind , i , iinfo , indrv1 ,        &
     &           indrv2 , indrv3 , indrv4 , indrv5 , its , j , j1 ,     &
     &           jblk , jmax , jr , nblk , nrmchk
      REAL(R8KIND) :: dtpcrt , eps , eps1 , nrm , onenrm , ortol ,      &
     &                pertol , scl , sep , tol , xj , xjm , ztr
      INTEGER , DIMENSION(4) :: iseed
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
! =====================================================================
!
!     .. Parameters ..
!     ..
!     .. Local Scalars ..
!     ..
!     .. Local Arrays ..
!     ..
!     .. External Functions ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      Info = 0
      DO i = 1 , M
         Ifail(i) = 0
      ENDDO
!
      IF ( N<0 ) THEN
         Info = -1
      ELSEIF ( M<0 .OR. M>N ) THEN
         Info = -4
      ELSEIF ( Ldz<MAX(1,N) ) THEN
         Info = -9
      ELSE
         DO j = 2 , M
            IF ( Iblock(j)<Iblock(j-1) ) THEN
               Info = -6
               EXIT
            ENDIF
            IF ( Iblock(j)==Iblock(j-1) .AND. W(j)<W(j-1) ) THEN
               Info = -5
               EXIT
            ENDIF
         ENDDO
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('ZSTEIN',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 .OR. M==0 ) THEN
         RETURN
      ELSEIF ( N==1 ) THEN
         Z(1,1) = CONE
         RETURN
      ENDIF
!
!     Get machine constants.
!
      eps = DLAMCH('Precision')
!
!     Initialize seed for random number generator DLARNV.
!
      DO i = 1 , 4
         iseed(i) = 1
      ENDDO
!
!     Initialize pointers.
!
      indrv1 = 0
      indrv2 = indrv1 + N
      indrv3 = indrv2 + N
      indrv4 = indrv3 + N
      indrv5 = indrv4 + N
!
!     Compute eigenvectors of matrix blocks.
!
      j1 = 1
      DO nblk = 1 , Iblock(M)
!
!        Find starting and ending indices of block nblk.
!
         IF ( nblk==1 ) THEN
            b1 = 1
         ELSE
            b1 = Isplit(nblk-1) + 1
         ENDIF
         bn = Isplit(nblk)
         blksiz = bn - b1 + 1
         IF ( blksiz/=1 ) THEN
            gpind = j1
!
!        Compute reorthogonalization criterion and stopping criterion.
!
            onenrm = ABS(D(b1)) + ABS(E(b1))
            onenrm = MAX(onenrm,ABS(D(bn))+ABS(E(bn-1)))
            DO i = b1 + 1 , bn - 1
               onenrm = MAX(onenrm,ABS(D(i))+ABS(E(i-1))+ABS(E(i)))
            ENDDO
            ortol = ODM3*onenrm
!
            dtpcrt = SQRT(ODM1/blksiz)
         ENDIF
!
!        Loop through eigenvalues of block nblk.
!
         jblk = 0
         DO j = j1 , M
            IF ( Iblock(j)/=nblk ) THEN
               j1 = j
               EXIT
            ENDIF
            jblk = jblk + 1
            xj = W(j)
!
!           Skip all the work if the block size is one.
!
            IF ( blksiz==1 ) THEN
               Work(indrv1+1) = ONE
               GOTO 20
            ENDIF
!
!           If eigenvalues j and j-1 are too close, add a relatively
!           small perturbation.
!
            IF ( jblk>1 ) THEN
               eps1 = ABS(eps*xj)
               pertol = TEN*eps1
               sep = xj - xjm
               IF ( sep<pertol ) xj = xjm + pertol
            ENDIF
!
            its = 0
            nrmchk = 0
!
!           Get random starting vector.
!
            CALL DLARNV(2,iseed,blksiz,Work(indrv1+1))
!
!           Copy the matrix T so it won't be destroyed in factorization.
!
            CALL DCOPY(blksiz,D(b1),1,Work(indrv4+1),1)
            CALL DCOPY(blksiz-1,E(b1),1,Work(indrv2+2),1)
            CALL DCOPY(blksiz-1,E(b1),1,Work(indrv3+1),1)
!
!           Compute LU factors with partial pivoting  ( PT = LU )
!
            tol = ZERO
            CALL DLAGTF(blksiz,Work(indrv4+1),xj,Work(indrv2+2),        &
     &                  Work(indrv3+1),tol,Work(indrv5+1),Iwork,iinfo)
            DO
!
!           Update iteration count.
!
               its = its + 1
               IF ( its>MAXITS ) THEN
!
!           If stopping criterion was not satisfied, update info and
!           store eigenvector number in array ifail.
!
                  Info = Info + 1
                  Ifail(Info) = j
                  EXIT
               ELSE
!
!           Normalize and scale the righthand side vector Pb.
!
                  jmax = IDAMAX(blksiz,Work(indrv1+1),1)
                  scl = blksiz*onenrm*MAX(eps,ABS(Work(indrv4+blksiz))) &
     &                  /ABS(Work(indrv1+jmax))
                  CALL DSCAL(blksiz,scl,Work(indrv1+1),1)
!
!           Solve the system LU = Pb.
!
                  CALL DLAGTS(-1,blksiz,Work(indrv4+1),Work(indrv2+2),  &
     &                        Work(indrv3+1),Work(indrv5+1),Iwork,      &
     &                        Work(indrv1+1),tol,iinfo)
!
!           Reorthogonalize by modified Gram-Schmidt if eigenvalues are
!           close enough.
!
                  IF ( jblk/=1 ) THEN
                     IF ( ABS(xj-xjm)>ortol ) gpind = j
                     IF ( gpind/=j ) THEN
                        DO i = gpind , j - 1
                           ztr = ZERO
                           DO jr = 1 , blksiz
                              ztr = ztr + Work(indrv1+jr)               &
     &                              *DBLE(Z(b1-1+jr,i))
                           ENDDO
                           DO jr = 1 , blksiz
                              Work(indrv1+jr) = Work(indrv1+jr)         &
     &                           - ztr*DBLE(Z(b1-1+jr,i))
                           ENDDO
                        ENDDO
                     ENDIF
                  ENDIF
!
!           Check the infinity norm of the iterate.
!
                  jmax = IDAMAX(blksiz,Work(indrv1+1),1)
                  nrm = ABS(Work(indrv1+jmax))
!
!           Continue for additional iterations after norm reaches
!           stopping criterion.
!
                  IF ( nrm>=dtpcrt ) THEN
                     nrmchk = nrmchk + 1
!
                     IF ( nrmchk>=EXTRA+1 ) EXIT
                  ENDIF
               ENDIF
            ENDDO
!
!           Accept iterate as jth eigenvector.
!
            scl = ONE/DNRM2(blksiz,Work(indrv1+1),1)
            jmax = IDAMAX(blksiz,Work(indrv1+1),1)
            IF ( Work(indrv1+jmax)<ZERO ) scl = -scl
            CALL DSCAL(blksiz,scl,Work(indrv1+1),1)
 20         DO i = 1 , N
               Z(i,j) = CZERO
            ENDDO
            DO i = 1 , blksiz
               Z(b1+i-1,j) = DCMPLX(Work(indrv1+i),ZERO)
            ENDDO
!
!           Save the shift to check eigenvalue spacing at next
!           iteration.
!
            xjm = xj
!
         ENDDO
      ENDDO
!
!
!     End of ZSTEIN
!
      END SUBROUTINE ZSTEIN
