!*==ilaenv.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b ilaenv
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       INTEGER FUNCTION ilaenv( ISPEC, NAME, OPTS, N1, N2, N3,
!                        N4 )
!
!       .. Scalar Arguments ..
!       CHARACTER*( * )    NAME, OPTS
!       INTEGER            ISPEC, N1, N2, N3, N4
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ilaenv returns problem-dependent parameters for the local
!> environment.  See ISPEC for a description of the parameters.
!>
!> In this version, the problem-dependent parameters are contained in
!> the integer array IPARMS in the common block CLAENV and the value
!> with index ISPEC is copied to ilaenv.  This version of ilaenv is
!> to be used in conjunction with xlaenv in TESTING and TIMING.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] ISPEC
!> \verbatim
!>          ISPEC is INTEGER
!>          Specifies the parameter to be returned as the value of
!>          ilaenv.
!>          = 1: the optimal blocksize; if this value is 1, an unblocked
!>               algorithm will give the best performance.
!>          = 2: the minimum block size for which the block routine
!>               should be used; if the usable block size is less than
!>               this value, an unblocked routine should be used.
!>          = 3: the crossover point (in a block routine, for N less
!>               than this value, an unblocked routine should be used)
!>          = 4: the number of shifts, used in the nonsymmetric
!>               eigenvalue routines
!>          = 5: the minimum column dimension for blocking to be used;
!>               rectangular blocks must have dimension at least k by m,
!>               where k is given by ilaenv(2,...) and m by ilaenv(5,...)
!>          = 6: the crossover point for the SVD (when reducing an m by n
!>               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds
!>               this value, a QR factorization is used first to reduce
!>               the matrix to a triangular form.)
!>          = 7: the number of processors
!>          = 8: the crossover point for the multishift QR and QZ methods
!>               for nonsymmetric eigenvalue problems.
!>          = 9: maximum size of the subproblems at the bottom of the
!>               computation tree in the divide-and-conquer algorithm
!>          =10: ieee NaN arithmetic can be trusted not to trap
!>          =11: infinity arithmetic can be trusted not to trap
!>
!>          Other specifications (up to 100) can be added later.
!> \endverbatim
!>
!> \param[in] NAME
!> \verbatim
!>          NAME is CHARACTER*(*)
!>          The name of the calling subroutine.
!> \endverbatim
!>
!> \param[in] OPTS
!> \verbatim
!>          OPTS is CHARACTER*(*)
!>          The character options to the subroutine NAME, concatenated
!>          into a single character string.  For example, UPLO = 'U',
!>          TRANS = 'T', and DIAG = 'N' for a triangular routine would
!>          be specified as OPTS = 'UTN'.
!> \endverbatim
!>
!> \param[in] N1
!> \verbatim
!>          N1 is INTEGER
!> \endverbatim
!>
!> \param[in] N2
!> \verbatim
!>          N2 is INTEGER
!> \endverbatim
!>
!> \param[in] N3
!> \verbatim
!>          N3 is INTEGER
!> \endverbatim
!>
!> \param[in] N4
!> \verbatim
!>          N4 is INTEGER
!>
!>          Problem dimensions for the subroutine NAME; these may not all
!>          be required.
!> \endverbatim
!>
!> \return ilaenv
!> \verbatim
!>          ilaenv is INTEGER
!>          >= 0: the value of the parameter specified by ISPEC
!>          < 0:  if ilaenv = -k, the k-th argument had an illegal value.
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
!> \ingroup aux_lin
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The following conventions have been used when calling ilaenv from the
!>  LAPACK routines:
!>  1)  OPTS is a concatenation of all of the character options to
!>      subroutine NAME, in the same order that they appear in the
!>      argument list for NAME, even if they are not used in determining
!>      the value of the parameter specified by ISPEC.
!>  2)  The problem dimensions N1, N2, N3, N4 are specified in the order
!>      that they appear in the argument list for NAME.  N1 is used
!>      first, N2 second, and so on, and unused problem dimensions are
!>      passed a value of -1.
!>  3)  The parameter value returned by ilaenv is checked for validity in
!>      the calling subroutine.  For example, ilaenv is used to retrieve
!>      the optimal blocksize for STRTRI as follows:
!>
!>      NB = ilaenv( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )
!>      IF( NB.LE.1 ) NB = MAX( 1, N )
!> \endverbatim
!>
!  =====================================================================
      INTEGER FUNCTION ilaenv(Ispec,Name,Opts,N1,N2,N3,N4)
      IMPLICIT NONE
!*--ILAENV153
!
!  -- LAPACK test routine (version 3.8.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2017
!
!     .. Scalar Arguments ..
      CHARACTER*(*) Name , Opts
      INTEGER Ispec , N1 , N2 , N3 , N4
!     ..
!
!  =====================================================================
!
!     .. Intrinsic Functions ..
      INTRINSIC INT , MIN , REAL
!     ..
!     .. External Functions ..
      INTEGER IEEECK
EXTERNAL    IEEECK
!     ..
!     .. Arrays in Common ..
      INTEGER IPArms(100)
!     ..
!     .. Common blocks ..
      COMMON /CLAENV/ IPArms
!     ..
!     .. Save statement ..
      SAVE /CLAENV/ 
!     ..
!     .. Executable Statements ..
!
      IF ( Ispec>=1 .AND. Ispec<=5 ) THEN
!
!        Return a value from the common block.
!
         IF ( Name(2:6)=='GEQR ' ) THEN
            IF ( N3==2 ) THEN
               ilaenv = IPArms(2)
            ELSE
               ilaenv = IPArms(1)
            ENDIF
         ELSEIF ( Name(2:6)=='GELQ ' ) THEN
            IF ( N3==2 ) THEN
               ilaenv = IPArms(2)
            ELSE
               ilaenv = IPArms(1)
            ENDIF
         ELSE
            ilaenv = IPArms(Ispec)
         ENDIF
!
      ELSEIF ( Ispec==6 ) THEN
!
!        Compute SVD crossover point.
!
         ilaenv = INT(REAL(MIN(N1,N2))*1.6E0)
!
      ELSEIF ( Ispec>=7 .AND. Ispec<=9 ) THEN
!
!        Return a value from the common block.
!
         ilaenv = IPArms(Ispec)
!
      ELSEIF ( Ispec==10 ) THEN
!
!        IEEE NaN arithmetic can be trusted not to trap
!
!        ilaenv = 0
         ilaenv = 1
         IF ( ilaenv==1 ) ilaenv = IEEECK(1,0.0,1.0)
!
      ELSEIF ( Ispec==11 ) THEN
!
!        Infinity arithmetic can be trusted not to trap
!
!        ilaenv = 0
         ilaenv = 1
         IF ( ilaenv==1 ) ilaenv = IEEECK(0,0.0,1.0)
!
      ELSE
!
!        Invalid value for ISPEC
!
         ilaenv = -1
      ENDIF
!
!
!     End of ilaenv
!
      END FUNCTION ilaenv
!*==ilaenv2stage.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
      INTEGER FUNCTION ILAENV2STAGE(Ispec,Name,Opts,N1,N2,N3,N4)
      IMPLICIT NONE
!*--ILAENV2STAGE247
!     .. Scalar Arguments ..
      CHARACTER*(*) Name , Opts
      INTEGER Ispec , N1 , N2 , N3 , N4
!     ..
!
!  =====================================================================
!
!     .. Local variables ..
      INTEGER iispec
!     .. External Functions ..
      INTEGER IPARAM2STAGE
EXTERNAL    IPARAM2STAGE
!     ..
!     .. Arrays in Common ..
      INTEGER IPArms(100)
!     ..
!     .. Common blocks ..
      COMMON /CLAENV/ IPArms
!     ..
!     .. Save statement ..
      SAVE /CLAENV/ 
!     ..
!     .. Executable Statements ..
!
      IF ( (Ispec<1) .OR. (Ispec>5) ) THEN
!
!        Invalid value for ISPEC
!
         ILAENV2STAGE = -1
!
!     1 <= ISPEC <= 5: 2stage eigenvalues SVD routines.
!
      ELSEIF ( Ispec==1 ) THEN
         ILAENV2STAGE = IPArms(1)
      ELSE
         iispec = 16 + Ispec
         ILAENV2STAGE = IPARAM2STAGE(iispec,Name,Opts,N1,N2,N3,N4)
!
      ENDIF
!
      END FUNCTION ILAENV2STAGE
