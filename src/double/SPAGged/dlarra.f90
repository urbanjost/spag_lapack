!*==dlarra.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DLARRA computes the splitting points with the specified threshold.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLARRA + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlarra.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlarra.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlarra.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLARRA( N, D, E, E2, SPLTOL, TNRM,
!                           NSPLIT, ISPLIT, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, N, NSPLIT
!       DOUBLE PRECISION    SPLTOL, TNRM
!       ..
!       .. Array Arguments ..
!       INTEGER            ISPLIT( * )
!       DOUBLE PRECISION   D( * ), E( * ), E2( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> Compute the splitting points with threshold SPLTOL.
!> DLARRA sets any "small" off-diagonal elements to zero.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix. N > 0.
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is DOUBLE PRECISION array, dimension (N)
!>          On entry, the N diagonal elements of the tridiagonal
!>          matrix T.
!> \endverbatim
!>
!> \param[in,out] E
!> \verbatim
!>          E is DOUBLE PRECISION array, dimension (N)
!>          On entry, the first (N-1) entries contain the subdiagonal
!>          elements of the tridiagonal matrix T; E(N) need not be set.
!>          On exit, the entries E( ISPLIT( I ) ), 1 <= I <= NSPLIT,
!>          are set to zero, the other entries of E are untouched.
!> \endverbatim
!>
!> \param[in,out] E2
!> \verbatim
!>          E2 is DOUBLE PRECISION array, dimension (N)
!>          On entry, the first (N-1) entries contain the SQUARES of the
!>          subdiagonal elements of the tridiagonal matrix T;
!>          E2(N) need not be set.
!>          On exit, the entries E2( ISPLIT( I ) ),
!>          1 <= I <= NSPLIT, have been set to zero
!> \endverbatim
!>
!> \param[in] SPLTOL
!> \verbatim
!>          SPLTOL is DOUBLE PRECISION
!>          The threshold for splitting. Two criteria can be used:
!>          SPLTOL<0 : criterion based on absolute off-diagonal value
!>          SPLTOL>0 : criterion that preserves relative accuracy
!> \endverbatim
!>
!> \param[in] TNRM
!> \verbatim
!>          TNRM is DOUBLE PRECISION
!>          The norm of the matrix.
!> \endverbatim
!>
!> \param[out] NSPLIT
!> \verbatim
!>          NSPLIT is INTEGER
!>          The number of blocks T splits into. 1 <= NSPLIT <= N.
!> \endverbatim
!>
!> \param[out] ISPLIT
!> \verbatim
!>          ISPLIT is INTEGER array, dimension (N)
!>          The splitting points, at which T breaks up into blocks.
!>          The first block consists of rows/columns 1 to ISPLIT(1),
!>          the second of rows/columns ISPLIT(1)+1 through ISPLIT(2),
!>          etc., and the NSPLIT-th consists of rows/columns
!>          ISPLIT(NSPLIT-1)+1 through ISPLIT(NSPLIT)=N.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
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
!> \date June 2017
!
!> \ingroup OTHERauxiliary
!
!> \par Contributors:
!  ==================
!>
!> Beresford Parlett, University of California, Berkeley, USA \n
!> Jim Demmel, University of California, Berkeley, USA \n
!> Inderjit Dhillon, University of Texas, Austin, USA \n
!> Osni Marques, LBNL/NERSC, USA \n
!> Christof Voemel, University of California, Berkeley, USA
!
!  =====================================================================
      SUBROUTINE DLARRA(N,D,E,E2,Spltol,Tnrm,Nsplit,Isplit,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!*--DLARRA140
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: D
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: E
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(*) :: E2
      REAL(R8KIND) , INTENT(IN) :: Spltol
      REAL(R8KIND) , INTENT(IN) :: Tnrm
      INTEGER , INTENT(INOUT) :: Nsplit
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Isplit
      INTEGER , INTENT(OUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: eabs , tmp1
      INTEGER :: i
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
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
      Info = 0
!
!     Quick return if possible
!
      IF ( N<=0 ) RETURN
!
!     Compute splitting points
      Nsplit = 1
      IF ( Spltol<ZERO ) THEN
!        Criterion based on absolute off-diagonal value
         tmp1 = ABS(Spltol)*Tnrm
         DO i = 1 , N - 1
            eabs = ABS(E(i))
            IF ( eabs<=tmp1 ) THEN
               E(i) = ZERO
               E2(i) = ZERO
               Isplit(Nsplit) = i
               Nsplit = Nsplit + 1
            ENDIF
         ENDDO
      ELSE
!        Criterion that guarantees relative accuracy
         DO i = 1 , N - 1
            eabs = ABS(E(i))
            IF ( eabs<=Spltol*SQRT(ABS(D(i)))*SQRT(ABS(D(i+1))) ) THEN
               E(i) = ZERO
               E2(i) = ZERO
               Isplit(Nsplit) = i
               Nsplit = Nsplit + 1
            ENDIF
         ENDDO
      ENDIF
      Isplit(Nsplit) = N
 
!
!     End of DLARRA
!
      END SUBROUTINE DLARRA
