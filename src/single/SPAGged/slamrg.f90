!*==slamrg.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b SLAMRG creates a permutation list to merge the entries of two independently sorted sets into a single set sorted in ascending order.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SLAMRG + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slamrg.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slamrg.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slamrg.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLAMRG( N1, N2, A, STRD1, STRD2, INDEX )
!
!       .. Scalar Arguments ..
!       INTEGER            N1, N2, STRD1, STRD2
!       ..
!       .. Array Arguments ..
!       INTEGER            INDEX( * )
!       REAL               A( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLAMRG will create a permutation list which will merge the elements
!> of A (which is composed of two independently sorted sets) into a
!> single set which is sorted in ascending order.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N1
!> \verbatim
!>          N1 is INTEGER
!> \endverbatim
!>
!> \param[in] N2
!> \verbatim
!>          N2 is INTEGER
!>         These arguments contain the respective lengths of the two
!>         sorted lists to be merged.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is REAL array, dimension (N1+N2)
!>         The first N1 elements of A contain a list of numbers which
!>         are sorted in either ascending or descending order.  Likewise
!>         for the final N2 elements.
!> \endverbatim
!>
!> \param[in] STRD1
!> \verbatim
!>          STRD1 is INTEGER
!> \endverbatim
!>
!> \param[in] STRD2
!> \verbatim
!>          STRD2 is INTEGER
!>         These are the strides to be taken through the array A.
!>         Allowable strides are 1 and -1.  They indicate whether a
!>         subset of A is sorted in ascending (STRDx = 1) or descending
!>         (STRDx = -1) order.
!> \endverbatim
!>
!> \param[out] INDEX
!> \verbatim
!>          INDEX is INTEGER array, dimension (N1+N2)
!>         On exit this array will contain a permutation such that
!>         if B( I ) = A( INDEX( I ) ) for I=1,N1+N2, then B will be
!>         sorted in ascending order.
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
!> \ingroup auxOTHERcomputational
!
!  =====================================================================
      SUBROUTINE SLAMRG(N1,N2,A,Strd1,Strd2,Index)
      IMPLICIT NONE
!*--SLAMRG103
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: N1
      INTEGER , INTENT(IN) :: N2
      REAL , INTENT(IN) , DIMENSION(*) :: A
      INTEGER , INTENT(IN) :: Strd1
      INTEGER , INTENT(IN) :: Strd2
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Index
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: i , ind1 , ind2 , n1sv , n2sv
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
!     ..
!     .. Executable Statements ..
!
      n1sv = N1
      n2sv = N2
      IF ( Strd1>0 ) THEN
         ind1 = 1
      ELSE
         ind1 = N1
      ENDIF
      IF ( Strd2>0 ) THEN
         ind2 = 1 + N1
      ELSE
         ind2 = N1 + N2
      ENDIF
      i = 1
      DO
!     while ( (N1SV > 0) & (N2SV > 0) )
         IF ( n1sv>0 .AND. n2sv>0 ) THEN
            IF ( A(ind1)<=A(ind2) ) THEN
               Index(i) = ind1
               i = i + 1
               ind1 = ind1 + Strd1
               n1sv = n1sv - 1
            ELSE
               Index(i) = ind2
               i = i + 1
               ind2 = ind2 + Strd2
               n2sv = n2sv - 1
            ENDIF
            CYCLE
         ENDIF
!     end while
         IF ( n1sv==0 ) THEN
            DO n1sv = 1 , n2sv
               Index(i) = ind2
               i = i + 1
               ind2 = ind2 + Strd2
            ENDDO
         ELSE
!     N2SV .EQ. 0
            DO n2sv = 1 , n1sv
               Index(i) = ind1
               i = i + 1
               ind1 = ind1 + Strd1
            ENDDO
         ENDIF
         EXIT
      ENDDO
!
!
!     End of SLAMRG
!
      END SUBROUTINE SLAMRG
