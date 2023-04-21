!*==dlamrg.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b DLAMRG creates a permutation list to merge the entries of two independently sorted sets into a single set sorted in ascending order.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLAMRG + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlamrg.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlamrg.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlamrg.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLAMRG( N1, N2, A, DTRD1, DTRD2, INDEX )
!
!       .. Scalar Arguments ..
!       INTEGER            DTRD1, DTRD2, N1, N2
!       ..
!       .. Array Arguments ..
!       INTEGER            INDEX( * )
!       DOUBLE PRECISION   A( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLAMRG will create a permutation list which will merge the elements
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
!>          A is DOUBLE PRECISION array, dimension (N1+N2)
!>         The first N1 elements of A contain a list of numbers which
!>         are sorted in either ascending or descending order.  Likewise
!>         for the final N2 elements.
!> \endverbatim
!>
!> \param[in] DTRD1
!> \verbatim
!>          DTRD1 is INTEGER
!> \endverbatim
!>
!> \param[in] DTRD2
!> \verbatim
!>          DTRD2 is INTEGER
!>         These are the strides to be taken through the array A.
!>         Allowable strides are 1 and -1.  They indicate whether a
!>         subset of A is sorted in ascending (DTRDx = 1) or descending
!>         (DTRDx = -1) order.
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
      SUBROUTINE DLAMRG(N1,N2,A,Dtrd1,Dtrd2,Index)
      IMPLICIT NONE
!*--DLAMRG103
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2016
!
!     .. Scalar Arguments ..
      INTEGER Dtrd1 , Dtrd2 , N1 , N2
!     ..
!     .. Array Arguments ..
      INTEGER Index(*)
      DOUBLE PRECISION A(*)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      INTEGER i , ind1 , ind2 , n1sv , n2sv
!     ..
!     .. Executable Statements ..
!
      n1sv = N1
      n2sv = N2
      IF ( Dtrd1>0 ) THEN
         ind1 = 1
      ELSE
         ind1 = N1
      ENDIF
      IF ( Dtrd2>0 ) THEN
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
               ind1 = ind1 + Dtrd1
               n1sv = n1sv - 1
            ELSE
               Index(i) = ind2
               i = i + 1
               ind2 = ind2 + Dtrd2
               n2sv = n2sv - 1
            ENDIF
            CYCLE
         ENDIF
!     end while
         IF ( n1sv==0 ) THEN
            DO n1sv = 1 , n2sv
               Index(i) = ind2
               i = i + 1
               ind2 = ind2 + Dtrd2
            ENDDO
         ELSE
!     N2SV .EQ. 0
            DO n2sv = 1 , n1sv
               Index(i) = ind1
               i = i + 1
               ind1 = ind1 + Dtrd1
            ENDDO
         ENDIF
         EXIT
      ENDDO
!
!
!     End of DLAMRG
!
      END SUBROUTINE DLAMRG
