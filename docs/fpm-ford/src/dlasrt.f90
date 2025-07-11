!*==dlasrt.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b DLASRT sorts numbers in increasing or decreasing order.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLASRT + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlasrt.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlasrt.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlasrt.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLASRT( ID, N, D, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          ID
!       INTEGER            INFO, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   D( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> Sort the numbers in D in increasing order (if ID = 'I') or
!> in decreasing order (if ID = 'D' ).
!>
!> Use Quick Sort, reverting to Insertion sort on arrays of
!> size <= 20. Dimension of STACK limits N to about 2**32.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] ID
!> \verbatim
!>          ID is CHARACTER*1
!>          = 'I': sort D in increasing order;
!>          = 'D': sort D in decreasing order.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The length of the array D.
!> \endverbatim
!>
!> \param[in,out] D
!> \verbatim
!>          D is DOUBLE PRECISION array, dimension (N)
!>          On entry, the array to be sorted.
!>          On exit, D has been sorted into increasing order
!>          (D(1) <= ... <= D(N) ) or into decreasing order
!>          (D(1) >= ... >= D(N) ), depending on ID.
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
!> \date June 2016
!
!> \ingroup auxOTHERcomputational
!
!  =====================================================================
      SUBROUTINE DLASRT(Id,N,D,Info)
      IMPLICIT NONE
!*--DLASRT92
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2016
!
!     .. Scalar Arguments ..
      CHARACTER Id
      INTEGER Info , N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION D(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER SELECT
      PARAMETER (SELECT=20)
!     ..
!     .. Local Scalars ..
      INTEGER dir , endd , i , j , start , stkpnt
      DOUBLE PRECISION d1 , d2 , d3 , dmnmx , tmp
!     ..
!     .. Local Arrays ..
      INTEGER stack(2,32)
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL XERBLA
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      Info = 0
      dir = -1
      IF ( LSAME(Id,'D') ) THEN
         dir = 0
      ELSEIF ( LSAME(Id,'I') ) THEN
         dir = 1
      ENDIF
      IF ( dir==-1 ) THEN
         Info = -1
      ELSEIF ( N<0 ) THEN
         Info = -2
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('DLASRT',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N<=1 ) RETURN
!
      stkpnt = 1
      stack(1,1) = 1
      stack(2,1) = N
      DO
         start = stack(1,stkpnt)
         endd = stack(2,stkpnt)
         stkpnt = stkpnt - 1
         IF ( endd-start<=SELECT .AND. endd>start ) THEN
!
!        Do Insertion sort on D( START:ENDD )
!
            IF ( dir==0 ) THEN
!
!           Sort into decreasing order
!
               DO i = start + 1 , endd
                  DO j = i , start + 1 , -1
                     IF ( D(j)<=D(j-1) ) EXIT
                     dmnmx = D(j)
                     D(j) = D(j-1)
                     D(j-1) = dmnmx
                  ENDDO
               ENDDO
!
            ELSE
!
!           Sort into increasing order
!
               DO i = start + 1 , endd
                  DO j = i , start + 1 , -1
                     IF ( D(j)>=D(j-1) ) EXIT
                     dmnmx = D(j)
                     D(j) = D(j-1)
                     D(j-1) = dmnmx
                  ENDDO
               ENDDO
!
            ENDIF
!
         ELSEIF ( endd-start>SELECT ) THEN
!
!        Partition D( START:ENDD ) and stack parts, largest one first
!
!        Choose partition entry as median of 3
!
            d1 = D(start)
            d2 = D(endd)
            i = (start+endd)/2
            d3 = D(i)
            IF ( d1<d2 ) THEN
               IF ( d3<d1 ) THEN
                  dmnmx = d1
               ELSEIF ( d3<d2 ) THEN
                  dmnmx = d3
               ELSE
                  dmnmx = d2
               ENDIF
            ELSEIF ( d3<d2 ) THEN
               dmnmx = d2
            ELSEIF ( d3<d1 ) THEN
               dmnmx = d3
            ELSE
               dmnmx = d1
            ENDIF
!
            IF ( dir==0 ) THEN
!
!           Sort into decreasing order
!
               i = start - 1
               j = endd + 1
               DO
                  j = j - 1
                  IF ( D(j)>=dmnmx ) THEN
                     DO
                        i = i + 1
                        IF ( D(i)<=dmnmx ) THEN
                           IF ( i<j ) THEN
                              tmp = D(i)
                              D(i) = D(j)
                              D(j) = tmp
                              EXIT
                           ENDIF
                           IF ( j-start>endd-j-1 ) THEN
                              stkpnt = stkpnt + 1
                              stack(1,stkpnt) = start
                              stack(2,stkpnt) = j
                              stkpnt = stkpnt + 1
                              stack(1,stkpnt) = j + 1
                              stack(2,stkpnt) = endd
                           ELSE
                              stkpnt = stkpnt + 1
                              stack(1,stkpnt) = j + 1
                              stack(2,stkpnt) = endd
                              stkpnt = stkpnt + 1
                              stack(1,stkpnt) = start
                              stack(2,stkpnt) = j
                           ENDIF
                           GOTO 50
                        ENDIF
                     ENDDO
                  ENDIF
               ENDDO
            ELSE
!
!           Sort into increasing order
!
               i = start - 1
               j = endd + 1
               DO
                  j = j - 1
                  IF ( D(j)<=dmnmx ) THEN
                     DO
                        i = i + 1
                        IF ( D(i)>=dmnmx ) THEN
                           IF ( i<j ) THEN
                              tmp = D(i)
                              D(i) = D(j)
                              D(j) = tmp
                              EXIT
                           ENDIF
                           IF ( j-start>endd-j-1 ) THEN
                              stkpnt = stkpnt + 1
                              stack(1,stkpnt) = start
                              stack(2,stkpnt) = j
                              stkpnt = stkpnt + 1
                              stack(1,stkpnt) = j + 1
                              stack(2,stkpnt) = endd
                           ELSE
                              stkpnt = stkpnt + 1
                              stack(1,stkpnt) = j + 1
                              stack(2,stkpnt) = endd
                              stkpnt = stkpnt + 1
                              stack(1,stkpnt) = start
                              stack(2,stkpnt) = j
                           ENDIF
                           GOTO 50
                        ENDIF
                     ENDDO
                  ENDIF
               ENDDO
            ENDIF
         ENDIF
 50      IF ( stkpnt<=0 ) EXIT
      ENDDO
!
!     End of DLASRT
!
      END SUBROUTINE DLASRT
