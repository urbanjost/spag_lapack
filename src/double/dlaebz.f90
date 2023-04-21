!*==dlaebz.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b DLAEBZ computes the number of eigenvalues of a real symmetric tridiagonal matrix which are less than or equal to a given value, and performs other tasks required by the routine sstebz.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLAEBZ + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaebz.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaebz.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaebz.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLAEBZ( IJOB, NITMAX, N, MMAX, MINP, NBMIN, ABSTOL,
!                          RELTOL, PIVMIN, D, E, E2, NVAL, AB, C, MOUT,
!                          NAB, WORK, IWORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            IJOB, INFO, MINP, MMAX, MOUT, N, NBMIN, NITMAX
!       DOUBLE PRECISION   ABSTOL, PIVMIN, RELTOL
!       ..
!       .. Array Arguments ..
!       INTEGER            IWORK( * ), NAB( MMAX, * ), NVAL( * )
!       DOUBLE PRECISION   AB( MMAX, * ), C( * ), D( * ), E( * ), E2( * ),
!      $                   WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLAEBZ contains the iteration loops which compute and use the
!> function N(w), which is the count of eigenvalues of a symmetric
!> tridiagonal matrix T less than or equal to its argument  w.  It
!> performs a choice of two types of loops:
!>
!> IJOB=1, followed by
!> IJOB=2: It takes as input a list of intervals and returns a list of
!>         sufficiently small intervals whose union contains the same
!>         eigenvalues as the union of the original intervals.
!>         The input intervals are (AB(j,1),AB(j,2)], j=1,...,MINP.
!>         The output interval (AB(j,1),AB(j,2)] will contain
!>         eigenvalues NAB(j,1)+1,...,NAB(j,2), where 1 <= j <= MOUT.
!>
!> IJOB=3: It performs a binary search in each input interval
!>         (AB(j,1),AB(j,2)] for a point  w(j)  such that
!>         N(w(j))=NVAL(j), and uses  C(j)  as the starting point of
!>         the search.  If such a w(j) is found, then on output
!>         AB(j,1)=AB(j,2)=w.  If no such w(j) is found, then on output
!>         (AB(j,1),AB(j,2)] will be a small interval containing the
!>         point where N(w) jumps through NVAL(j), unless that point
!>         lies outside the initial interval.
!>
!> Note that the intervals are in all cases half-open intervals,
!> i.e., of the form  (a,b] , which includes  b  but not  a .
!>
!> To avoid underflow, the matrix should be scaled so that its largest
!> element is no greater than  overflow**(1/2) * underflow**(1/4)
!> in absolute value.  To assure the most accurate computation
!> of small eigenvalues, the matrix should be scaled to be
!> not much smaller than that, either.
!>
!> See W. Kahan "Accurate Eigenvalues of a Symmetric Tridiagonal
!> Matrix", Report CS41, Computer Science Dept., Stanford
!> University, July 21, 1966
!>
!> Note: the arguments are, in general, *not* checked for unreasonable
!> values.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] IJOB
!> \verbatim
!>          IJOB is INTEGER
!>          Specifies what is to be done:
!>          = 1:  Compute NAB for the initial intervals.
!>          = 2:  Perform bisection iteration to find eigenvalues of T.
!>          = 3:  Perform bisection iteration to invert N(w), i.e.,
!>                to find a point which has a specified number of
!>                eigenvalues of T to its left.
!>          Other values will cause DLAEBZ to return with INFO=-1.
!> \endverbatim
!>
!> \param[in] NITMAX
!> \verbatim
!>          NITMAX is INTEGER
!>          The maximum number of "levels" of bisection to be
!>          performed, i.e., an interval of width W will not be made
!>          smaller than 2^(-NITMAX) * W.  If not all intervals
!>          have converged after NITMAX iterations, then INFO is set
!>          to the number of non-converged intervals.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The dimension n of the tridiagonal matrix T.  It must be at
!>          least 1.
!> \endverbatim
!>
!> \param[in] MMAX
!> \verbatim
!>          MMAX is INTEGER
!>          The maximum number of intervals.  If more than MMAX intervals
!>          are generated, then DLAEBZ will quit with INFO=MMAX+1.
!> \endverbatim
!>
!> \param[in] MINP
!> \verbatim
!>          MINP is INTEGER
!>          The initial number of intervals.  It may not be greater than
!>          MMAX.
!> \endverbatim
!>
!> \param[in] NBMIN
!> \verbatim
!>          NBMIN is INTEGER
!>          The smallest number of intervals that should be processed
!>          using a vector loop.  If zero, then only the scalar loop
!>          will be used.
!> \endverbatim
!>
!> \param[in] ABSTOL
!> \verbatim
!>          ABSTOL is DOUBLE PRECISION
!>          The minimum (absolute) width of an interval.  When an
!>          interval is narrower than ABSTOL, or than RELTOL times the
!>          larger (in magnitude) endpoint, then it is considered to be
!>          sufficiently small, i.e., converged.  This must be at least
!>          zero.
!> \endverbatim
!>
!> \param[in] RELTOL
!> \verbatim
!>          RELTOL is DOUBLE PRECISION
!>          The minimum relative width of an interval.  When an interval
!>          is narrower than ABSTOL, or than RELTOL times the larger (in
!>          magnitude) endpoint, then it is considered to be
!>          sufficiently small, i.e., converged.  Note: this should
!>          always be at least radix*machine epsilon.
!> \endverbatim
!>
!> \param[in] PIVMIN
!> \verbatim
!>          PIVMIN is DOUBLE PRECISION
!>          The minimum absolute value of a "pivot" in the Sturm
!>          sequence loop.
!>          This must be at least  max |e(j)**2|*safe_min  and at
!>          least safe_min, where safe_min is at least
!>          the smallest number that can divide one without overflow.
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is DOUBLE PRECISION array, dimension (N)
!>          The diagonal elements of the tridiagonal matrix T.
!> \endverbatim
!>
!> \param[in] E
!> \verbatim
!>          E is DOUBLE PRECISION array, dimension (N)
!>          The offdiagonal elements of the tridiagonal matrix T in
!>          positions 1 through N-1.  E(N) is arbitrary.
!> \endverbatim
!>
!> \param[in] E2
!> \verbatim
!>          E2 is DOUBLE PRECISION array, dimension (N)
!>          The squares of the offdiagonal elements of the tridiagonal
!>          matrix T.  E2(N) is ignored.
!> \endverbatim
!>
!> \param[in,out] NVAL
!> \verbatim
!>          NVAL is INTEGER array, dimension (MINP)
!>          If IJOB=1 or 2, not referenced.
!>          If IJOB=3, the desired values of N(w).  The elements of NVAL
!>          will be reordered to correspond with the intervals in AB.
!>          Thus, NVAL(j) on output will not, in general be the same as
!>          NVAL(j) on input, but it will correspond with the interval
!>          (AB(j,1),AB(j,2)] on output.
!> \endverbatim
!>
!> \param[in,out] AB
!> \verbatim
!>          AB is DOUBLE PRECISION array, dimension (MMAX,2)
!>          The endpoints of the intervals.  AB(j,1) is  a(j), the left
!>          endpoint of the j-th interval, and AB(j,2) is b(j), the
!>          right endpoint of the j-th interval.  The input intervals
!>          will, in general, be modified, split, and reordered by the
!>          calculation.
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is DOUBLE PRECISION array, dimension (MMAX)
!>          If IJOB=1, ignored.
!>          If IJOB=2, workspace.
!>          If IJOB=3, then on input C(j) should be initialized to the
!>          first search point in the binary search.
!> \endverbatim
!>
!> \param[out] MOUT
!> \verbatim
!>          MOUT is INTEGER
!>          If IJOB=1, the number of eigenvalues in the intervals.
!>          If IJOB=2 or 3, the number of intervals output.
!>          If IJOB=3, MOUT will equal MINP.
!> \endverbatim
!>
!> \param[in,out] NAB
!> \verbatim
!>          NAB is INTEGER array, dimension (MMAX,2)
!>          If IJOB=1, then on output NAB(i,j) will be set to N(AB(i,j)).
!>          If IJOB=2, then on input, NAB(i,j) should be set.  It must
!>             satisfy the condition:
!>             N(AB(i,1)) <= NAB(i,1) <= NAB(i,2) <= N(AB(i,2)),
!>             which means that in interval i only eigenvalues
!>             NAB(i,1)+1,...,NAB(i,2) will be considered.  Usually,
!>             NAB(i,j)=N(AB(i,j)), from a previous call to DLAEBZ with
!>             IJOB=1.
!>             On output, NAB(i,j) will contain
!>             max(na(k),min(nb(k),N(AB(i,j)))), where k is the index of
!>             the input interval that the output interval
!>             (AB(j,1),AB(j,2)] came from, and na(k) and nb(k) are the
!>             the input values of NAB(k,1) and NAB(k,2).
!>          If IJOB=3, then on output, NAB(i,j) contains N(AB(i,j)),
!>             unless N(w) > NVAL(i) for all search points  w , in which
!>             case NAB(i,1) will not be modified, i.e., the output
!>             value will be the same as the input value (modulo
!>             reorderings -- see NVAL and AB), or unless N(w) < NVAL(i)
!>             for all search points  w , in which case NAB(i,2) will
!>             not be modified.  Normally, NAB should be set to some
!>             distinctive value(s) before DLAEBZ is called.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (MMAX)
!>          Workspace.
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (MMAX)
!>          Workspace.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:       All intervals converged.
!>          = 1--MMAX: The last INFO intervals did not converge.
!>          = MMAX+1:  More than MMAX intervals were generated.
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
!> \ingroup OTHERauxiliary
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>      This routine is intended to be called only by other LAPACK
!>  routines, thus the interface is less user-friendly.  It is intended
!>  for two purposes:
!>
!>  (a) finding eigenvalues.  In this case, DLAEBZ should have one or
!>      more initial intervals set up in AB, and DLAEBZ should be called
!>      with IJOB=1.  This sets up NAB, and also counts the eigenvalues.
!>      Intervals with no eigenvalues would usually be thrown out at
!>      this point.  Also, if not all the eigenvalues in an interval i
!>      are desired, NAB(i,1) can be increased or NAB(i,2) decreased.
!>      For example, set NAB(i,1)=NAB(i,2)-1 to get the largest
!>      eigenvalue.  DLAEBZ is then called with IJOB=2 and MMAX
!>      no smaller than the value of MOUT returned by the call with
!>      IJOB=1.  After this (IJOB=2) call, eigenvalues NAB(i,1)+1
!>      through NAB(i,2) are approximately AB(i,1) (or AB(i,2)) to the
!>      tolerance specified by ABSTOL and RELTOL.
!>
!>  (b) finding an interval (a',b'] containing eigenvalues w(f),...,w(l).
!>      In this case, start with a Gershgorin interval  (a,b).  Set up
!>      AB to contain 2 search intervals, both initially (a,b).  One
!>      NVAL element should contain  f-1  and the other should contain  l
!>      , while C should contain a and b, resp.  NAB(i,1) should be -1
!>      and NAB(i,2) should be N+1, to flag an error if the desired
!>      interval does not lie in (a,b).  DLAEBZ is then called with
!>      IJOB=3.  On exit, if w(f-1) < w(f), then one of the intervals --
!>      j -- will have AB(j,1)=AB(j,2) and NAB(j,1)=NAB(j,2)=f-1, while
!>      if, to the specified tolerance, w(f-k)=...=w(f+r), k > 0 and r
!>      >= 0, then the interval will have  N(AB(j,1))=NAB(j,1)=f-k and
!>      N(AB(j,2))=NAB(j,2)=f+r.  The cases w(l) < w(l+1) and
!>      w(l-r)=...=w(l+k) are handled similarly.
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE DLAEBZ(Ijob,Nitmax,N,Mmax,Minp,Nbmin,Abstol,Reltol,    &
     &                  Pivmin,D,E,E2,Nval,Ab,C,Mout,Nab,Work,Iwork,    &
     &                  Info)
      IMPLICIT NONE
!*--DLAEBZ323
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Ijob , Info , Minp , Mmax , Mout , N , Nbmin , Nitmax
      DOUBLE PRECISION Abstol , Pivmin , Reltol
!     ..
!     .. Array Arguments ..
      INTEGER Iwork(*) , Nab(Mmax,*) , Nval(*)
      DOUBLE PRECISION Ab(Mmax,*) , C(*) , D(*) , E(*) , E2(*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , TWO , HALF
      PARAMETER (ZERO=0.0D0,TWO=2.0D0,HALF=1.0D0/TWO)
!     ..
!     .. Local Scalars ..
      INTEGER itmp1 , itmp2 , j , ji , jit , jp , kf , kfnew , kl ,     &
     &        klnew
      DOUBLE PRECISION tmp1 , tmp2
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX , MIN
!     ..
!     .. Executable Statements ..
!
!     Check for Errors
!
      Info = 0
      IF ( Ijob<1 .OR. Ijob>3 ) THEN
         Info = -1
         RETURN
      ENDIF
!
!     Initialize NAB
!
      IF ( Ijob==1 ) THEN
!
!        Compute the number of eigenvalues in the initial intervals.
!
         Mout = 0
         DO ji = 1 , Minp
            DO jp = 1 , 2
               tmp1 = D(1) - Ab(ji,jp)
               IF ( ABS(tmp1)<Pivmin ) tmp1 = -Pivmin
               Nab(ji,jp) = 0
               IF ( tmp1<=ZERO ) Nab(ji,jp) = 1
!
               DO j = 2 , N
                  tmp1 = D(j) - E2(j-1)/tmp1 - Ab(ji,jp)
                  IF ( ABS(tmp1)<Pivmin ) tmp1 = -Pivmin
                  IF ( tmp1<=ZERO ) Nab(ji,jp) = Nab(ji,jp) + 1
               ENDDO
            ENDDO
            Mout = Mout + Nab(ji,2) - Nab(ji,1)
         ENDDO
         RETURN
      ENDIF
!
!     Initialize for loop
!
!     KF and KL have the following meaning:
!        Intervals 1,...,KF-1 have converged.
!        Intervals KF,...,KL  still need to be refined.
!
      kf = 1
      kl = Minp
!
!     If IJOB=2, initialize C.
!     If IJOB=3, use the user-supplied starting point.
!
      IF ( Ijob==2 ) THEN
         DO ji = 1 , Minp
            C(ji) = HALF*(Ab(ji,1)+Ab(ji,2))
         ENDDO
      ENDIF
!
!     Iteration loop
!
      DO jit = 1 , Nitmax
!
!        Loop over intervals
!
         IF ( kl-kf+1>=Nbmin .AND. Nbmin>0 ) THEN
!
!           Begin of Parallel Version of the loop
!
            DO ji = kf , kl
!
!              Compute N(c), the number of eigenvalues less than c
!
               Work(ji) = D(1) - C(ji)
               Iwork(ji) = 0
               IF ( Work(ji)<=Pivmin ) THEN
                  Iwork(ji) = 1
                  Work(ji) = MIN(Work(ji),-Pivmin)
               ENDIF
!
               DO j = 2 , N
                  Work(ji) = D(j) - E2(j-1)/Work(ji) - C(ji)
                  IF ( Work(ji)<=Pivmin ) THEN
                     Iwork(ji) = Iwork(ji) + 1
                     Work(ji) = MIN(Work(ji),-Pivmin)
                  ENDIF
               ENDDO
            ENDDO
!
            IF ( Ijob<=2 ) THEN
!
!              IJOB=2: Choose all intervals containing eigenvalues.
!
               klnew = kl
               DO ji = kf , kl
!
!                 Insure that N(w) is monotone
!
                  Iwork(ji) = MIN(Nab(ji,2),MAX(Nab(ji,1),Iwork(ji)))
!
!                 Update the Queue -- add intervals if both halves
!                 contain eigenvalues.
!
                  IF ( Iwork(ji)==Nab(ji,2) ) THEN
!
!                    No eigenvalue in the upper interval:
!                    just use the lower interval.
!
                     Ab(ji,2) = C(ji)
!
                  ELSEIF ( Iwork(ji)==Nab(ji,1) ) THEN
!
!                    No eigenvalue in the lower interval:
!                    just use the upper interval.
!
                     Ab(ji,1) = C(ji)
                  ELSE
                     klnew = klnew + 1
                     IF ( klnew<=Mmax ) THEN
!
!                       Eigenvalue in both intervals -- add upper to
!                       queue.
!
                        Ab(klnew,2) = Ab(ji,2)
                        Nab(klnew,2) = Nab(ji,2)
                        Ab(klnew,1) = C(ji)
                        Nab(klnew,1) = Iwork(ji)
                        Ab(ji,2) = C(ji)
                        Nab(ji,2) = Iwork(ji)
                     ELSE
                        Info = Mmax + 1
                     ENDIF
                  ENDIF
               ENDDO
               IF ( Info/=0 ) RETURN
               kl = klnew
            ELSE
!
!              IJOB=3: Binary search.  Keep only the interval containing
!                      w   s.t. N(w) = NVAL
!
               DO ji = kf , kl
                  IF ( Iwork(ji)<=Nval(ji) ) THEN
                     Ab(ji,1) = C(ji)
                     Nab(ji,1) = Iwork(ji)
                  ENDIF
                  IF ( Iwork(ji)>=Nval(ji) ) THEN
                     Ab(ji,2) = C(ji)
                     Nab(ji,2) = Iwork(ji)
                  ENDIF
               ENDDO
            ENDIF
!
         ELSE
!
!           End of Parallel Version of the loop
!
!           Begin of Serial Version of the loop
!
            klnew = kl
            DO ji = kf , kl
!
!              Compute N(w), the number of eigenvalues less than w
!
               tmp1 = C(ji)
               tmp2 = D(1) - tmp1
               itmp1 = 0
               IF ( tmp2<=Pivmin ) THEN
                  itmp1 = 1
                  tmp2 = MIN(tmp2,-Pivmin)
               ENDIF
!
               DO j = 2 , N
                  tmp2 = D(j) - E2(j-1)/tmp2 - tmp1
                  IF ( tmp2<=Pivmin ) THEN
                     itmp1 = itmp1 + 1
                     tmp2 = MIN(tmp2,-Pivmin)
                  ENDIF
               ENDDO
!
               IF ( Ijob<=2 ) THEN
!
!                 IJOB=2: Choose all intervals containing eigenvalues.
!
!                 Insure that N(w) is monotone
!
                  itmp1 = MIN(Nab(ji,2),MAX(Nab(ji,1),itmp1))
!
!                 Update the Queue -- add intervals if both halves
!                 contain eigenvalues.
!
                  IF ( itmp1==Nab(ji,2) ) THEN
!
!                    No eigenvalue in the upper interval:
!                    just use the lower interval.
!
                     Ab(ji,2) = tmp1
!
                  ELSEIF ( itmp1==Nab(ji,1) ) THEN
!
!                    No eigenvalue in the lower interval:
!                    just use the upper interval.
!
                     Ab(ji,1) = tmp1
                  ELSEIF ( klnew<Mmax ) THEN
!
!                    Eigenvalue in both intervals -- add upper to queue.
!
                     klnew = klnew + 1
                     Ab(klnew,2) = Ab(ji,2)
                     Nab(klnew,2) = Nab(ji,2)
                     Ab(klnew,1) = tmp1
                     Nab(klnew,1) = itmp1
                     Ab(ji,2) = tmp1
                     Nab(ji,2) = itmp1
                  ELSE
                     Info = Mmax + 1
                     RETURN
                  ENDIF
               ELSE
!
!                 IJOB=3: Binary search.  Keep only the interval
!                         containing  w  s.t. N(w) = NVAL
!
                  IF ( itmp1<=Nval(ji) ) THEN
                     Ab(ji,1) = tmp1
                     Nab(ji,1) = itmp1
                  ENDIF
                  IF ( itmp1>=Nval(ji) ) THEN
                     Ab(ji,2) = tmp1
                     Nab(ji,2) = itmp1
                  ENDIF
               ENDIF
            ENDDO
            kl = klnew
!
         ENDIF
!
!        Check for convergence
!
         kfnew = kf
         DO ji = kf , kl
            tmp1 = ABS(Ab(ji,2)-Ab(ji,1))
            tmp2 = MAX(ABS(Ab(ji,2)),ABS(Ab(ji,1)))
            IF ( tmp1<MAX(Abstol,Pivmin,Reltol*tmp2) .OR. Nab(ji,1)     &
     &           >=Nab(ji,2) ) THEN
!
!              Converged -- Swap with position KFNEW,
!                           then increment KFNEW
!
               IF ( ji>kfnew ) THEN
                  tmp1 = Ab(ji,1)
                  tmp2 = Ab(ji,2)
                  itmp1 = Nab(ji,1)
                  itmp2 = Nab(ji,2)
                  Ab(ji,1) = Ab(kfnew,1)
                  Ab(ji,2) = Ab(kfnew,2)
                  Nab(ji,1) = Nab(kfnew,1)
                  Nab(ji,2) = Nab(kfnew,2)
                  Ab(kfnew,1) = tmp1
                  Ab(kfnew,2) = tmp2
                  Nab(kfnew,1) = itmp1
                  Nab(kfnew,2) = itmp2
                  IF ( Ijob==3 ) THEN
                     itmp1 = Nval(ji)
                     Nval(ji) = Nval(kfnew)
                     Nval(kfnew) = itmp1
                  ENDIF
               ENDIF
               kfnew = kfnew + 1
            ENDIF
         ENDDO
         kf = kfnew
!
!        Choose Midpoints
!
         DO ji = kf , kl
            C(ji) = HALF*(Ab(ji,1)+Ab(ji,2))
         ENDDO
!
!        If no more intervals to refine, quit.
!
         IF ( kf>kl ) EXIT
      ENDDO
!
!     Converged
!
      Info = MAX(kl+1-kf,0)
      Mout = kl
!
!
!     End of DLAEBZ
!
      END SUBROUTINE DLAEBZ
