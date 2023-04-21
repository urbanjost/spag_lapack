!> \brief \b CCHKEE
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       PROGRAM CCHKEE
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CCHKEE tests the COMPLEX LAPACK subroutines for the matrix
!> eigenvalue problem.  The test paths in this version are
!>
!> NEP (Nonsymmetric Eigenvalue Problem):
!>     Test CGEHRD, CUNGHR, CHSEQR, CTREVC, CHSEIN, and CUNMHR
!>
!> SEP (Hermitian Eigenvalue Problem):
!>     Test CHETRD, CUNGTR, CSTEQR, CSTERF, CSTEIN, CSTEDC,
!>     and drivers CHEEV(X), CHBEV(X), CHPEV(X),
!>                 CHEEVD,   CHBEVD,   CHPEVD
!>
!> SVD (Singular Value Decomposition):
!>     Test CGEBRD, CUNGBR, and CBDSQR
!>     and the drivers CGESVD, CGESDD
!>
!> CEV (Nonsymmetric Eigenvalue/eigenvector Driver):
!>     Test CGEEV
!>
!> CES (Nonsymmetric Schur form Driver):
!>     Test CGEES
!>
!> CVX (Nonsymmetric Eigenvalue/eigenvector Expert Driver):
!>     Test CGEEVX
!>
!> CSX (Nonsymmetric Schur form Expert Driver):
!>     Test CGEESX
!>
!> CGG (Generalized Nonsymmetric Eigenvalue Problem):
!>     Test CGGHD3, CGGBAL, CGGBAK, CHGEQZ, and CTGEVC
!>
!> CGS (Generalized Nonsymmetric Schur form Driver):
!>     Test CGGES
!>
!> CGV (Generalized Nonsymmetric Eigenvalue/eigenvector Driver):
!>     Test CGGEV
!>
!> CGX (Generalized Nonsymmetric Schur form Expert Driver):
!>     Test CGGESX
!>
!> CXV (Generalized Nonsymmetric Eigenvalue/eigenvector Expert Driver):
!>     Test CGGEVX
!>
!> CSG (Hermitian Generalized Eigenvalue Problem):
!>     Test CHEGST, CHEGV, CHEGVD, CHEGVX, CHPGST, CHPGV, CHPGVD,
!>     CHPGVX, CHBGST, CHBGV, CHBGVD, and CHBGVX
!>
!> CHB (Hermitian Band Eigenvalue Problem):
!>     Test CHBTRD
!>
!> CBB (Band Singular Value Decomposition):
!>     Test CGBBRD
!>
!> CEC (Eigencondition estimation):
!>     Test CTRSYL, CTREXC, CTRSNA, and CTRSEN
!>
!> CBL (Balancing a general matrix)
!>     Test CGEBAL
!>
!> CBK (Back transformation on a balanced matrix)
!>     Test CGEBAK
!>
!> CGL (Balancing a matrix pair)
!>     Test CGGBAL
!>
!> CGK (Back transformation on a matrix pair)
!>     Test CGGBAK
!>
!> GLM (Generalized Linear Regression Model):
!>     Tests CGGGLM
!>
!> GQR (Generalized QR and RQ factorizations):
!>     Tests CGGQRF and CGGRQF
!>
!> GSV (Generalized Singular Value Decomposition):
!>     Tests CGGSVD, CGGSVP, CTGSJA, CLAGS2, CLAPLL, and CLAPMT
!>
!> CSD (CS decomposition):
!>     Tests CUNCSD
!>
!> LSE (Constrained Linear Least Squares):
!>     Tests CGGLSE
!>
!> Each test path has a different set of inputs, but the data sets for
!> the driver routines xEV, xES, xVX, and xSX can be concatenated in a
!> single input file.  The first line of input should contain one of the
!> 3-character path names in columns 1-3.  The number of remaining lines
!> depends on what is found on the first line.
!>
!> The number of matrix types used in testing is often controllable from
!> the input file.  The number of matrix types for each path, and the
!> test routine that describes them, is as follows:
!>
!> Path name(s)  Types    Test routine
!>
!> CHS or NEP      21     CCHKHS
!> CST or SEP      21     CCHKST (routines)
!>                 18     CDRVST (drivers)
!> CBD or SVD      16     CCHKBD (routines)
!>                  5     CDRVBD (drivers)
!> CEV             21     CDRVEV
!> CES             21     CDRVES
!> CVX             21     CDRVVX
!> CSX             21     CDRVSX
!> CGG             26     CCHKGG (routines)
!> CGS             26     CDRGES
!> CGX              5     CDRGSX
!> CGV             26     CDRGEV
!> CXV              2     CDRGVX
!> CSG             21     CDRVSG
!> CHB             15     CCHKHB
!> CBB             15     CCHKBB
!> CEC              -     CCHKEC
!> CBL              -     CCHKBL
!> CBK              -     CCHKBK
!> CGL              -     CCHKGL
!> CGK              -     CCHKGK
!> GLM              8     CCKGLM
!> GQR              8     CCKGQR
!> GSV              8     CCKGSV
!> CSD              3     CCKCSD
!> LSE              8     CCKLSE
!>
!>-----------------------------------------------------------------------
!>
!> NEP input file:
!>
!> line 2:  NN, INTEGER
!>          Number of values of N.
!>
!> line 3:  NVAL, INTEGER array, dimension (NN)
!>          The values for the matrix dimension N.
!>
!> line 4:  NPARMS, INTEGER
!>          Number of values of the parameters NB, NBMIN, NX, NS, and
!>          MAXB.
!>
!> line 5:  NBVAL, INTEGER array, dimension (NPARMS)
!>          The values for the blocksize NB.
!>
!> line 6:  NBMIN, INTEGER array, dimension (NPARMS)
!>          The values for the minimum blocksize NBMIN.
!>
!> line 7:  NXVAL, INTEGER array, dimension (NPARMS)
!>          The values for the crossover point NX.
!>
!> line 8:  INMIN, INTEGER array, dimension (NPARMS)
!>          LAHQR vs TTQRE crossover point, >= 11
!>
!> line 9:  INWIN, INTEGER array, dimension (NPARMS)
!>          recommended deflation window size
!>
!> line 10: INIBL, INTEGER array, dimension (NPARMS)
!>          nibble crossover point
!>
!> line 11:  ISHFTS, INTEGER array, dimension (NPARMS)
!>          number of simultaneous shifts)
!>
!> line 12:  IACC22, INTEGER array, dimension (NPARMS)
!>          select structured matrix multiply: 0, 1 or 2)
!>
!> line 13: THRESH
!>          Threshold value for the test ratios.  Information will be
!>          printed about each test for which the test ratio is greater
!>          than or equal to the threshold.  To have all of the test
!>          ratios printed, use THRESH = 0.0 .
!>
!> line 14: NEWSD, INTEGER
!>          A code indicating how to set the random number seed.
!>          = 0:  Set the seed to a default value before each run
!>          = 1:  Initialize the seed to a default value only before the
!>                first run
!>          = 2:  Like 1, but use the seed values on the next line
!>
!> If line 14 was 2:
!>
!> line 15: INTEGER array, dimension (4)
!>          Four integer values for the random number seed.
!>
!> lines 15-EOF:  The remaining lines occur in sets of 1 or 2 and allow
!>          the user to specify the matrix types.  Each line contains
!>          a 3-character path name in columns 1-3, and the number
!>          of matrix types must be the first nonblank item in columns
!>          4-80.  If the number of matrix types is at least 1 but is
!>          less than the maximum number of possible types, a second
!>          line will be read to get the numbers of the matrix types to
!>          be used.  For example,
!> NEP 21
!>          requests all of the matrix types for the nonsymmetric
!>          eigenvalue problem, while
!> NEP  4
!> 9 10 11 12
!>          requests only matrices of type 9, 10, 11, and 12.
!>
!>          The valid 3-character path names are 'NEP' or 'CHS' for the
!>          nonsymmetric eigenvalue routines.
!>
!>-----------------------------------------------------------------------
!>
!> SEP or CSG input file:
!>
!> line 2:  NN, INTEGER
!>          Number of values of N.
!>
!> line 3:  NVAL, INTEGER array, dimension (NN)
!>          The values for the matrix dimension N.
!>
!> line 4:  NPARMS, INTEGER
!>          Number of values of the parameters NB, NBMIN, and NX.
!>
!> line 5:  NBVAL, INTEGER array, dimension (NPARMS)
!>          The values for the blocksize NB.
!>
!> line 6:  NBMIN, INTEGER array, dimension (NPARMS)
!>          The values for the minimum blocksize NBMIN.
!>
!> line 7:  NXVAL, INTEGER array, dimension (NPARMS)
!>          The values for the crossover point NX.
!>
!> line 8:  THRESH
!>          Threshold value for the test ratios.  Information will be
!>          printed about each test for which the test ratio is greater
!>          than or equal to the threshold.
!>
!> line 9:  TSTCHK, LOGICAL
!>          Flag indicating whether or not to test the LAPACK routines.
!>
!> line 10: TSTDRV, LOGICAL
!>          Flag indicating whether or not to test the driver routines.
!>
!> line 11: TSTERR, LOGICAL
!>          Flag indicating whether or not to test the error exits for
!>          the LAPACK routines and driver routines.
!>
!> line 12: NEWSD, INTEGER
!>          A code indicating how to set the random number seed.
!>          = 0:  Set the seed to a default value before each run
!>          = 1:  Initialize the seed to a default value only before the
!>                first run
!>          = 2:  Like 1, but use the seed values on the next line
!>
!> If line 12 was 2:
!>
!> line 13: INTEGER array, dimension (4)
!>          Four integer values for the random number seed.
!>
!> lines 13-EOF:  Lines specifying matrix types, as for NEP.
!>          The valid 3-character path names are 'SEP' or 'CST' for the
!>          Hermitian eigenvalue routines and driver routines, and
!>          'CSG' for the routines for the Hermitian generalized
!>          eigenvalue problem.
!>
!>-----------------------------------------------------------------------
!>
!> SVD input file:
!>
!> line 2:  NN, INTEGER
!>          Number of values of M and N.
!>
!> line 3:  MVAL, INTEGER array, dimension (NN)
!>          The values for the matrix row dimension M.
!>
!> line 4:  NVAL, INTEGER array, dimension (NN)
!>          The values for the matrix column dimension N.
!>
!> line 5:  NPARMS, INTEGER
!>          Number of values of the parameter NB, NBMIN, NX, and NRHS.
!>
!> line 6:  NBVAL, INTEGER array, dimension (NPARMS)
!>          The values for the blocksize NB.
!>
!> line 7:  NBMIN, INTEGER array, dimension (NPARMS)
!>          The values for the minimum blocksize NBMIN.
!>
!> line 8:  NXVAL, INTEGER array, dimension (NPARMS)
!>          The values for the crossover point NX.
!>
!> line 9:  NSVAL, INTEGER array, dimension (NPARMS)
!>          The values for the number of right hand sides NRHS.
!>
!> line 10: THRESH
!>          Threshold value for the test ratios.  Information will be
!>          printed about each test for which the test ratio is greater
!>          than or equal to the threshold.
!>
!> line 11: TSTCHK, LOGICAL
!>          Flag indicating whether or not to test the LAPACK routines.
!>
!> line 12: TSTDRV, LOGICAL
!>          Flag indicating whether or not to test the driver routines.
!>
!> line 13: TSTERR, LOGICAL
!>          Flag indicating whether or not to test the error exits for
!>          the LAPACK routines and driver routines.
!>
!> line 14: NEWSD, INTEGER
!>          A code indicating how to set the random number seed.
!>          = 0:  Set the seed to a default value before each run
!>          = 1:  Initialize the seed to a default value only before the
!>                first run
!>          = 2:  Like 1, but use the seed values on the next line
!>
!> If line 14 was 2:
!>
!> line 15: INTEGER array, dimension (4)
!>          Four integer values for the random number seed.
!>
!> lines 15-EOF:  Lines specifying matrix types, as for NEP.
!>          The 3-character path names are 'SVD' or 'CBD' for both the
!>          SVD routines and the SVD driver routines.
!>
!>-----------------------------------------------------------------------
!>
!> CEV and CES data files:
!>
!> line 1:  'CEV' or 'CES' in columns 1 to 3.
!>
!> line 2:  NSIZES, INTEGER
!>          Number of sizes of matrices to use. Should be at least 0
!>          and at most 20. If NSIZES = 0, no testing is done
!>          (although the remaining  3 lines are still read).
!>
!> line 3:  NN, INTEGER array, dimension(NSIZES)
!>          Dimensions of matrices to be tested.
!>
!> line 4:  NB, NBMIN, NX, NS, NBCOL, INTEGERs
!>          These integer parameters determine how blocking is done
!>          (see ILAENV for details)
!>          NB     : block size
!>          NBMIN  : minimum block size
!>          NX     : minimum dimension for blocking
!>          NS     : number of shifts in xHSEQR
!>          NBCOL  : minimum column dimension for blocking
!>
!> line 5:  THRESH, REAL
!>          The test threshold against which computed residuals are
!>          compared. Should generally be in the range from 10. to 20.
!>          If it is 0., all test case data will be printed.
!>
!> line 6:  NEWSD, INTEGER
!>          A code indicating how to set the random number seed.
!>          = 0:  Set the seed to a default value before each run
!>          = 1:  Initialize the seed to a default value only before the
!>                first run
!>          = 2:  Like 1, but use the seed values on the next line
!>
!> If line 6 was 2:
!>
!> line 7:  INTEGER array, dimension (4)
!>          Four integer values for the random number seed.
!>
!> lines 8 and following:  Lines specifying matrix types, as for NEP.
!>          The 3-character path name is 'CEV' to test CGEEV, or
!>          'CES' to test CGEES.
!>
!>-----------------------------------------------------------------------
!>
!> The CVX data has two parts. The first part is identical to CEV,
!> and the second part consists of test matrices with precomputed
!> solutions.
!>
!> line 1:  'CVX' in columns 1-3.
!>
!> line 2:  NSIZES, INTEGER
!>          If NSIZES = 0, no testing of randomly generated examples
!>          is done, but any precomputed examples are tested.
!>
!> line 3:  NN, INTEGER array, dimension(NSIZES)
!>
!> line 4:  NB, NBMIN, NX, NS, NBCOL, INTEGERs
!>
!> line 5:  THRESH, REAL
!>
!> line 6:  NEWSD, INTEGER
!>
!> If line 6 was 2:
!>
!> line 7:  INTEGER array, dimension (4)
!>
!> lines 8 and following: The first line contains 'CVX' in columns 1-3
!>          followed by the number of matrix types, possibly with
!>          a second line to specify certain matrix types.
!>          If the number of matrix types = 0, no testing of randomly
!>          generated examples is done, but any precomputed examples
!>          are tested.
!>
!> remaining lines : Each matrix is stored on 1+N+N**2 lines, where N is
!>          its dimension. The first line contains the dimension N and
!>          ISRT (two integers). ISRT indicates whether the last N lines
!>          are sorted by increasing real part of the eigenvalue
!>          (ISRT=0) or by increasing imaginary part (ISRT=1). The next
!>          N**2 lines contain the matrix rowwise, one entry per line.
!>          The last N lines correspond to each eigenvalue. Each of
!>          these last N lines contains 4 real values: the real part of
!>          the eigenvalues, the imaginary part of the eigenvalue, the
!>          reciprocal condition number of the eigenvalues, and the
!>          reciprocal condition number of the vector eigenvector. The
!>          end of data is indicated by dimension N=0. Even if no data
!>          is to be tested, there must be at least one line containing
!>          N=0.
!>
!>-----------------------------------------------------------------------
!>
!> The CSX data is like CVX. The first part is identical to CEV, and the
!> second part consists of test matrices with precomputed solutions.
!>
!> line 1:  'CSX' in columns 1-3.
!>
!> line 2:  NSIZES, INTEGER
!>          If NSIZES = 0, no testing of randomly generated examples
!>          is done, but any precomputed examples are tested.
!>
!> line 3:  NN, INTEGER array, dimension(NSIZES)
!>
!> line 4:  NB, NBMIN, NX, NS, NBCOL, INTEGERs
!>
!> line 5:  THRESH, REAL
!>
!> line 6:  NEWSD, INTEGER
!>
!> If line 6 was 2:
!>
!> line 7:  INTEGER array, dimension (4)
!>
!> lines 8 and following: The first line contains 'CSX' in columns 1-3
!>          followed by the number of matrix types, possibly with
!>          a second line to specify certain matrix types.
!>          If the number of matrix types = 0, no testing of randomly
!>          generated examples is done, but any precomputed examples
!>          are tested.
!>
!> remaining lines : Each matrix is stored on 3+N**2 lines, where N is
!>          its dimension. The first line contains the dimension N, the
!>          dimension M of an invariant subspace, and ISRT. The second
!>          line contains M integers, identifying the eigenvalues in the
!>          invariant subspace (by their position in a list of
!>          eigenvalues ordered by increasing real part (if ISRT=0) or
!>          by increasing imaginary part (if ISRT=1)). The next N**2
!>          lines contain the matrix rowwise. The last line contains the
!>          reciprocal condition number for the average of the selected
!>          eigenvalues, and the reciprocal condition number for the
!>          corresponding right invariant subspace. The end of data in
!>          indicated by a line containing N=0, M=0, and ISRT = 0.  Even
!>          if no data is to be tested, there must be at least one line
!>          containing N=0, M=0 and ISRT=0.
!>
!>-----------------------------------------------------------------------
!>
!> CGG input file:
!>
!> line 2:  NN, INTEGER
!>          Number of values of N.
!>
!> line 3:  NVAL, INTEGER array, dimension (NN)
!>          The values for the matrix dimension N.
!>
!> line 4:  NPARMS, INTEGER
!>          Number of values of the parameters NB, NBMIN, NBCOL, NS, and
!>          MAXB.
!>
!> line 5:  NBVAL, INTEGER array, dimension (NPARMS)
!>          The values for the blocksize NB.
!>
!> line 6:  NBMIN, INTEGER array, dimension (NPARMS)
!>          The values for NBMIN, the minimum row dimension for blocks.
!>
!> line 7:  NSVAL, INTEGER array, dimension (NPARMS)
!>          The values for the number of shifts.
!>
!> line 8:  MXBVAL, INTEGER array, dimension (NPARMS)
!>          The values for MAXB, used in determining minimum blocksize.
!>
!> line 9:  IACC22, INTEGER array, dimension (NPARMS)
!>          select structured matrix multiply: 1 or 2)
!>
!> line 10: NBCOL, INTEGER array, dimension (NPARMS)
!>          The values for NBCOL, the minimum column dimension for
!>          blocks.
!>
!> line 11: THRESH
!>          Threshold value for the test ratios.  Information will be
!>          printed about each test for which the test ratio is greater
!>          than or equal to the threshold.
!>
!> line 12: TSTCHK, LOGICAL
!>          Flag indicating whether or not to test the LAPACK routines.
!>
!> line 13: TSTDRV, LOGICAL
!>          Flag indicating whether or not to test the driver routines.
!>
!> line 14: TSTERR, LOGICAL
!>          Flag indicating whether or not to test the error exits for
!>          the LAPACK routines and driver routines.
!>
!> line 15: NEWSD, INTEGER
!>          A code indicating how to set the random number seed.
!>          = 0:  Set the seed to a default value before each run
!>          = 1:  Initialize the seed to a default value only before the
!>                first run
!>          = 2:  Like 1, but use the seed values on the next line
!>
!> If line 15 was 2:
!>
!> line 16: INTEGER array, dimension (4)
!>          Four integer values for the random number seed.
!>
!> lines 17-EOF:  Lines specifying matrix types, as for NEP.
!>          The 3-character path name is 'CGG' for the generalized
!>          eigenvalue problem routines and driver routines.
!>
!>-----------------------------------------------------------------------
!>
!> CGS and CGV input files:
!>
!> line 1:  'CGS' or 'CGV' in columns 1 to 3.
!>
!> line 2:  NN, INTEGER
!>          Number of values of N.
!>
!> line 3:  NVAL, INTEGER array, dimension(NN)
!>          Dimensions of matrices to be tested.
!>
!> line 4:  NB, NBMIN, NX, NS, NBCOL, INTEGERs
!>          These integer parameters determine how blocking is done
!>          (see ILAENV for details)
!>          NB     : block size
!>          NBMIN  : minimum block size
!>          NX     : minimum dimension for blocking
!>          NS     : number of shifts in xHGEQR
!>          NBCOL  : minimum column dimension for blocking
!>
!> line 5:  THRESH, REAL
!>          The test threshold against which computed residuals are
!>          compared. Should generally be in the range from 10. to 20.
!>          If it is 0., all test case data will be printed.
!>
!> line 6:  TSTERR, LOGICAL
!>          Flag indicating whether or not to test the error exits.
!>
!> line 7:  NEWSD, INTEGER
!>          A code indicating how to set the random number seed.
!>          = 0:  Set the seed to a default value before each run
!>          = 1:  Initialize the seed to a default value only before the
!>                first run
!>          = 2:  Like 1, but use the seed values on the next line
!>
!> If line 17 was 2:
!>
!> line 7:  INTEGER array, dimension (4)
!>          Four integer values for the random number seed.
!>
!> lines 7-EOF:  Lines specifying matrix types, as for NEP.
!>          The 3-character path name is 'CGS' for the generalized
!>          eigenvalue problem routines and driver routines.
!>
!>-----------------------------------------------------------------------
!>
!> CGX input file:
!> line 1:  'CGX' in columns 1 to 3.
!>
!> line 2:  N, INTEGER
!>          Value of N.
!>
!> line 3:  NB, NBMIN, NX, NS, NBCOL, INTEGERs
!>          These integer parameters determine how blocking is done
!>          (see ILAENV for details)
!>          NB     : block size
!>          NBMIN  : minimum block size
!>          NX     : minimum dimension for blocking
!>          NS     : number of shifts in xHGEQR
!>          NBCOL  : minimum column dimension for blocking
!>
!> line 4:  THRESH, REAL
!>          The test threshold against which computed residuals are
!>          compared. Should generally be in the range from 10. to 20.
!>          Information will be printed about each test for which the
!>          test ratio is greater than or equal to the threshold.
!>
!> line 5:  TSTERR, LOGICAL
!>          Flag indicating whether or not to test the error exits for
!>          the LAPACK routines and driver routines.
!>
!> line 6:  NEWSD, INTEGER
!>          A code indicating how to set the random number seed.
!>          = 0:  Set the seed to a default value before each run
!>          = 1:  Initialize the seed to a default value only before the
!>                first run
!>          = 2:  Like 1, but use the seed values on the next line
!>
!> If line 6 was 2:
!>
!> line 7: INTEGER array, dimension (4)
!>          Four integer values for the random number seed.
!>
!> If line 2 was 0:
!>
!> line 7-EOF: Precomputed examples are tested.
!>
!> remaining lines : Each example is stored on 3+2*N*N lines, where N is
!>          its dimension. The first line contains the dimension (a
!>          single integer).  The next line contains an integer k such
!>          that only the last k eigenvalues will be selected and appear
!>          in the leading diagonal blocks of $A$ and $B$. The next N*N
!>          lines contain the matrix A, one element per line. The next N*N
!>          lines contain the matrix B. The last line contains the
!>          reciprocal of the eigenvalue cluster condition number and the
!>          reciprocal of the deflating subspace (associated with the
!>          selected eigencluster) condition number.  The end of data is
!>          indicated by dimension N=0.  Even if no data is to be tested,
!>          there must be at least one line containing N=0.
!>
!>-----------------------------------------------------------------------
!>
!> CXV input files:
!> line 1:  'CXV' in columns 1 to 3.
!>
!> line 2:  N, INTEGER
!>          Value of N.
!>
!> line 3:  NB, NBMIN, NX, NS, NBCOL, INTEGERs
!>          These integer parameters determine how blocking is done
!>          (see ILAENV for details)
!>          NB     : block size
!>          NBMIN  : minimum block size
!>          NX     : minimum dimension for blocking
!>          NS     : number of shifts in xHGEQR
!>          NBCOL  : minimum column dimension for blocking
!>
!> line 4:  THRESH, REAL
!>          The test threshold against which computed residuals are
!>          compared. Should generally be in the range from 10. to 20.
!>          Information will be printed about each test for which the
!>          test ratio is greater than or equal to the threshold.
!>
!> line 5:  TSTERR, LOGICAL
!>          Flag indicating whether or not to test the error exits for
!>          the LAPACK routines and driver routines.
!>
!> line 6:  NEWSD, INTEGER
!>          A code indicating how to set the random number seed.
!>          = 0:  Set the seed to a default value before each run
!>          = 1:  Initialize the seed to a default value only before the
!>                first run
!>          = 2:  Like 1, but use the seed values on the next line
!>
!> If line 6 was 2:
!>
!> line 7: INTEGER array, dimension (4)
!>          Four integer values for the random number seed.
!>
!> If line 2 was 0:
!>
!> line 7-EOF: Precomputed examples are tested.
!>
!> remaining lines : Each example is stored on 3+2*N*N lines, where N is
!>          its dimension. The first line contains the dimension (a
!>          single integer). The next N*N lines contain the matrix A, one
!>          element per line. The next N*N lines contain the matrix B.
!>          The next line contains the reciprocals of the eigenvalue
!>          condition numbers.  The last line contains the reciprocals of
!>          the eigenvector condition numbers.  The end of data is
!>          indicated by dimension N=0.  Even if no data is to be tested,
!>          there must be at least one line containing N=0.
!>
!>-----------------------------------------------------------------------
!>
!> CHB input file:
!>
!> line 2:  NN, INTEGER
!>          Number of values of N.
!>
!> line 3:  NVAL, INTEGER array, dimension (NN)
!>          The values for the matrix dimension N.
!>
!> line 4:  NK, INTEGER
!>          Number of values of K.
!>
!> line 5:  KVAL, INTEGER array, dimension (NK)
!>          The values for the matrix dimension K.
!>
!> line 6:  THRESH
!>          Threshold value for the test ratios.  Information will be
!>          printed about each test for which the test ratio is greater
!>          than or equal to the threshold.
!>
!> line 7:  NEWSD, INTEGER
!>          A code indicating how to set the random number seed.
!>          = 0:  Set the seed to a default value before each run
!>          = 1:  Initialize the seed to a default value only before the
!>                first run
!>          = 2:  Like 1, but use the seed values on the next line
!>
!> If line 7 was 2:
!>
!> line 8:  INTEGER array, dimension (4)
!>          Four integer values for the random number seed.
!>
!> lines 8-EOF:  Lines specifying matrix types, as for NEP.
!>          The 3-character path name is 'CHB'.
!>
!>-----------------------------------------------------------------------
!>
!> CBB input file:
!>
!> line 2:  NN, INTEGER
!>          Number of values of M and N.
!>
!> line 3:  MVAL, INTEGER array, dimension (NN)
!>          The values for the matrix row dimension M.
!>
!> line 4:  NVAL, INTEGER array, dimension (NN)
!>          The values for the matrix column dimension N.
!>
!> line 4:  NK, INTEGER
!>          Number of values of K.
!>
!> line 5:  KVAL, INTEGER array, dimension (NK)
!>          The values for the matrix bandwidth K.
!>
!> line 6:  NPARMS, INTEGER
!>          Number of values of the parameter NRHS
!>
!> line 7:  NSVAL, INTEGER array, dimension (NPARMS)
!>          The values for the number of right hand sides NRHS.
!>
!> line 8:  THRESH
!>          Threshold value for the test ratios.  Information will be
!>          printed about each test for which the test ratio is greater
!>          than or equal to the threshold.
!>
!> line 9:  NEWSD, INTEGER
!>          A code indicating how to set the random number seed.
!>          = 0:  Set the seed to a default value before each run
!>          = 1:  Initialize the seed to a default value only before the
!>                first run
!>          = 2:  Like 1, but use the seed values on the next line
!>
!> If line 9 was 2:
!>
!> line 10: INTEGER array, dimension (4)
!>          Four integer values for the random number seed.
!>
!> lines 10-EOF:  Lines specifying matrix types, as for SVD.
!>          The 3-character path name is 'CBB'.
!>
!>-----------------------------------------------------------------------
!>
!> CEC input file:
!>
!> line  2: THRESH, REAL
!>          Threshold value for the test ratios.  Information will be
!>          printed about each test for which the test ratio is greater
!>          than or equal to the threshold.
!>
!> lines  3-EOF:
!>
!> Input for testing the eigencondition routines consists of a set of
!> specially constructed test cases and their solutions.  The data
!> format is not intended to be modified by the user.
!>
!>-----------------------------------------------------------------------
!>
!> CBL and CBK input files:
!>
!> line 1:  'CBL' in columns 1-3 to test CGEBAL, or 'CBK' in
!>          columns 1-3 to test CGEBAK.
!>
!> The remaining lines consist of specially constructed test cases.
!>
!>-----------------------------------------------------------------------
!>
!> CGL and CGK input files:
!>
!> line 1:  'CGL' in columns 1-3 to test CGGBAL, or 'CGK' in
!>          columns 1-3 to test CGGBAK.
!>
!> The remaining lines consist of specially constructed test cases.
!>
!>-----------------------------------------------------------------------
!>
!> GLM data file:
!>
!> line 1:  'GLM' in columns 1 to 3.
!>
!> line 2:  NN, INTEGER
!>          Number of values of M, P, and N.
!>
!> line 3:  MVAL, INTEGER array, dimension(NN)
!>          Values of M (row dimension).
!>
!> line 4:  PVAL, INTEGER array, dimension(NN)
!>          Values of P (row dimension).
!>
!> line 5:  NVAL, INTEGER array, dimension(NN)
!>          Values of N (column dimension), note M <= N <= M+P.
!>
!> line 6:  THRESH, REAL
!>          Threshold value for the test ratios.  Information will be
!>          printed about each test for which the test ratio is greater
!>          than or equal to the threshold.
!>
!> line 7:  TSTERR, LOGICAL
!>          Flag indicating whether or not to test the error exits for
!>          the LAPACK routines and driver routines.
!>
!> line 8:  NEWSD, INTEGER
!>          A code indicating how to set the random number seed.
!>          = 0:  Set the seed to a default value before each run
!>          = 1:  Initialize the seed to a default value only before the
!>                first run
!>          = 2:  Like 1, but use the seed values on the next line
!>
!> If line 8 was 2:
!>
!> line 9:  INTEGER array, dimension (4)
!>          Four integer values for the random number seed.
!>
!> lines 9-EOF:  Lines specifying matrix types, as for NEP.
!>          The 3-character path name is 'GLM' for the generalized
!>          linear regression model routines.
!>
!>-----------------------------------------------------------------------
!>
!> GQR data file:
!>
!> line 1:  'GQR' in columns 1 to 3.
!>
!> line 2:  NN, INTEGER
!>          Number of values of M, P, and N.
!>
!> line 3:  MVAL, INTEGER array, dimension(NN)
!>          Values of M.
!>
!> line 4:  PVAL, INTEGER array, dimension(NN)
!>          Values of P.
!>
!> line 5:  NVAL, INTEGER array, dimension(NN)
!>          Values of N.
!>
!> line 6:  THRESH, REAL
!>          Threshold value for the test ratios.  Information will be
!>          printed about each test for which the test ratio is greater
!>          than or equal to the threshold.
!>
!> line 7:  TSTERR, LOGICAL
!>          Flag indicating whether or not to test the error exits for
!>          the LAPACK routines and driver routines.
!>
!> line 8:  NEWSD, INTEGER
!>          A code indicating how to set the random number seed.
!>          = 0:  Set the seed to a default value before each run
!>          = 1:  Initialize the seed to a default value only before the
!>                first run
!>          = 2:  Like 1, but use the seed values on the next line
!>
!> If line 8 was 2:
!>
!> line 9:  INTEGER array, dimension (4)
!>          Four integer values for the random number seed.
!>
!> lines 9-EOF:  Lines specifying matrix types, as for NEP.
!>          The 3-character path name is 'GQR' for the generalized
!>          QR and RQ routines.
!>
!>-----------------------------------------------------------------------
!>
!> GSV data file:
!>
!> line 1:  'GSV' in columns 1 to 3.
!>
!> line 2:  NN, INTEGER
!>          Number of values of M, P, and N.
!>
!> line 3:  MVAL, INTEGER array, dimension(NN)
!>          Values of M (row dimension).
!>
!> line 4:  PVAL, INTEGER array, dimension(NN)
!>          Values of P (row dimension).
!>
!> line 5:  NVAL, INTEGER array, dimension(NN)
!>          Values of N (column dimension).
!>
!> line 6:  THRESH, REAL
!>          Threshold value for the test ratios.  Information will be
!>          printed about each test for which the test ratio is greater
!>          than or equal to the threshold.
!>
!> line 7:  TSTERR, LOGICAL
!>          Flag indicating whether or not to test the error exits for
!>          the LAPACK routines and driver routines.
!>
!> line 8:  NEWSD, INTEGER
!>          A code indicating how to set the random number seed.
!>          = 0:  Set the seed to a default value before each run
!>          = 1:  Initialize the seed to a default value only before the
!>                first run
!>          = 2:  Like 1, but use the seed values on the next line
!>
!> If line 8 was 2:
!>
!> line 9:  INTEGER array, dimension (4)
!>          Four integer values for the random number seed.
!>
!> lines 9-EOF:  Lines specifying matrix types, as for NEP.
!>          The 3-character path name is 'GSV' for the generalized
!>          SVD routines.
!>
!>-----------------------------------------------------------------------
!>
!> CSD data file:
!>
!> line 1:  'CSD' in columns 1 to 3.
!>
!> line 2:  NM, INTEGER
!>          Number of values of M, P, and N.
!>
!> line 3:  MVAL, INTEGER array, dimension(NM)
!>          Values of M (row and column dimension of orthogonal matrix).
!>
!> line 4:  PVAL, INTEGER array, dimension(NM)
!>          Values of P (row dimension of top-left block).
!>
!> line 5:  NVAL, INTEGER array, dimension(NM)
!>          Values of N (column dimension of top-left block).
!>
!> line 6:  THRESH, REAL
!>          Threshold value for the test ratios.  Information will be
!>          printed about each test for which the test ratio is greater
!>          than or equal to the threshold.
!>
!> line 7:  TSTERR, LOGICAL
!>          Flag indicating whether or not to test the error exits for
!>          the LAPACK routines and driver routines.
!>
!> line 8:  NEWSD, INTEGER
!>          A code indicating how to set the random number seed.
!>          = 0:  Set the seed to a default value before each run
!>          = 1:  Initialize the seed to a default value only before the
!>                first run
!>          = 2:  Like 1, but use the seed values on the next line
!>
!> If line 8 was 2:
!>
!> line 9:  INTEGER array, dimension (4)
!>          Four integer values for the random number seed.
!>
!> lines 9-EOF:  Lines specifying matrix types, as for NEP.
!>          The 3-character path name is 'CSD' for the CSD routine.
!>
!>-----------------------------------------------------------------------
!>
!> LSE data file:
!>
!> line 1:  'LSE' in columns 1 to 3.
!>
!> line 2:  NN, INTEGER
!>          Number of values of M, P, and N.
!>
!> line 3:  MVAL, INTEGER array, dimension(NN)
!>          Values of M.
!>
!> line 4:  PVAL, INTEGER array, dimension(NN)
!>          Values of P.
!>
!> line 5:  NVAL, INTEGER array, dimension(NN)
!>          Values of N, note P <= N <= P+M.
!>
!> line 6:  THRESH, REAL
!>          Threshold value for the test ratios.  Information will be
!>          printed about each test for which the test ratio is greater
!>          than or equal to the threshold.
!>
!> line 7:  TSTERR, LOGICAL
!>          Flag indicating whether or not to test the error exits for
!>          the LAPACK routines and driver routines.
!>
!> line 8:  NEWSD, INTEGER
!>          A code indicating how to set the random number seed.
!>          = 0:  Set the seed to a default value before each run
!>          = 1:  Initialize the seed to a default value only before the
!>                first run
!>          = 2:  Like 1, but use the seed values on the next line
!>
!> If line 8 was 2:
!>
!> line 9:  INTEGER array, dimension (4)
!>          Four integer values for the random number seed.
!>
!> lines 9-EOF:  Lines specifying matrix types, as for NEP.
!>          The 3-character path name is 'GSV' for the generalized
!>          SVD routines.
!>
!>-----------------------------------------------------------------------
!>
!> NMAX is currently set to 132 and must be at least 12 for some of the
!> precomputed examples, and LWORK = NMAX*(5*NMAX+20) in the parameter
!> statements below.  For SVD, we assume NRHS may be as big as N.  The
!> parameter NEED is set to 14 to allow for 14 N-by-N matrices for CGG.
!> \endverbatim
!
!  Arguments:
!  ==========
!
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
!> \ingroup complex_eig
!
!  =====================================================================
PROGRAM CCHKEE
   USE M_TST_EIG , ONLY:ALAREQ  , CCHKBB     , CCHKBD     , CCHKBK
   USE M_TST_EIG , ONLY:CCHKBL  , CCHKEC     , CCHKGG     , CCHKGK
   USE M_TST_EIG , ONLY:CCHKGL  , CCHKHB     , CCHKHB2STG , CCHKHS
   USE M_TST_EIG , ONLY:CCHKST  , CCHKST2STG , CCKCSD     , CCKGLM
   USE M_TST_EIG , ONLY:CCKGQR  , CCKGSV     , CCKLSE     , CDRGES
   USE M_TST_EIG , ONLY:CDRGES3 , CDRGEV     , CDRGEV3    , CDRGSX
   USE M_TST_EIG , ONLY:CDRGVX  , CDRVBD     , CDRVES     , CDRVEV
   USE M_TST_EIG , ONLY:CDRVSG  , CDRVST     , CDRVST2STG , CDRVSX
   USE M_TST_EIG , ONLY:CDRVVX  , CERRBD     , CERRED     , CERRGG
   USE M_TST_EIG , ONLY:CERRHS  , CERRST     , XLAENV
   USE M_TST_EIG , ONLY:CDRVSG2STG
#if defined(_OPENMP)
   USE OMP_LIB
#endif
      IMPLICIT NONE
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2016
!  =====================================================================
!     .. Parameters ..
      INTEGER,PARAMETER :: NMAX=132
      INTEGER,PARAMETER :: NCMAX=20
      INTEGER,PARAMETER :: NEED=14
      INTEGER,PARAMETER :: LWORK=NMAX*(5*NMAX+20)
      INTEGER,PARAMETER :: LIWORK=NMAX*(NMAX+20)
      INTEGER,PARAMETER :: MAXIN=20
      INTEGER,PARAMETER :: MAXT=30
      INTEGER,PARAMETER :: NIN=5,NOUT=6
!     ..
!     .. Local Scalars ..
      LOGICAL cbb , cbk , cbl , ces , cev , cgg , cgk , cgl , cgs ,     &
     &        cgv , cgx , chb , csd , csx , cvx , cxv , fatal , glm ,   &
     &        gqr , gsv , lse , nep , sep , svd , tstchk , tstdif ,     &
     &        tstdrv , tsterr
      CHARACTER c1
      CHARACTER*3 c3 , path
      CHARACTER*32 vname
      CHARACTER*10 intstr
      CHARACTER*80 line
      INTEGER i , i1 , ic , info , itmp , k , lenp , maxtyp , newsd ,   &
     &        nk , nn , nparms , nrhs , ntypes , vers_major ,           &
     &        vers_minor , vers_patch , n_threads
      REAL eps , s1 , s2 , thresh , thrshn
!     ..
!     .. Local Arrays ..
      LOGICAL dotype(MAXT) , logwrk(NMAX)
      INTEGER ioldsd(4) , iseed(4) , iwork(LIWORK) , kval(MAXIN) ,      &
     &        mval(MAXIN) , mxbval(MAXIN) , nbcol(MAXIN) , nbmin(MAXIN) &
     &        , nbval(MAXIN) , nsval(MAXIN) , nval(MAXIN) , nxval(MAXIN)&
     &        , pval(MAXIN)
      INTEGER inmin(MAXIN) , inwin(MAXIN) , inibl(MAXIN) , ishfts(MAXIN)&
     &        , iacc22(MAXIN)
      REAL alpha(NMAX) , beta(NMAX) , dr(NMAX,12) , result(500)
      COMPLEX dc(NMAX,6) , taua(NMAX) , taub(NMAX) , x(5*NMAX)
!     ..
!     .. Allocatable Arrays ..
      INTEGER allocatestatus
      REAL , DIMENSION(:) , ALLOCATABLE  ::  rwork , s
      COMPLEX , DIMENSION(:) , ALLOCATABLE  ::  work
      COMPLEX , DIMENSION(:,:) , ALLOCATABLE  ::  a , b , c
!     ..
!     .. External Functions ..
      LOGICAL LSAMEN
      REAL SECOND , SLAMCH
      EXTERNAL LSAMEN , SECOND , SLAMCH
!     ..
!     .. External Subroutines ..
!      EXTERNAL           ALAREQ, CCHKBB, CCHKBD, CCHKBK, CCHKBL, CCHKEC,
!     $                   CCHKGG, CCHKGK, CCHKGL, CCHKHB, CCHKHS, CCHKST,
!     $                   CCKCSD, CCKGLM, CCKGQR, CCKGSV, CCKLSE, CDRGES,
!     $                   CDRGEV, CDRGSX, CDRGVX, CDRVBD, CDRVES, CDRVEV,
!     $                   CDRVSG, CDRVST, CDRVSX, CDRVVX, CERRBD,
!     $                   CERRED, CERRGG, CERRHS, CERRST, ILAVER, XLAENV,
!     $                   CDRGES3, CDRGEV3,
!     $                   CCHKST2STG, CDRVST2STG, CCHKHB2STG
      EXTERNAL ILAVER
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC LEN , MIN
!     ..
!     .. Scalars in Common ..
      LOGICAL LERr , OK
      CHARACTER*32 SRNamt
      INTEGER INFot , MAXb , NPRoc , NSHift , NUNit , SELdim , SELopt
!     ..
!     .. Arrays in Common ..
      LOGICAL SELval(20)
      INTEGER IPArms(100)
      REAL SELwi(20) , SELwr(20)
!     ..
!     .. Common blocks ..
      COMMON /CENVIR/ NPRoc , NSHift , MAXb
      COMMON /CLAENV/ IPArms
      COMMON /INFOC / INFot , NUNit , OK , LERr
      COMMON /SRNAMC/ SRNamt
      COMMON /SSLCT / SELopt , SELdim , SELval , SELwr , SELwi
!     ..
!     .. Data statements ..
      DATA intstr/'0123456789'/
      DATA ioldsd/0 , 0 , 0 , 1/
!     ..
!     .. Allocate memory dynamically ..
!
      ALLOCATE (s(NMAX*NMAX),STAT=allocatestatus)
      IF ( allocatestatus/=0 ) STOP "*** Not enough memory ***"
      ALLOCATE (a(NMAX*NMAX,NEED),STAT=allocatestatus)
      IF ( allocatestatus/=0 ) STOP "*** Not enough memory ***"
      ALLOCATE (b(NMAX*NMAX,5),STAT=allocatestatus)
      IF ( allocatestatus/=0 ) STOP "*** Not enough memory ***"
      ALLOCATE (c(NCMAX*NCMAX,NCMAX*NCMAX),STAT=allocatestatus)
      IF ( allocatestatus/=0 ) STOP "*** Not enough memory ***"
      ALLOCATE (rwork(LWORK),STAT=allocatestatus)
      IF ( allocatestatus/=0 ) STOP "*** Not enough memory ***"
      ALLOCATE (work(LWORK),STAT=allocatestatus)
      IF ( allocatestatus/=0 ) STOP "*** Not enough memory ***"
!     ..
!     .. Executable Statements ..
!
      a = 0.0
      b = 0.0
      c = 0.0
      dc = 0.0
      s1 = SECOND()
      fatal = .FALSE.
      NUNit = NOUT
 100  DO
!
!     Return to here to read multiple sets of data
!
!
!     Read the first line and set the 3-character test path
!
         READ (NIN,FMT='(A80)',END=200) line
         path = line(1:3)
         nep = LSAMEN(3,path,'NEP') .OR. LSAMEN(3,path,'CHS')
         sep = LSAMEN(3,path,'SEP') .OR. LSAMEN(3,path,'CST') .OR.      &
     &         LSAMEN(3,path,'CSG') .OR. LSAMEN(3,path,'SE2')
         svd = LSAMEN(3,path,'SVD') .OR. LSAMEN(3,path,'CBD')
         cev = LSAMEN(3,path,'CEV')
         ces = LSAMEN(3,path,'CES')
         cvx = LSAMEN(3,path,'CVX')
         csx = LSAMEN(3,path,'CSX')
         cgg = LSAMEN(3,path,'CGG')
         cgs = LSAMEN(3,path,'CGS')
         cgx = LSAMEN(3,path,'CGX')
         cgv = LSAMEN(3,path,'CGV')
         cxv = LSAMEN(3,path,'CXV')
         chb = LSAMEN(3,path,'CHB')
         cbb = LSAMEN(3,path,'CBB')
         glm = LSAMEN(3,path,'GLM')
         gqr = LSAMEN(3,path,'GQR') .OR. LSAMEN(3,path,'GRQ')
         gsv = LSAMEN(3,path,'GSV')
         csd = LSAMEN(3,path,'CSD')
         lse = LSAMEN(3,path,'LSE')
         cbl = LSAMEN(3,path,'CBL')
         cbk = LSAMEN(3,path,'CBK')
         cgl = LSAMEN(3,path,'CGL')
         cgk = LSAMEN(3,path,'CGK')
!
!     Report values of parameters.
!
         IF ( path=='   ' ) CYCLE
         IF ( nep ) THEN
            WRITE (NOUT,FMT=99012)
         ELSEIF ( sep ) THEN
            WRITE (NOUT,FMT=99013)
         ELSEIF ( svd ) THEN
            WRITE (NOUT,FMT=99014)
         ELSEIF ( cev ) THEN
            WRITE (NOUT,FMT=99020)
         ELSEIF ( ces ) THEN
            WRITE (NOUT,FMT=99021)
         ELSEIF ( cvx ) THEN
            WRITE (NOUT,FMT=99022)
         ELSEIF ( csx ) THEN
            WRITE (NOUT,FMT=99023)
         ELSEIF ( cgg ) THEN
            WRITE (NOUT,FMT=99024)
         ELSEIF ( cgs ) THEN
            WRITE (NOUT,FMT=99035)
         ELSEIF ( cgx ) THEN
            WRITE (NOUT,FMT=99034)
         ELSEIF ( cgv ) THEN
            WRITE (NOUT,FMT=99036)
         ELSEIF ( cxv ) THEN
            WRITE (NOUT,FMT=99037)
         ELSEIF ( chb ) THEN
            WRITE (NOUT,FMT=99025)
         ELSEIF ( cbb ) THEN
            WRITE (NOUT,FMT=99032)
         ELSEIF ( glm ) THEN
            WRITE (NOUT,FMT=99028)
         ELSEIF ( gqr ) THEN
            WRITE (NOUT,FMT=99029)
         ELSEIF ( gsv ) THEN
            WRITE (NOUT,FMT=99030)
         ELSEIF ( csd ) THEN
            WRITE (NOUT,FMT=99039)
         ELSEIF ( lse ) THEN
            WRITE (NOUT,FMT=99031)
         ELSEIF ( cbl ) THEN
!
!        CGEBAL:  Balancing
!
            CALL CCHKBL(NIN,NOUT)
            GOTO 200
         ELSEIF ( cbk ) THEN
!
!        CGEBAK:  Back transformation
!
            CALL CCHKBK(NIN,NOUT)
            GOTO 200
         ELSEIF ( cgl ) THEN
!
!        CGGBAL:  Balancing
!
            CALL CCHKGL(NIN,NOUT)
            GOTO 200
         ELSEIF ( cgk ) THEN
!
!        CGGBAK:  Back transformation
!
            CALL CCHKGK(NIN,NOUT)
            GOTO 200
         ELSEIF ( LSAMEN(3,path,'CEC') ) THEN
!
!        CEC:  Eigencondition estimation
!
            READ (NIN,FMT=*) thresh
            CALL XLAENV(1,1)
            CALL XLAENV(12,1)
            tsterr = .TRUE.
            CALL CCHKEC(thresh,tsterr,NIN,NOUT)
            GOTO 200
         ELSE
            WRITE (NOUT,FMT=99007) path
            GOTO 200
         ENDIF
         CALL ILAVER(vers_major,vers_minor,vers_patch)
         WRITE (NOUT,FMT=99027) vers_major , vers_minor , vers_patch
         WRITE (NOUT,FMT=99015)
!
!     Read the number of values of M, P, and N.
!
         READ (NIN,FMT=*) nn
         IF ( nn<0 ) THEN
            WRITE (NOUT,FMT=99010) '   NN ' , nn , 1
            nn = 0
            fatal = .TRUE.
         ELSEIF ( nn>MAXIN ) THEN
            WRITE (NOUT,FMT=99011) '   NN ' , nn , MAXIN
            nn = 0
            fatal = .TRUE.
         ENDIF
!
!     Read the values of M
!
         IF ( .NOT.(cgx .OR. cxv) ) THEN
            READ (NIN,FMT=*) (mval(i),i=1,nn)
            IF ( svd ) THEN
               vname = '    M '
            ELSE
               vname = '    N '
            ENDIF
            DO i = 1 , nn
               IF ( mval(i)<0 ) THEN
                  WRITE (NOUT,FMT=99010) vname , mval(i) , 0
                  fatal = .TRUE.
               ELSEIF ( mval(i)>NMAX ) THEN
                  WRITE (NOUT,FMT=99011) vname , mval(i) , NMAX
                  fatal = .TRUE.
               ENDIF
            ENDDO
            WRITE (NOUT,FMT=99016) 'M:    ' , (mval(i),i=1,nn)
         ENDIF
!
!     Read the values of P
!
         IF ( glm .OR. gqr .OR. gsv .OR. csd .OR. lse ) THEN
            READ (NIN,FMT=*) (pval(i),i=1,nn)
            DO i = 1 , nn
               IF ( pval(i)<0 ) THEN
                  WRITE (NOUT,FMT=99010) ' P  ' , pval(i) , 0
                  fatal = .TRUE.
               ELSEIF ( pval(i)>NMAX ) THEN
                  WRITE (NOUT,FMT=99011) ' P  ' , pval(i) , NMAX
                  fatal = .TRUE.
               ENDIF
            ENDDO
            WRITE (NOUT,FMT=99016) 'P:    ' , (pval(i),i=1,nn)
         ENDIF
!
!     Read the values of N
!
         IF ( svd .OR. cbb .OR. glm .OR. gqr .OR. gsv .OR. csd .OR.     &
     &        lse ) THEN
            READ (NIN,FMT=*) (nval(i),i=1,nn)
            DO i = 1 , nn
               IF ( nval(i)<0 ) THEN
                  WRITE (NOUT,FMT=99010) '    N ' , nval(i) , 0
                  fatal = .TRUE.
               ELSEIF ( nval(i)>NMAX ) THEN
                  WRITE (NOUT,FMT=99011) '    N ' , nval(i) , NMAX
                  fatal = .TRUE.
               ENDIF
            ENDDO
         ELSE
            DO i = 1 , nn
               nval(i) = mval(i)
            ENDDO
         ENDIF
         IF ( .NOT.(cgx .OR. cxv) ) THEN
            WRITE (NOUT,FMT=99016) 'N:    ' , (nval(i),i=1,nn)
         ELSE
            WRITE (NOUT,FMT=99016) 'N:    ' , nn
         ENDIF
!
!     Read the number of values of K, followed by the values of K
!
         IF ( chb .OR. cbb ) THEN
            READ (NIN,FMT=*) nk
            READ (NIN,FMT=*) (kval(i),i=1,nk)
            DO i = 1 , nk
               IF ( kval(i)<0 ) THEN
                  WRITE (NOUT,FMT=99010) '    K ' , kval(i) , 0
                  fatal = .TRUE.
               ELSEIF ( kval(i)>NMAX ) THEN
                  WRITE (NOUT,FMT=99011) '    K ' , kval(i) , NMAX
                  fatal = .TRUE.
               ENDIF
            ENDDO
            WRITE (NOUT,FMT=99016) 'K:    ' , (kval(i),i=1,nk)
         ENDIF
!
         IF ( cev .OR. ces .OR. cvx .OR. csx ) THEN
!
!        For the nonsymmetric QR driver routines, only one set of
!        parameters is allowed.
!
            READ (NIN,FMT=*) nbval(1) , nbmin(1) , nxval(1) , inmin(1) ,&
     &                       inwin(1) , inibl(1) , ishfts(1) , iacc22(1)
            IF ( nbval(1)<1 ) THEN
               WRITE (NOUT,FMT=99010) '   NB ' , nbval(1) , 1
               fatal = .TRUE.
            ELSEIF ( nbmin(1)<1 ) THEN
               WRITE (NOUT,FMT=99010) 'NBMIN ' , nbmin(1) , 1
               fatal = .TRUE.
            ELSEIF ( nxval(1)<1 ) THEN
               WRITE (NOUT,FMT=99010) '   NX ' , nxval(1) , 1
               fatal = .TRUE.
            ELSEIF ( inmin(1)<1 ) THEN
               WRITE (NOUT,FMT=99010) '   INMIN ' , inmin(1) , 1
               fatal = .TRUE.
            ELSEIF ( inwin(1)<1 ) THEN
               WRITE (NOUT,FMT=99010) '   INWIN ' , inwin(1) , 1
               fatal = .TRUE.
            ELSEIF ( inibl(1)<1 ) THEN
               WRITE (NOUT,FMT=99010) '   INIBL ' , inibl(1) , 1
               fatal = .TRUE.
            ELSEIF ( ishfts(1)<1 ) THEN
               WRITE (NOUT,FMT=99010) '   ISHFTS ' , ishfts(1) , 1
               fatal = .TRUE.
            ELSEIF ( iacc22(1)<0 ) THEN
               WRITE (NOUT,FMT=99010) '   IACC22 ' , iacc22(1) , 0
               fatal = .TRUE.
            ENDIF
            CALL XLAENV(1,nbval(1))
            CALL XLAENV(2,nbmin(1))
            CALL XLAENV(3,nxval(1))
            CALL XLAENV(12,MAX(11,inmin(1)))
            CALL XLAENV(13,inwin(1))
            CALL XLAENV(14,inibl(1))
            CALL XLAENV(15,ishfts(1))
            CALL XLAENV(16,iacc22(1))
            WRITE (NOUT,FMT=99016) 'NB:   ' , nbval(1)
            WRITE (NOUT,FMT=99016) 'NBMIN:' , nbmin(1)
            WRITE (NOUT,FMT=99016) 'NX:   ' , nxval(1)
            WRITE (NOUT,FMT=99016) 'INMIN:   ' , inmin(1)
            WRITE (NOUT,FMT=99016) 'INWIN: ' , inwin(1)
            WRITE (NOUT,FMT=99016) 'INIBL: ' , inibl(1)
            WRITE (NOUT,FMT=99016) 'ISHFTS: ' , ishfts(1)
            WRITE (NOUT,FMT=99016) 'IACC22: ' , iacc22(1)
!
         ELSEIF ( cgs .OR. cgx .OR. cgv .OR. cxv ) THEN
!
!        For the nonsymmetric generalized driver routines, only one set of
!        parameters is allowed.
!
            READ (NIN,FMT=*) nbval(1) , nbmin(1) , nxval(1) , nsval(1) ,&
     &                       mxbval(1)
            IF ( nbval(1)<1 ) THEN
               WRITE (NOUT,FMT=99010) '   NB ' , nbval(1) , 1
               fatal = .TRUE.
            ELSEIF ( nbmin(1)<1 ) THEN
               WRITE (NOUT,FMT=99010) 'NBMIN ' , nbmin(1) , 1
               fatal = .TRUE.
            ELSEIF ( nxval(1)<1 ) THEN
               WRITE (NOUT,FMT=99010) '   NX ' , nxval(1) , 1
               fatal = .TRUE.
            ELSEIF ( nsval(1)<2 ) THEN
               WRITE (NOUT,FMT=99010) '   NS ' , nsval(1) , 2
               fatal = .TRUE.
            ELSEIF ( mxbval(1)<1 ) THEN
               WRITE (NOUT,FMT=99010) ' MAXB ' , mxbval(1) , 1
               fatal = .TRUE.
            ENDIF
            CALL XLAENV(1,nbval(1))
            CALL XLAENV(2,nbmin(1))
            CALL XLAENV(3,nxval(1))
            CALL XLAENV(4,nsval(1))
            CALL XLAENV(8,mxbval(1))
            WRITE (NOUT,FMT=99016) 'NB:   ' , nbval(1)
            WRITE (NOUT,FMT=99016) 'NBMIN:' , nbmin(1)
            WRITE (NOUT,FMT=99016) 'NX:   ' , nxval(1)
            WRITE (NOUT,FMT=99016) 'NS:   ' , nsval(1)
            WRITE (NOUT,FMT=99016) 'MAXB: ' , mxbval(1)
         ELSEIF ( .NOT.chb .AND. .NOT.glm .AND. .NOT.gqr .AND.          &
     &            .NOT.gsv .AND. .NOT.csd .AND. .NOT.lse ) THEN
!
!        For the other paths, the number of parameters can be varied
!        from the input file.  Read the number of parameter values.
!
            READ (NIN,FMT=*) nparms
            IF ( nparms<1 ) THEN
               WRITE (NOUT,FMT=99010) 'NPARMS' , nparms , 1
               nparms = 0
               fatal = .TRUE.
            ELSEIF ( nparms>MAXIN ) THEN
               WRITE (NOUT,FMT=99011) 'NPARMS' , nparms , MAXIN
               nparms = 0
               fatal = .TRUE.
            ENDIF
!
!        Read the values of NB
!
            IF ( .NOT.cbb ) THEN
               READ (NIN,FMT=*) (nbval(i),i=1,nparms)
               DO i = 1 , nparms
                  IF ( nbval(i)<0 ) THEN
                     WRITE (NOUT,FMT=99010) '   NB ' , nbval(i) , 0
                     fatal = .TRUE.
                  ELSEIF ( nbval(i)>NMAX ) THEN
                     WRITE (NOUT,FMT=99011) '   NB ' , nbval(i) , NMAX
                     fatal = .TRUE.
                  ENDIF
               ENDDO
               WRITE (NOUT,FMT=99016) 'NB:   ' , (nbval(i),i=1,nparms)
            ENDIF
!
!        Read the values of NBMIN
!
            IF ( nep .OR. sep .OR. svd .OR. cgg ) THEN
               READ (NIN,FMT=*) (nbmin(i),i=1,nparms)
               DO i = 1 , nparms
                  IF ( nbmin(i)<0 ) THEN
                     WRITE (NOUT,FMT=99010) 'NBMIN ' , nbmin(i) , 0
                     fatal = .TRUE.
                  ELSEIF ( nbmin(i)>NMAX ) THEN
                     WRITE (NOUT,FMT=99011) 'NBMIN ' , nbmin(i) , NMAX
                     fatal = .TRUE.
                  ENDIF
               ENDDO
               WRITE (NOUT,FMT=99016) 'NBMIN:' , (nbmin(i),i=1,nparms)
            ELSE
               DO i = 1 , nparms
                  nbmin(i) = 1
               ENDDO
            ENDIF
!
!        Read the values of NX
!
            IF ( nep .OR. sep .OR. svd ) THEN
               READ (NIN,FMT=*) (nxval(i),i=1,nparms)
               DO i = 1 , nparms
                  IF ( nxval(i)<0 ) THEN
                     WRITE (NOUT,FMT=99010) '   NX ' , nxval(i) , 0
                     fatal = .TRUE.
                  ELSEIF ( nxval(i)>NMAX ) THEN
                     WRITE (NOUT,FMT=99011) '   NX ' , nxval(i) , NMAX
                     fatal = .TRUE.
                  ENDIF
               ENDDO
               WRITE (NOUT,FMT=99016) 'NX:   ' , (nxval(i),i=1,nparms)
            ELSE
               DO i = 1 , nparms
                  nxval(i) = 1
               ENDDO
            ENDIF
!
!        Read the values of NSHIFT (if CGG) or NRHS (if SVD
!        or CBB).
!
            IF ( svd .OR. cbb .OR. cgg ) THEN
               READ (NIN,FMT=*) (nsval(i),i=1,nparms)
               DO i = 1 , nparms
                  IF ( nsval(i)<0 ) THEN
                     WRITE (NOUT,FMT=99010) '   NS ' , nsval(i) , 0
                     fatal = .TRUE.
                  ELSEIF ( nsval(i)>NMAX ) THEN
                     WRITE (NOUT,FMT=99011) '   NS ' , nsval(i) , NMAX
                     fatal = .TRUE.
                  ENDIF
               ENDDO
               WRITE (NOUT,FMT=99016) 'NS:   ' , (nsval(i),i=1,nparms)
            ELSE
               DO i = 1 , nparms
                  nsval(i) = 1
               ENDDO
            ENDIF
!
!        Read the values for MAXB.
!
            IF ( cgg ) THEN
               READ (NIN,FMT=*) (mxbval(i),i=1,nparms)
               DO i = 1 , nparms
                  IF ( mxbval(i)<0 ) THEN
                     WRITE (NOUT,FMT=99010) ' MAXB ' , mxbval(i) , 0
                     fatal = .TRUE.
                  ELSEIF ( mxbval(i)>NMAX ) THEN
                     WRITE (NOUT,FMT=99011) ' MAXB ' , mxbval(i) , NMAX
                     fatal = .TRUE.
                  ENDIF
               ENDDO
               WRITE (NOUT,FMT=99016) 'MAXB: ' , (mxbval(i),i=1,nparms)
            ELSE
               DO i = 1 , nparms
                  mxbval(i) = 1
               ENDDO
            ENDIF
!
!        Read the values for INMIN.
!
            IF ( nep ) THEN
               READ (NIN,FMT=*) (inmin(i),i=1,nparms)
               DO i = 1 , nparms
                  IF ( inmin(i)<0 ) THEN
                     WRITE (NOUT,FMT=99010) ' INMIN ' , inmin(i) , 0
                     fatal = .TRUE.
                  ENDIF
               ENDDO
               WRITE (NOUT,FMT=99016) 'INMIN: ' , (inmin(i),i=1,nparms)
            ELSE
               DO i = 1 , nparms
                  inmin(i) = 1
               ENDDO
            ENDIF
!
!        Read the values for INWIN.
!
            IF ( nep ) THEN
               READ (NIN,FMT=*) (inwin(i),i=1,nparms)
               DO i = 1 , nparms
                  IF ( inwin(i)<0 ) THEN
                     WRITE (NOUT,FMT=99010) ' INWIN ' , inwin(i) , 0
                     fatal = .TRUE.
                  ENDIF
               ENDDO
               WRITE (NOUT,FMT=99016) 'INWIN: ' , (inwin(i),i=1,nparms)
            ELSE
               DO i = 1 , nparms
                  inwin(i) = 1
               ENDDO
            ENDIF
!
!        Read the values for INIBL.
!
            IF ( nep ) THEN
               READ (NIN,FMT=*) (inibl(i),i=1,nparms)
               DO i = 1 , nparms
                  IF ( inibl(i)<0 ) THEN
                     WRITE (NOUT,FMT=99010) ' INIBL ' , inibl(i) , 0
                     fatal = .TRUE.
                  ENDIF
               ENDDO
               WRITE (NOUT,FMT=99016) 'INIBL: ' , (inibl(i),i=1,nparms)
            ELSE
               DO i = 1 , nparms
                  inibl(i) = 1
               ENDDO
            ENDIF
!
!        Read the values for ISHFTS.
!
            IF ( nep ) THEN
               READ (NIN,FMT=*) (ishfts(i),i=1,nparms)
               DO i = 1 , nparms
                  IF ( ishfts(i)<0 ) THEN
                     WRITE (NOUT,FMT=99010) ' ISHFTS ' , ishfts(i) , 0
                     fatal = .TRUE.
                  ENDIF
               ENDDO
               WRITE (NOUT,FMT=99016) 'ISHFTS: ' ,                      &
     &                                (ishfts(i),i=1,nparms)
            ELSE
               DO i = 1 , nparms
                  ishfts(i) = 1
               ENDDO
            ENDIF
!
!        Read the values for IACC22.
!
            IF ( nep .OR. cgg ) THEN
               READ (NIN,FMT=*) (iacc22(i),i=1,nparms)
               DO i = 1 , nparms
                  IF ( iacc22(i)<0 ) THEN
                     WRITE (NOUT,FMT=99010) ' IACC22 ' , iacc22(i) , 0
                     fatal = .TRUE.
                  ENDIF
               ENDDO
               WRITE (NOUT,FMT=99016) 'IACC22: ' ,                      &
     &                                (iacc22(i),i=1,nparms)
            ELSE
               DO i = 1 , nparms
                  iacc22(i) = 1
               ENDDO
            ENDIF
!
!        Read the values for NBCOL.
!
            IF ( cgg ) THEN
               READ (NIN,FMT=*) (nbcol(i),i=1,nparms)
               DO i = 1 , nparms
                  IF ( nbcol(i)<0 ) THEN
                     WRITE (NOUT,FMT=99010) 'NBCOL ' , nbcol(i) , 0
                     fatal = .TRUE.
                  ELSEIF ( nbcol(i)>NMAX ) THEN
                     WRITE (NOUT,FMT=99011) 'NBCOL ' , nbcol(i) , NMAX
                     fatal = .TRUE.
                  ENDIF
               ENDDO
               WRITE (NOUT,FMT=99016) 'NBCOL:' , (nbcol(i),i=1,nparms)
            ELSE
               DO i = 1 , nparms
                  nbcol(i) = 1
               ENDDO
            ENDIF
         ENDIF
!
!     Calculate and print the machine dependent constants.
!
         WRITE (NOUT,FMT=*)
         eps = SLAMCH('Underflow threshold')
         WRITE (NOUT,FMT=99018) 'underflow' , eps
         eps = SLAMCH('Overflow threshold')
         WRITE (NOUT,FMT=99018) 'overflow ' , eps
         eps = SLAMCH('Epsilon')
         WRITE (NOUT,FMT=99018) 'precision' , eps
!
!     Read the threshold value for the test ratios.
!
         READ (NIN,FMT=*) thresh
         WRITE (NOUT,FMT=99017) thresh
         IF ( sep .OR. svd .OR. cgg ) THEN
!
!        Read the flag that indicates whether to test LAPACK routines.
!
            READ (NIN,FMT=*) tstchk
!
!        Read the flag that indicates whether to test driver routines.
!
            READ (NIN,FMT=*) tstdrv
         ENDIF
!
!     Read the flag that indicates whether to test the error exits.
!
         READ (NIN,FMT=*) tsterr
!
!     Read the code describing how to set the random number seed.
!
         READ (NIN,FMT=*) newsd
!
!     If NEWSD = 2, read another line with 4 integers for the seed.
!
         IF ( newsd==2 ) READ (NIN,FMT=*) (ioldsd(i),i=1,4)
!
         DO i = 1 , 4
            iseed(i) = ioldsd(i)
         ENDDO
!
         IF ( fatal ) THEN
            WRITE (NOUT,FMT=99001)
            STOP
         ENDIF
         EXIT
      ENDDO
      DO
!
!     Read the input lines indicating the test path and its parameters.
!     The first three characters indicate the test path, and the number
!     of test matrix types must be the first nonblank item in columns
!     4-80.
!
!
         IF ( .NOT.(cgx .OR. cxv) ) THEN
!
 120        READ (NIN,FMT='(A80)',END=200) line
            c3 = line(1:3)
            lenp = LEN(line)
            i = 3
            itmp = 0
            i1 = 0
            DO
               i = i + 1
               IF ( i>lenp ) THEN
                  IF ( i1<=0 ) ntypes = MAXT
                  EXIT
               ENDIF
               IF ( line(i:i)/=' ' .AND. line(i:i)/=',' ) THEN
                  i1 = i
                  c1 = line(i1:i1)
!
!        Check that a valid integer was read
!
                  DO k = 1 , 10
                     IF ( c1==intstr(k:k) ) THEN
                        ic = k - 1
                        GOTO 125
                     ENDIF
                  ENDDO
                  WRITE (NOUT,FMT=99008) i , line
                  GOTO 120
 125              itmp = 10*itmp + ic
               ELSEIF ( i1>0 ) THEN
                  EXIT
               ENDIF
            ENDDO
            ntypes = itmp
!
!     Skip the tests if NTYPES is <= 0.
!
            IF ( .NOT.(cev .OR. ces .OR. cvx .OR. csx .OR. cgv .OR. cgs)&
     &           .AND. ntypes<=0 ) THEN
               WRITE (NOUT,FMT=99009) c3
               GOTO 120
            ENDIF
!
         ELSE
            IF ( cgx ) c3 = 'CGX'
            IF ( cxv ) c3 = 'CXV'
         ENDIF
!
!     Reset the random number seed.
!
         IF ( newsd==0 ) THEN
            DO k = 1 , 4
               iseed(k) = ioldsd(k)
            ENDDO
         ENDIF
!
         IF ( LSAMEN(3,c3,'CHS') .OR. LSAMEN(3,c3,'NEP') ) THEN
!
!        -------------------------------------
!        NEP:  Nonsymmetric Eigenvalue Problem
!        -------------------------------------
!        Vary the parameters
!           NB    = block size
!           NBMIN = minimum block size
!           NX    = crossover point
!           NS    = number of shifts
!           MAXB  = minimum submatrix size
!
            maxtyp = 21
            ntypes = MIN(maxtyp,ntypes)
            CALL ALAREQ(c3,ntypes,dotype,maxtyp,NIN,NOUT)
            CALL XLAENV(1,1)
            IF ( tsterr ) CALL CERRHS('CHSEQR',NOUT)
            DO i = 1 , nparms
               CALL XLAENV(1,nbval(i))
               CALL XLAENV(2,nbmin(i))
               CALL XLAENV(3,nxval(i))
               CALL XLAENV(12,MAX(11,inmin(i)))
               CALL XLAENV(13,inwin(i))
               CALL XLAENV(14,inibl(i))
               CALL XLAENV(15,ishfts(i))
               CALL XLAENV(16,iacc22(i))
!
               IF ( newsd==0 ) THEN
                  DO k = 1 , 4
                     iseed(k) = ioldsd(k)
                  ENDDO
               ENDIF
               WRITE (NOUT,FMT=99038) c3 , nbval(i) , nbmin(i) ,        &
     &                                nxval(i) , MAX(11,inmin(i)) ,     &
     &                                inwin(i) , inibl(i) , ishfts(i) , &
     &                                iacc22(i)
               CALL CCHKHS(nn,nval,maxtyp,dotype,iseed,thresh,NOUT,     &
     &                     a(1,1),NMAX,a(1,2),a(1,3),a(1,4),a(1,5),NMAX,&
     &                     a(1,6),a(1,7),dc(1,1),dc(1,2),a(1,8),a(1,9), &
     &                     a(1,10),a(1,11),a(1,12),dc(1,3),work,LWORK,  &
     &                     rwork,iwork,logwrk,result,info)
               IF ( info/=0 ) WRITE (NOUT,FMT=99019) 'CCHKHS' , info
            ENDDO
!
         ELSEIF ( LSAMEN(3,c3,'CST') .OR. LSAMEN(3,c3,'SEP') .OR.       &
     &            LSAMEN(3,c3,'SE2') ) THEN
!
!        ----------------------------------
!        SEP:  Symmetric Eigenvalue Problem
!        ----------------------------------
!        Vary the parameters
!           NB    = block size
!           NBMIN = minimum block size
!           NX    = crossover point
!
            maxtyp = 21
            ntypes = MIN(maxtyp,ntypes)
            CALL ALAREQ(c3,ntypes,dotype,maxtyp,NIN,NOUT)
            CALL XLAENV(1,1)
            CALL XLAENV(9,25)
            IF ( tsterr ) THEN
#if defined(_OPENMP)
               n_threads = OMP_GET_NUM_THREADS()
               CALL OMP_SET_NUM_THREADS(1)
#endif
               CALL CERRST('CST',NOUT)
#if defined(_OPENMP)
               CALL OMP_SET_NUM_THREADS(n_threads)
#endif
            ENDIF
            DO i = 1 , nparms
               CALL XLAENV(1,nbval(i))
               CALL XLAENV(2,nbmin(i))
               CALL XLAENV(3,nxval(i))
!
               IF ( newsd==0 ) THEN
                  DO k = 1 , 4
                     iseed(k) = ioldsd(k)
                  ENDDO
               ENDIF
               WRITE (NOUT,FMT=99002) c3 , nbval(i) , nbmin(i) ,        &
     &                                nxval(i)
               IF ( tstchk ) THEN
                  IF ( LSAMEN(3,c3,'SE2') ) THEN
                     CALL CCHKST2STG(nn,nval,maxtyp,dotype,iseed,thresh,&
     &                               NOUT,a(1,1),NMAX,a(1,2),dr(1,1),   &
     &                               dr(1,2),dr(1,3),dr(1,4),dr(1,5),   &
     &                               dr(1,6),dr(1,7),dr(1,8),dr(1,9),   &
     &                               dr(1,10),dr(1,11),a(1,3),NMAX,     &
     &                               a(1,4),a(1,5),dc(1,1),a(1,6),work, &
     &                               LWORK,rwork,LWORK,iwork,LIWORK,    &
     &                               result,info)
                  ELSE
                     CALL CCHKST(nn,nval,maxtyp,dotype,iseed,thresh,    &
     &                           NOUT,a(1,1),NMAX,a(1,2),dr(1,1),dr(1,2)&
     &                           ,dr(1,3),dr(1,4),dr(1,5),dr(1,6),      &
     &                           dr(1,7),dr(1,8),dr(1,9),dr(1,10),      &
     &                           dr(1,11),a(1,3),NMAX,a(1,4),a(1,5),    &
     &                           dc(1,1),a(1,6),work,LWORK,rwork,LWORK, &
     &                           iwork,LIWORK,result,info)
                  ENDIF
                  IF ( info/=0 ) WRITE (NOUT,FMT=99019) 'CCHKST' , info
               ENDIF
               IF ( tstdrv ) THEN
                  IF ( LSAMEN(3,c3,'SE2') ) THEN
                     CALL CDRVST2STG(nn,nval,18,dotype,iseed,thresh,    &
     &                               NOUT,a(1,1),NMAX,dr(1,3),dr(1,4),  &
     &                               dr(1,5),dr(1,8),dr(1,9),dr(1,10),  &
     &                               a(1,2),NMAX,a(1,3),dc(1,1),a(1,4), &
     &                               work,LWORK,rwork,LWORK,iwork,      &
     &                               LIWORK,result,info)
                  ELSE
                     CALL CDRVST(nn,nval,18,dotype,iseed,thresh,NOUT,   &
     &                           a(1,1),NMAX,dr(1,3),dr(1,4),dr(1,5),   &
     &                           dr(1,8),dr(1,9),dr(1,10),a(1,2),NMAX,  &
     &                           a(1,3),dc(1,1),a(1,4),work,LWORK,rwork,&
     &                           LWORK,iwork,LIWORK,result,info)
                  ENDIF
                  IF ( info/=0 ) WRITE (NOUT,FMT=99019) 'CDRVST' , info
               ENDIF
            ENDDO
!
         ELSEIF ( LSAMEN(3,c3,'CSG') ) THEN
!
!        ----------------------------------------------
!        CSG:  Hermitian Generalized Eigenvalue Problem
!        ----------------------------------------------
!        Vary the parameters
!           NB    = block size
!           NBMIN = minimum block size
!           NX    = crossover point
!
            maxtyp = 21
            ntypes = MIN(maxtyp,ntypes)
            CALL ALAREQ(c3,ntypes,dotype,maxtyp,NIN,NOUT)
            CALL XLAENV(9,25)
            DO i = 1 , nparms
               CALL XLAENV(1,nbval(i))
               CALL XLAENV(2,nbmin(i))
               CALL XLAENV(3,nxval(i))
!
               IF ( newsd==0 ) THEN
                  DO k = 1 , 4
                     iseed(k) = ioldsd(k)
                  ENDDO
               ENDIF
               WRITE (NOUT,FMT=99002) c3 , nbval(i) , nbmin(i) ,        &
     &                                nxval(i)
               IF ( tstchk ) THEN
!               CALL CDRVSG( NN, NVAL, MAXTYP, DOTYPE, ISEED, THRESH,
!     $                      NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), NMAX,
!     $                      DR( 1, 3 ), A( 1, 3 ), NMAX, A( 1, 4 ),
!     $                      A( 1, 5 ), A( 1, 6 ), A( 1, 7 ), WORK,
!     $                      LWORK, RWORK, LWORK, IWORK, LIWORK, RESULT,
!     $                      INFO )
                  CALL CDRVSG2STG(nn,nval,maxtyp,dotype,iseed,thresh,   &
     &                            NOUT,a(1,1),NMAX,a(1,2),NMAX,dr(1,3), &
     &                            dr(1,4),a(1,3),NMAX,a(1,4),a(1,5),    &
     &                            a(1,6),a(1,7),work,LWORK,rwork,LWORK, &
     &                            iwork,LIWORK,result,info)
                  IF ( info/=0 ) WRITE (NOUT,FMT=99019) 'CDRVSG' , info
               ENDIF
            ENDDO
!
         ELSEIF ( LSAMEN(3,c3,'CBD') .OR. LSAMEN(3,c3,'SVD') ) THEN
!
!        ----------------------------------
!        SVD:  Singular Value Decomposition
!        ----------------------------------
!        Vary the parameters
!           NB    = block size
!           NBMIN = minimum block size
!           NX    = crossover point
!           NRHS  = number of right hand sides
!
            maxtyp = 16
            ntypes = MIN(maxtyp,ntypes)
            CALL ALAREQ(c3,ntypes,dotype,maxtyp,NIN,NOUT)
            CALL XLAENV(9,25)
!
!        Test the error exits
!
            CALL XLAENV(1,1)
            IF ( tsterr .AND. tstchk ) CALL CERRBD('CBD',NOUT)
            IF ( tsterr .AND. tstdrv ) CALL CERRED('CBD',NOUT)
!
            DO i = 1 , nparms
               nrhs = nsval(i)
               CALL XLAENV(1,nbval(i))
               CALL XLAENV(2,nbmin(i))
               CALL XLAENV(3,nxval(i))
               IF ( newsd==0 ) THEN
                  DO k = 1 , 4
                     iseed(k) = ioldsd(k)
                  ENDDO
               ENDIF
               WRITE (NOUT,FMT=99004) c3 , nbval(i) , nbmin(i) ,        &
     &                                nxval(i) , nrhs
               IF ( tstchk ) THEN
                  CALL CCHKBD(nn,mval,nval,maxtyp,dotype,nrhs,iseed,    &
     &                        thresh,a(1,1),NMAX,dr(1,1),dr(1,2),dr(1,3)&
     &                        ,dr(1,4),a(1,2),NMAX,a(1,3),a(1,4),a(1,5),&
     &                        NMAX,a(1,6),NMAX,a(1,7),a(1,8),work,LWORK,&
     &                        rwork,NOUT,info)
                  IF ( info/=0 ) WRITE (NOUT,FMT=99019) 'CCHKBD' , info
               ENDIF
               IF ( tstdrv ) CALL CDRVBD(nn,mval,nval,maxtyp,dotype,    &
     &              iseed,thresh,a(1,1),NMAX,a(1,2),NMAX,a(1,3),NMAX,   &
     &              a(1,4),a(1,5),a(1,6),dr(1,1),dr(1,2),dr(1,3),work,  &
     &              LWORK,rwork,iwork,NOUT,info)
            ENDDO
!
         ELSEIF ( LSAMEN(3,c3,'CEV') ) THEN
!
!        --------------------------------------------
!        CEV:  Nonsymmetric Eigenvalue Problem Driver
!              CGEEV (eigenvalues and eigenvectors)
!        --------------------------------------------
!
            maxtyp = 21
            ntypes = MIN(maxtyp,ntypes)
            IF ( ntypes<=0 ) THEN
               WRITE (NOUT,FMT=99009) c3
            ELSE
               IF ( tsterr ) CALL CERRED(c3,NOUT)
               CALL ALAREQ(c3,ntypes,dotype,maxtyp,NIN,NOUT)
               CALL CDRVEV(nn,nval,ntypes,dotype,iseed,thresh,NOUT,     &
     &                     a(1,1),NMAX,a(1,2),dc(1,1),dc(1,2),a(1,3),   &
     &                     NMAX,a(1,4),NMAX,a(1,5),NMAX,result,work,    &
     &                     LWORK,rwork,iwork,info)
               IF ( info/=0 ) WRITE (NOUT,FMT=99019) 'CGEEV' , info
            ENDIF
            WRITE (NOUT,FMT=99026)
            GOTO 100
!
         ELSEIF ( LSAMEN(3,c3,'CES') ) THEN
!
!        --------------------------------------------
!        CES:  Nonsymmetric Eigenvalue Problem Driver
!              CGEES (Schur form)
!        --------------------------------------------
!
            maxtyp = 21
            ntypes = MIN(maxtyp,ntypes)
            IF ( ntypes<=0 ) THEN
               WRITE (NOUT,FMT=99009) c3
            ELSE
               IF ( tsterr ) CALL CERRED(c3,NOUT)
               CALL ALAREQ(c3,ntypes,dotype,maxtyp,NIN,NOUT)
               CALL CDRVES(nn,nval,ntypes,dotype,iseed,thresh,NOUT,     &
     &                     a(1,1),NMAX,a(1,2),a(1,3),dc(1,1),dc(1,2),   &
     &                     a(1,4),NMAX,result,work,LWORK,rwork,iwork,   &
     &                     logwrk,info)
               IF ( info/=0 ) WRITE (NOUT,FMT=99019) 'CGEES' , info
            ENDIF
            WRITE (NOUT,FMT=99026)
            GOTO 100
!
         ELSEIF ( LSAMEN(3,c3,'CVX') ) THEN
!
!        --------------------------------------------------------------
!        CVX:  Nonsymmetric Eigenvalue Problem Expert Driver
!              CGEEVX (eigenvalues, eigenvectors and condition numbers)
!        --------------------------------------------------------------
!
            maxtyp = 21
            ntypes = MIN(maxtyp,ntypes)
            IF ( ntypes<0 ) THEN
               WRITE (NOUT,FMT=99009) c3
            ELSE
               IF ( tsterr ) CALL CERRED(c3,NOUT)
               CALL ALAREQ(c3,ntypes,dotype,maxtyp,NIN,NOUT)
               CALL CDRVVX(nn,nval,ntypes,dotype,iseed,thresh,NIN,NOUT, &
     &                     a(1,1),NMAX,a(1,2),dc(1,1),dc(1,2),a(1,3),   &
     &                     NMAX,a(1,4),NMAX,a(1,5),NMAX,dr(1,1),dr(1,2),&
     &                     dr(1,3),dr(1,4),dr(1,5),dr(1,6),dr(1,7),     &
     &                     dr(1,8),result,work,LWORK,rwork,info)
               IF ( info/=0 ) WRITE (NOUT,FMT=99019) 'CGEEVX' , info
            ENDIF
            WRITE (NOUT,FMT=99026)
            GOTO 100
!
         ELSEIF ( LSAMEN(3,c3,'CSX') ) THEN
!
!        ---------------------------------------------------
!        CSX:  Nonsymmetric Eigenvalue Problem Expert Driver
!              CGEESX (Schur form and condition numbers)
!        ---------------------------------------------------
!
            maxtyp = 21
            ntypes = MIN(maxtyp,ntypes)
            IF ( ntypes<0 ) THEN
               WRITE (NOUT,FMT=99009) c3
            ELSE
               IF ( tsterr ) CALL CERRED(c3,NOUT)
               CALL ALAREQ(c3,ntypes,dotype,maxtyp,NIN,NOUT)
               CALL CDRVSX(nn,nval,ntypes,dotype,iseed,thresh,NIN,NOUT, &
     &                     a(1,1),NMAX,a(1,2),a(1,3),dc(1,1),dc(1,2),   &
     &                     dc(1,3),a(1,4),NMAX,a(1,5),result,work,LWORK,&
     &                     rwork,logwrk,info)
               IF ( info/=0 ) WRITE (NOUT,FMT=99019) 'CGEESX' , info
            ENDIF
            WRITE (NOUT,FMT=99026)
            GOTO 100
!
         ELSEIF ( LSAMEN(3,c3,'CGG') ) THEN
!
!        -------------------------------------------------
!        CGG:  Generalized Nonsymmetric Eigenvalue Problem
!        -------------------------------------------------
!        Vary the parameters
!           NB    = block size
!           NBMIN = minimum block size
!           NS    = number of shifts
!           MAXB  = minimum submatrix size
!           IACC22: structured matrix multiply
!           NBCOL = minimum column dimension for blocks
!
            maxtyp = 26
            ntypes = MIN(maxtyp,ntypes)
            CALL ALAREQ(c3,ntypes,dotype,maxtyp,NIN,NOUT)
            CALL XLAENV(1,1)
            IF ( tstchk .AND. tsterr ) CALL CERRGG(c3,NOUT)
            DO i = 1 , nparms
               CALL XLAENV(1,nbval(i))
               CALL XLAENV(2,nbmin(i))
               CALL XLAENV(4,nsval(i))
               CALL XLAENV(8,mxbval(i))
               CALL XLAENV(16,iacc22(i))
               CALL XLAENV(5,nbcol(i))
!
               IF ( newsd==0 ) THEN
                  DO k = 1 , 4
                     iseed(k) = ioldsd(k)
                  ENDDO
               ENDIF
               WRITE (NOUT,FMT=99003) c3 , nbval(i) , nbmin(i) ,        &
     &                                nsval(i) , mxbval(i) , iacc22(i) ,&
     &                                nbcol(i)
               tstdif = .FALSE.
               thrshn = 10.
               IF ( tstchk ) THEN
                  CALL CCHKGG(nn,nval,maxtyp,dotype,iseed,thresh,tstdif,&
     &                        thrshn,NOUT,a(1,1),NMAX,a(1,2),a(1,3),    &
     &                        a(1,4),a(1,5),a(1,6),a(1,7),a(1,8),a(1,9),&
     &                        NMAX,a(1,10),a(1,11),a(1,12),dc(1,1),     &
     &                        dc(1,2),dc(1,3),dc(1,4),a(1,13),a(1,14),  &
     &                        work,LWORK,rwork,logwrk,result,info)
                  IF ( info/=0 ) WRITE (NOUT,FMT=99019) 'CCHKGG' , info
               ENDIF
            ENDDO
!
         ELSEIF ( LSAMEN(3,c3,'CGS') ) THEN
!
!        -------------------------------------------------
!        CGS:  Generalized Nonsymmetric Eigenvalue Problem
!              CGGES (Schur form)
!        -------------------------------------------------
!
            maxtyp = 26
            ntypes = MIN(maxtyp,ntypes)
            IF ( ntypes<=0 ) THEN
               WRITE (NOUT,FMT=99009) c3
            ELSE
               IF ( tsterr ) CALL CERRGG(c3,NOUT)
               CALL ALAREQ(c3,ntypes,dotype,maxtyp,NIN,NOUT)
               CALL CDRGES(nn,nval,maxtyp,dotype,iseed,thresh,NOUT,     &
     &                     a(1,1),NMAX,a(1,2),a(1,3),a(1,4),a(1,7),NMAX,&
     &                     a(1,8),dc(1,1),dc(1,2),work,LWORK,rwork,     &
     &                     result,logwrk,info)
!
               IF ( info/=0 ) WRITE (NOUT,FMT=99019) 'CDRGES' , info
!
! Blocked version
!
               CALL XLAENV(16,2)
               CALL CDRGES3(nn,nval,maxtyp,dotype,iseed,thresh,NOUT,    &
     &                      a(1,1),NMAX,a(1,2),a(1,3),a(1,4),a(1,7),    &
     &                      NMAX,a(1,8),dc(1,1),dc(1,2),work,LWORK,     &
     &                      rwork,result,logwrk,info)
!
               IF ( info/=0 ) WRITE (NOUT,FMT=99019) 'CDRGES3' , info
            ENDIF
            WRITE (NOUT,FMT=99026)
 
            GOTO 100
!
         ELSEIF ( cgx ) THEN
!
!        -------------------------------------------------
!        CGX  Generalized Nonsymmetric Eigenvalue Problem
!              CGGESX (Schur form and condition numbers)
!        -------------------------------------------------
!
            maxtyp = 5
            ntypes = maxtyp
            IF ( nn<0 ) THEN
               WRITE (NOUT,FMT=99009) c3
            ELSE
               IF ( tsterr ) CALL CERRGG(c3,NOUT)
               CALL ALAREQ(c3,ntypes,dotype,maxtyp,NIN,NOUT)
               CALL XLAENV(5,2)
               CALL CDRGSX(nn,NCMAX,thresh,NIN,NOUT,a(1,1),NMAX,a(1,2), &
     &                     a(1,3),a(1,4),a(1,5),a(1,6),dc(1,1),dc(1,2), &
     &                     c,NCMAX*NCMAX,s,work,LWORK,rwork,iwork,      &
     &                     LIWORK,logwrk,info)
               IF ( info/=0 ) WRITE (NOUT,FMT=99019) 'CDRGSX' , info
            ENDIF
            WRITE (NOUT,FMT=99026)
            GOTO 100
!
         ELSEIF ( LSAMEN(3,c3,'CGV') ) THEN
!
!        -------------------------------------------------
!        CGV:  Generalized Nonsymmetric Eigenvalue Problem
!              CGGEV (Eigenvalue/vector form)
!        -------------------------------------------------
!
            maxtyp = 26
            ntypes = MIN(maxtyp,ntypes)
            IF ( ntypes<=0 ) THEN
               WRITE (NOUT,FMT=99009) c3
            ELSE
               IF ( tsterr ) CALL CERRGG(c3,NOUT)
               CALL ALAREQ(c3,ntypes,dotype,maxtyp,NIN,NOUT)
               CALL CDRGEV(nn,nval,maxtyp,dotype,iseed,thresh,NOUT,     &
     &                     a(1,1),NMAX,a(1,2),a(1,3),a(1,4),a(1,7),NMAX,&
     &                     a(1,8),a(1,9),NMAX,dc(1,1),dc(1,2),dc(1,3),  &
     &                     dc(1,4),work,LWORK,rwork,result,info)
               IF ( info/=0 ) WRITE (NOUT,FMT=99019) 'CDRGEV' , info
!
! Blocked version
!
               CALL XLAENV(16,2)
               CALL CDRGEV3(nn,nval,maxtyp,dotype,iseed,thresh,NOUT,    &
     &                      a(1,1),NMAX,a(1,2),a(1,3),a(1,4),a(1,7),    &
     &                      NMAX,a(1,8),a(1,9),NMAX,dc(1,1),dc(1,2),    &
     &                      dc(1,3),dc(1,4),work,LWORK,rwork,result,    &
     &                      info)
               IF ( info/=0 ) WRITE (NOUT,FMT=99019) 'CDRGEV3' , info
            ENDIF
            WRITE (NOUT,FMT=99026)
            GOTO 100
!
         ELSEIF ( cxv ) THEN
!
!        -------------------------------------------------
!        CXV:  Generalized Nonsymmetric Eigenvalue Problem
!              CGGEVX (eigenvalue/vector with condition numbers)
!        -------------------------------------------------
!
            maxtyp = 2
            ntypes = maxtyp
            IF ( nn<0 ) THEN
               WRITE (NOUT,FMT=99009) c3
            ELSE
               IF ( tsterr ) CALL CERRGG(c3,NOUT)
               CALL ALAREQ(c3,ntypes,dotype,maxtyp,NIN,NOUT)
               CALL CDRGVX(nn,thresh,NIN,NOUT,a(1,1),NMAX,a(1,2),a(1,3),&
     &                     a(1,4),dc(1,1),dc(1,2),a(1,5),a(1,6),iwork(1)&
     &                     ,iwork(2),dr(1,1),dr(1,2),dr(1,3),dr(1,4),   &
     &                     dr(1,5),dr(1,6),work,LWORK,rwork,iwork(3),   &
     &                     LIWORK-2,result,logwrk,info)
!
               IF ( info/=0 ) WRITE (NOUT,FMT=99019) 'CDRGVX' , info
            ENDIF
            WRITE (NOUT,FMT=99026)
            GOTO 100
!
         ELSEIF ( LSAMEN(3,c3,'CHB') ) THEN
!
!        ------------------------------
!        CHB:  Hermitian Band Reduction
!        ------------------------------
!
            maxtyp = 15
            ntypes = MIN(maxtyp,ntypes)
            CALL ALAREQ(c3,ntypes,dotype,maxtyp,NIN,NOUT)
            IF ( tsterr ) THEN
#if defined(_OPENMP)
               n_threads = OMP_GET_NUM_THREADS()
               CALL OMP_SET_NUM_THREADS(1)
#endif
               CALL CERRST('CHB',NOUT)
#if defined(_OPENMP)
               CALL OMP_SET_NUM_THREADS(n_threads)
#endif
            ENDIF
!         CALL CCHKHB( NN, NVAL, NK, KVAL, MAXTYP, DOTYPE, ISEED, THRESH,
!     $                NOUT, A( 1, 1 ), NMAX, DR( 1, 1 ), DR( 1, 2 ),
!     $                A( 1, 2 ), NMAX, WORK, LWORK, RWORK, RESULT,
!     $                INFO )
            CALL CCHKHB2STG(nn,nval,nk,kval,maxtyp,dotype,iseed,thresh, &
     &                      NOUT,a(1,1),NMAX,dr(1,1),dr(1,2),dr(1,3),   &
     &                      dr(1,4),dr(1,5),a(1,2),NMAX,work,LWORK,     &
     &                      rwork,result,info)
            IF ( info/=0 ) WRITE (NOUT,FMT=99019) 'CCHKHB' , info
!
         ELSEIF ( LSAMEN(3,c3,'CBB') ) THEN
!
!        ------------------------------
!        CBB:  General Band Reduction
!        ------------------------------
!
            maxtyp = 15
            ntypes = MIN(maxtyp,ntypes)
            CALL ALAREQ(c3,ntypes,dotype,maxtyp,NIN,NOUT)
            DO i = 1 , nparms
               nrhs = nsval(i)
!
               IF ( newsd==0 ) THEN
                  DO k = 1 , 4
                     iseed(k) = ioldsd(k)
                  ENDDO
               ENDIF
               WRITE (NOUT,FMT=99033) c3 , nrhs
               CALL CCHKBB(nn,mval,nval,nk,kval,maxtyp,dotype,nrhs,     &
     &                     iseed,thresh,NOUT,a(1,1),NMAX,a(1,2),2*NMAX, &
     &                     dr(1,1),dr(1,2),a(1,4),NMAX,a(1,5),NMAX,     &
     &                     a(1,6),NMAX,a(1,7),work,LWORK,rwork,result,  &
     &                     info)
               IF ( info/=0 ) WRITE (NOUT,FMT=99019) 'CCHKBB' , info
            ENDDO
!
         ELSEIF ( LSAMEN(3,c3,'GLM') ) THEN
!
!        -----------------------------------------
!        GLM:  Generalized Linear Regression Model
!        -----------------------------------------
!
            CALL XLAENV(1,1)
            IF ( tsterr ) CALL CERRGG('GLM',NOUT)
            CALL CCKGLM(nn,nval,mval,pval,ntypes,iseed,thresh,NMAX,     &
     &                  a(1,1),a(1,2),b(1,1),b(1,2),x,work,dr(1,1),NIN, &
     &                  NOUT,info)
            IF ( info/=0 ) WRITE (NOUT,FMT=99019) 'CCKGLM' , info
!
         ELSEIF ( LSAMEN(3,c3,'GQR') ) THEN
!
!        ------------------------------------------
!        GQR:  Generalized QR and RQ factorizations
!        ------------------------------------------
!
            CALL XLAENV(1,1)
            IF ( tsterr ) CALL CERRGG('GQR',NOUT)
            CALL CCKGQR(nn,mval,nn,pval,nn,nval,ntypes,iseed,thresh,    &
     &                  NMAX,a(1,1),a(1,2),a(1,3),a(1,4),taua,b(1,1),   &
     &                  b(1,2),b(1,3),b(1,4),b(1,5),taub,work,dr(1,1),  &
     &                  NIN,NOUT,info)
            IF ( info/=0 ) WRITE (NOUT,FMT=99019) 'CCKGQR' , info
!
         ELSEIF ( LSAMEN(3,c3,'GSV') ) THEN
!
!        ----------------------------------------------
!        GSV:  Generalized Singular Value Decomposition
!        ----------------------------------------------
!
            CALL XLAENV(1,1)
            IF ( tsterr ) CALL CERRGG('GSV',NOUT)
            CALL CCKGSV(nn,mval,pval,nval,ntypes,iseed,thresh,NMAX,     &
     &                  a(1,1),a(1,2),b(1,1),b(1,2),a(1,3),b(1,3),a(1,4)&
     &                  ,alpha,beta,b(1,4),iwork,work,dr(1,1),NIN,NOUT, &
     &                  info)
            IF ( info/=0 ) WRITE (NOUT,FMT=99019) 'CCKGSV' , info
!
         ELSEIF ( LSAMEN(3,c3,'CSD') ) THEN
!
!        ----------------------------------------------
!        CSD:  CS Decomposition
!        ----------------------------------------------
!
            CALL XLAENV(1,1)
            IF ( tsterr ) CALL CERRGG('CSD',NOUT)
            CALL CCKCSD(nn,mval,pval,nval,ntypes,iseed,thresh,NMAX,     &
     &                  a(1,1),a(1,2),a(1,3),a(1,4),a(1,5),a(1,6),rwork,&
     &                  iwork,work,dr(1,1),NIN,NOUT,info)
            IF ( info/=0 ) WRITE (NOUT,FMT=99019) 'CCKCSD' , info
!
         ELSEIF ( LSAMEN(3,c3,'LSE') ) THEN
!
!        --------------------------------------
!        LSE:  Constrained Linear Least Squares
!        --------------------------------------
!
            CALL XLAENV(1,1)
            IF ( tsterr ) CALL CERRGG('LSE',NOUT)
            CALL CCKLSE(nn,mval,pval,nval,ntypes,iseed,thresh,NMAX,     &
     &                  a(1,1),a(1,2),b(1,1),b(1,2),x,work,dr(1,1),NIN, &
     &                  NOUT,info)
            IF ( info/=0 ) WRITE (NOUT,FMT=99019) 'CCKLSE' , info
         ELSE
            WRITE (NOUT,FMT=*)
            WRITE (NOUT,FMT=*)
            WRITE (NOUT,FMT=99007) c3
         ENDIF
         IF ( cgx .OR. cxv ) EXIT
      ENDDO
 200  WRITE (NOUT,FMT=99005)
      s2 = SECOND()
      WRITE (NOUT,FMT=99006) s2 - s1
!
      DEALLOCATE (s,STAT=allocatestatus)
      DEALLOCATE (a,STAT=allocatestatus)
      DEALLOCATE (b,STAT=allocatestatus)
      DEALLOCATE (c,STAT=allocatestatus)
      DEALLOCATE (rwork,STAT=allocatestatus)
      DEALLOCATE (work,STAT=allocatestatus)
!
99001 FORMAT (/' Execution not attempted due to input errors')
99002 FORMAT (//1X,A3,':  NB =',I4,', NBMIN =',I4,', NX =',I4)
99003 FORMAT (//1X,A3,':  NB =',I4,', NBMIN =',I4,', NS =',I4,          &
     &        ', MAXB =',I4,', IACC22 =',I4,', NBCOL =',I4)
99004 FORMAT (//1X,A3,':  NB =',I4,', NBMIN =',I4,', NX =',I4,          &
     &        ', NRHS =',I4)
99005 FORMAT (//' End of tests')
99006 FORMAT (' Total time used = ',F12.2,' seconds',/)
99007 FORMAT (1X,A3,':  Unrecognized path name')
99008 FORMAT (//' *** Invalid integer value in column ',I2,' of input', &
     &        ' line:',/A79)
99009 FORMAT (//1X,A3,' routines were not tested')
99010 FORMAT (' Invalid input value: ',A,'=',I6,'; must be >=',I6)
99011 FORMAT (' Invalid input value: ',A,'=',I6,'; must be <=',I6)
99012 FORMAT (' Tests of the Nonsymmetric Eigenvalue Problem routines')
99013 FORMAT (' Tests of the Hermitian Eigenvalue Problem routines')
99014 FORMAT (' Tests of the Singular Value Decomposition routines')
99015 FORMAT (/' The following parameter values will be used:')
99016 FORMAT (4X,A,10I6,/10X,10I6)
99017 FORMAT (/' Routines pass computational tests if test ratio is ',  &
     &        'less than',F8.2,/)
99018 FORMAT (' Relative machine ',A,' is taken to be',E16.6)
99019 FORMAT (' *** Error code from ',A,' = ',I4)
99020 FORMAT (/' Tests of the Nonsymmetric Eigenvalue Problem Driver',  &
     &        /'    CGEEV (eigenvalues and eigevectors)')
99021 FORMAT (/' Tests of the Nonsymmetric Eigenvalue Problem Driver',  &
     &        /'    CGEES (Schur form)')
99022 FORMAT (/' Tests of the Nonsymmetric Eigenvalue Problem Expert',  &
     &        ' Driver',/'    CGEEVX (eigenvalues, eigenvectors and',   &
     &        ' condition numbers)')
99023 FORMAT (/' Tests of the Nonsymmetric Eigenvalue Problem Expert',  &
     &        ' Driver',/'    CGEESX (Schur form and condition',        &
     &        ' numbers)')
99024 FORMAT (/' Tests of the Generalized Nonsymmetric Eigenvalue ',    &
     &        'Problem routines')
99025 FORMAT (' Tests of CHBTRD',/' (reduction of a Hermitian band ',   &
     &        'matrix to real tridiagonal form)')
99026 FORMAT (/1X,71('-'))
99027 FORMAT (/' LAPACK VERSION ',I1,'.',I1,'.',I1)
99028 FORMAT (/' Tests of the Generalized Linear Regression Model ',    &
     &        'routines')
99029 FORMAT (/' Tests of the Generalized QR and RQ routines')
99030 FORMAT (/' Tests of the Generalized Singular Value',              &
     &        ' Decomposition routines')
99031 FORMAT (/' Tests of the Linear Least Squares routines')
99032 FORMAT (' Tests of CGBBRD',/' (reduction of a general band ',     &
     &        'matrix to real bidiagonal form)')
99033 FORMAT (//1X,A3,':  NRHS =',I4)
99034 FORMAT (/' Tests of the Generalized Nonsymmetric Eigenvalue ',    &
     &        'Problem Expert Driver CGGESX')
99035 FORMAT (/' Tests of the Generalized Nonsymmetric Eigenvalue ',    &
     &        'Problem Driver CGGES')
99036 FORMAT (/' Tests of the Generalized Nonsymmetric Eigenvalue ',    &
     &        'Problem Driver CGGEV')
99037 FORMAT (/' Tests of the Generalized Nonsymmetric Eigenvalue ',    &
     &        'Problem Expert Driver CGGEVX')
99038 FORMAT (//1X,A3,':  NB =',I4,', NBMIN =',I4,', NX =',I4,          &
     &        ', INMIN=',I4,', INWIN =',I4,', INIBL =',I4,', ISHFTS =', &
     &        I4,', IACC22 =',I4)
99039 FORMAT (/' Tests of the CS Decomposition routines')
!
!     End of CCHKEE
!
END PROGRAM CCHKEE
