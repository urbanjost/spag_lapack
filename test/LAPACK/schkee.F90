!*==schkee.f90  processed by SPAG 7.51RB at 23:29 on  6 Mar 2022
!> \brief \b SCHKEE
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       PROGRAM SCHKEE
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SCHKEE tests the REAL LAPACK subroutines for the matrix
!> eigenvalue problem.  The test paths in this version are
!>
!> NEP (Nonsymmetric Eigenvalue Problem):
!>     Test SGEHRD, SORGHR, SHSEQR, STREVC, SHSEIN, and SORMHR
!>
!> SEP (Symmetric Eigenvalue Problem):
!>     Test SSYTRD, SORGTR, SSTEQR, SSTERF, SSTEIN, SSTEDC,
!>     and drivers SSYEV(X), SSBEV(X), SSPEV(X), SSTEV(X),
!>                 SSYEVD,   SSBEVD,   SSPEVD,   SSTEVD
!>
!> SVD (Singular Value Decomposition):
!>     Test SGEBRD, SORGBR, SBDSQR, SBDSDC
!>     and the drivers SGESVD, SGESDD
!>
!> SEV (Nonsymmetric Eigenvalue/eigenvector Driver):
!>     Test SGEEV
!>
!> SES (Nonsymmetric Schur form Driver):
!>     Test SGEES
!>
!> SVX (Nonsymmetric Eigenvalue/eigenvector Expert Driver):
!>     Test SGEEVX
!>
!> SSX (Nonsymmetric Schur form Expert Driver):
!>     Test SGEESX
!>
!> SGG (Generalized Nonsymmetric Eigenvalue Problem):
!>     Test SGGHD3, SGGBAL, SGGBAK, SHGEQZ, and STGEVC
!>
!> SGS (Generalized Nonsymmetric Schur form Driver):
!>     Test SGGES
!>
!> SGV (Generalized Nonsymmetric Eigenvalue/eigenvector Driver):
!>     Test SGGEV
!>
!> SGX (Generalized Nonsymmetric Schur form Expert Driver):
!>     Test SGGESX
!>
!> SXV (Generalized Nonsymmetric Eigenvalue/eigenvector Expert Driver):
!>     Test SGGEVX
!>
!> SSG (Symmetric Generalized Eigenvalue Problem):
!>     Test SSYGST, SSYGV, SSYGVD, SSYGVX, SSPGST, SSPGV, SSPGVD,
!>     SSPGVX, SSBGST, SSBGV, SSBGVD, and SSBGVX
!>
!> SSB (Symmetric Band Eigenvalue Problem):
!>     Test SSBTRD
!>
!> SBB (Band Singular Value Decomposition):
!>     Test SGBBRD
!>
!> SEC (Eigencondition estimation):
!>     Test SLALN2, SLASY2, SLAEQU, SLAEXC, STRSYL, STREXC, STRSNA,
!>     STRSEN, and SLAQTR
!>
!> SBL (Balancing a general matrix)
!>     Test SGEBAL
!>
!> SBK (Back transformation on a balanced matrix)
!>     Test SGEBAK
!>
!> SGL (Balancing a matrix pair)
!>     Test SGGBAL
!>
!> SGK (Back transformation on a matrix pair)
!>     Test SGGBAK
!>
!> GLM (Generalized Linear Regression Model):
!>     Tests SGGGLM
!>
!> GQR (Generalized QR and RQ factorizations):
!>     Tests SGGQRF and SGGRQF
!>
!> GSV (Generalized Singular Value Decomposition):
!>     Tests SGGSVD, SGGSVP, STGSJA, SLAGS2, SLAPLL, and SLAPMT
!>
!> CSD (CS decomposition):
!>     Tests SORCSD
!>
!> LSE (Constrained Linear Least Squares):
!>     Tests SGGLSE
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
!> SHS or NEP      21     SCHKHS
!> SST or SEP      21     SCHKST (routines)
!>                 18     SDRVST (drivers)
!> SBD or SVD      16     SCHKBD (routines)
!>                  5     SDRVBD (drivers)
!> SEV             21     SDRVEV
!> SES             21     SDRVES
!> SVX             21     SDRVVX
!> SSX             21     SDRVSX
!> SGG             26     SCHKGG (routines)
!> SGS             26     SDRGES
!> SGX              5     SDRGSX
!> SGV             26     SDRGEV
!> SXV              2     SDRGVX
!> SSG             21     SDRVSG
!> SSB             15     SCHKSB
!> SBB             15     SCHKBB
!> SEC              -     SCHKEC
!> SBL              -     SCHKBL
!> SBK              -     SCHKBK
!> SGL              -     SCHKGL
!> SGK              -     SCHKGK
!> GLM              8     SCKGLM
!> GQR              8     SCKGQR
!> GSV              8     SCKGSV
!> CSD              3     SCKCSD
!> LSE              8     SCKLSE
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
!>          The valid 3-character path names are 'NEP' or 'SHS' for the
!>          nonsymmetric eigenvalue routines.
!>
!>-----------------------------------------------------------------------
!>
!> SEP or SSG input file:
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
!>          The 3-character path names are 'SEP' or 'SST' for the
!>          symmetric eigenvalue routines and driver routines, and
!>          'SSG' for the routines for the symmetric generalized
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
!>          The 3-character path names are 'SVD' or 'SBD' for both the
!>          SVD routines and the SVD driver routines.
!>
!>-----------------------------------------------------------------------
!>
!> SEV and SES data files:
!>
!> line 1:  'SEV' or 'SES' in columns 1 to 3.
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
!> If line 7 was 2:
!>
!> line 8:  INTEGER array, dimension (4)
!>          Four integer values for the random number seed.
!>
!> lines 9 and following:  Lines specifying matrix types, as for NEP.
!>          The 3-character path name is 'SEV' to test SGEEV, or
!>          'SES' to test SGEES.
!>
!>-----------------------------------------------------------------------
!>
!> The SVX data has two parts. The first part is identical to SEV,
!> and the second part consists of test matrices with precomputed
!> solutions.
!>
!> line 1:  'SVX' in columns 1-3.
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
!> line 6:  TSTERR, LOGICAL
!>
!> line 7:  NEWSD, INTEGER
!>
!> If line 7 was 2:
!>
!> line 8:  INTEGER array, dimension (4)
!>
!> lines 9 and following: The first line contains 'SVX' in columns 1-3
!>          followed by the number of matrix types, possibly with
!>          a second line to specify certain matrix types.
!>          If the number of matrix types = 0, no testing of randomly
!>          generated examples is done, but any precomputed examples
!>          are tested.
!>
!> remaining lines : Each matrix is stored on 1+2*N lines, where N is
!>          its dimension. The first line contains the dimension (a
!>          single integer). The next N lines contain the matrix, one
!>          row per line. The last N lines correspond to each
!>          eigenvalue. Each of these last N lines contains 4 real
!>          values: the real part of the eigenvalue, the imaginary
!>          part of the eigenvalue, the reciprocal condition number of
!>          the eigenvalues, and the reciprocal condition number of the
!>          eigenvector.  The end of data is indicated by dimension N=0.
!>          Even if no data is to be tested, there must be at least one
!>          line containing N=0.
!>
!>-----------------------------------------------------------------------
!>
!> The SSX data is like SVX. The first part is identical to SEV, and the
!> second part consists of test matrices with precomputed solutions.
!>
!> line 1:  'SSX' in columns 1-3.
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
!> line 6:  TSTERR, LOGICAL
!>
!> line 7:  NEWSD, INTEGER
!>
!> If line 7 was 2:
!>
!> line 8:  INTEGER array, dimension (4)
!>
!> lines 9 and following: The first line contains 'SSX' in columns 1-3
!>          followed by the number of matrix types, possibly with
!>          a second line to specify certain matrix types.
!>          If the number of matrix types = 0, no testing of randomly
!>          generated examples is done, but any precomputed examples
!>          are tested.
!>
!> remaining lines : Each matrix is stored on 3+N lines, where N is its
!>          dimension. The first line contains the dimension N and the
!>          dimension M of an invariant subspace. The second line
!>          contains M integers, identifying the eigenvalues in the
!>          invariant subspace (by their position in a list of
!>          eigenvalues ordered by increasing real part). The next N
!>          lines contain the matrix. The last line contains the
!>          reciprocal condition number for the average of the selected
!>          eigenvalues, and the reciprocal condition number for the
!>          corresponding right invariant subspace. The end of data is
!>          indicated by a line containing N=0 and M=0. Even if no data
!>          is to be tested, there must be at least one line containing
!>          N=0 and M=0.
!>
!>-----------------------------------------------------------------------
!>
!> SGG input file:
!>
!> line 2:  NN, INTEGER
!>          Number of values of N.
!>
!> line 3:  NVAL, INTEGER array, dimension (NN)
!>          The values for the matrix dimension N.
!>
!> line 4:  NPARMS, INTEGER
!>          Number of values of the parameters NB, NBMIN, NS, MAXB, and
!>          NBCOL.
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
!>          The 3-character path name is 'SGG' for the generalized
!>          eigenvalue problem routines and driver routines.
!>
!>-----------------------------------------------------------------------
!>
!> SGS and SGV input files:
!>
!> line 1:  'SGS' or 'SGV' in columns 1 to 3.
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
!>          The 3-character path name is 'SGS' for the generalized
!>          eigenvalue problem routines and driver routines.
!>
!>-----------------------------------------------------------------------
!>
!> SXV input files:
!>
!> line 1:  'SXV' in columns 1 to 3.
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
!> remaining lines : Each example is stored on 3+2*N lines, where N is
!>          its dimension. The first line contains the dimension (a
!>          single integer). The next N lines contain the matrix A, one
!>          row per line. The next N lines contain the matrix B.  The
!>          next line contains the reciprocals of the eigenvalue
!>          condition numbers.  The last line contains the reciprocals of
!>          the eigenvector condition numbers.  The end of data is
!>          indicated by dimension N=0.  Even if no data is to be tested,
!>          there must be at least one line containing N=0.
!>
!>-----------------------------------------------------------------------
!>
!> SGX input files:
!>
!> line 1:  'SGX' in columns 1 to 3.
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
!> remaining lines : Each example is stored on 3+2*N lines, where N is
!>          its dimension. The first line contains the dimension (a
!>          single integer).  The next line contains an integer k such
!>          that only the last k eigenvalues will be selected and appear
!>          in the leading diagonal blocks of $A$ and $B$. The next N
!>          lines contain the matrix A, one row per line.  The next N
!>          lines contain the matrix B.  The last line contains the
!>          reciprocal of the eigenvalue cluster condition number and the
!>          reciprocal of the deflating subspace (associated with the
!>          selected eigencluster) condition number.  The end of data is
!>          indicated by dimension N=0.  Even if no data is to be tested,
!>          there must be at least one line containing N=0.
!>
!>-----------------------------------------------------------------------
!>
!> SSB input file:
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
!>          The 3-character path name is 'SSB'.
!>
!>-----------------------------------------------------------------------
!>
!> SBB input file:
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
!>          The 3-character path name is 'SBB'.
!>
!>-----------------------------------------------------------------------
!>
!> SEC input file:
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
!> SBL and SBK input files:
!>
!> line 1:  'SBL' in columns 1-3 to test SGEBAL, or 'SBK' in
!>          columns 1-3 to test SGEBAK.
!>
!> The remaining lines consist of specially constructed test cases.
!>
!>-----------------------------------------------------------------------
!>
!> SGL and SGK input files:
!>
!> line 1:  'SGL' in columns 1-3 to test SGGBAL, or 'SGK' in
!>          columns 1-3 to test SGGBAK.
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
!> precomputed examples, and LWORK = NMAX*(5*NMAX+5)+1 in the parameter
!> statements below.  For SVD, we assume NRHS may be as big as N.  The
!> parameter NEED is set to 14 to allow for 14 N-by-N matrices for SGG.
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
!> \ingroup single_eig
!
!  =====================================================================
      PROGRAM SCHKEE
      USE M_TST__EIG , ONLY:ALAREQ , XLAENV , SCHKBB , SCHKBD
      USE M_TST__EIG , ONLY:SCHKBK , SCHKBL , SCHKEC , SCHKGG
      USE M_TST__EIG , ONLY:SCHKGK , SCHKGL , SCHKHS , SCHKSB2STG
      USE M_TST__EIG , ONLY:SCHKST2STG , SCHKST , SCKCSD , SCKGLM
      USE M_TST__EIG , ONLY:SCKGQR , SCKGSV , SCKLSE , SDRGES3
      USE M_TST__EIG , ONLY:SDRGES , SDRGEV3 , SDRGEV , SDRGSX
      USE M_TST__EIG , ONLY:SDRGVX , SDRVBD , SDRVES , SDRVEV
      USE M_TST__EIG , ONLY:SDRVSG2STG , SDRVST2STG , SDRVST , SDRVSX
      USE M_TST__EIG , ONLY:SDRVVX , SERRBD , SERRED , SERRGG
      USE M_TST__EIG , ONLY:SERRHS , SERRST
#if defined(_OPENMP)
      USE OMP_LIB
#endif
      IMPLICIT NONE
!*--SCHKEE1055
!*** Start of declarations inserted by SPAG
!*** End of declarations inserted by SPAG
!
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2016
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER NMAX
      PARAMETER (NMAX=132)
      INTEGER NCMAX
      PARAMETER (NCMAX=20)
      INTEGER NEED
      PARAMETER (NEED=14)
      INTEGER LWORK
      PARAMETER (LWORK=NMAX*(5*NMAX+5)+1)
      INTEGER LIWORK
      PARAMETER (LIWORK=NMAX*(5*NMAX+20))
      INTEGER MAXIN
      PARAMETER (MAXIN=20)
      INTEGER MAXT
      PARAMETER (MAXT=30)
      INTEGER NIN , NOUT
      PARAMETER (NIN=5,NOUT=6)
!     ..
!     .. Local Scalars ..
      LOGICAL csd , fatal , glm , gqr , gsv , lse , nep , sbb , sbk ,   &
     &        sbl , sep , ses , sev , sgg , sgk , sgl , sgs , sgv ,     &
     &        sgx , ssb , ssx , svd , svx , sxv , tstchk , tstdif ,     &
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
      REAL d(NMAX,12) , result(500) , taua(NMAX) , taub(NMAX) ,         &
     &     x(5*NMAX)
!     ..
!     .. Allocatable Arrays ..
      INTEGER allocatestatus
      REAL , DIMENSION(:) , ALLOCATABLE  ::  work
      REAL , DIMENSION(:,:) , ALLOCATABLE  ::  a , b , c
!     ..
!     .. External Functions ..
      LOGICAL LSAMEN
      REAL SECOND , SLAMCH
      EXTERNAL LSAMEN , SECOND , SLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL SCHKSB , SDRVSG , ILAVER
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
      ALLOCATE (a(NMAX*NMAX,NEED),STAT=allocatestatus)
      IF ( allocatestatus/=0 ) STOP "*** Not enough memory ***"
      ALLOCATE (b(NMAX*NMAX,5),STAT=allocatestatus)
      IF ( allocatestatus/=0 ) STOP "*** Not enough memory ***"
      ALLOCATE (c(NCMAX*NCMAX,NCMAX*NCMAX),STAT=allocatestatus)
      IF ( allocatestatus/=0 ) STOP "*** Not enough memory ***"
      ALLOCATE (work(LWORK),STAT=allocatestatus)
      IF ( allocatestatus/=0 ) STOP "*** Not enough memory ***"
!     ..
!     .. Executable Statements ..
!
      a = 0.0
      b = 0.0
      c = 0.0
      d = 0.0
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
         nep = LSAMEN(3,path,'NEP') .OR. LSAMEN(3,path,'SHS')
         sep = LSAMEN(3,path,'SEP') .OR. LSAMEN(3,path,'SST') .OR.      &
     &         LSAMEN(3,path,'SSG') .OR. LSAMEN(3,path,'SE2')
         svd = LSAMEN(3,path,'SVD') .OR. LSAMEN(3,path,'DBD')
         svd = LSAMEN(3,path,'SVD') .OR. LSAMEN(3,path,'SBD')
         sev = LSAMEN(3,path,'SEV')
         ses = LSAMEN(3,path,'SES')
         svx = LSAMEN(3,path,'SVX')
         ssx = LSAMEN(3,path,'SSX')
         sgg = LSAMEN(3,path,'SGG')
         sgs = LSAMEN(3,path,'SGS')
         sgx = LSAMEN(3,path,'SGX')
         sgv = LSAMEN(3,path,'SGV')
         sxv = LSAMEN(3,path,'SXV')
         ssb = LSAMEN(3,path,'SSB')
         sbb = LSAMEN(3,path,'SBB')
         glm = LSAMEN(3,path,'GLM')
         gqr = LSAMEN(3,path,'GQR') .OR. LSAMEN(3,path,'GRQ')
         gsv = LSAMEN(3,path,'GSV')
         csd = LSAMEN(3,path,'CSD')
         lse = LSAMEN(3,path,'LSE')
         sbl = LSAMEN(3,path,'SBL')
         sbk = LSAMEN(3,path,'SBK')
         sgl = LSAMEN(3,path,'SGL')
         sgk = LSAMEN(3,path,'SGK')
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
         ELSEIF ( sev ) THEN
            WRITE (NOUT,FMT=99020)
         ELSEIF ( ses ) THEN
            WRITE (NOUT,FMT=99021)
         ELSEIF ( svx ) THEN
            WRITE (NOUT,FMT=99022)
         ELSEIF ( ssx ) THEN
            WRITE (NOUT,FMT=99023)
         ELSEIF ( sgg ) THEN
            WRITE (NOUT,FMT=99024)
         ELSEIF ( sgs ) THEN
            WRITE (NOUT,FMT=99035)
         ELSEIF ( sgx ) THEN
            WRITE (NOUT,FMT=99034)
         ELSEIF ( sgv ) THEN
            WRITE (NOUT,FMT=99036)
         ELSEIF ( sxv ) THEN
            WRITE (NOUT,FMT=99037)
         ELSEIF ( ssb ) THEN
            WRITE (NOUT,FMT=99025)
         ELSEIF ( sbb ) THEN
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
         ELSEIF ( sbl ) THEN
!
!        SGEBAL:  Balancing
!
            CALL SCHKBL(NIN,NOUT)
            CYCLE
         ELSEIF ( sbk ) THEN
!
!        SGEBAK:  Back transformation
!
            CALL SCHKBK(NIN,NOUT)
            CYCLE
         ELSEIF ( sgl ) THEN
!
!        SGGBAL:  Balancing
!
            CALL SCHKGL(NIN,NOUT)
            CYCLE
         ELSEIF ( sgk ) THEN
!
!        SGGBAK:  Back transformation
!
            CALL SCHKGK(NIN,NOUT)
            CYCLE
         ELSEIF ( LSAMEN(3,path,'SEC') ) THEN
!
!        SEC:  Eigencondition estimation
!
            READ (NIN,FMT=*) thresh
            CALL XLAENV(1,1)
            CALL XLAENV(12,11)
            CALL XLAENV(13,2)
            CALL XLAENV(14,0)
            CALL XLAENV(15,2)
            CALL XLAENV(16,2)
            tsterr = .TRUE.
            CALL SCHKEC(thresh,tsterr,NIN,NOUT)
            CYCLE
         ELSE
            WRITE (NOUT,FMT=99007) path
            CYCLE
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
         IF ( .NOT.(sgx .OR. sxv) ) THEN
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
         IF ( svd .OR. sbb .OR. glm .OR. gqr .OR. gsv .OR. csd .OR.     &
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
         IF ( .NOT.(sgx .OR. sxv) ) THEN
            WRITE (NOUT,FMT=99016) 'N:    ' , (nval(i),i=1,nn)
         ELSE
            WRITE (NOUT,FMT=99016) 'N:    ' , nn
         ENDIF
!
!     Read the number of values of K, followed by the values of K
!
         IF ( ssb .OR. sbb ) THEN
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
         IF ( sev .OR. ses .OR. svx .OR. ssx ) THEN
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
         ELSEIF ( sgs .OR. sgx .OR. sgv .OR. sxv ) THEN
!
!        For the nonsymmetric generalized driver routines, only one set
!        of parameters is allowed.
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
!
         ELSEIF ( .NOT.ssb .AND. .NOT.glm .AND. .NOT.gqr .AND.          &
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
            IF ( .NOT.sbb ) THEN
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
            IF ( nep .OR. sep .OR. svd .OR. sgg ) THEN
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
!        Read the values of NSHIFT (if SGG) or NRHS (if SVD
!        or SBB).
!
            IF ( svd .OR. sbb .OR. sgg ) THEN
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
            IF ( sgg ) THEN
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
            IF ( nep .OR. sgg ) THEN
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
            IF ( sgg ) THEN
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
         IF ( sep .OR. svd .OR. sgg ) THEN
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
         IF ( .NOT.(sgx .OR. sxv) ) THEN
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
            IF ( .NOT.(sev .OR. ses .OR. svx .OR. ssx .OR. sgv .OR. sgs)&
     &           .AND. ntypes<=0 ) THEN
               WRITE (NOUT,FMT=99009) c3
               GOTO 120
            ENDIF
!
         ELSE
            IF ( sxv ) c3 = 'SXV'
            IF ( sgx ) c3 = 'SGX'
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
         IF ( LSAMEN(3,c3,'SHS') .OR. LSAMEN(3,c3,'NEP') ) THEN
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
            IF ( tsterr ) CALL SERRHS('SHSEQR',NOUT)
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
               CALL SCHKHS(nn,nval,maxtyp,dotype,iseed,thresh,NOUT,     &
     &                     a(1,1),NMAX,a(1,2),a(1,3),a(1,4),a(1,5),NMAX,&
     &                     a(1,6),a(1,7),d(1,1),d(1,2),d(1,3),d(1,4),   &
     &                     d(1,5),d(1,6),a(1,8),a(1,9),a(1,10),a(1,11), &
     &                     a(1,12),d(1,7),work,LWORK,iwork,logwrk,      &
     &                     result,info)
               IF ( info/=0 ) WRITE (NOUT,FMT=99019) 'SCHKHS' , info
            ENDDO
!
         ELSEIF ( LSAMEN(3,c3,'SST') .OR. LSAMEN(3,c3,'SEP') .OR.       &
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
               CALL SERRST('SST',NOUT)
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
                     CALL SCHKST2STG(nn,nval,maxtyp,dotype,iseed,thresh,&
     &                               NOUT,a(1,1),NMAX,a(1,2),d(1,1),    &
     &                               d(1,2),d(1,3),d(1,4),d(1,5),d(1,6),&
     &                               d(1,7),d(1,8),d(1,9),d(1,10),      &
     &                               d(1,11),a(1,3),NMAX,a(1,4),a(1,5), &
     &                               d(1,12),a(1,6),work,LWORK,iwork,   &
     &                               LIWORK,result,info)
                  ELSE
                     CALL SCHKST(nn,nval,maxtyp,dotype,iseed,thresh,    &
     &                           NOUT,a(1,1),NMAX,a(1,2),d(1,1),d(1,2), &
     &                           d(1,3),d(1,4),d(1,5),d(1,6),d(1,7),    &
     &                           d(1,8),d(1,9),d(1,10),d(1,11),a(1,3),  &
     &                           NMAX,a(1,4),a(1,5),d(1,12),a(1,6),work,&
     &                           LWORK,iwork,LIWORK,result,info)
                  ENDIF
                  IF ( info/=0 ) WRITE (NOUT,FMT=99019) 'SCHKST' , info
               ENDIF
               IF ( tstdrv ) THEN
                  IF ( LSAMEN(3,c3,'SE2') ) THEN
                     CALL SDRVST2STG(nn,nval,18,dotype,iseed,thresh,    &
     &                               NOUT,a(1,1),NMAX,d(1,3),d(1,4),    &
     &                               d(1,5),d(1,6),d(1,8),d(1,9),d(1,10)&
     &                               ,d(1,11),a(1,2),NMAX,a(1,3),d(1,12)&
     &                               ,a(1,4),work,LWORK,iwork,LIWORK,   &
     &                               result,info)
                  ELSE
                     CALL SDRVST(nn,nval,18,dotype,iseed,thresh,NOUT,   &
     &                           a(1,1),NMAX,d(1,3),d(1,4),d(1,5),d(1,6)&
     &                           ,d(1,8),d(1,9),d(1,10),d(1,11),a(1,2), &
     &                           NMAX,a(1,3),d(1,12),a(1,4),work,LWORK, &
     &                           iwork,LIWORK,result,info)
                  ENDIF
                  IF ( info/=0 ) WRITE (NOUT,FMT=99019) 'SDRVST' , info
               ENDIF
            ENDDO
!
         ELSEIF ( LSAMEN(3,c3,'SSG') ) THEN
!
!        ----------------------------------------------
!        SSG:  Symmetric Generalized Eigenvalue Problem
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
!               CALL SDRVSG( NN, NVAL, MAXTYP, DOTYPE, ISEED, THRESH,
!     $                      NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), NMAX,
!     $                      D( 1, 3 ), A( 1, 3 ), NMAX, A( 1, 4 ),
!     $                      A( 1, 5 ), A( 1, 6 ), A( 1, 7 ), WORK,
!     $                      LWORK, IWORK, LIWORK, RESULT, INFO )
                  CALL SDRVSG2STG(nn,nval,maxtyp,dotype,iseed,thresh,   &
     &                            NOUT,a(1,1),NMAX,a(1,2),NMAX,d(1,3),  &
     &                            d(1,3),a(1,3),NMAX,a(1,4),a(1,5),     &
     &                            a(1,6),a(1,7),work,LWORK,iwork,LIWORK,&
     &                            result,info)
                  IF ( info/=0 ) WRITE (NOUT,FMT=99019) 'SDRVSG' , info
               ENDIF
            ENDDO
!
         ELSEIF ( LSAMEN(3,c3,'SBD') .OR. LSAMEN(3,c3,'SVD') ) THEN
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
            CALL XLAENV(1,1)
            CALL XLAENV(9,25)
!
!        Test the error exits
!
            IF ( tsterr .AND. tstchk ) CALL SERRBD('SBD',NOUT)
            IF ( tsterr .AND. tstdrv ) CALL SERRED('SBD',NOUT)
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
                  CALL SCHKBD(nn,mval,nval,maxtyp,dotype,nrhs,iseed,    &
     &                        thresh,a(1,1),NMAX,d(1,1),d(1,2),d(1,3),  &
     &                        d(1,4),a(1,2),NMAX,a(1,3),a(1,4),a(1,5),  &
     &                        NMAX,a(1,6),NMAX,a(1,7),a(1,8),work,LWORK,&
     &                        iwork,NOUT,info)
                  IF ( info/=0 ) WRITE (NOUT,FMT=99019) 'SCHKBD' , info
               ENDIF
               IF ( tstdrv ) CALL SDRVBD(nn,mval,nval,maxtyp,dotype,    &
     &              iseed,thresh,a(1,1),NMAX,a(1,2),NMAX,a(1,3),NMAX,   &
     &              a(1,4),a(1,5),a(1,6),d(1,1),d(1,2),d(1,3),work,     &
     &              LWORK,iwork,NOUT,info)
            ENDDO
!
         ELSEIF ( LSAMEN(3,c3,'SEV') ) THEN
!
!        --------------------------------------------
!        SEV:  Nonsymmetric Eigenvalue Problem Driver
!              SGEEV (eigenvalues and eigenvectors)
!        --------------------------------------------
!
            maxtyp = 21
            ntypes = MIN(maxtyp,ntypes)
            IF ( ntypes<=0 ) THEN
               WRITE (NOUT,FMT=99009) c3
            ELSE
               IF ( tsterr ) CALL SERRED(c3,NOUT)
               CALL ALAREQ(c3,ntypes,dotype,maxtyp,NIN,NOUT)
               CALL SDRVEV(nn,nval,ntypes,dotype,iseed,thresh,NOUT,     &
     &                     a(1,1),NMAX,a(1,2),d(1,1),d(1,2),d(1,3),     &
     &                     d(1,4),a(1,3),NMAX,a(1,4),NMAX,a(1,5),NMAX,  &
     &                     result,work,LWORK,iwork,info)
               IF ( info/=0 ) WRITE (NOUT,FMT=99019) 'SGEEV' , info
            ENDIF
            WRITE (NOUT,FMT=99026)
            GOTO 100
!
         ELSEIF ( LSAMEN(3,c3,'SES') ) THEN
!
!        --------------------------------------------
!        SES:  Nonsymmetric Eigenvalue Problem Driver
!              SGEES (Schur form)
!        --------------------------------------------
!
            maxtyp = 21
            ntypes = MIN(maxtyp,ntypes)
            IF ( ntypes<=0 ) THEN
               WRITE (NOUT,FMT=99009) c3
            ELSE
               IF ( tsterr ) CALL SERRED(c3,NOUT)
               CALL ALAREQ(c3,ntypes,dotype,maxtyp,NIN,NOUT)
               CALL SDRVES(nn,nval,ntypes,dotype,iseed,thresh,NOUT,     &
     &                     a(1,1),NMAX,a(1,2),a(1,3),d(1,1),d(1,2),     &
     &                     d(1,3),d(1,4),a(1,4),NMAX,result,work,LWORK, &
     &                     iwork,logwrk,info)
               IF ( info/=0 ) WRITE (NOUT,FMT=99019) 'SGEES' , info
            ENDIF
            WRITE (NOUT,FMT=99026)
            GOTO 100
!
         ELSEIF ( LSAMEN(3,c3,'SVX') ) THEN
!
!        --------------------------------------------------------------
!        SVX:  Nonsymmetric Eigenvalue Problem Expert Driver
!              SGEEVX (eigenvalues, eigenvectors and condition numbers)
!        --------------------------------------------------------------
!
            maxtyp = 21
            ntypes = MIN(maxtyp,ntypes)
            IF ( ntypes<0 ) THEN
               WRITE (NOUT,FMT=99009) c3
            ELSE
               IF ( tsterr ) CALL SERRED(c3,NOUT)
               CALL ALAREQ(c3,ntypes,dotype,maxtyp,NIN,NOUT)
               CALL SDRVVX(nn,nval,ntypes,dotype,iseed,thresh,NIN,NOUT, &
     &                     a(1,1),NMAX,a(1,2),d(1,1),d(1,2),d(1,3),     &
     &                     d(1,4),a(1,3),NMAX,a(1,4),NMAX,a(1,5),NMAX,  &
     &                     d(1,5),d(1,6),d(1,7),d(1,8),d(1,9),d(1,10),  &
     &                     d(1,11),d(1,12),result,work,LWORK,iwork,info)
               IF ( info/=0 ) WRITE (NOUT,FMT=99019) 'SGEEVX' , info
            ENDIF
            WRITE (NOUT,FMT=99026)
            GOTO 100
!
         ELSEIF ( LSAMEN(3,c3,'SSX') ) THEN
!
!        ---------------------------------------------------
!        SSX:  Nonsymmetric Eigenvalue Problem Expert Driver
!              SGEESX (Schur form and condition numbers)
!        ---------------------------------------------------
!
            maxtyp = 21
            ntypes = MIN(maxtyp,ntypes)
            IF ( ntypes<0 ) THEN
               WRITE (NOUT,FMT=99009) c3
            ELSE
               IF ( tsterr ) CALL SERRED(c3,NOUT)
               CALL ALAREQ(c3,ntypes,dotype,maxtyp,NIN,NOUT)
               CALL SDRVSX(nn,nval,ntypes,dotype,iseed,thresh,NIN,NOUT, &
     &                     a(1,1),NMAX,a(1,2),a(1,3),d(1,1),d(1,2),     &
     &                     d(1,3),d(1,4),d(1,5),d(1,6),a(1,4),NMAX,     &
     &                     a(1,5),result,work,LWORK,iwork,logwrk,info)
               IF ( info/=0 ) WRITE (NOUT,FMT=99019) 'SGEESX' , info
            ENDIF
            WRITE (NOUT,FMT=99026)
            GOTO 100
!
         ELSEIF ( LSAMEN(3,c3,'SGG') ) THEN
!
!        -------------------------------------------------
!        SGG:  Generalized Nonsymmetric Eigenvalue Problem
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
            IF ( tstchk .AND. tsterr ) CALL SERRGG(c3,NOUT)
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
                  CALL SCHKGG(nn,nval,maxtyp,dotype,iseed,thresh,tstdif,&
     &                        thrshn,NOUT,a(1,1),NMAX,a(1,2),a(1,3),    &
     &                        a(1,4),a(1,5),a(1,6),a(1,7),a(1,8),a(1,9),&
     &                        NMAX,a(1,10),a(1,11),a(1,12),d(1,1),d(1,2)&
     &                        ,d(1,3),d(1,4),d(1,5),d(1,6),a(1,13),     &
     &                        a(1,14),work,LWORK,logwrk,result,info)
                  IF ( info/=0 ) WRITE (NOUT,FMT=99019) 'SCHKGG' , info
               ENDIF
            ENDDO
!
         ELSEIF ( LSAMEN(3,c3,'SGS') ) THEN
!
!        -------------------------------------------------
!        SGS:  Generalized Nonsymmetric Eigenvalue Problem
!              SGGES (Schur form)
!        -------------------------------------------------
!
            maxtyp = 26
            ntypes = MIN(maxtyp,ntypes)
            IF ( ntypes<=0 ) THEN
               WRITE (NOUT,FMT=99009) c3
            ELSE
               IF ( tsterr ) CALL SERRGG(c3,NOUT)
               CALL ALAREQ(c3,ntypes,dotype,maxtyp,NIN,NOUT)
               CALL SDRGES(nn,nval,maxtyp,dotype,iseed,thresh,NOUT,     &
     &                     a(1,1),NMAX,a(1,2),a(1,3),a(1,4),a(1,7),NMAX,&
     &                     a(1,8),d(1,1),d(1,2),d(1,3),work,LWORK,      &
     &                     result,logwrk,info)
!
               IF ( info/=0 ) WRITE (NOUT,FMT=99019) 'SDRGES' , info
!
!     Blocked version
!
               CALL XLAENV(16,1)
               CALL SDRGES3(nn,nval,maxtyp,dotype,iseed,thresh,NOUT,    &
     &                      a(1,1),NMAX,a(1,2),a(1,3),a(1,4),a(1,7),    &
     &                      NMAX,a(1,8),d(1,1),d(1,2),d(1,3),work,LWORK,&
     &                      result,logwrk,info)
!
               IF ( info/=0 ) WRITE (NOUT,FMT=99019) 'SDRGES3' , info
            ENDIF
            WRITE (NOUT,FMT=99026)
            GOTO 100
!
         ELSEIF ( sgx ) THEN
!
!        -------------------------------------------------
!        SGX:  Generalized Nonsymmetric Eigenvalue Problem
!              SGGESX (Schur form and condition numbers)
!        -------------------------------------------------
!
            maxtyp = 5
            ntypes = maxtyp
            IF ( nn<0 ) THEN
               WRITE (NOUT,FMT=99009) c3
            ELSE
               IF ( tsterr ) CALL SERRGG(c3,NOUT)
               CALL ALAREQ(c3,ntypes,dotype,maxtyp,NIN,NOUT)
               CALL XLAENV(5,2)
               CALL SDRGSX(nn,NCMAX,thresh,NIN,NOUT,a(1,1),NMAX,a(1,2), &
     &                     a(1,3),a(1,4),a(1,5),a(1,6),d(1,1),d(1,2),   &
     &                     d(1,3),c(1,1),NCMAX*NCMAX,a(1,12),work,LWORK,&
     &                     iwork,LIWORK,logwrk,info)
               IF ( info/=0 ) WRITE (NOUT,FMT=99019) 'SDRGSX' , info
            ENDIF
            WRITE (NOUT,FMT=99026)
            GOTO 100
!
         ELSEIF ( LSAMEN(3,c3,'SGV') ) THEN
!
!        -------------------------------------------------
!        SGV:  Generalized Nonsymmetric Eigenvalue Problem
!              SGGEV (Eigenvalue/vector form)
!        -------------------------------------------------
!
            maxtyp = 26
            ntypes = MIN(maxtyp,ntypes)
            IF ( ntypes<=0 ) THEN
               WRITE (NOUT,FMT=99009) c3
            ELSE
               IF ( tsterr ) CALL SERRGG(c3,NOUT)
               CALL ALAREQ(c3,ntypes,dotype,maxtyp,NIN,NOUT)
               CALL SDRGEV(nn,nval,maxtyp,dotype,iseed,thresh,NOUT,     &
     &                     a(1,1),NMAX,a(1,2),a(1,3),a(1,4),a(1,7),NMAX,&
     &                     a(1,8),a(1,9),NMAX,d(1,1),d(1,2),d(1,3),     &
     &                     d(1,4),d(1,5),d(1,6),work,LWORK,result,info)
               IF ( info/=0 ) WRITE (NOUT,FMT=99019) 'SDRGEV' , info
!
! Blocked version
!
               CALL SDRGEV3(nn,nval,maxtyp,dotype,iseed,thresh,NOUT,    &
     &                      a(1,1),NMAX,a(1,2),a(1,3),a(1,4),a(1,7),    &
     &                      NMAX,a(1,8),a(1,9),NMAX,d(1,1),d(1,2),d(1,3)&
     &                      ,d(1,4),d(1,5),d(1,6),work,LWORK,result,    &
     &                      info)
               IF ( info/=0 ) WRITE (NOUT,FMT=99019) 'SDRGEV3' , info
            ENDIF
            WRITE (NOUT,FMT=99026)
            GOTO 100
!
         ELSEIF ( sxv ) THEN
!
!        -------------------------------------------------
!        SXV:  Generalized Nonsymmetric Eigenvalue Problem
!              SGGEVX (eigenvalue/vector with condition numbers)
!        -------------------------------------------------
!
            maxtyp = 2
            ntypes = maxtyp
            IF ( nn<0 ) THEN
               WRITE (NOUT,FMT=99009) c3
            ELSE
               IF ( tsterr ) CALL SERRGG(c3,NOUT)
               CALL ALAREQ(c3,ntypes,dotype,maxtyp,NIN,NOUT)
               CALL SDRGVX(nn,thresh,NIN,NOUT,a(1,1),NMAX,a(1,2),a(1,3),&
     &                     a(1,4),d(1,1),d(1,2),d(1,3),a(1,5),a(1,6),   &
     &                     iwork(1),iwork(2),d(1,4),d(1,5),d(1,6),d(1,7)&
     &                     ,d(1,8),d(1,9),work,LWORK,iwork(3),LIWORK-2, &
     &                     result,logwrk,info)
!
               IF ( info/=0 ) WRITE (NOUT,FMT=99019) 'SDRGVX' , info
            ENDIF
            WRITE (NOUT,FMT=99026)
            GOTO 100
!
         ELSEIF ( LSAMEN(3,c3,'SSB') ) THEN
!
!        ------------------------------
!        SSB:  Symmetric Band Reduction
!        ------------------------------
!
            maxtyp = 15
            ntypes = MIN(maxtyp,ntypes)
            CALL ALAREQ(c3,ntypes,dotype,maxtyp,NIN,NOUT)
            IF ( tsterr ) CALL SERRST('SSB',NOUT)
!         CALL SCHKSB( NN, NVAL, NK, KVAL, MAXTYP, DOTYPE, ISEED, THRESH,
!     $                NOUT, A( 1, 1 ), NMAX, D( 1, 1 ), D( 1, 2 ),
!     $                A( 1, 2 ), NMAX, WORK, LWORK, RESULT, INFO )
            CALL SCHKSB2STG(nn,nval,nk,kval,maxtyp,dotype,iseed,thresh, &
     &                      NOUT,a(1,1),NMAX,d(1,1),d(1,2),d(1,3),d(1,4)&
     &                      ,d(1,5),a(1,2),NMAX,work,LWORK,result,info)
            IF ( info/=0 ) WRITE (NOUT,FMT=99019) 'SCHKSB' , info
!
         ELSEIF ( LSAMEN(3,c3,'SBB') ) THEN
!
!        ------------------------------
!        SBB:  General Band Reduction
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
               CALL SCHKBB(nn,mval,nval,nk,kval,maxtyp,dotype,nrhs,     &
     &                     iseed,thresh,NOUT,a(1,1),NMAX,a(1,2),2*NMAX, &
     &                     d(1,1),d(1,2),a(1,4),NMAX,a(1,5),NMAX,a(1,6),&
     &                     NMAX,a(1,7),work,LWORK,result,info)
               IF ( info/=0 ) WRITE (NOUT,FMT=99019) 'SCHKBB' , info
            ENDDO
!
         ELSEIF ( LSAMEN(3,c3,'GLM') ) THEN
!
!        -----------------------------------------
!        GLM:  Generalized Linear Regression Model
!        -----------------------------------------
!
            CALL XLAENV(1,1)
            IF ( tsterr ) CALL SERRGG('GLM',NOUT)
            CALL SCKGLM(nn,mval,pval,nval,ntypes,iseed,thresh,NMAX,     &
     &                  a(1,1),a(1,2),b(1,1),b(1,2),x,work,d(1,1),NIN,  &
     &                  NOUT,info)
            IF ( info/=0 ) WRITE (NOUT,FMT=99019) 'SCKGLM' , info
!
         ELSEIF ( LSAMEN(3,c3,'GQR') ) THEN
!
!        ------------------------------------------
!        GQR:  Generalized QR and RQ factorizations
!        ------------------------------------------
!
            CALL XLAENV(1,1)
            IF ( tsterr ) CALL SERRGG('GQR',NOUT)
            CALL SCKGQR(nn,mval,nn,pval,nn,nval,ntypes,iseed,thresh,    &
     &                  NMAX,a(1,1),a(1,2),a(1,3),a(1,4),taua,b(1,1),   &
     &                  b(1,2),b(1,3),b(1,4),b(1,5),taub,work,d(1,1),   &
     &                  NIN,NOUT,info)
            IF ( info/=0 ) WRITE (NOUT,FMT=99019) 'SCKGQR' , info
!
         ELSEIF ( LSAMEN(3,c3,'GSV') ) THEN
!
!        ----------------------------------------------
!        GSV:  Generalized Singular Value Decomposition
!        ----------------------------------------------
!
            CALL XLAENV(1,1)
            IF ( tsterr ) CALL SERRGG('GSV',NOUT)
            CALL SCKGSV(nn,mval,pval,nval,ntypes,iseed,thresh,NMAX,     &
     &                  a(1,1),a(1,2),b(1,1),b(1,2),a(1,3),b(1,3),a(1,4)&
     &                  ,taua,taub,b(1,4),iwork,work,d(1,1),NIN,NOUT,   &
     &                  info)
            IF ( info/=0 ) WRITE (NOUT,FMT=99019) 'SCKGSV' , info
!
         ELSEIF ( LSAMEN(3,c3,'CSD') ) THEN
!
!        ----------------------------------------------
!        CSD:  CS Decomposition
!        ----------------------------------------------
!
            CALL XLAENV(1,1)
            IF ( tsterr ) CALL SERRGG('CSD',NOUT)
            CALL SCKCSD(nn,mval,pval,nval,ntypes,iseed,thresh,NMAX,     &
     &                  a(1,1),a(1,2),a(1,3),a(1,4),a(1,5),a(1,6),a(1,7)&
     &                  ,iwork,work,d(1,1),NIN,NOUT,info)
            IF ( info/=0 ) WRITE (NOUT,FMT=99019) 'SCKCSD' , info
!
         ELSEIF ( LSAMEN(3,c3,'LSE') ) THEN
!
!        --------------------------------------
!        LSE:  Constrained Linear Least Squares
!        --------------------------------------
!
            CALL XLAENV(1,1)
            IF ( tsterr ) CALL SERRGG('LSE',NOUT)
            CALL SCKLSE(nn,mval,pval,nval,ntypes,iseed,thresh,NMAX,     &
     &                  a(1,1),a(1,2),b(1,1),b(1,2),x,work,d(1,1),NIN,  &
     &                  NOUT,info)
            IF ( info/=0 ) WRITE (NOUT,FMT=99019) 'SCKLSE' , info
!
         ELSE
            WRITE (NOUT,FMT=*)
            WRITE (NOUT,FMT=*)
            WRITE (NOUT,FMT=99007) c3
         ENDIF
         IF ( sgx .OR. sxv ) EXIT
      ENDDO
 200  WRITE (NOUT,FMT=99005)
      s2 = SECOND()
      WRITE (NOUT,FMT=99006) s2 - s1
!
      DEALLOCATE (a,STAT=allocatestatus)
      DEALLOCATE (b,STAT=allocatestatus)
      DEALLOCATE (c,STAT=allocatestatus)
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
99013 FORMAT (' Tests of the Symmetric Eigenvalue Problem routines')
99014 FORMAT (' Tests of the Singular Value Decomposition routines')
99015 FORMAT (/' The following parameter values will be used:')
99016 FORMAT (4X,A,10I6,/10X,10I6)
99017 FORMAT (/' Routines pass computational tests if test ratio is ',  &
     &        'less than',F8.2,/)
99018 FORMAT (' Relative machine ',A,' is taken to be',E16.6)
99019 FORMAT (' *** Error code from ',A,' = ',I4)
99020 FORMAT (/' Tests of the Nonsymmetric Eigenvalue Problem Driver',  &
     &        /'    SGEEV (eigenvalues and eigevectors)')
99021 FORMAT (/' Tests of the Nonsymmetric Eigenvalue Problem Driver',  &
     &        /'    SGEES (Schur form)')
99022 FORMAT (/' Tests of the Nonsymmetric Eigenvalue Problem Expert',  &
     &        ' Driver',/'    SGEEVX (eigenvalues, eigenvectors and',   &
     &        ' condition numbers)')
99023 FORMAT (/' Tests of the Nonsymmetric Eigenvalue Problem Expert',  &
     &        ' Driver',/'    SGEESX (Schur form and condition',        &
     &        ' numbers)')
99024 FORMAT (/' Tests of the Generalized Nonsymmetric Eigenvalue ',    &
     &        'Problem routines')
99025 FORMAT (' Tests of SSBTRD',/' (reduction of a symmetric band ',   &
     &        'matrix to tridiagonal form)')
99026 FORMAT (/1X,71('-'))
99027 FORMAT (/' LAPACK VERSION ',I1,'.',I1,'.',I1)
99028 FORMAT (/' Tests of the Generalized Linear Regression Model ',    &
     &        'routines')
99029 FORMAT (/' Tests of the Generalized QR and RQ routines')
99030 FORMAT (/' Tests of the Generalized Singular Value',              &
     &        ' Decomposition routines')
99031 FORMAT (/' Tests of the Linear Least Squares routines')
99032 FORMAT (' Tests of SGBBRD',/' (reduction of a general band ',     &
     &        'matrix to real bidiagonal form)')
99033 FORMAT (//1X,A3,':  NRHS =',I4)
99034 FORMAT (/' Tests of the Generalized Nonsymmetric Eigenvalue ',    &
     &        'Problem Expert Driver SGGESX')
99035 FORMAT (/' Tests of the Generalized Nonsymmetric Eigenvalue ',    &
     &        'Problem Driver SGGES')
99036 FORMAT (/' Tests of the Generalized Nonsymmetric Eigenvalue ',    &
     &        'Problem Driver SGGEV')
99037 FORMAT (/' Tests of the Generalized Nonsymmetric Eigenvalue ',    &
     &        'Problem Expert Driver SGGEVX')
99038 FORMAT (//1X,A3,':  NB =',I4,', NBMIN =',I4,', NX =',I4,          &
     &        ', INMIN=',I4,', INWIN =',I4,', INIBL =',I4,', ISHFTS =', &
     &        I4,', IACC22 =',I4)
99039 FORMAT (/' Tests of the CS Decomposition routines')
!
!     End of SCHKEE
!
      END PROGRAM SCHKEE
