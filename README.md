# EXPERIMENT
  clone of LAPACK being run through plusFORT spag utility

  Not in production or tested

  See the many other optimized versions of LAPACK or the reference
  sites for a production version of LAPACK.

  This has other uses.

  * package is so big it hits limits in EXECUTE_COMMAND_LINE(3f) on some platforms when using fpm(1).
    using a response file or breaking the ar(1) command are possible solutions.

  * had trouble getting the test programs to load, left the main programs in test/ but made a package
    out of the routines and modules and placed them in test_aux and used them as a local dependency 
    without a separate .git history.

    Cannot have routines outside of a module in test/ directory be found?

  * even with auto discovery on had to put one test program into the fpm.toml file or it would not
    allow a test.dependency line; might be a syntax that works without doing that but oould not come
    up with one yet. Need to check the code.

  * to run the codes with just a simple "fpm test" need to change the file opens. Would be nice if
    there were a way to indicate the test commands.

  * there were originally a lot of duplicates which caused load problems, but changed test procedures
    into modules, leaving the files as individual files but using INCLUDE to make them easier to 
    manage; build does not currently detect changes in include files unless remove build directory
    or make a trivial change to the files that do the INCLUDE to rebuild after changing the included
    files.

  * the test programs had a lot of real and complex values passed where there should have been complex
    and double complex, double values passed where real should have been passed and so on. I think I
    got all of them in the test/ and test_aux/ code.

  * because the actual library routines do non-standard things like passing arrays without regard to
    shape might take changes in the interface beyond just a 'USE M_LAPACK' and 'USE_BLAS', might be
    more practical to make a wrapper module and make names generic via that and so on.

  * test/LAPACK/EIG/schkec.f90(192): error #8284: If the actual argument
  is scalar, the dummy argument shall be scalar unless the actual argument
  is of type character or is an element of an array that is not assumed
  shape, pointer, or polymorphic.   [NINFO]

      CALL SGET40(rtgexc,ltgexc,ntgexc,ktgexc,Nin)
      -----------^

## TODO
  * need to verify all tests and make a script to run them properly to get a unit test running that
    does not require additional intrastructure
  * for now, no XBLAS, no C/C++ interfaces
  * need to make tests for OpenMP version
  * set up to generate Doxygen documentation and generate man-pages and HTML documents
  * perhaps an fman(1) for the LAPACK routines; seperate test routines from standard library routines
  * make optional build for built-in timing

# LAPACK
  From the original documentation ...

[![Build Status](https://travis-ci.org/Reference-LAPACK/lapack.svg?branch=master)](https://travis-ci.org/Reference-LAPACK/lapack)
[![Appveyor](https://ci.appveyor.com/api/projects/status/bh38iin398msrbtr?svg=true)](https://ci.appveyor.com/project/langou/lapack/)
[![codecov](https://codecov.io/gh/Reference-LAPACK/lapack/branch/master/graph/badge.svg)](https://codecov.io/gh/Reference-LAPACK/lapack)
[![Packaging status](https://repology.org/badge/tiny-repos/lapack.svg)](https://repology.org/metapackage/lapack/versions)


* VERSION 1.0   :  February 29, 1992
* VERSION 1.0a  :  June 30, 1992
* VERSION 1.0b  :  October 31, 1992
* VERSION 1.1   :  March 31, 1993
* VERSION 2.0   :  September 30, 1994
* VERSION 3.0   :  June 30, 1999
* VERSION 3.0 + update :  October 31, 1999
* VERSION 3.0 + update :  May 31, 2000
* VERSION 3.1   : November 2006
* VERSION 3.1.1 : February 2007
* VERSION 3.2   : November 2008
* VERSION 3.2.1 : April 2009
* VERSION 3.2.2 : June 2010
* VERSION 3.3.0 : November 2010
* VERSION 3.3.1 : April 2011
* VERSION 3.4.0 : November 2011
* VERSION 3.4.1 : April 2012
* VERSION 3.4.2 : September 2012
* VERSION 3.5.0 : November 2013
* VERSION 3.6.0 : November 2015
* VERSION 3.6.1 : June 2016
* VERSION 3.7.0 : December 2016
* VERSION 3.7.1 : June 2017
* VERSION 3.8.0 : November 2017
* VERSION 3.9.0 : November 2019

LAPACK is a library of Fortran subroutines for solving the most commonly
occurring problems in numerical linear algebra.

LAPACK is a freely-available software package. It can be included in commercial
software packages (and has been). We only ask that that proper credit be given
to the authors, for example by citing the LAPACK Users' Guide. The license used
for the software is the [modified BSD license](https://github.com/Reference-LAPACK/lapack/blob/master/LICENSE).

Like all software, it is copyrighted. It is not trademarked, but we do ask the
following: if you modify the source for these routines we ask that you change
the name of the routine and comment the changes made to the original.

We will gladly answer any questions regarding the software. If a modification
is done, however, it is the responsibility of the person who modified the
routine to provide support.

LAPACK is [available from GitHub](https://github.com/Reference-LAPACK/lapack).
LAPACK releases are also [available on netlib](http://www.netlib.org/lapack/).

The distribution contains (1) the Fortran source for LAPACK, and (2) its
testing programs.  It also contains (3) the Fortran reference implementation of
the Basic Linear Algebra Subprograms (the Level 1, 2, and 3 BLAS) needed by
LAPACK.  However this code is intended for use only if there is no other
implementation of the BLAS already available on your machine; the efficiency of
LAPACK depends very much on the efficiency of the BLAS.  It also contains (4)
CBLAS, a C interface to the BLAS, and (5) LAPACKE, a C interface to LAPACK.

## Installation

 - LAPACK can be installed with `make`. The configuration must be set in the
   `make.inc` file. A `make.inc.example` for a Linux machine running GNU compilers
   is given in the main directory. Some specific `make.inc` are also available in
   the `INSTALL` directory.
 - LAPACK includes also the [CMake](https://cmake.org/) build.  You will need
   to have CMake installed on your machine.  CMake will allow an easy
   installation on a Windows Machine.  An example CMake build to install the
   LAPACK library under `$HOME/.local/lapack/` is:
   ```sh
   mkdir build
   cd build
   cmake -DCMAKE_INSTALL_LIBDIR=$HOME/.local/lapack ..
   cmake --build . -j --target install
   ```


## User Support

LAPACK has been thoroughly tested, on many different types of computers. The
LAPACK project supports the package in the sense that reports of errors or poor
performance will gain immediate attention from the developers. Such reports,
descriptions of interesting applications, and other comments should be sent by
email to [the LAPACK team](mailto:lapack@icl.utk.edu).

A list of known problems, bugs, and compiler errors for LAPACK is
[maintained on netlib](http://www.netlib.org/lapack/release_notes.html).
Please see as well the [GitHub issue tracker](https://github.com/Reference-LAPACK/lapack/issues).

For further information on LAPACK please read our [FAQ](http://www.netlib.org/lapack/faq.html)
and [Users' Guide](http://www.netlib.org/lapack/lug/lapack_lug.html).
A [user forum](http://icl.cs.utk.edu/lapack-forum/) and specific information for
[running LAPACK under Windows](http://icl.cs.utk.edu/lapack-for-windows/lapack/).
is also available to help you with the LAPACK library.


## Testing

LAPACK includes a thorough test suite. We recommend that, after compilation,
you run the test suite.

For complete information on the LAPACK Testing please consult LAPACK Working
Note 41 "Installation Guide for LAPACK".


## LAPACKE

LAPACK now includes the [LAPACKE](http://www.netlib.org/lapack/lapacke.html)
package.  LAPACKE is a Standard C language API for LAPACK that was born from a
collaboration of the LAPACK and INTEL Math Kernel Library teams.
