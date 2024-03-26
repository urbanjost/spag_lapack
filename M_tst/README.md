# test_aux

Modules and procedures needed by spag_lapack test
programs built as a package as such files placed
in the test/ directory will not be loaded by the
test programs at this time.

Building the test dependencies as an internal
package resolves the loader issues.
