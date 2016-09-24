# Process this file with autoconf to produce a configure script.
AC_INIT([mpibwa], [1.0])
AM_INIT_AUTOMAKE([foreign -Wall])

# Checks for programs.
AC_PROG_CC([mpicc])
AC_PROG_RANLIB
AM_PROG_AR

# Checks for libraries.
AC_CHECK_LIB([m], [sqrt])
AC_CHECK_LIB([z], [gzdopen])

# Check for large-file support.
AC_SYS_LARGEFILE

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
