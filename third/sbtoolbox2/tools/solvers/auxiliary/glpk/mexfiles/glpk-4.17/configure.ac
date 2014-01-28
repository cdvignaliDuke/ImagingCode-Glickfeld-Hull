dnl Process this file with autoconf to produce a configure script.

AC_INIT([GLPK],[4.17],[bug-glpk@gnu.org])
AC_CONFIG_SRCDIR([include/glpk.h])
AM_INIT_AUTOMAKE
AC_CONFIG_HEADERS([config.h])

AC_ARG_ENABLE(gmp,
AC_HELP_STRING([--enable-gmp],[enable using GNU MP library [[default=yes]]]),
[case $enableval in
yes|no) ;;
*) AC_MSG_ERROR([bad value $enableval for --enable-gmp, need yes or no])
esac],
[enable_gmp=yes])

dnl Check for programs.
AC_PROG_CC
AC_PROG_INSTALL

dnl Check for libtool and disable unnecessary libtool tests.
define([AC_LIBTOOL_LANG_CXX_CONFIG],[:])
define([AC_LIBTOOL_LANG_F77_CONFIG],[:])
define([AC_LIBTOOL_LANG_GCJ_CONFIG],[:])
AC_PROG_LIBTOOL

dnl Check for libraries.
AC_CHECK_LIB([m],[sqrt])

if test "$enable_gmp" = "yes"; then
AC_CHECK_HEADERS([gmp.h])
if test "$ac_cv_header_gmp_h" = "yes"; then
LIBS="-lgmp $LIBS"
else
AC_MSG_WARN([GNU MP not installed.])
AC_MSG_NOTICE([*** To attain a much better performance on using exact arithmetic])
AC_MSG_NOTICE([*** it is highly recommended to install the GNU MP library and then])
AC_MSG_NOTICE([*** repeat configuring and compiling GLPK.])
fi
fi

AC_CONFIG_FILES([include/Makefile src/Makefile examples/Makefile Makefile])
AC_OUTPUT
