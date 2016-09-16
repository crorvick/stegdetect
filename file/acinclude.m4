dnl from autoconf 2.13 acgeneral.m4, with patch:
dnl Date: Fri, 15 Jan 1999 05:52:41 -0800
dnl Message-ID: <199901151352.FAA18237@shade.twinsun.com>
dnl From: eggert@twinsun.com (Paul Eggert)
dnl Subject: autoconf 2.13 AC_CHECK_TYPE doesn't allow shell vars
dnl Newsgroups: gnu.utils.bug

dnl from autoconf 2.13 acspecific.m4, with changes to check for daylight

AC_DEFUN(AC_STRUCT_TIMEZONE_DAYLIGHT,
[AC_REQUIRE([AC_STRUCT_TM])dnl
AC_CACHE_CHECK([for tm_zone in struct tm], ac_cv_struct_tm_zone,
[AC_TRY_COMPILE([#include <sys/types.h>
#include <$ac_cv_struct_tm>], [struct tm tm; tm.tm_zone;],
  ac_cv_struct_tm_zone=yes, ac_cv_struct_tm_zone=no)])
if test "$ac_cv_struct_tm_zone" = yes; then
  AC_DEFINE(HAVE_TM_ZONE, 1, [Define if you have tm_size])
fi
AC_CACHE_CHECK(for tzname, ac_cv_var_tzname,
[AC_TRY_LINK(
changequote(<<, >>)dnl
<<#include <time.h>
#ifndef tzname /* For SGI.  */
extern char *tzname[]; /* RS6000 and others reject char **tzname.  */
#endif>>,
changequote([, ])dnl
[atoi(*tzname);], ac_cv_var_tzname=yes, ac_cv_var_tzname=no)])
  if test $ac_cv_var_tzname = yes; then
    AC_DEFINE(HAVE_TZNAME, 1, [Define if you have tzname])
  fi

AC_CACHE_CHECK([for tm_isdst in struct tm], ac_cv_struct_tm_isdst,
[AC_TRY_COMPILE([#include <sys/types.h>
#include <$ac_cv_struct_tm>], [struct tm tm; tm.tm_isdst;],
  ac_cv_struct_tm_isdst=yes, ac_cv_struct_tm_isdst=no)])
if test "$ac_cv_struct_tm_isdst" = yes; then
  AC_DEFINE(HAVE_TM_ISDST, 1, [Define if you have tm_isdst])
fi
AC_CACHE_CHECK(for daylight, ac_cv_var_daylight,
[AC_TRY_LINK(
changequote(<<, >>)dnl
<<#include <time.h>
#ifndef daylight /* In case IRIX #defines this, too  */
extern int daylight;
#endif>>,
changequote([, ])dnl
[atoi(daylight);], ac_cv_var_daylight=yes, ac_cv_var_daylight=no)])
  if test $ac_cv_var_daylight = yes; then
    AC_DEFINE(HAVE_DAYLIGHT, 1, [Define if you have daylight])
  fi
])

dnl AC_CHECK_TYPE2(TYPE, DEFAULT)
AC_DEFUN(AC_CHECK_TYPE2,
[AC_REQUIRE([AC_HEADER_STDC])dnl
AC_MSG_CHECKING(for $1)
AC_CACHE_VAL(ac_cv_type_$1,
[AC_EGREP_CPP(dnl
changequote(<<,>>)dnl
<<(^|[^a-zA-Z_0-9])$1[^a-zA-Z_0-9]>>dnl
changequote([,]), [#include <sys/types.h>
#if STDC_HEADERS
#include <stdlib.h>
#include <stddef.h>
#endif], eval "ac_cv_type_$1=yes", eval "ac_cv_type_$1=no")])dnl
if eval "test \"`echo '$ac_cv_type_'$1`\" = yes"; then
  AC_MSG_RESULT(yes)
else
  AC_MSG_RESULT(no)
  AC_DEFINE_UNQUOTED($1, $2)
fi
])
