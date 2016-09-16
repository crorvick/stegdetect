dnl Process this file with autoconf to produce a configure script.
AC_INIT([file.c], [3.37])
AM_INIT_AUTOMAKE
AC_CONFIG_HEADER(config.h)

AC_MSG_CHECKING(for builtin ELF support)
AC_ARG_ENABLE(elf,
[  --enable-elf            disable builtin ELF support],
[if test "${enableval}" = yes; then
  AC_MSG_RESULT(yes)
  AC_DEFINE(BUILTIN_ELF, 1, [Disable builtin ELF support])
else
  AC_MSG_RESULT(no)
fi], [
  # enable by default
  AC_MSG_RESULT(NO)
])

AC_MSG_CHECKING(for ELF core file support)
AC_ARG_ENABLE(elf-core,
[  --disable-elf-core       disable ELF core file support],
[if test "${enableval}" = yes; then
  AC_MSG_RESULT(yes)
  AC_DEFINE(ELFCORE, 1, [Disable ELF core file support])
else
  AC_MSG_RESULT(no)
fi], [
  # enable by default
  AC_MSG_RESULT(yes)
  AC_DEFINE(ELFCORE, 1, [Disable ELF core file support])
])

AC_MSG_CHECKING(for file formats in man section 5)
AC_ARG_ENABLE(fsect-man5,
[  --enable-fsect-man5      enable file formats in man section 5],
[if test "${enableval}" = yes; then
  AC_MSG_RESULT(yes)
  fsect=5
else
  AC_MSG_RESULT(no)
  fsect=4
fi], [
  # disable by default
  AC_MSG_RESULT(no)
  fsect=4
])
AC_SUBST(fsect)
AM_CONDITIONAL(FSECT5, test x$fsect = x5)

dnl Checks for programs.
AC_PROG_CC
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_RANLIB

dnl Checks for headers
AC_HEADER_STDC
AC_HEADER_MAJOR
AC_HEADER_SYS_WAIT
AC_CHECK_HEADERS(unistd.h)
AC_CHECK_HEADERS(locale.h)
AC_CHECK_HEADERS(sys/mman.h)

dnl Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_TYPE_OFF_T
AC_TYPE_SIZE_T
AC_TYPE_UINT8_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T
AC_STRUCT_ST_RDEV
AC_STRUCT_TIMEZONE_DAYLIGHT

dnl Checks for functions
AC_CHECK_FUNCS(mmap strerror strtoul)

AC_OUTPUT(Makefile)
