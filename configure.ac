# Configure.ac file
AC_PREREQ([2.69])
AC_INIT([FULL-PACKAGE-NAME], [VERSION], [BUG-REPORT-ADDRESS])
AM_INIT_AUTOMAKE(foreign)
AC_CONFIG_SRCDIR([src/main.cpp])
AC_CONFIG_HEADERS([config.h])
AM_MAINTAINER_MODE([disable])
AC_PREFIX_DEFAULT(/usr/local/)

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS([fcntl.h inttypes.h limits.h stdlib.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_SUBST(AM_CXXFLAGS, "-std=c++11 -g -Wall -O3 -Wextra -Wno-unknown-pragmas -Wno-unused-parameter -Wno-sign-compare")
AC_SUBST(LDFLAGS, "$LDFLAGS -pthread -lz")
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_MODE_T
AC_TYPE_SIZE_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_MMAP
AC_CHECK_FUNCS([mkdir munmap realpath])

AC_CONFIG_FILES([Makefile
                 src/Makefile])
AC_OUTPUT
