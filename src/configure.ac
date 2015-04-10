AC_INIT([s2-geometry-library], [1.0], [opensource@google.com])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([dist-bzip2 foreign silent-rules subdir-objects])

AM_SILENT_RULES([yes])

AC_PROG_CXX

LT_INIT([disable-shared])
PKG_PROG_PKG_CONFIG

AC_LANG([C++])

# https://code.google.com/p/gflags/
PKG_CHECK_MODULES([GFLAGS], [libgflags], [dnl
  CPPFLAGS="$CPPFLAGS $GFLAGS_CFLAGS"
  LIBS="$LIBS $GFLAGS_LIBS"
])

# https://code.google.com/p/google-glog/
PKG_CHECK_MODULES([GLOG], [libglog], [dnl
  CPPFLAGS="$CPPFLAGS $GLOG_CFLAGS"
  LIBS="$LIBS $GLOG_LIBS"
])

# http://www.openssl.org/
PKG_CHECK_MODULES([LIBCRYPTO], [libcrypto], [dnl
  CPPFLAGS="$CPPFLAGS $LIBCRYPTO_CFLAGS"
  LIBS="$LIBS $LIBCRYPTO_LIBS"
])

# http://swtch.com/plan9port/unix/
# Does not use pkg-config.
AC_CHECK_HEADERS([utf.h], [],
                 AC_MSG_ERROR([libutf is required]))

AX_C___ATTRIBUTE__
if test x"$ac_cv___attribute__" = x"yes"; then
   AC_SUBST(ac_google_attribute, 1)
else
   AC_SUBST(ac_google_attribute, 0)
fi

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
