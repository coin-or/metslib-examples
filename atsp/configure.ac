dnl --------------------------------
dnl Initialization macros.
dnl --------------------------------

AC_INIT(tsp, 0.4.4, mirko.maischberger@gmail.com)
AM_CONFIG_HEADER(src/config.h)
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_AUX_DIR([config])

dnl -----------------------------------------------
dnl Package name and version number (user defined)
dnl -----------------------------------------------

AM_INIT_AUTOMAKE(tsp, 0.4.4, mirko.maischberger@gmail.com)

dnl -----------------------------------------------
dnl Checks for programs.
dnl -----------------------------------------------

AC_PROG_CXX
AM_PROG_LIBTOOL
AM_SANITY_CHECK

dnl -----------------------------------------------
dnl Checks for libraries.
dnl -----------------------------------------------

AC_LANG_CPLUSPLUS

PKG_CHECK_MODULES(metslib, metslib >= 0.4.4)
AC_SUBST(metslib_CFLAGS)
AC_SUBST(metslib_LIBS)


dnl ---------------------------------------------
dnl Turn off optimizations on demand
dnl ---------------------------------------------
release_version=yes
AC_ARG_ENABLE(release, [  --enable-debug        Turn off optimizations.  ], \
 if test "$enableval" = "yes"; then \
	release_version=no
fi;)	

if test $release_version = no; then 
  CFLAGS=`echo $CFLAGS | sed s/-O2//` -g -fno-default-inline
  CXXFLAGS=`echo $CXXFLAGS | sed s/-O2//` -g -fno-default-inline
  CPPFLAGS=`echo $CPPFLAGS | sed s/-O2//` -g -fno-default-inline
  CFLAGS=`echo $CFLAGS | sed s/-O3//`
  CXXFLAGS=`echo $CXXFLAGS | sed s/-O3//`
  CPPFLAGS=`echo $CPPFLAGS | sed s/-O3//`
fi

dnl ---------------------------------------------
dnl g++ specific options
dnl ---------------------------------------------

if test $CXX = g++; then
  CXXFLAGS="$CXXFLAGS -Wall -Weffc++"
  if test $release_version = yes; then
    CXXFLAGS="$CXXFLAGS -O3"
  fi
fi

dnl -----------------------------------------------
dnl Generates Makefiles.
dnl -----------------------------------------------

AC_OUTPUT(Makefile \
	  src/Makefile 
)

