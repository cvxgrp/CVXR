dnl R interface to ecos Embedded Conic Solver routines
dnl 
dnl This configure inspired and hacked from Rcsdp
dnl
dnl By default the included ecos code is built and linked, otherwise the user can supply
dnl the location of the ecos library and includes, or the location of the ecos installation directory
dnl
dnl 1) If the user sets PKG_CFLAGS, PKG_CPPFLAGS and PKG_LIBS (all three) then these are used
dnl
dnl 2) If user provides location of the csdp installation we use the lib and include directories there
dnl     --with-ecos=<location of ecos DIR>
dnl 3) If none of these are given the included ecos code is built and linked
dnl
dnl Compiler flags for ecos can be given using --with-ecos-cflags

AC_INIT(DESCRIPTION)

FOUND_ALL=""

dnl Check if global variables are already set, if so we are done
if test [-n "$PKG_CFLAGS" -a -n "$PKG_CPPFLAGS" -a -n "$PKG_LIBS"] ; then
    FOUND_ALL="yes"
else
    if test [-n "$PKG_CFLAGS" -o -n "$PKG_CPPFLAGS" -o -n "$PKG_LIBS"] ; then
	AC_MSG_ERROR([All three of PKG_CFLAGS PKG_CPPFLAGS and PKG_LIBS must be given together])
    fi
fi

dnl get arguments
AC_ARG_WITH([ecos-cflags],
    [AS_HELP_STRING([--with-ecos-cflags],[compiler flags for ecos])],
    ECOS_CFLAGS=$withval,
    ECOS_CFLAGS="")

AC_ARG_WITH([ecos-dir],
    [AS_HELP_STRING([--with-ecos-dir],[specify location of ecos installation])],
    ECOS_DIR=$withval,
    ECOS_DIR="")

if test [-z "$FOUND_ALL" -a -n "$ECOS_DIR"] ; then
    if test [! -d "$ECOS_DIR"] ; then
	AC_MSG_ERROR([ECOS installation directory $ECOS_DIR does not exist])
    fi

    PKG_CPPFLAGS="-I$ECOS_DIR/include -I$ECOS_DIR/external/ldl/include -I$ECOS_DIR/external/amd/include -I$ECOS_DIR/external/SuiteSparse_config"
    PKG_LIBS="-L$ECOS_DIR"
    FOUND_ALL="yes"
fi

dnl No location given so build and link provided ecos
if test [-z "$FOUND_ALL"] ; then
    echo "Using bundled ECOS"
    PKG_CPPFLAGS="-Iecos/include -Iecos/external/ldl/include -Iecos/external/amd/include -Iecos/external/SuiteSparse_config"
    PKG_LIBS="-Lecos -lecos -lecos_bb"
    ECOS_TS="ecos.ts"
else
    ECOS_TS=""
fi

dnl Find the R home directory.
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "Could not determine R_HOME."
  exit 1
fi

CC=`"${R_HOME}/bin/R" CMD config CC`
CPP=`"${R_HOME}/bin/R" CMD config CPP`
CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`
CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`

AC_PROG_CC
AC_PROG_CPP

PKG_LIBS="${PKG_LIBS} -lecos -lecos_bb"
PKG_CFLAGS="${ECOS_CFLAGS}"

dnl if given by user, check header and library
if test [-n "$FOUND_ALL"] ; then
    LIBS="${PKG_LIBS}"
    CFLAGS="${CFLAGS} ${PKG_CFLAGS}"
    CPPFLAGS="${CPPFLAGS} ${PKG_CPPFLAGS}"
    
    AC_CHECK_HEADER(ecos.h)
    if test "${ac_cv_header_ecos_h}" == no; then
	AC_MSG_ERROR([ECOS header ecos.h not found in provided location])
    fi
    
    AC_SEARCH_LIBS(ECOS_setup,ecos,,AC_MSG_ERROR([Could not link ECOS library using $PKG_LIBS]))
fi

dnl substitute in src/Makevars.in
AC_SUBST(PKG_CFLAGS)
AC_SUBST(PKG_CPPFLAGS)
AC_SUBST(PKG_LIBS)
AC_SUBST(ECOS_TS)
AC_OUTPUT(src/Makevars)

exit 0

