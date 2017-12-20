AC_PREREQ(2.62)
AC_INIT([DESCRIPTION])

: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi

LTO=-flto
CXX=$(${R_HOME}/bin/R CMD config CXX) $(LTO)

AC_LANG(C++)
AC_REQUIRE_CPP
AC_PROG_CXX

AC_ARG_ENABLE([lto],
[AS_HELP_STRING([--enable-lto],[enable link-time optimization @<:@yes@:>@])],
[want_lto="${enableval}"], [want_lto=yes])
if test "x${want_lto}" != xno; then
  AC_LINK_IFELSE([AC_LANG_PROGRAM([[]], [[int main (void) { return 0; }]])],
  [AC_MSG_RESULT(yes)], [AC_MSG_RESULT(no) LTO=""])
fi

if test "${GXX}" = yes && test -z "${LTO}"; then
  gxx_version=`${CXX} -v 2>&1 | grep "^.*g.. version" | sed -e 's/^.*g.. version *//'`
  case ${gxx_version} in
    1.*|2.*|3.*|4.0.*|4.1.*|4.2.*|4.3.*|4.4.*|4.5.*|4.6.*|4.7.*|4.8.*)
      AC_MSG_WARN([It is preferable to use g++ version 4.9 or greater with rstanarm.])
   ;;
   *)
      LTO=-flto=jobserver
  esac
fi

if [[ $CXX == *"clang"* ]] && test -z "${LTO}"; then
  clang_version=`${CXX} --version 2>&1 | grep "version"`
  if [[ $clang_version == *"Apple"* ]]; then
    clang_version=`${CXX} --version 2>&1 | grep "version"` | \
     sed -e 's/^Apple LLVM version *//'
     case ${clang_version} in 2.*|3.*|4.*5.*|6.*|7.*|8.*)
      AC_MSG_WARN([It is preferable to use Apple LLVM version 9 or greater with rstanarm.])
     ;;
     *)
       LTO="-flto=thin"
    esac
  fi
  elif
    case ${clang_version} in 1.*|2.*|3.1.*|3.2.*|3.3.*|3.4.*|3.5.*|3.6.*|3.7.*|3.8.*)
      AC_MSG_WARN([It is preferable to use clang++ version 3.9 or greater with rstanarm.])
     ;;
     *)
       LTO="-flto=thin"
     esac
  fi
fi

AC_SUBST(LTO)
AM_CONDITIONAL(BUILD_LTO, [test "x${want_lto}" != xno])

AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
