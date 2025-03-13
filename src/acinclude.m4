
dnl check for CUDA
dnl
AC_DEFUN([TC_CUDA],
[AC_CACHE_CHECK([whether CUDA is available], tc_cv_cuda_dir,
 [AC_ARG_WITH([cuda], AS_HELP_STRING([--with-cuda=DIR],
 	[specify the CUDA toolkit path]),
	[tc_cv_cuda_dir="$with_cuda"])
  HAVE_CUDA="${tc_cv_cuda_dir:-no}"
  if test "$HAVE_CUDA" != no
  then
    HAVE_CUDA="yes"
    CUDA_DIR=$tc_cv_cuda_dir
    CUDA_INCLUDES="-I$tc_cv_cuda_dir/include -I$tc_cv_cuda_dir/src"
    case $host in
      x86_64-*-*) CUDA_LIBS="-L$tc_cv_cuda_dir/lib64 -Wl,-rpath,$tc_cv_cuda_dir/lib64 -lcublas -lcusparse -lcurand -lcudart -lrt" ;;
      i?86-*-*) CUDA_LIBS="-L$tc_cv_cuda_dir/lib32 -Wl,-rpath,$tc_cv_cuda_dir/lib32 -lcublas -lcusparse -lcurand -lcudart -lrt" ;;
      *) tc_cv_cuda_dir="no"; HAVE_CUDA="no"; CUDA_LIBDIR= ; CUDA_INCLUDEDIR= ;;
    esac
    AC_SUBST([CUDA_DIR])
    AC_SUBST([CUDA_LIBS])
    AC_SUBST([CUDA_INCLUDES])
    AC_SUBST([HAVE_CUDA])
  fi
 ])dnl
])dnl

dnl decide thread model
dnl
AC_DEFUN([TC_THREAD_LIBRARY],
[AC_CACHE_CHECK([thread runtime library to be used (intel or gnu)], tc_cv_thread_library,
 [AC_ARG_WITH([thread-library], AS_HELP_STRING([--with-thread-library=intel|gnu],
   [specify thread runtime library to be used]),
   [tc_cv_thread_library="$with_thread_library"])
 THREAD_LIBRARY="${tc_cv_thread_library:-gnu}"
 ])dnl
 AC_SUBST([THREAD_LIBRARY])
])dnl




dnl check for MKL
dnl
AC_DEFUN([TC_MKL],
[AC_REQUIRE([TC_THREAD_LIBRARY])dnl
 AC_CACHE_CHECK([whether MKL is available], tc_cv_mkl_dir,
 [AC_ARG_WITH([mkl], AS_HELP_STRING([--with-mkl=DIR],
 	[specify the MKL installation path]),
	[tc_cv_mkl_dir="$with_mkl"])
  HAVE_MKL="${tc_cv_mkl_dir:-no}"
  if test "$HAVE_MKL" != no
  then
    HAVE_MKL="yes"
    MKL_DIR=$tc_cv_mkl_dir
    MKL_INCLUDES="-I$tc_cv_mkl_dir/include"
    case $host in
      x86_64-*-*)
        if test "$THREAD_LIBRARY" == intel
        then
          MKL_LIBS="-L$tc_cv_mkl_dir/lib/intel64 -Wl,-rpath,$tc_cv_mkl_dir/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5  -lm -lpthread"
        else
          MKL_LIBS="-L$tc_cv_mkl_dir/lib/intel64 -Wl,-rpath,$tc_cv_mkl_dir/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp  -lm -lpthread"
        fi ;;
      *)
       AC_MSG_NOTICE(["Cannot use MKL for this architecture"])
       tc_cv_mkl_dir="no"; HAVE_MKL="no"; MKL_LIBDIR= ; MKL_INCLUDEDIR= ;
      ;;
    esac
    AC_SUBST([MKL_LIBS])
    AC_SUBST([MKL_DIR])
    AC_SUBST([MKL_INCLUDES])
    AC_SUBST([HAVE_MKL])
  fi
 ])dnl
])dnl


dnl check for lapack
dnl
AC_DEFUN([TC_LAPACK],
[AC_CACHE_CHECK([for lapack libraries], tc_cv_lapack_libs,
 [AC_ARG_WITH([lapack], AS_HELP_STRING([--with-lapack=libs],
 	[specify the lapack libraries]),
	[tc_cv_lapack_libs="$with_lapack"])
    LAPACK_LIBS=$tc_cv_lapack_libs
    AC_SUBST([LAPACK_LIBS])
 ])dnl
])dnl

dnl which SVN version to use
dnl
AC_DEFUN([TC_SUBVERSION],
[tc_subversion="svn"
 AC_ARG_WITH([subversion], AS_HELP_STRING([--with-subversion=PROG],
 	[specify the subversion executable]),
	[tc_subversion=$with_subversion])
 if test ! -x $tc_subversion
 then
  AC_PATH_PROG([SVN], [$tc_subversion])
 else
  SVN=$tc_subversion
  AC_SUBST([SVN])
 fi
])dnl


