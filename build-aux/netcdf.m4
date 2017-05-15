dnl
dnl AC_PATH_NETCDF()
dnl
dnl Defines NETCDF_INCDIR, NETCDF_LIBDIR, and NETCDF_LIBNAME, or exits with an error message.
dnl I.e., the intent is to have this used when the application REQUIRES the netcdf library.
dnl Note: the NETCDF_LIBNAME has the leading 'lib' and trailing '.a' stripped from it!  I.e.,
dnl you would use it like this in Makefile.in:
dnl	$(CC) -o progfile progfile.c -I@NETCDF_INCDIR@ -L@NETCDF_LIBDIR@ -l@NETCDF_LIBNAME@
dnl
dnl Easiest way to use this file is to add its contents to /usr/share/autoconf/aclocal.m4
dnl
dnl version 1.0
dnl 1 Nov 2004
dnl David W. Pierce
dnl Climate Research Division
dnl Scripps Institution of Oceanography
dnl dpierce@ucsd.edu
dnl
dnl *******************************************************************************
dnl This code is in the public domain, and can be used for any purposes whatsoever.
dnl *******************************************************************************
dnl
dnl
AC_DEFUN([AC_PATH_NETCDF],[
AC_ARG_WITH( netcdf_incdir, [ --with-netcdf_incdir=dir directory containing netCDF includes],  NETCDF_INCDIR=$withval)
dnl
dnl
dnl =================================================================================
dnl check for netcdf include directory
dnl
err=0
if test x$NETCDF_INCDIR != x; then
        AC_CHECK_HEADER( $NETCDF_INCDIR/netcdf.h,
          echo "Using user-specified netCDF include dir=$NETCDF_INCDIR",
          err=1 )
fi
if test $err -eq 1; then
        echo "Error: user specified netCDF include directory does not have netcdf.h!"
        exit -1
fi
if test x$NETCDF_INCDIR = x; then
        AC_CHECK_HEADER( /usr/local/include/netcdf.h, NETCDF_INCDIR=/usr/local/include )
fi
if test x$NETCDF_INCDIR = x; then
        AC_CHECK_HEADER( /usr/include/netcdf.h, NETCDF_INCDIR=/usr/include )
fi
if test x$NETCDF_INCDIR = x; then
        AC_CHECK_HEADER( $HOME/include/netcdf.h, NETCDF_INCDIR=$HOME/include )
fi
if test x$NETCDF_INCDIR = x; then
        AC_CHECK_HEADER( $HOME/src/netcdf/netcdf.h, NETCDF_INCDIR=$HOME/src/netcdf )
fi
if test x$NETCDF_INCDIR = x; then
        AC_CHECK_HEADER( $HOME/src/netcdf/netcdf-3.4/netcdf.h, NETCDF_INCDIR=$HOME/src/netcdfnetcdf-3.4 )
fi
if test x$NETCDF_INCDIR = x; then
        AC_CHECK_HEADER( $HOME/src/netcdf/netcdf-3.4/src/libsrc/netcdf.h, NETCDF_INCDIR=$HOME/src/netcdfnetcdf-3.4/src/libsrc )
fi
if test x$NETCDF_INCDIR = x; then
        AC_CHECK_HEADER( /sw/include/netcdf.h, NETCDF_INCDIR=/sw/include )
fi
if test x$NETCDF_INCDIR = x; then
	echo " "
        echo "Fatal error: I cannot find the directory that holds the netcdf include file netcdf.h!"
	echo "You can specify it as follows:"
	echo "      ./configure --with-netcdf_incdir=directory_with_file_netcdf.h"
	echo " "
	echo " *** Special note for R CMD INSTALL users: *********************************"
	echo "     The syntax for specifying multiple --configure-args does not seem to be"
	echo "     well documented in R.  If you have installed the netcdf include and library"
	echo "     directories in some non-standard location, you can specify BOTH these"
	echo "     during the R CMD INSTALL process using the following syntax:"
	echo " "
	echo "   R CMD INSTALL --configure-args=\"-with-netcdf_incdir=/path/to/netcdf/incdir -with-netcdf_libdir=/path/to/netcdf/libdir\" ncdf_1.1.tar.gz"
	echo " "
	echo "     where you should, of course, specify your own netcdf include and library"
	echo "     directories, and the actual package name."
	echo " ***************************************************************************"
	echo " "
        exit -1
fi
echo "Found netcdf.h in: $NETCDF_INCDIR"
dnl
dnl
dnl =================================================================================
dnl check for name of netcdf library
dnl
NETCDF_LIBNAME=libnetcdf.a
AC_ARG_WITH( netcdf_libname, [--with-netcdf_libname=fname name of netcdf library file [libnetcdf.a]], NETCDF_LIBNAME=$withval)
dnl
dnl
dnl =================================================================================
dnl check for netcdf lib directory
dnl
AC_ARG_WITH( netcdf_libdir, [ --with-netcdf_libdir=dir directory containing netCDF library],  NETCDF_LIBDIR=$withval)
err=0
if test x$NETCDF_LIBDIR != x; then
        AC_CHECK_FILE( $NETCDF_LIBDIR/$NETCDF_LIBNAME,
          echo "Using user-specified netCDF library dir=$NETCDF_LIBDIR",
          err=1 )
fi
if test $err -eq 1; then
        echo "Error: user specified netCDF library directory does not have $NETCDF_LIBNAME !"
        exit -1
fi
if test x$NETCDF_LIBDIR = x; then
        AC_CHECK_FILE( /usr/local/lib/$NETCDF_LIBNAME, NETCDF_LIBDIR=/usr/local/lib )
fi
if test x$NETCDF_LIBDIR = x; then
        AC_CHECK_FILE( /usr/local/lib32/r4i4/$NETCDF_LIBNAME, NETCDF_LIBDIR=/usr/local/lib32/r4i4 )
fi
if test x$NETCDF_LIBDIR = x; then
        AC_CHECK_FILE( /usr/lib/$NETCDF_LIBNAME, NETCDF_LIBDIR=/usr/lib )
fi
if test x$NETCDF_LIBDIR = x; then
        AC_CHECK_FILE( /lib/$NETCDF_LIBNAME, NETCDF_LIBDIR=/lib )
fi
if test x$NETCDF_LIBDIR = x; then
        AC_CHECK_FILE( $HOME/lib/$NETCDF_LIBNAME, NETCDF_LIBDIR=$HOME/lib )
fi
if test x$NETCDF_LIBDIR = x; then
        AC_CHECK_FILE( $HOME/src/netcdf/$NETCDF_LIBNAME, NETCDF_LIBDIR=$HOME/src/netcdf )
fi
if test x$NETCDF_LIBDIR = x; then
        AC_CHECK_FILE( $HOME/src/netcdf/netcdf-3.4/src/libsrc/$NETCDF_LIBNAME, NETCDF_LIBDIR=$HOME/src/netcdf/netcdf-3.4/src/libsrc/ )
fi
if test x$NETCDF_LIBDIR = x; then
        AC_CHECK_FILE( /sw/lib/$NETCDF_LIBNAME, NETCDF_LIBDIR=/sw/lib/ )
fi
if test x$NETCDF_LIBDIR = x; then
        echo "Fatal error: I cannot find the directory that holds the netcdf library file $NETCDF_LIBNAME !"
        echo "The default library file is named libnetcdf.a."
	echo "If I should look for a netcdf library file with a different name than the"
	echo "default, you can specify it as follows:"
	echo "      ./configure --with-netcdf_libname=lib_file_name.a"
	echo "If I should look in a different directory for the library file,"
        echo "you can specify it as follows:"
	echo "      ./configure --with-netcdf_libdir=directory_with_library_file"
	echo " "
	echo " *** Special note for R CMD INSTALL users: *********************************"
	echo "     The syntax for specifying multiple --configure-args does not seem to be"
	echo "     well documented in R.  If you have installed the netcdf include and library"
	echo "     directories in some non-standard location, you can specify BOTH these"
	echo "     during the R CMD INSTALL process using the following syntax:"
	echo " "
	echo "   R CMD INSTALL --configure-args=\"-with-netcdf_incdir=/path/to/netcdf/incdir -with-netcdf_libdir=/path/to/netcdf/libdir\" ncdf_1.1.tar.gz"
	echo " "
	echo "     where you should, of course, specify your own netcdf include and library"
	echo "     directories, and the actual package name."
	echo " ***************************************************************************"
	echo " "
        exit -1
fi
echo "Found netcdf library file $NETCDF_LIBNAME in directory $NETCDF_LIBDIR"
NETCDF_LIBNAME=`echo $NETCDF_LIBNAME | sed s/lib// | sed s/\.a//`
dnl
dnl Export our variables
dnl
AC_SUBST(NETCDF_INCDIR)
AC_SUBST(NETCDF_LIBDIR)
AC_SUBST(NETCDF_LIBNAME)
dnl
])
