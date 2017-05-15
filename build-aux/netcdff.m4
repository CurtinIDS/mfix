dnl
dnl AC_PATH_NETCDFF()
dnl
dnl Defines NETCDFF_INCDIR, NETCDFF_LIBDIR, and NETCDFF_LIBNAME, or exits with an error message.
dnl I.e., the intent is to have this used when the application REQUIRES the netcdff library.
dnl Note: the NETCDFF_LIBNAME has the leading 'lib' and trailing '.a' stripped from it!  I.e.,
dnl you would use it like this in Makefile.in:
dnl	$(CC) -o progfile progfile.c -I@NETCDFF_INCDIR@ -L@NETCDFF_LIBDIR@ -l@NETCDFF_LIBNAME@
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
AC_DEFUN([AC_PATH_NETCDFF],[
AC_ARG_WITH( netcdff_incdir, [ --with-netcdff_incdir=dir directory containing netCDF includes],  NETCDFF_INCDIR=$withval)
dnl
dnl
dnl =================================================================================
dnl check for netcdff include directory
dnl
err=0
if test x$NETCDFF_INCDIR != x; then
        AC_CHECK_HEADER( $NETCDFF_INCDIR/netcdff.h,
          echo "Using user-specified netCDF include dir=$NETCDFF_INCDIR",
          err=1 )
fi
if test $err -eq 1; then
        echo "Error: user specified netCDF include directory does not have netcdff.h!"
        exit -1
fi
if test x$NETCDFF_INCDIR = x; then
        AC_CHECK_HEADER( /usr/local/include/netcdff.h, NETCDFF_INCDIR=/usr/local/include )
fi
if test x$NETCDFF_INCDIR = x; then
        AC_CHECK_HEADER( /usr/include/netcdff.h, NETCDFF_INCDIR=/usr/include )
fi
if test x$NETCDFF_INCDIR = x; then
        AC_CHECK_HEADER( $HOME/include/netcdff.h, NETCDFF_INCDIR=$HOME/include )
fi
if test x$NETCDFF_INCDIR = x; then
        AC_CHECK_HEADER( $HOME/src/netcdff/netcdff.h, NETCDFF_INCDIR=$HOME/src/netcdff )
fi
if test x$NETCDFF_INCDIR = x; then
        AC_CHECK_HEADER( $HOME/src/netcdff/netcdff-3.4/netcdff.h, NETCDFF_INCDIR=$HOME/src/netcdffnetcdff-3.4 )
fi
if test x$NETCDFF_INCDIR = x; then
        AC_CHECK_HEADER( $HOME/src/netcdff/netcdff-3.4/src/libsrc/netcdff.h, NETCDFF_INCDIR=$HOME/src/netcdffnetcdff-3.4/src/libsrc )
fi
if test x$NETCDFF_INCDIR = x; then
        AC_CHECK_HEADER( /sw/include/netcdff.h, NETCDFF_INCDIR=/sw/include )
fi
if test x$NETCDFF_INCDIR = x; then
	echo " "
        echo "Fatal error: I cannot find the directory that holds the netcdff include file netcdff.h!"
	echo "You can specify it as follows:"
	echo "      ./configure --with-netcdff_incdir=directory_with_file_netcdff.h"
	echo " "
	echo " *** Special note for R CMD INSTALL users: *********************************"
	echo "     The syntax for specifying multiple --configure-args does not seem to be"
	echo "     well documented in R.  If you have installed the netcdff include and library"
	echo "     directories in some non-standard location, you can specify BOTH these"
	echo "     during the R CMD INSTALL process using the following syntax:"
	echo " "
	echo "   R CMD INSTALL --configure-args=\"-with-netcdff_incdir=/path/to/netcdff/incdir -with-netcdff_libdir=/path/to/netcdff/libdir\" ncdf_1.1.tar.gz"
	echo " "
	echo "     where you should, of course, specify your own netcdff include and library"
	echo "     directories, and the actual package name."
	echo " ***************************************************************************"
	echo " "
        exit -1
fi
echo "Found netcdff.h in: $NETCDFF_INCDIR"
dnl
dnl
dnl =================================================================================
dnl check for name of netcdff library
dnl
NETCDFF_LIBNAME=libnetcdff.a
AC_ARG_WITH( netcdff_libname, [--with-netcdff_libname=fname name of netcdff library file [libnetcdff.a]], NETCDFF_LIBNAME=$withval)
dnl
dnl
dnl =================================================================================
dnl check for netcdff lib directory
dnl
AC_ARG_WITH( netcdff_libdir, [ --with-netcdff_libdir=dir directory containing netCDF library],  NETCDFF_LIBDIR=$withval)
err=0
if test x$NETCDFF_LIBDIR != x; then
        AC_CHECK_FILE( $NETCDFF_LIBDIR/$NETCDFF_LIBNAME,
          echo "Using user-specified netCDF library dir=$NETCDFF_LIBDIR",
          err=1 )
fi
if test $err -eq 1; then
        echo "Error: user specified netCDF library directory does not have $NETCDFF_LIBNAME !"
        exit -1
fi
if test x$NETCDFF_LIBDIR = x; then
        AC_CHECK_FILE( /usr/local/lib/$NETCDFF_LIBNAME, NETCDFF_LIBDIR=/usr/local/lib )
fi
if test x$NETCDFF_LIBDIR = x; then
        AC_CHECK_FILE( /usr/local/lib32/r4i4/$NETCDFF_LIBNAME, NETCDFF_LIBDIR=/usr/local/lib32/r4i4 )
fi
if test x$NETCDFF_LIBDIR = x; then
        AC_CHECK_FILE( /usr/lib/$NETCDFF_LIBNAME, NETCDFF_LIBDIR=/usr/lib )
fi
if test x$NETCDFF_LIBDIR = x; then
        AC_CHECK_FILE( /lib/$NETCDFF_LIBNAME, NETCDFF_LIBDIR=/lib )
fi
if test x$NETCDFF_LIBDIR = x; then
        AC_CHECK_FILE( $HOME/lib/$NETCDFF_LIBNAME, NETCDFF_LIBDIR=$HOME/lib )
fi
if test x$NETCDFF_LIBDIR = x; then
        AC_CHECK_FILE( $HOME/src/netcdff/$NETCDFF_LIBNAME, NETCDFF_LIBDIR=$HOME/src/netcdff )
fi
if test x$NETCDFF_LIBDIR = x; then
        AC_CHECK_FILE( $HOME/src/netcdff/netcdff-3.4/src/libsrc/$NETCDFF_LIBNAME, NETCDFF_LIBDIR=$HOME/src/netcdff/netcdff-3.4/src/libsrc/ )
fi
if test x$NETCDFF_LIBDIR = x; then
        AC_CHECK_FILE( /sw/lib/$NETCDFF_LIBNAME, NETCDFF_LIBDIR=/sw/lib/ )
fi
if test x$NETCDFF_LIBDIR = x; then
        echo "Fatal error: I cannot find the directory that holds the netcdff library file $NETCDFF_LIBNAME !"
        echo "The default library file is named libnetcdff.a."
	echo "If I should look for a netcdff library file with a different name than the"
	echo "default, you can specify it as follows:"
	echo "      ./configure --with-netcdff_libname=lib_file_name.a"
	echo "If I should look in a different directory for the library file,"
        echo "you can specify it as follows:"
	echo "      ./configure --with-netcdff_libdir=directory_with_library_file"
	echo " "
	echo " *** Special note for R CMD INSTALL users: *********************************"
	echo "     The syntax for specifying multiple --configure-args does not seem to be"
	echo "     well documented in R.  If you have installed the netcdff include and library"
	echo "     directories in some non-standard location, you can specify BOTH these"
	echo "     during the R CMD INSTALL process using the following syntax:"
	echo " "
	echo "   R CMD INSTALL --configure-args=\"-with-netcdff_incdir=/path/to/netcdff/incdir -with-netcdff_libdir=/path/to/netcdff/libdir\" ncdf_1.1.tar.gz"
	echo " "
	echo "     where you should, of course, specify your own netcdff include and library"
	echo "     directories, and the actual package name."
	echo " ***************************************************************************"
	echo " "
        exit -1
fi
echo "Found netcdff library file $NETCDFF_LIBNAME in directory $NETCDFF_LIBDIR"
NETCDFF_LIBNAME=`echo $NETCDFF_LIBNAME | sed s/lib// | sed s/\.a//`
dnl
dnl Export our variables
dnl
AC_SUBST(NETCDFF_INCDIR)
AC_SUBST(NETCDFF_LIBDIR)
AC_SUBST(NETCDFF_LIBNAME)
dnl
])
