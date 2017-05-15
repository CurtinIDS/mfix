##
## Programs and options...
##
include /root/applications/visit/1.5.3/linux-intel/include/make-variables
TOPDIR=/root/applications/visit/1.5.3/linux-intel


##
## Libraries and includes...
##
VTK_INCLUDE= \
  -I/root/applications/visit/1.5.3/linux-intel/include/vtk \
  -I/root/applications/visit/1.5.3/linux-intel/include/vtk/Common \
  -I/root/applications/visit/1.5.3/linux-intel/include/vtk/Filtering \
  -I/root/applications/visit/1.5.3/linux-intel/include/vtk/Graphics \
  -I/root/applications/visit/1.5.3/linux-intel/include/vtk/Hybrid \
  -I/root/applications/visit/1.5.3/linux-intel/include/vtk/IO \
  -I/root/applications/visit/1.5.3/linux-intel/include/vtk/Imaging \
  -I/root/applications/visit/1.5.3/linux-intel/include/vtk/Rendering \
  -I/root/applications/visit/1.5.3/linux-intel/include/vtk/Utilities
MOC=/root/applications/visit/1.5.3/linux-intel/bin/moc
CXXFLAGS=$(CXXFLAGSORIG)
CPPFLAGS=$(CPPFLAGSORIG) $(VTK_INCLUDE) -I. -I/root/applications/visit/1.5.3/linux-intel/include -I/root/applications/visit/1.5.3/linux-intel/include/visit
LDFLAGS=$(LDFLAGSORIG)  -L/root/.visit/linux-intel/plugins/databases

##
## Files...
##
PLUGINDIR=databases
PLUGINNAME=MFIXDatabase

WIDGETS=
ISRC=MFIXPluginInfo.C
COMMONSRC=MFIXPluginInfo.C MFIXCommonPluginInfo.C
MSRC=MFIXMDServerPluginInfo.C avtMFIXFileFormat.C
ESRC=MFIXEnginePluginInfo.C avtMFIXFileFormat.C


SRC=$(ISRC) $(COMMONSRC) $(MSRC) $(ESRC)

SHLIB_FORCED=

ELIBS_FOR_MACOSX_PREBINDING=$(BZIP2_LIBS)
MLIBS_FOR_MACOSX_PREBINDING=$(BZIP2_LIBS)
ILIBS=
MLIBS=-lpipeline_ser -ldatabase_ser -lmir_ser -lavtmath_ser -ldbatts -lavtexceptions -lstate -lcomm -lmisc -lexpr -lparser -lplugin -lutility -lvisit_vtk -llightweight_visit_vtk  -L/root/.visit/linux-intel/plugins/databases -lI$(PLUGINNAME) $(MLIBS_FOR_MACOSX_PREBINDING) $(VTK_LIBS)
ESERLIBS=-lpipeline_ser -ldatabase_ser -lmir_ser -lavtmath_ser -ldbatts -lavtexceptions -lstate -lcomm -lmisc -lexpr -lparser -lplugin -lutility -lvisit_vtk -llightweight_visit_vtk  -L/root/.visit/linux-intel/plugins/databases -lI$(PLUGINNAME) $(ELIBS_FOR_MACOSX_PREBINDING) $(VTK_LIBS)
EPARLIBS=-lpipeline_par -ldatabase_par -lmir_par -lavtmath_par -ldbatts -lavtexceptions -lstate -lcomm -lmisc -lexpr -lparser -lplugin -lutility -lvisit_vtk -llightweight_visit_vtk  -L/root/.visit/linux-intel/plugins/databases -lI$(PLUGINNAME) $(ELIBS_FOR_MACOSX_PREBINDING) $(VTK_LIBS) $(SHLIB_MPI_LIBS)

IDSO=/root/.visit/linux-intel/plugins/databases/libIMFIXDatabase.so
MDSO=/root/.visit/linux-intel/plugins/databases/libMMFIXDatabase.so
ESERDSO=/root/.visit/linux-intel/plugins/databases/libEMFIXDatabase_ser.so
EPARDSO=/root/.visit/linux-intel/plugins/databases/libEMFIXDatabase_par.so


DISTRIB=

##
## Derived objects
##
IOBJ=$(ISRC:.C=.o)
GOBJ=$(COMMONSRC:.C=.o) $(GSRC:.C=.o)
SOBJ=$(COMMONSRC:.C=.o) $(SSRC:.C=.o)
VOBJ=$(COMMONSRC:.C=.o) $(VSRC:.C=.o)
MOBJ=$(COMMONSRC:.C=.o) $(MSRC:.C=.o) $(MSPECIFICSRC:.C=_mds.o)
ESEROBJ=$(COMMONSRC:.C=.o) $(ESRC:.C=.o)
EPAROBJ=$(COMMONSRC:.C=.o) $(ESRC:.C=_par.o)

MOCSRC = $(WIDGETS:.h=_moc.C)
MOCOBJ = $(MOCSRC:.C=.o)

##
## Standard targets...
##
all: message $(IDSO) $(MDSERVERLIB) $(ENGINELIBSER) $(ENGINELIBPAR)

clean:
	$(RM) $(IOBJ) $(COMMONOBJ)
	$(RM) $(GOBJ) $(SOBJ) $(VOBJ) $(MOBJ) $(ESEROBJ) $(EPAROBJ)
	$(RM) $(MOCSRC) $(MOCOBJ)
	$(RM) $(IDSO) $(GDSO) $(SDSO) $(VDSO) $(MDSO) $(ESERDSO) $(EPARDSO)
	$(RM) $(JAVAOBJ)

##
## Other targets...
##
message:
	@echo
	@echo "****************************************************************************"
	@echo "*** Building MFIX Database Plugin"
	@echo "****************************************************************************"

##
## moc
##
$(MOCSRC) or_no_widgets: $(WIDGETS)
	@rm -f $@
	$(MOC) $(@:_moc.C=.h) > $@


##
## Automatic dependency stuff
##
include /root/applications/visit/1.5.3/linux-intel/include/make-targets
