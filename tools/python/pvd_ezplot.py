###########################################################
# Author : Jeff Dietiker
# Date :   January, 20th, 2015
#
# Purpose: This script reads a PVD files and generates
#          a series of data files (extension .dat),
#          and a series of XY plots (displayed on screen
#          and saved in a .png file), one file per variable.
#
#          Developed and tested with python 2.7.9
#
# usage: pvd_ezplot.py [-h] [--tmin TMIN] [--tmax TMAX] [--maxcells MAXCELLS]
#                      [--no-mean] [--no-fig] [--no-png] [--no-dat]
#                      [--legend LEGEND] [--frame_alpha FRAME_ALPHA]
#                      [--font-family FONT_FAMILY] [--font-weight FONT_WEIGHT]
#                      [--font-size FONT_SIZE] [--linewidth LINEWIDTH]
#                      filename
#
# positional arguments:
#   filename              PVD file to be parsed (include path).By default, the
#                         data is plotted on the screen, saved in .png file, and
#                         exported in .dat file (one file per variable, one
#                         column per cell). The mean value of all cells is also
#                         computed, plotted and saved.
#
# optional arguments:
#   -h, --help            show this help message and exit
#   --tmin TMIN           Lower limit of time range (seconds).
#   --tmax TMAX           Upper limit of time range (seconds).
#   --maxcells MAXCELLS   Maximum number of cells. When the number of cells is
#                         larger than maxcells, cell data is not stored in dat
#                         files nor plotted (default =20). This limit is used to
#                         prevent unreasonable data processing for large cell
#                         count. Setting maxcells to 0 will result in only the
#                         mean value to be stored and plotted unless --no-mean
#                         is used (nothing is stored or plotted in that case).
#   --no-mean             Do not compute mean.
#   --no-fig              Do not show Figure.
#   --no-png              Do not save Figure.
#   --no-dat              Do not save Data.
#   --legend LEGEND       Plot legend location (default='best', Acceptable
#                         values = 'right', 'center left', 'upper right', 'lower
#                         right', 'best','center', 'lower left', 'center right',
#                         'upper left', 'upper center', 'lower center'.
#   --frame_alpha FRAME_ALPHA
#                         Plot legend frame transparency (0.0=tranparent,
#                         1.0=opaque, default=1.0)
#   --font-family FONT_FAMILY
#                         Plot font family (default=sans-serif)
#   --font-weight FONT_WEIGHT
#                         Plot font weight (default=normal)
#   --font-size FONT_SIZE
#                         Plot font size (default=16)
#   --linewidth LINEWIDTH
#                         Plot line width (default=3)
#
###########################################################
# Import libraries
###########################################################
import numpy as np
import pylab
import matplotlib.font_manager
from vtk import *
from vtk.util.numpy_support import vtk_to_numpy
import argparse
import os

###########################################################
# Argument parser
###########################################################
UNDEFINED = 9.87654321E31
parser = argparse.ArgumentParser()
# File name
parser.add_argument(    "filename", \
                        help="PVD file to be parsed (include path)."
                        "By default, the data is plotted on the screen, saved in .png file, and exported in .dat file"
                        " (one file per variable, one column per cell)."
                        " The mean value of all cells is also computed, plotted and saved.", \
                        type=str)

# Time range
parser.add_argument(    "--tmin", \
                        help="Lower limit of time range (seconds).", \
                        dest="tmin", \
                        type=float, \
                        default=0.0)

parser.add_argument(    "--tmax", \
                        help="Upper limit of time range (seconds).", \
                        dest="tmax", \
                        type=float, \
                        default=UNDEFINED)

parser.add_argument(    "--maxcells", \
                        help="Maximum number of cells. When the number of cells is larger than maxcells, "
                        "cell data is not stored in dat files nor plotted (default =20). "
                        "This limit is used to prevent unreasonable data processing for large cell count. "
                        "Setting maxcells to 0 will result in only the mean value to be stored and plotted "
                        "unless --no-mean is used (nothing is stored or plotted in that case).", \
                        dest="maxcells", \
                        type=float, \
                        default=20)

# Option to skip mean computation
parser.add_argument(    '--no-mean', \
                        help="Do not compute mean.", \
                        dest='mean', \
                        action='store_false')

# Option to skip showing Figure
parser.add_argument(    '--no-fig', \
                        help="Do not show Figure.", \
                        dest='figure', \
                        action='store_false')

# Option to skip saving Figure in png file
parser.add_argument(    '--no-png', \
                        help="Do not save Figure.", \
                        dest='png', \
                        action='store_false')

# Option to save data in .dat file
parser.add_argument(    '--no-dat', \
                        help="Do not save Data.", \
                        dest='dat', \
                        action='store_false')

# Show legend and legend location
parser.add_argument(    "--legend", \
                        help="Plot legend location (default='best', "
                        "Acceptable values = 'right', 'center left', 'upper right', 'lower right', 'best',"
                        "'center', 'lower left', 'center right', 'upper left', 'upper center', 'lower center'.", \
                        dest="legend", \
                        type=str, \
                        default='best')

# Legend frame transparency (used in plots)
parser.add_argument(    "--frame_alpha", \
                        help="Plot legend frame transparency (0.0=tranparent, 1.0=opaque, default=1.0)", \
                        dest="frame_alpha", \
                        type=float, \
                        default=1.0)


# Font familty (used in plots)
parser.add_argument(    "--font-family", \
                        help="Plot font family (default=sans-serif)", \
                        dest="font_family", \
                        type=str, \
                        default='sans-serif')

# Font weight (used in plots)
parser.add_argument(    "--font-weight", \
                        help="Plot font weight (default=normal)", \
                        dest="font_weight", \
                        type=str, \
                        default='normal')

# Font size (used in plots)
parser.add_argument(    "--font-size", \
                        help="Plot font size (default=16)", \
                        dest="font_size", \
                        type=int, \
                        default=16)

# Plot line width
parser.add_argument(    "--linewidth", \
                        help="Plot line width (default=3)", \
                        dest="linewidth", \
                        type=int, \
                        default=3)

args = parser.parse_args()

###########################################################
# Parse arguments
###########################################################
# Parse file name
PVD_filename = args.filename
PWD      = os.path.dirname(PVD_filename)
PVD_File = os.path.basename(PVD_filename)

# Parse time range
undefined_time_range=(0.0,UNDEFINED)
time_range=(args.tmin,args.tmax)

# Option to compute mean
compute_mean = args.mean

# Max number of cells
maxcells = args.maxcells

# Show legend and legend location
if args.legend>0:
    show_legend = True
    legend_loc = args.legend
else:
    show_legend = False

# Option to show Figure
show_fig = args.figure

# Option to save Figure in png file
save_fig = args.png

# Option to save data in .dat file
save_dat = args.dat


###########################################################
# Functions
###########################################################
def save_data(header,data,extra):
    for n in range(0,n_arrays):
        array_name = arrayname[n]
        outputfile = os.path.join(PWD,PVD_base+"_"+array_name+extra+".dat")
        f=open(outputfile,'w')
        f.write("# {0}\n".format(array_name))
        f.write("# {0}".format(header[n][0]))
        for s in header[n][1:]:
            f.write(" , {0}".format(s))
        f.write("\n")
        for row in range(0,nsteps):
            f.write("{0:14.8E}".format(data[n][0][row]))
            ncols=len(data[n])
            for col in range(1,ncols):
                f.write(" , {0:14.8E}".format(data[n][col][row]))
            f.write("\n")
        f.close()


def Generate_Figures(header,data,extra1,extra2):
    for n in range(0,n_arrays):
        array_name = arrayname[n]
        pylab.figure()
        nplots = len(data[n])
        for p in range(1,nplots):
            pylab.plot( data[n][0], data[n][p],lw=linewidth, label=header[n][p])
        pylab.title(extra1+arrayname[n])
        pylab.xlabel("Time [sec]")
        pylab.xlim(time_range[0],time_range[1])
        if show_legend:
            pylab.legend(loc = legend_loc,framealpha=args.frame_alpha)
            # pylab.legend(loc='best', fancybox=True, framealpha=0.5)
            # pylab.legend(loc='best', fancybox=False, framealpha=0.5)
        if save_fig:
            pngfile = os.path.join(PWD,PVD_base+"_"+array_name+extra2+".png")
            pylab.savefig(pngfile)

###########################################################
# Parse PVD file
###########################################################
# Extract base and extension
PVD_base, PVD_ext = PVD_File.split(".")

file=open(PVD_filename,'r')

print "Processing {0} ...".format(PVD_File)

row = file.readlines()

# Extract time info and corresponding vtu file
vtu_filelist = []
steps=[]

for line in row:
    if line.find("timestep") > -1:
        s = line.split('"')
        vtu_filelist.append((float(s[1]),os.path.join(PWD, s[7])))
        steps.append(float(s[1]))
nsteps = len(vtu_filelist)
print "Found {0} time steps in original file.".format(nsteps)
print "PVD file time range = [{0};{1}] seconds.".format(steps[0],steps[-1])

# Clip time range
if time_range != undefined_time_range:
    steps=[]
    clipped_vtu_filelist=[]
    for time,vtu_file in vtu_filelist:
        if time_range[0]<=time and time<=time_range[1]:
            clipped_vtu_filelist.append((time,vtu_file))
            steps.append(time)
    nsteps = len(clipped_vtu_filelist)
    print "After clipping to user-specified time range, {0} time steps remain.".format(nsteps)
    print "New Time range = [{0};{1}] seconds.".format(steps[0],steps[-1])
else:
    clipped_vtu_filelist=vtu_filelist
time_range = (steps[0],steps[-1])

# List initializations
arrayname=[]
data=[]
data_mean=[]
header=[]
header_mean=[]

velcomp=["U","V","W"]
velcompmean=["U_mean","V_mean","W_mean"]

###########################################################
# Parse first vtu file
###########################################################
# Open first vtu file to get file info and
# set up output

time = clipped_vtu_filelist[0][0]
vtu_file = clipped_vtu_filelist[0][1]
reader = vtk.vtkXMLUnstructuredGridReader()
reader.SetFileName(vtu_file)
reader.Update()

# Get number of arrays and number of cells
n_arrays = reader.GetNumberOfCellArrays()
n_cells  = reader.GetNumberOfCells()

print "Found {0} arrays and {1} cells.".format(n_arrays,n_cells)

# If the number of cells ia larger that maxcells, then only compute the mean
# This is done to avoid plotting to many plots
process_cells = True
if n_cells > maxcells:
    compute_mean = True
    process_cells = False
if n_cells ==1:
    compute_mean = False
    process_cells = True

# Inspect arrays to generate the data header
for n in range(0,n_arrays):
    header.append([])
    array_name = reader.GetCellArrayName(n)
    arrayname.append(array_name)

#Get the data from the vtk file, and extract number of components
    vtk_array = reader.GetOutput().GetCellData().GetArray(array_name)
    ncomp=vtk_array.GetNumberOfComponents()

# Prepare data labels:
# If the variable is a scalar (1 component) then use the variable name.
# If the variable is a vector (3 components), use U,V and W.
# Then append the cell index (only if there are more than one cell).

    header[n].append("Time")
#mean data header
    if compute_mean:
        header_mean.append([])
        header_mean[n].append("Time")
        if ncomp==1:
            header_mean[n].append(array_name+"_mean")
        else:
            for nc in range(0,ncomp):header_mean[n].append(velcompmean[nc])

        data_mean.append([])
        for col in range(0,len(header_mean[n])):
            data_mean[n].append([])

#data header
    if process_cells:
        extra=""
        for cell in range(0,n_cells):
            if n_cells>1: extra = "_"+str(cell)
            if ncomp==1:
                header[n].append(array_name+extra)
            else:
                for nc in range(0,ncomp):header[n].append(velcomp[nc]+extra)

    data.append([])
    for col in range(0,len(header[n])):
        data[n].append([])


###########################################################
# Process all vtu files in the clipped time range
###########################################################
# Loop through the vtu files
for time,vtu_file in clipped_vtu_filelist:

# Load each vtu file
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(vtu_file)
    reader.Update()


    for n in range(0,n_arrays):

#Get the data from the vtk file, and extract number of components
        array_name = reader.GetCellArrayName(n)
        vtk_array = reader.GetOutput().GetCellData().GetArray(array_name)
        ncomp=vtk_array.GetNumberOfComponents()
        numpy_array = vtk_to_numpy(vtk_array )

        mean = np.mean(numpy_array,axis=0)

# First column is the time
        data[n][0].append(time)

# Next is the mean (same as cell data if there is only one cell)
        if compute_mean:
            data_mean[n][0].append(time)
            col=0
            if ncomp==1:
                col+=1
                data_mean[n][col].append(mean)
            else:
                for e in mean.tolist():
                    col+=1
                    data_mean[n][col].append(e)
#data
        if process_cells:
            col=0
            for e in numpy_array.tolist():
                if ncomp==1:
                    col+=1
                    data[n][col].append(e)
                else:
                    for e2 in e:
                        col+=1
                        data[n][col].append(e2)


# Save data in text files
if process_cells and save_dat:
    save_data(header,data,"")
    print "Saving Data in .dat file(s)..."

# Save mean data_in text files
if compute_mean  and save_dat:
    save_data(header_mean,data_mean,"_mean")
    print "Saving Mean data in .dat file(s)..."

# Plot data
font = {'family' : args.font_family,
        'weight' : args.font_weight,
        'size'   : args.font_size}

matplotlib.rc('font', **font)

linewidth = args.linewidth


# Plot data
if process_cells and show_fig:
    print "Generating Figure(s) for each cell..."
    Generate_Figures(header,data,"","")


# Plot mean data
if compute_mean and show_fig:
    print "Generating Figure(s) for mean data..."
    Generate_Figures(header_mean,data_mean,"Mean ","_mean")


if show_fig:
    pylab.show(block=False)
    raw_input('Press Enter to close all Figures.')
    pylab.close()

print "Done."
