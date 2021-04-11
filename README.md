###### COMPARING VARIABLES BETWEEN TWO NETCDF4 FILES
    
    This is a simple Python3 command-line program that iterates through
    the matchines variables of two NetCDF4 files, and computes some 
    basic statistics (mean, min., max., and standard deviation), writing
    those statistics to a text file. For each variable, an array of differences 
    is computed (e.g. the absolute value of the first array minus the second array), 
    and the same statistics are computed and written to the output text file. 
    Variables that do not have the same x,y NetCDF dimensions are skipped, and 
    only variable keys that exist in BOTH NetCDF4 files, passed-in via command-line, 
    are considered.
    
    Note that when the number of elements is displayed for each array,
    however, that these may be differ between the two NetCDFs for that
    variable. This is because the number of NoData (fill values) may
    differ between the two files for that variable. So the number of valid
    pixel values (grid points) USED to compute the statistics may differ.
    
    So the statistics for some variable may differ, because the number of 
    its "NoData" (fill) values differ; these are set to -9999.0 in the script.
    
    The --outdir (or -o) flag is required in which the text file is written to.
    The --plot (-p) flag is optional to create output plots. These PNG plots, 
    if created, also go to the output directory specified by --outdir (or -o flag).
   
###### INSTALLATION:

    To use this software, it is expected you will have docker (https://www.docker.com)
    and git (https://git-scm.com/) installed. If not, you may download the code from
    this repository directly as a .zip. To use git and docker, try these commands:

    $ git clone https://github.com/gerasimosmichalitsianos/compare_netcdf
    $ cd compare_netcdf/
    $ docker build -t compare_netcdf .
     
###### USAGE:
 
    To use this code, pass in the filenames of the two NetCDF4 files you wish to compare.
    
    The first way to use this program is to use Docker e.g.:

    $ ls -l /home/someuser/data
    20210312150000-STAR-L2P_GHRSST-SSTsubskin-AHI_H08-ACSPO_V2.70-v02.0-fv01.0.nc
    20210312150000-STAR-L2P_GHRSST-SSTsubskin-AHI_H08-ACSPO_V2.71-v02.0-fv01.0.nc
    $ DIR=/home/someuser/data
    $ docker run -v $DIR:$DIR compare_netcdf 
      --netcdf1 $DIR/20210312150000-STAR-L2P_GHRSST-SSTsubskin-AHI_H08-ACSPO_V2.71-v02.0-fv01.0.nc 
      --netcdf2 $DIR/20210312150000-STAR-L2P_GHRSST-SSTsubskin-AHI_H08-ACSPO_V2.70-v02.0-fv01.0.nc
      --outdir $DIR --plot

    The second way is to use the Python source file directly (bin/compare_netcdf.py):
    
    $ cd bin/
    $ python3 compare_netcdf.py
      --netcdf1 $DIR/20210312150000-STAR-L2P_GHRSST-SSTsubskin-AHI_H08-ACSPO_V2.71-v02.0-fv01.0.nc 
      --netcdf2 $DIR/20210312150000-STAR-L2P_GHRSST-SSTsubskin-AHI_H08-ACSPO_V2.70-v02.0-fv01.0.nc
      --outdir $DIR
      
    Remember the --plot (-p) is optional for PNG creation.

    Note that when using the Python script directory, you will need to have the following Python3
    modules installed: 
    
      (1) netCDF4 
      (2) NumPy 
      (3) Matplotlib
      
    They can be installed like so:
      pip3 install numpy
      pip3 install netCDF4 
      pip3 install matplotlib
    
    Matplotlib is used here to create the plots.

###### ALL COMMAND-LINE OPTIONS

    --version, -v  : display version info.
    --usage,   -h  : display this usage messsage
    --netcdf1, -n1 : pass in name of 1-band Geotiff holding 1-band panchromatic Geotiff image (high resolution, required)
    --netcdf2, -n2 : pass in name of 3 or 4 band multispectral Geotiff image file (low-resolution, required)
    --plot   , -p  : create output PNG plots for each variable (one plot for variable from first NetCDF1, one from second, then a plot differences)
    --outdir , -o  : output directory (required!)
      
###### PYTHON VERSION:
     
    Supports Python 2.7.x, 3.x
       
###### SAMPLE OUTPUT TO TEXT FILE


    Running comparison at:
    2021-04-11-03:01:59

    netcdf file 1: ../sample_data/20210312150000-STAR-L2P_GHRSST-SSTsubskin-AHI_H08-ACSPO_V2.70-v02.0-fv01.0.nc
    netcdf file 2: ../sample_data/20210312150000-STAR-L2P_GHRSST-SSTsubskin-AHI_H08-ACSPO_V2.71-v02.0-fv01.0.nc

    *Note: comparison plots will be created in: /home/gmichali/sample_data/

    1. brightness_temperature_08um6

      ../sample_data/20210312150000-STAR-L2P_GHRSST-SSTsubskin-AHI_H08-ACSPO_V2.70-v02.0-fv01.0.nc (brightness_temperature_08um6)
        Mean      : 290.02698
        Min,      : 263.3
        Max.      : 296.15
        Std. Dev. : 5.2306266
      ../sample_data/20210312150000-STAR-L2P_GHRSST-SSTsubskin-AHI_H08-ACSPO_V2.71-v02.0-fv01.0.nc (brightness_temperature_08um6)
        Mean      : 290.13742
        Min,      : 262.97998
        Max.      : 296.13998
        Std. Dev. : 5.2312717

      ../sample_data/20210312150000-STAR-L2P_GHRSST-SSTsubskin-AHI_H08-ACSPO_V2.70-v02.0-fv01.0.nc (brightness_temperature_08um6)
      ../sample_data/20210312150000-STAR-L2P_GHRSST-SSTsubskin-AHI_H08-ACSPO_V2.71-v02.0-fv01.0.nc (brightness_temperature_08um6)
        Difference number of elements     : 3556082
        Difference Mean,Min.,Max.,St. dev : 454.91278, 0.0, 10294.86, 2114.8115

      *Note: ../sample_data/20210312150000-STAR-L2P_GHRSST-SSTsubskin-AHI_H08-ACSPO_V2.70-v02.0-fv01.0_BRIGHTNESS_TEMPERATURE_08UM6._1.png created.
      *Note: ../sample_data/20210312150000-STAR-L2P_GHRSST-SSTsubskin-AHI_H08-ACSPO_V2.71-v02.0-fv01.0_BRIGHTNESS_TEMPERATURE_08UM6._2.png created.
      *Note: ../sample_data/20210312150000-STAR-L2P_GHRSST-SSTsubskin-AHI_H08-ACSPO_V2.71-v02.0-fv01.0_BRIGHTNESS_TEMPERATURE_08UM6._DIFFERENCE.png created.
      
###### @author: 

    Gerasimos Michalitsianos
    gerasimosmichalitsianos@gmail.com
    11 April 2021
