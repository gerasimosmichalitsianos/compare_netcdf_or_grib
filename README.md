###### COMPARING VARIABLES BETWEEN TWO NETCDF4 FILES
    
    This is a simple Python3 command-line program that iterates through
    the matchines variables of two NetCDF4 files, and computes some 
    basic statistics (mean, min., max., and standard deviation), writing
    those statistics to the console (e.g. UNIX/Linux terminal). For each
    variable, an array of differences is computed (e.g. the absolute value
    of the first array minus the second array), and the same statistics 
    are computed and written to the console. Variables that do not have
    the same x,y NetCDF dimensions are skipped, and only variable keys that exist 
    in BOTH NetCDF4 files, passed-in via command-line, are considered.
    
    Note that when the number of elements is displayed for each array,
    however, that these may be differ between the two NetCDFs for that
    variable. This is because the number of NoData (fill values) may
    differ between the two files for that variable. So the number of valid
    pixel values (grid points) USED to compute the statistics may differ.
   
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

    The second way is to use the Python source file directly (bin/compare_netcdf.py):
    
    $ python3 compare_netcdf.py
      --netcdf1 $DIR/20210312150000-STAR-L2P_GHRSST-SSTsubskin-AHI_H08-ACSPO_V2.71-v02.0-fv01.0.nc 
      --netcdf2 $DIR/20210312150000-STAR-L2P_GHRSST-SSTsubskin-AHI_H08-ACSPO_V2.70-v02.0-fv01.0.nc

    Note that when using the Python script directory, you will need to have the following Python3
    modules installed: 
    
      (1) netCDF4 
      (2) NumPy 
      (3) Matplotlib
    
    Matplotlib will be used in future versions of the code for making plots.

###### ALL COMMAND-LINE OPTIONS

    --version, -v  : display version info.
    --usage,   -h  : display this usage messsage
    --netcdf1, -n1 : pass in name of 1-band Geotiff holding 1-band panchromatic Geotiff image (high resolution, required)
    --netcdf2, -n2 : pass in name of 3 or 4 band multispectral Geotiff image file (low-resolution, required)
      
###### PYTHON VERSION:
     
    Supports Python 3.x
       
###### Sample Console/Terminal Output

    Running comparison at: 
      2021-03-13-07:47:42
 
    netcdf file 1: /home/someuser/data/20210312150000-STAR-L2P_GHRSST-SSTsubskin-AHI_H08-ACSPO_V2.71-v02.0-fv01.0.nc
    netcdf file 2: /home/someuser/data/20210312150000-STAR-L2P_GHRSST-SSTsubskin-AHI_H08-ACSPO_V2.70-v02.0-fv01.0.nc
   
    ----------------------------------------------------------------------------- 

    1. brightness_temperature_08um6

      NetCDF1 Mean,Min.,Max.,St. dev    : 290.13742, 262.97998, 296.13998, 5.2312717
      NetCDF2 Mean,Min.,Max.,St. dev    : 290.02698, 263.3, 296.15, 5.2306266
      NetCDF1 number of elements        : 3690379
      NetCDF2 number of elements        : 3556082
      Difference number of elements     : 3690379
      Difference Mean,Min.,Max.,St. dev : 812.9029, 0.0, 10294.75, 2775.4597
   
    ----------------------------------------------------------------------------- 

    2. brightness_temperature_10um4

      NetCDF1 Mean,Min.,Max.,St. dev    : 292.83353, 266.28, 298.91998, 5.2316866
      NetCDF2 Mean,Min.,Max.,St. dev    : 292.7556, 266.6, 298.93, 5.228918
      NetCDF1 number of elements        : 3690379
      NetCDF2 number of elements        : 3556082
      Difference number of elements     : 3690379
      Difference Mean,Min.,Max.,St. dev : 813.1033, 0.0, 10297.26, 2776.1316
      
###### @author: 

    Gerasimos Michalitsianos
    gerasimosmichalitsianos@gmail.com
    13 March 2021
