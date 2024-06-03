import os
import glob
import sys
import numpy as np
import datetime
# also do a pip3 install cfgrib
import xarray
import argparse
import matplotlib
import matplotlib.pyplot as plt
import warnings as warn
from PIL import Image
from netCDF4 import Dataset
from matplotlib.pylab import *
from scipy.stats import pearsonr
matplotlib.use('Agg')
matplotlib.rcParams['agg.path.chunksize'] = 1000000

# -------------------------
# turn off interactive mode
# -------------------------
ioff()

# create custom exception
# -----------------------
class EmptyByteStrException(Exception):
  pass

def compareNetcdfsOrGribs(netcdf_or_grib1, netcdf_or_grib2, toPlot, outDir):

  # get output directory for plotting if necessary
  # also make sure we use full absolute paths for input files
  # ---------------------------------------------------------
  outputDir = outDir
  netcdf_or_grib1 = os.path.abspath(netcdf_or_grib1)
  netcdf_or_grib2 = os.path.abspath(netcdf_or_grib2)

  # --------------------------------------
  # open up NetCDFs as new dataset objects
  # --------------------------------------
  if netcdf_or_grib1.lower().endswith('.nc') and netcdf_or_grib2.lower().endswith('.nc'):
    # read the variable names using netcdf4 library, if we are dealing with netcdfs
    # -----------------------------------------------------------------------------
    dataset1 = Dataset(netcdf_or_grib1, 'r')
    dataset2 = Dataset(netcdf_or_grib2, 'r')

    variablesFile1 = [v.strip() for v in list(dataset1.variables)]
    variablesFile2 = [v.strip() for v in list(dataset2.variables)]
  else:
    dataset1 = xarray.open_dataset(netcdf_or_grib1, engine = 'cfgrib')
    dataset2 = xarray.open_dataset(netcdf_or_grib2, engine = 'cfgrib')
    variablesFile1 = [v.strip() for v in list(dataset1.keys())] 
    variablesFile2 = [v.strip() for v in list(dataset2.keys())] 

  # --------------------------
  # write out name of variable
  # --------------------------
  variablesFile1.sort()
  variablesFile2.sort()
  
  os.chdir(outputDir)
  timeStamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
  baseNameStats = 'comparison_statistics_' + timeStamp + '.txt'
  outNameStats = os.path.join(outputDir, baseNameStats)
  f = open(outNameStats, 'w')

  # -------------------------------------------
  # get common variables from both NetCDF files
  # -------------------------------------------
  startDate = datetime.datetime.now().strftime('%Y-%m-%d-%H:%M:%S')
  commonVariableKeys = list(set(variablesFile1).intersection(variablesFile2))
  if len(commonVariableKeys) < 1:
    f.write(' \n No common variables keys between NetCDF files. Exiting ... \n ')
    sys.exit(1)

  # open up text file for writing pixel statistics
  # ----------------------------------------------
  f.write('%s\n'%' ')
  f.write('%s\n'%'  Running comparison at: ')
  f.write('%s\n'%('    '+startDate))
  f.write('%s\n'%' ')
  f.write('%s\n'%('  netcdf file 1: '+os.path.basename(netcdf_or_grib1)))
  f.write('%s\n'%('  netcdf file 2: '+os.path.basename(netcdf_or_grib2)))
  f.write('%s\n'%' ')
 
  if toPlot:
    f.write(  '  *Note: comparison plots will be created in: '+outputDir+'\n\n' )
  noData = -9999.0
  variableCounter = 1
  allVariables = variablesFile1 + variablesFile2
  allVariables = list(set(allVariables))

  for var in allVariables:

    # ------------------------------------------
    # make sure both variables are in the common
    # variable keys list
    # ------------------------------------------
    var1 = var
    var2 = var
    if var1 not in commonVariableKeys or var2 not in commonVariableKeys:
      continue

    # ---------------------------------------------
    # initialize empty lists to hold MISSING values
    # ---------------------------------------------
    missingValues1 = []
    missingValues2 = []

    # ---------------------------------------------
    # try to get the MISSING value(s) from the
    # NetCDF datasets
    # ---------------------------------------------
    try:
      missingValues1 = dataset1[var1]._FillValue
      missingValues2 = dataset2[var2]._FillValue
    except: pass

    # --------------------------------------
    # convert MISSING values to regular list
    # --------------------------------------
    if type(missingValues1) is not list:
      missingValues1 = [missingValues1]
    if type(missingValues2) is not list:
      missingValues2 = [missingValues2]

    # using filename, get xarray engine format to read the file
    # ---------------------------------------------------------
    if not netcdf_or_grib1.lower().endswith('.nc'):
      engine1 = 'cfgrib'
    else:
      engine1 = 'netcdf4'

    if not netcdf_or_grib2.lower().endswith('.nc'):
      engine2 = 'cfgrib'   
    else:
      engine2 = 'netcdf4'

    # check for empty binary string for variable
    # ------------------------------------------
    ds1 = xarray.open_dataset(netcdf_or_grib1, engine = engine1)
    ds2 = xarray.open_dataset(netcdf_or_grib2, engine = engine2)

    # get arrays of data
    # ------------------
    try:
      ds1 = xarray.open_dataset(netcdf_or_grib1, engine = engine1)
      ds2 = xarray.open_dataset(netcdf_or_grib2, engine = engine2)
      if np.all(ds1[var1].values == b'') or np.all(ds2[var2].values == b''):
        raise EmptyByteStrException('WARNING: for var: '+str(var1) + \
          ', var. is empty byte string literal. Continuing ...') 
      netCDF_array1 = np.array(ds1[var1].values).astype(float)
      netCDF_array2 = np.array(ds2[var2].values).astype(float)
    except EmptyByteStrException:
      continue
    except Exception as e:
      print('ERROR (fatal): Unable to read one of the input file(s). Error message is: ')
      print('   => ' + str(e) + ' . Exiting ... ')
      sys.exit(1)
 
    # convert from masked array to regular numpy array (if necessary)
    # ---------------------------------------------------------------
    try:
      netCDF_array1 = np.array(netCDF_array1, dtype = np.float32)
      netCDF_array2 = np.array(netCDF_array2, dtype = np.float32)
    except:
      print("WARNING: for variable '{}' unable to convert array to float: ".format(var))
      print('  first array: ')
      print('    '+str(netCDF_array1))
      print('  second array: ')
      print('    '+str(netCDF_array2))
      continue

    # move onto next variable if either variable
    # has zero elements
    # ------------------------------------------
    if netCDF_array1.size == 0 or netCDF_array2.size == 0:
      continue

    # if variables have different number of elements,
    # then move on
    # -----------------------------------------------
    if netCDF_array1.size != netCDF_array2.size:
      continue

    # set noData values to -9999 for purpose of validation
    # ----------------------------------------------------
    if len(missingValues1) > 0:
      for missingValue in missingValues1:
        if netCDF_array1.size > 1:
          netCDF_array1[(netCDF_array1 == missingValue)] = noData

    if len(missingValues2) > 0:
      for missingValue in missingValues2:
        if netCDF_array1.size > 1:
          netCDF_array2[(netCDF_array2 == missingValue)] = noData

    # get valid  pixel values (e.g. that are -9999.0)
    # -----------------------------------------------
    if netCDF_array1.size>1 and netCDF_array2.size>1:
      validPixels_array1 = netCDF_array1[(netCDF_array1 != noData)]
      validPixels_array2 = netCDF_array2[(netCDF_array2 != noData)]
      validPixels_array1 = validPixels_array1[~np.isnan(validPixels_array1)]
      validPixels_array2 = validPixels_array2[~np.isnan(validPixels_array2)]
    else:
      validPixels_array1 = np.copy(netCDF_array1)
      validPixels_array2 = np.copy(netCDF_array2)
      validPixels_array1 = validPixels_array1[~np.isnan(validPixels_array1)]
      validPixels_array2 = validPixels_array2[~np.isnan(validPixels_array2)]
   
    f.write('%s\n'%('  '+str(variableCounter)+'. '+var1+'\n'))
    if( validPixels_array1.size>0 ):
      f.write( '    '+os.path.basename(netcdf_or_grib1)+' ('+var1+') '+'\n' )
      f.write( '%s\n'%( '      Mean      : '+str(validPixels_array1.mean())))
      f.write( '%s\n'%( '      Min,      : '+str(validPixels_array1.min())))
      f.write( '%s\n'%( '      Max.      : '+str(validPixels_array1.max())))
      f.write( '%s\n'%( '      Std. Dev. : '+str(validPixels_array1.std())))
    else:
      f.write( '%s\n'%( '      Mean      : '+' (no valid pixels, all filler or noData) '))
      f.write( '%s\n'%( '      Min,      : '+' (no valid pixels, all filler or noData) '))
      f.write( '%s\n'%( '      Max.      : '+' (no valid pixels, all filler or noData) '))
      f.write( '%s\n'%( '      Std. Dev. : '+' (no valid pixels, all filler or noData) '))

    if( validPixels_array2.size>0 ):
      f.write( '    '+os.path.basename(netcdf_or_grib2)+' ('+var2+') '+'\n' )
      f.write( '%s\n'%( '      Mean      : '+str(validPixels_array2.mean())))
      f.write( '%s\n'%( '      Min,      : '+str(validPixels_array2.min())))
      f.write( '%s\n'%( '      Max.      : '+str(validPixels_array2.max())))
      f.write( '%s\n'%( '      Std. Dev. : '+str(validPixels_array2.std())))
    else:
      f.write( '%s\n'%( '      Mean      : '+' (no valid pixels, all filler or noData) '))
      f.write( '%s\n'%( '      Min,      : '+' (no valid pixels, all filler or noData) '))
      f.write( '%s\n'%( '      Max.      : '+' (no valid pixels, all filler or noData) '))
      f.write( '%s\n'%( '      Std. Dev. : '+' (no valid pixels, all filler or noData) '))

    # get absolute value of differences between two arrays
    # ----------------------------------------------------
    f.write( '%s\n' % ' ')
    f.write( '%s\n' %('    '+os.path.basename(netcdf_or_grib1)+' ('+var1+') ' ))
    f.write( '%s\n' %('    '+os.path.basename(netcdf_or_grib2)+' ('+var2+') ' ))
    differences = np.absolute(netCDF_array1 - netCDF_array2)
   
    '''
    Below we are trying to make sure that after creating a new array
    that contains the differences between the two arrays, that
    anywhere there is a noData value in either array, this is set
    in the differences array
    '''
    # if first netcdf array is scalar value, set
    # it to 1 element array that is the noData value
    # -----------------------------------------------
    if netCDF_array1.size == 1 or netCDF_array2.size == 1:
      if netCDF_array1 == noData or netCDF_array2 == noData:
        differences = np.array([noData],dtype = float)
    else:
      differences[(netCDF_array1 == noData)] = noData
   
    if differences.size > 1:
      differences = differences[~np.isnan(differences)]

      validDifferences = differences[(differences != noData)]
      if validDifferences.size > 0:
        f.write('%s\n'%('      difference number of elements     : '+\
                str(differences[(differences != noData)].size)))
        f.write('%s\n'%('      difference Mean,Min.,Max.,St. dev : '+str(
            differences[(differences != noData)].mean())+', '+\
                    str(differences[(differences!=noData)].min())+', '+\
                    str(differences[(differences!=noData)].max())+', '+\
                    str(differences[(differences!=noData)].std())))
      else:
        f.write('%s\n'%('      difference number of elements     : '+\
                str(differences[(differences != noData)].size)))
        f.write('%s\n'%('      difference Mean,Min.,Max.,St. dev : '+str(
            noData)+', '+\
                    str(noData)+', '+\
                    str(noData)+', '+\
                    str(noData)))
      f.write('%s\n'%' ')
    else:
      f.write( '%s\n'%('      differences number of elements  : {}'.format(str(netCDF_array1.size))))
      f.write( '%s\n'%('      difference Mean,Min.,Max.,St. dev : noData,noData,noData,noData'))
      f.write( '%s\n'% ' ')

    # --------------------------
    # increment variable counter
    # --------------------------
    variableCounter += 1

    if toPlot:
     
      dims1 = list(netCDF_array1.shape)
      dims2 = list(netCDF_array2.shape)

      if not (dims1 == dims2):
        f.write(  '%s\n' % '   \n   *WARNING: variables ' + var1.upper() + ' does not have same dimensions for NetCDFs. Unable to plot.\n')
      else:

        # ----------------------------------------
        # output names to contain comparison plots
        # ----------------------------------------
        outName1 = os.path.join( outputDir,
          os.path.basename(netcdf_or_grib1).replace('.nc','') + '_' + var1.upper() + '_1.png' )
        outName2 = os.path.join( outputDir,
          os.path.basename(netcdf_or_grib2).replace('.nc','') + '_' + var2.upper() + '_2.png' )
        outNameDiff = os.path.join( outputDir,
          os.path.basename(netcdf_or_grib2).replace('.nc', '') + '_' + var1.upper() + '_DIFFERENCE.png' )

        # ----------------------------------------------
        # attempt to extract 2D dimensions from array(s)
        # ----------------------------------------------
        dims = np.array(dims1)
        dims = dims[(dims > 1)]

        if(dims.size == 2):

          arr1 = np.reshape(netCDF_array1, (dims[0],dims[1]))
          arr2 = np.reshape(netCDF_array2, (dims[0],dims[1]))
         
          # cast to "double precision" as python sees it
          # --------------------------------------------
          InvalidLocs = np.where( (arr1 == noData) | (arr2 == noData))
          arr1 = np.array(arr1,dtype = np.float64)
          arr2 = np.array(arr2,dtype = np.float64)
          diff = np.absolute(arr1 - arr2)

          arr1[(arr1 == noData)] = np.nan
          arr2[(arr2 == noData)] = np.nan
          diff[InvalidLocs] = np.nan

          matshow(arr1)
          plt.title(var1,fontsize=8)
          plt.colorbar()
          plt.grid()
          plt.savefig(outName1,dpi=150)
          plt.close()

          matshow(arr2)
          plt.title(var2,fontsize=8)
          plt.colorbar()
          plt.grid()
          plt.savefig(outName2,dpi=150)
          plt.close()

          matshow(diff)
          plt.title(var1 + ' DIFFERENCE ', fontsize=8)
          plt.colorbar()
          plt.grid()
          plt.savefig( outNameDiff,dpi=150 )
          plt.close()
          
          # if variables are same dimension, do scatter plot
          # ------------------------------------------------
          if arr1.size == arr2.size:
            arr1_1d = arr1.flatten()
            arr2_1d = arr2.flatten()

            valid_locs = np.where((~np.isnan(arr1_1d))&(~np.isnan(arr2_1d)))
            if valid_locs[0].size > 0:
              arr1_1d = arr1_1d[valid_locs]
              arr2_1d = arr2_1d[valid_locs] 
              corr, _ = pearsonr(arr1_1d, arr2_1d)
              corr = np.around(corr, decimals = 3)
              plt.title(var1 + ' CORRELATION ' + ' - Pearson Corr. Coef: (' + str(corr) + ')', fontsize = 8)
            else:            
              plt.title(var1 + ' CORRELATION', fontsize=8)

            plt.scatter(arr1_1d, arr2_1d, s = 0.1)
            plt.grid()
            plt.savefig(outName1.replace('_1.png', '_correlation.png'))
            plt.close() 

          del arr1, arr2, diff

        elif(dims.size == 1):

          arr1 = netCDF_array1.flatten()
          arr2 = netCDF_array2.flatten()

          # --------------------------------------------
          # cast to "double precision" as python sees it
          # --------------------------------------------
          arr1 = np.array(arr1,dtype = np.float64)
          arr2 = np.array(arr2,dtype = np.float64)
          diff = np.absolute(arr1 - arr2)

          arr1[(arr1 == noData)] = np.nan
          arr2[(arr1 == noData)] = np.nan
          diff[(arr1 == noData)|(arr2 == noData)] = np.nan

          plt.title(var1,fontsize=5)
          plt.plot(arr1, 'r' )
          plt.grid()
          plt.savefig(outName1, dpi=150)
          plt.close()

          plt.title(var2, fontsize=5)
          plt.plot(arr2, 'b')
          plt.grid()
          plt.savefig(outName2,dpi=150)
          plt.close()

          plt.title(var1+' DIFFERENCE ', fontsize=5)
          plt.plot(diff, 'k')
          plt.grid()
          plt.savefig(outNameDiff,dpi=150)
          plt.close()
          
          # if variables are same dimension, do scatter plot
          # ------------------------------------------------
          if arr1.size == arr2.size:
            valid_locs = np.where((~np.isnan(arr1_1d))&(~np.isnan(arr2_1d)))
            if valid_locs[0].size > 0:
              arr1_1d = arr1_1d[valid_locs]
              arr2_1d = arr2_1d[valid_locs] 
              corr, _ = pearsonr(arr1_1d, arr2_1d)
              corr = np.around(corr, decimals = 3)
              plt.title(var1 + ' CORRELATION ' + ' - Pearson Corr. Coef: (' + str(corr) + ')', fontsize=8)
            else:            
              plt.title(var1 + ' CORRELATION', fontsize=8)

            plt.scatter(arr1_1d, arr2_1d, s=0.1)
            plt.grid()
            plt.savefig(outName1.replace('_1.png', '_correlation.png'))
            plt.close() 

          del arr1,arr2,diff

        elif( dims.size == 1 ):

          arr1 = netCDF_array1.flatten()
          arr2 = netCDF_array2.flatten()

          # --------------------------------------------
          # cast to "double precision" as python sees it
          del arr1,arr2,diff
       
        if not False in [os.path.isfile(f) for f in [outName1,outName2,outNameDiff]]:
          # merge images into single 3-panel PNG
          # ------------------------------------
          imgs = [Image.open(i) for i in [outName1,outName2,outNameDiff]]
          widths, heights = zip(*(img.size for img in imgs))

          sum_width = sum(widths)
          max_height = max(heights)
          new_im = Image.new('RGB',(sum_width, max_height))
          xoffset = 0

          for im in imgs:
            new_im.paste(im, (xoffset,0))
            xoffset += im.size[0]
          outNameMerged = os.path.join( outputDir,
            os.path.basename(netcdf_or_grib2).replace('.nc','')+'_'+var1.upper()+'_MERGED.png' )
          new_im.save(outNameMerged)

  # -----------------
  # close stats. file
  # -----------------
  f.close()

  # --------------------------
  # close both NetCDF datasets
  # --------------------------
  dataset1.close()
  dataset2.close()

def show_usage():
  print('''
    USAGE: $ python3 compare_netcdf_or_grib.py
      --netcdf_or_grib_1, -f1 <netcdf1|grib1> (required)
      --netcdf_or_grib_2, -f2 <netcdf2|grib2> (required)
      --version, -v (get version info.)
      --usage,   -u (get usage info e.g. print this message)
      --plot,    -p (plot the pairs of variables and their differences to PNGs.
      --outdir,  -o (specify output directory, required)
  ''')
  sys.exit(1)

def show_version():
  print('''
    \n
      CompareNetcdf.py
      Version 1.0.1

      @author  : Gerasimos Michalitsianos
      @updated : 20 August 2023 
      @contact : gerasimos.michalitsianos@noaa.gov 
    \n
  ''')
  sys.exit(1)

def main():

  # parse input arguments
  # ---------------------
  netcdf_or_grib_1 = ''
  netcdf_or_grib_2 = ''
  
  plot = None
  outDir = None

  parser = argparse.ArgumentParser(description='python3 NetCDF/GRIB comparison tool.')
  parser.add_argument('-f1','--netcdf_or_grib_1', required = False, type = str,
    dest='netcdf_or_grib_1',help = 'first netcdf (or grib) file')
  parser.add_argument('-f2','--netcdf_or_grib_2', required = False, type = str,
    dest='netcdf_or_grib_2',help = 'second netcdf (or grib) file')
  parser.add_argument('-v','--version', required = False,
    dest='version',help = 'show version.', action = 'store_true')
  parser.add_argument('-u','--usage', required=False,
    dest='showhelp',action = 'store_true', help = 'show help message.')
  parser.add_argument('-p','--plot', required = False,
    dest='createPlot',help = 'create plots.', action = 'store_true')
  parser.add_argument('-o','--outdir', required = True,type = str,
    dest='outDir',help = 'output directory (optional)')
  args = parser.parse_args()

  version  = args.version
  showHelp = args.showhelp
  netcdf_or_grib_1 = args.netcdf_or_grib_1
  netcdf_or_grib_2 = args.netcdf_or_grib_2
  plot     = args.createPlot
  outDir   = args.outDir

  # show version if desired
  # -----------------------
  if version:
    show_version()
 
  # show help message if desired
  # ----------------------------
  if showHelp:
    show_usage()

  # make sure user passed-in an output directory.
  # ---------------------------------------------
  if outDir is None:
    print( '  \n Pass in output directory with -o or --outdir flag.' )
    show_usage()

  # check to see if user passed-in "plot" flag
  # ------------------------------------------
  if plot is not None:
    plot = True
  else:
    plot = False

  # make sure both netcdf files were passed-in
  # via the command-line
  # ------------------------------------------
  if netcdf_or_grib_1 == '' or netcdf_or_grib_2 == '':
    print( '  \n Please pass in two netcdf files.' )
    print( '    Use the --netcdf1 (-n1) and --netcdf2 (-n2) command-line flags. Exiting ... ' )
    sys.exit(1)

  # make sure user passed the command-line arguments.
  # -------------------------------------------------
  if netcdf_or_grib_1 is None or netcdf_or_grib_2 is None:
    print( '  \n Please pass-in filenames of two netcdf files to compare using --netcdf1,--netcd2 flags.')
    show_usage()

  # check to make sure both netcdf files passed-in
  # are indeeed existing files
  # ----------------------------------------------
  if not os.path.isfile(netcdf_or_grib_1):
    print( '  \n Not an existing file: '+netcdf_or_grib_1 + '. ')
    sys.exit(1)
  elif not os.path.isfile(netcdf_or_grib_2):
    print( '  \n Not an existing file: '+netcdf_or_grib_2 + '. ')
    sys.exit(1)
  else: pass

  # make sure both flags have appropriate extensions
  # ------------------------------------------------
  if not netcdf_or_grib_1.endswith('.nc') and not netcdf_or_grib_1.lower().endswith('.grb') \
    and not netcdf_or_grib_1.endswith('.grib') and not netcdf_or_grib_1.lower().endswith('.grib2'):
    print( '  Not a NetCDF or GRIB: ', netcdf_or_grib_1 )
    sys.exit(1)
  elif not netcdf_or_grib_2.endswith('.nc') and not netcdf_or_grib_2.lower().endswith('.grb') \
    and not netcdf_or_grib_2.endswith('.grib') and not netcdf_or_grib_2.lower().endswith('.grib2'):
    print( ' Not a NetCDF or GRIB: ' , netcdf_or_grib_2 )
    sys.exit(1)
  else: pass
  outDir = os.path.abspath(outDir)
  compareNetcdfsOrGribs(netcdf_or_grib_1,netcdf_or_grib_2, plot , outDir)
  
  # clean up pngs not needed anymore
  # --------------------------------
  for png in glob.glob(outDir + '/*.png'):
    if png.lower().endswith('_merged.png'):
      continue
    os.remove(png)

if __name__ == '__main__':
  with warn.catch_warnings():
    warn.filterwarnings("ignore",category=DeprecationWarning)
    main()
