import os
import sys
import numpy as np
import datetime
import xarray
import argparse
import matplotlib
import matplotlib.pyplot as plt
import warnings as warn
from PIL import Image
from netCDF4 import Dataset
from matplotlib.pylab import *
matplotlib.use('Agg')
matplotlib.rcParams['agg.path.chunksize'] = 1000000

# -------------------------
# turn off interactive mode
# -------------------------
ioff()

def CompareNetcdfs( NetCDF1,NetCDF2,ToPlot,OutDir ):

  # ----------------------------------------------
  # get output directory for plotting if necessary
  # ----------------------------------------------
  OutputDir = OutDir 

  # --------------------------------------
  # open up NetCDFs as new dataset objects
  # --------------------------------------
  np.set_printoptions(precision=5)
  Dataset1 = Dataset( NetCDF1, 'r' )
  Dataset2 = Dataset( NetCDF2, 'r' ) 

  # -----------------------------
  # get the names of the variable 
  # keys from both NetCDF files
  # -----------------------------
  VariablesNetCDF1 = [ v.strip() for v in list(Dataset1.variables) ]
  VariablesNetCDF2 = [ v.strip() for v in list(Dataset2.variables) ]
  VariablesNetCDF1.sort()
  VariablesNetCDF2.sort()

  # --------------------------
  # write out name of variable
  # --------------------------
  os.chdir( OutputDir )
  TimeStamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
  BaseNameStats = 'comparison_statistics_'+TimeStamp+'.txt'
  OutNameStats = os.path.join( OutputDir,BaseNameStats )
  f=open( OutNameStats,'w' )

  # -------------------------------------------
  # get common variables from both NetCDF files
  # -------------------------------------------
  StartDate = datetime.datetime.now().strftime('%Y-%m-%d-%H:%M:%S')
  CommonVariableKeys = list( set(VariablesNetCDF1).intersection(VariablesNetCDF2))
  if len( CommonVariableKeys )<1:
    f.write( ' \n No common variables keys between NetCDF files. Exiting ... \n ')
    sys.exit(1)

  # open up text file for writing pixel statistics
  # ----------------------------------------------
  f.write('%s\n'%' ')
  f.write('%s\n'%'  Running comparison at: ')
  f.write('%s\n'%('    '+StartDate))
  f.write('%s\n'%' ')
  f.write('%s\n'%('  netcdf file 1: '+NetCDF1))
  f.write('%s\n'%('  netcdf file 2: '+NetCDF2))
  f.write('%s\n'%' ')
  
  if ToPlot:
    f.write(  '  *Note: comparison plots will be created in: '+OutputDir+'\n\n' )

  # set a nodata value 
  NoData = -9999.0
  VariableCounter=1

  for Var1,Var2 in zip( VariablesNetCDF1,VariablesNetCDF2 ):

    # ------------------------------------------
    # make sure both variables are in the common 
    # variable keys list
    # ------------------------------------------
    if Var1 not in CommonVariableKeys or Var2 not in CommonVariableKeys:
      continue

    # ---------------------------------------------
    # initialize empty lists to hold MISSING values
    # ---------------------------------------------
    MissingValues1=[]
    MissingValues2=[]

    # ---------------------------------------------
    # try to get the MISSING value(s) from the 
    # NetCDF datasets 
    # ---------------------------------------------
    try:
      MissingValues1=Dataset1[Var1]._FillValue
      MissingValues2=Dataset2[Var2]._FillValue
    except: pass

    # --------------------------------------
    # convert MISSING values to regular list
    # --------------------------------------
    if type( MissingValues1 ) is not list:
      MissingValues1 = [ MissingValues1 ]
    if type( MissingValues2 ) is not list:
      MissingValues2 = [ MissingValues2 ]

    # ------------------
    # get arrays of data
    # ------------------
    try:
      ds1 = xarray.open_dataset( NetCDF1, engine='netcdf4' )
      ds2 = xarray.open_dataset( NetCDF2, engine='netcdf4' )
      NetCDF_array1 = np.array(ds1[Var1].values).astype(float)
      NetCDF_array2 = np.array(ds2[Var2].values).astype(float)
    except Exception as e:
      NetCDF_array1 = Dataset1[Var1][:]
      NetCDF_array2 = Dataset2[Var2][:] 
    
    # convert from masked array to regular numpy array (if necessary)
    # ---------------------------------------------------------------
    NetCDF_array1 = np.array(NetCDF_array1,dtype=np.float32)
    NetCDF_array2 = np.array(NetCDF_array2,dtype=np.float32)

    # ------------------------------------------
    # move onto next variable if either variable
    # has zero elements
    # ------------------------------------------
    if NetCDF_array1.size == 0 or NetCDF_array2.size == 0:
      continue

    # -----------------------------------------------
    # if variables have different number of elements,
    # then move on
    # -----------------------------------------------
    if NetCDF_array1.size!=NetCDF_array2.size:
      continue

    # ----------------------------------------------------
    # set NoData values to -9999 for purpose of validation
    # ----------------------------------------------------
    if len(MissingValues1)>0:
      for MissingValue in MissingValues1:
        if NetCDF_array1.size>1:
          NetCDF_array1[(NetCDF_array1==MissingValue)]=NoData

    if len(MissingValues2)>0:
      for MissingValue in MissingValues2:
        if NetCDF_array1.size>1:
          NetCDF_array2[(NetCDF_array2==MissingValue)]=NoData

    # -----------------------------------------------
    # get valid  pixel values (e.g. that are -9999.0)
    # -----------------------------------------------
    if NetCDF_array1.size>1 and NetCDF_array2.size>1:
      ValidPixels_array1 = NetCDF_array1[( NetCDF_array1!=NoData)]
      ValidPixels_array2 = NetCDF_array2[( NetCDF_array2!=NoData)]
      ValidPixels_array1 = ValidPixels_array1[~np.isnan(ValidPixels_array1)]
      ValidPixels_array2 = ValidPixels_array2[~np.isnan(ValidPixels_array2)]
    else:
      ValidPixels_array1 = np.copy(NetCDF_array1)
      ValidPixels_array2 = np.copy(NetCDF_array2)
      ValidPixels_array1 = ValidPixels_array1[~np.isnan(ValidPixels_array1)]
      ValidPixels_array2 = ValidPixels_array2[~np.isnan(ValidPixels_array2)]
    
    f.write('%s\n'%('  '+str(VariableCounter)+'. '+Var1+'\n'))
    if( ValidPixels_array1.size>0 ):
      f.write( '    '+NetCDF1+' ('+Var1+') '+'\n' ) 
      f.write( '%s\n'%( '      Mean      : '+str(ValidPixels_array1.mean())))
      f.write( '%s\n'%( '      Min,      : '+str(ValidPixels_array1.min())))
      f.write( '%s\n'%( '      Max.      : '+str(ValidPixels_array1.max())))
      f.write( '%s\n'%( '      Std. Dev. : '+str(ValidPixels_array1.std())))
    else:
      f.write( '%s\n'%( '      Mean      : '+' (no valid pixels, all filler or NoData) '))
      f.write( '%s\n'%( '      Min,      : '+' (no valid pixels, all filler or NoData) '))
      f.write( '%s\n'%( '      Max.      : '+' (no valid pixels, all filler or NoData) '))
      f.write( '%s\n'%( '      Std. Dev. : '+' (no valid pixels, all filler or NoData) '))

    if( ValidPixels_array2.size>0 ):
      f.write( '    '+NetCDF2+' ('+Var2+') '+'\n' ) 
      f.write( '%s\n'%( '      Mean      : '+str(ValidPixels_array2.mean())))
      f.write( '%s\n'%( '      Min,      : '+str(ValidPixels_array2.min())))
      f.write( '%s\n'%( '      Max.      : '+str(ValidPixels_array2.max())))
      f.write( '%s\n'%( '      Std. Dev. : '+str(ValidPixels_array2.std())))
    else:
      f.write( '%s\n'%( '      Mean      : '+' (no valid pixels, all filler or NoData) '))
      f.write( '%s\n'%( '      Min,      : '+' (no valid pixels, all filler or NoData) '))
      f.write( '%s\n'%( '      Max.      : '+' (no valid pixels, all filler or NoData) '))
      f.write( '%s\n'%( '      Std. Dev. : '+' (no valid pixels, all filler or NoData) '))

    # get absolute value of differences between two arrays
    # ----------------------------------------------------
    f.write( '%s\n' % ' ')
    f.write( '%s\n' %('    '+NetCDF1+' ('+Var1+') ' ))
    f.write( '%s\n' %('    '+NetCDF2+' ('+Var2+') ' ))
    Differences = np.absolute( NetCDF_array1 - NetCDF_array2 )
   
    '''
    Below we are trying to make sure that after creating a new array
    that contains the differences between the two arrays, that 
    anywhere there is a NoData value in either array, this is set
    in the differences array
    '''
    # if first netcdf array is scalar value, set 
    # it to 1 element array that is the NoData value
    # -----------------------------------------------
    if NetCDF_array1.size == 1 or NetCDF_array2.size == 1:
      if NetCDF_array1 == NoData or NetCDF_array2 == NoData:
        Differences = np.array([NoData],dtype=float)
    else:
      Differences[(NetCDF_array1==NoData)]=NoData
    
    if Differences.size>1:
      Differences = Differences[~np.isnan(Differences)]

      ValidDifferences = Differences[(Differences!=NoData)]
      if ValidDifferences.size>0:
        f.write('%s\n'%('      Difference number of elements     : '+\
                str(Differences[(Differences!=NoData)].size)))
        f.write('%s\n'%('      Difference Mean,Min.,Max.,St. dev : '+str(
            Differences[(Differences!=NoData)].mean())+', '+\
                    str(Differences[(Differences!=NoData)].min())+', '+\
                    str(Differences[(Differences!=NoData)].max())+', '+\
                    str(Differences[(Differences!=NoData)].std())))
      else:
        f.write('%s\n'%('      Difference number of elements     : '+\
                str(Differences[(Differences!=NoData)].size)))
        f.write('%s\n'%('      Difference Mean,Min.,Max.,St. dev : '+str(
            NoData)+', '+\
                    str(NoData)+', '+\
                    str(NoData)+', '+\
                    str(NoData)))
      f.write('%s\n'%' ')
    else:
      f.write( '%s\n'%('      Differences number of elements  : NoData'))
      f.write( '%s\n'%('      Difference Mean,Min.,Max.,St. dev : NoData,NoData,NoData,NoData'))
      f.write( '%s\n'% ' ')

    # --------------------------
    # increment variable counter
    # --------------------------
    VariableCounter+=1

    if ToPlot:
      
      Dims1 = list(  NetCDF_array1.shape )
      Dims2 = list(  NetCDF_array2.shape )

      if not ( Dims1 == Dims2 ):
        f.write(  '%s\n' % '   \n   *WARNING: variables '+Var1.upper()+' does not have same dimensions for NetCDFs. Unable to plot.\n')
      else:

        # ----------------------------------------
        # output names to contain comparison plots
        # ----------------------------------------
        OutName1 = os.path.join( OutputDir,
          os.path.basename(NetCDF1).replace('.nc','')+'_'+Var1.upper()+'_1.png' )
        OutName2 = os.path.join( OutputDir,
          os.path.basename(NetCDF2).replace('.nc','')+'_'+Var2.upper()+'_2.png' )
        OutNameDiff = os.path.join( OutputDir,
          os.path.basename(NetCDF2).replace('.nc','')+'_'+Var1.upper()+'_DIFFERENCE.png' )

        # ----------------------------------------------
        # attempt to extract 2D dimensions from array(s)
        # ----------------------------------------------
        Dims=np.array(Dims1)
        Dims=Dims[(Dims>1)]
        if( Dims.size == 2 ):

          Arr1 = np.reshape( NetCDF_array1, (Dims[0],Dims[1]))
          Arr2 = np.reshape( NetCDF_array2, (Dims[0],Dims[1]))
          
          # --------------------------------------------
          # cast to "double precision" as python sees it
          # --------------------------------------------
          InvalidLocs = np.where( (Arr1==NoData)|(Arr2==NoData))
          Arr1 = np.array(Arr1,dtype=np.float64)
          Arr2 = np.array(Arr2,dtype=np.float64)
          Diff = np.absolute( Arr1 - Arr2 )

          Arr1[(Arr1==NoData)]=np.nan
          Arr2[(Arr1==NoData)]=np.nan
          Diff[InvalidLocs]=np.nan

          matshow( Arr1 )
          plt.title( Var1,fontsize=8 )
          plt.colorbar()
          plt.grid()
          plt.savefig( OutName1,dpi=150 )
          plt.close()
          f.write( '%s\n' %('    *Note: '+OutName1+' created.' ))

          matshow( Arr2 )
          plt.title( Var2,fontsize=8 )
          plt.colorbar()
          plt.grid()
          plt.savefig( OutName2,dpi=150 )
          plt.close()
          f.write( '%s\n' %('    *Note: '+OutName2+' created.' ))

          matshow( Diff )
          plt.title( Var1+' DIFFERENCE ',fontsize=8 )
          plt.colorbar()
          plt.grid()
          plt.savefig( OutNameDiff,dpi=150 )
          plt.close()
          f.write( '%s\n' %('    *Note: '+OutNameDiff+' created.\n' ))
          del Arr1,Arr2,Diff

        elif( Dims.size == 1 ):

          Arr1 = NetCDF_array1.flatten() 
          Arr2 = NetCDF_array2.flatten() 

          # --------------------------------------------
          # cast to "double precision" as python sees it
          # --------------------------------------------
          Arr1 = np.array(Arr1,dtype=np.float64)
          Arr2 = np.array(Arr2,dtype=np.float64)
          Diff = np.absolute( Arr1 - Arr2 )

          Arr1[(Arr1==NoData)]=np.nan
          Arr2[(Arr1==NoData)]=np.nan
          Diff[(Arr1==NoData)|(Arr2==NoData)]=np.nan

          plt.title( Var1,fontsize=5 )
          plt.plot( Arr1 ,'r' )
          plt.grid()
          plt.savefig( OutName1,dpi=150 )
          plt.close()
          f.write(  '%s\n' %('    *Note: '+OutName1+' created.' ))

          plt.title( Var2,fontsize=5 )
          plt.plot( Arr2, 'b' )
          plt.grid()
          plt.savefig( OutName2,dpi=150 )
          plt.close()
          f.write(  '%s\n' %('    *Note: '+OutName2+' created.' ))

          plt.title( Var1+' DIFFERENCE ',fontsize=5 )
          plt.plot( Diff, 'k')
          plt.grid()
          plt.savefig( OutNameDiff,dpi=150 )
          plt.close()
          f.write(  '%s\n' %('    *Note: '+OutNameDiff+' created.\n' ))
          del Arr1,Arr2,Diff
       
        if not False in [os.path.isfile(f) for f in [OutName1,OutName2,OutNameDiff]]:
          # merge images into single 3-panel PNG
          # ------------------------------------ 
          imgs = [Image.open(i) for i in [OutName1,OutName2,OutNameDiff]]
          widths, heights = zip(*(img.size for img in imgs))

          sum_width = sum(widths)
          max_height = max(heights)
          new_im = Image.new('RGB',(sum_width, max_height))
          xoffset = 0

          for im in imgs:
            new_im.paste(im, (xoffset,0))
            xoffset += im.size[0]
          OutNameMerged = os.path.join( OutputDir,
            os.path.basename(NetCDF2).replace('.nc','')+'_'+Var1.upper()+'_MERGED.png' )
          new_im.save(OutNameMerged)

  # -----------------
  # close stats. file
  # -----------------
  f.write('%s\n'%'   \n  ----------------------------------------------------------------------------- \n')
  f.close()

  # --------------------------
  # close both NetCDF datasets
  # --------------------------
  Dataset1.close()
  Dataset2.close()

def show_usage():
  print('''
    USAGE: $ python3 CompareNetcdf.py 
      --netcdf1, -n1 <netcdf1> (required)
      --netcdf2, -n2 <netcdf2> (required)
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
      Version 1.0.0

      @author  : Gerasimos Michalitsianos
      @updated : 24 January 2023 
      @contact : Lakithra@protonmail.com 
    \n
  ''')
  sys.exit(1)

def main():

  # ---------------------
  # parse input arguments
  # ---------------------
  netcdf1=''
  netcdf2=''
  
  plot=None
  outDir=None

  parser = argparse.ArgumentParser( description='python3 NetCDF comparison tool.' )
  parser.add_argument( '-n1','--netcdf1',required=False,type=str,
    dest='netcdf1',help='first netcdf file')
  parser.add_argument( '-n2','--netcdf2',required=False,type=str,
    dest='netcdf2',help='second netcdf file')
  parser.add_argument('-v','--version',required=False,
    dest='version',help='show version.',action='store_true')
  parser.add_argument('-u','--usage',required=False,
    dest='showhelp',action='store_true',help='show help message.')
  parser.add_argument('-p','--plot',required=False,
    dest='createPlot',help='create plots.',action='store_true')
  parser.add_argument( '-o','--outdir',required=True,type=str,
    dest='outDir',help='output directory (optional)')
  args = parser.parse_args()

  version  = args.version
  showHelp = args.showhelp
  netcdf1  = args.netcdf1
  netcdf2  = args.netcdf2
  plot     = args.createPlot
  outDir   = args.outDir

  # -----------------------
  # show version if desired 
  # -----------------------
  if version:
    show_version()
 
  # ----------------------------
  # show help message if desired
  # ----------------------------
  if showHelp:
    show_usage()

  # ---------------------------------------------
  # make sure user passed-in an output directory.
  # ---------------------------------------------
  if outDir is None:
    print( '  \n Pass in output directory with -o or --outdir flag.' )
    show_usage()

  # ------------------------------------------
  # check to see if user passed-in "plot" flag 
  # ------------------------------------------
  if plot is not None:
    plot=True
  else:
    plot=False

  # ------------------------------------------
  # make sure both netcdf files were passed-in
  # via the command-line
  # ------------------------------------------
  if netcdf1 == '' or netcdf2 == '':
    print( '  \n Please pass in two netcdf files.' )
    print( '    Use the --netcdf1 (-n1) and --netcdf2 (-n2) command-line flags. Exiting ... ' )
    sys.exit(1)

  # -------------------------------------------------
  # make sure user passed the command-line arguments.
  # -------------------------------------------------
  if netcdf1 is None or netcdf2 is None:
    print( '  \n Please pass-in filenames of two netcdf files to compare using --netcdf1,--netcd2 flags.')
    show_usage()

  # ----------------------------------------------
  # check to make sure both netcdf files passed-in
  # are indeeed existing files
  # ----------------------------------------------
  if not os.path.isfile( netcdf1 ):
    print( '  \n Not an existing file: '+netcdf1+'. ')
    sys.exit(1)
  elif not os.path.isfile( netcdf2 ):
    print( '  \n Not an existing file: '+netcdf2+'. ')
    sys.exit(1)
  else: pass

  # ------------------------------------------------
  # make sure both flags have appropriate extensions
  # ------------------------------------------------
  if not netcdf1.endswith('.nc'):
    print( '  Not a NetCDF: ', netcdf1 )
    sys.exit(1)
  elif not netcdf2.endswith('.nc'):
    print( ' Not a NetCDF: ' , netcdf2 )
    sys.exit(1)
  else: pass
  outDir = os.path.abspath(outDir)
  CompareNetcdfs( netcdf1,netcdf2, plot , outDir )

if __name__ == '__main__':
  with warn.catch_warnings():
    warn.filterwarnings("ignore",category=DeprecationWarning)
    main()
