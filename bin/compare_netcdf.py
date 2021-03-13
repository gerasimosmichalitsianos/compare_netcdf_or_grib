import os
import sys
import numpy as np
import datetime
import argparse
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from matplotlib.pylab import *

# -------------------------
# turn off interactive mode
# -------------------------
ioff()

def CompareNetcdfs( NetCDF1,NetCDF2 ):

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

  # -------------------------------------------
  # get common variables from both NetCDF files
  # -------------------------------------------

  StartDate = datetime.datetime.now().strftime('%Y-%m-%d-%H:%M:%S')
  CommonVariableKeys = list( set(VariablesNetCDF1).intersection(VariablesNetCDF2))
  if len( CommonVariableKeys )<1:
    sys.stdout.write( ' \n No common variables keys between NetCDF files. Exiting ... \n ')
    sys.exit(1)

  # open up text file for writing pixel statistics
  sys.stdout.write('%s\n'%' ')
  sys.stdout.write('%s\n'%'  Running comparison at: ')
  sys.stdout.write('%s\n'%('    '+StartDate))
  sys.stdout.write('%s\n'%' ')
  sys.stdout.write('%s\n'%('  netcdf file 1: '+NetCDF1))
  sys.stdout.write('%s\n'%('  netcdf file 2: '+NetCDF2))
  sys.stdout.write('%s\n'%' ')
  sys.stdout.write('%s\n'%'   \n  ----------------------------------------------------------------------------- \n')

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
    NetCDF_array1 = np.array( Dataset1[Var1][:] )
    NetCDF_array2 = np.array( Dataset2[Var2][:] )

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
        NetCDF_array1[(NetCDF_array1==MissingValue)]=NoData

    if len(MissingValues2)>0:
      for MissingValue in MissingValues2:
        NetCDF_array2[(NetCDF_array2==MissingValue)]=NoData

    sys.stdout.write('%s\n'%('  '+str(VariableCounter)+'. '+Var1+'\n'))
    sys.stdout.write('%s\n'%('    NetCDF1 Mean,Min.,Max.,St. dev    : '+str(
        NetCDF_array1[(NetCDF_array1!=NoData)].mean())+', '+\
                str(NetCDF_array1[(NetCDF_array1!=NoData)].min())+', '+\
                str(NetCDF_array1[(NetCDF_array1!=NoData)].max())+', '+\
                str(NetCDF_array1[(NetCDF_array1!=NoData)].std())))
    sys.stdout.write('%s\n'%('    NetCDF2 Mean,Min.,Max.,St. dev    : '+str(
        NetCDF_array2[(NetCDF_array2!=NoData)].mean())+', '+\
                str(NetCDF_array2[(NetCDF_array2!=NoData)].min())+', '+\
                str(NetCDF_array2[(NetCDF_array2!=NoData)].max())+', '+\
                str(NetCDF_array2[(NetCDF_array2!=NoData)].std())))
    sys.stdout.write('%s\n'%('    NetCDF1 number of elements        : '+\
            str(NetCDF_array1[(NetCDF_array1!=NoData)].size)))
    sys.stdout.write('%s\n'%('    NetCDF2 number of elements        : '+\
            str(NetCDF_array2[(NetCDF_array2!=NoData)].size)))

    # ----------------------------------------------------
    # get absolute value of differences between two arrays
    # ----------------------------------------------------
    Differences = np.absolute( NetCDF_array1 - NetCDF_array2 )
   
    '''
    Below we are trying to make sure that after creating a new array
    that contains the differences between the two arrays, that 
    anywhere there is a NoData value in either array, this is set
    in the differences array
    '''
    # -----------------------------------------------
    # if first netcdf array is scalar value, set 
    # it to 1 element array that is the NoData value
    # -----------------------------------------------
    if NetCDF_array1.size == 1 or NetCDF_array2.size == 1:
      if NetCDF_array1 == NoData or NetCDF_array2 == NoData:
        Differences = np.array([NoData],dtype=float)
    else:
      Differences[(NetCDF_array1==NoData)]=NoData
    
    if Differences.size>1:
      sys.stdout.write('%s\n'%('    Difference number of elements     : '+\
              str(Differences[(Differences!=NoData)].size)))
      sys.stdout.write('%s\n'%('    Difference Mean,Min.,Max.,St. dev : '+str(
          Differences[(Differences!=NoData)].mean())+', '+\
                  str(Differences[(Differences!=NoData)].min())+', '+\
                  str(Differences[(Differences!=NoData)].max())+', '+\
                  str(Differences[(Differences!=NoData)].std())))
    else:
      sys.stdout.write( '%s\n'%('    Differences number of elements  : NoData'))
      sys.stdout.write( '%s\n'%('    Difference Mean,Min.,Max.,St. dev : NoData,NoData,NoData,NoData'))
    sys.stdout.write('%s\n'%'   \n  ----------------------------------------------------------------------------- \n')

    # --------------------------
    # increment variable counter
    # --------------------------
    VariableCounter+=1

  # --------------------------
  # close both NetCDF datasets
  # --------------------------
  Dataset1.close()
  Dataset2.close()

def show_usage():
  print('''
    USAGE: $ python3 CompareNetcdf.py 
      --netcdf1 <netcdf1> 
      --netcdf2 <netcdf2>
  ''')
  sys.exit(1)

def show_version():
  print('''
    \n 
      CompareNetcdf.py
      Version 1.0.0

      @author  : Gerasimos Michalitsianos
      @updated : 12 March 2021
      @contact : gerasimosmichalitsianos@gmail.com
    \n
  ''')
  sys.exit(1)

def main():

  # ---------------------
  # parse input arguments
  # ---------------------
  netcdf1=''
  netcdf2=''
  parser = argparse.ArgumentParser( description='python3 NetCDF comparison tool.' )
  parser.add_argument( '-n1','--netcdf1',required=False,type=str,
    dest='netcdf1',help='first netcdf file')
  parser.add_argument( '-n2','--netcdf2',required=False,type=str,
    dest='netcdf2',help='second netcdf file')
  parser.add_argument('-v','--version',required=False,
    dest='version',help='show version.',action='store_true')
  parser.add_argument('-u','--usage',required=False,
    dest='showhelp',action='store_true',help='show help message.')
  args = parser.parse_args()

  version  = args.version
  showHelp = args.showhelp
  netcdf1  = args.netcdf1
  netcdf2  = args.netcdf2

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

  # ------------------------------------------
  # make sure both netcdf files were passed-in
  # via the command-line
  # ------------------------------------------
  if netcdf1 == '' or netcdf2 == '':
    print( '  \n Please pass in two netcdf files.' )
    print( '    Use the --netcdf1 (-n1) and --netcdf2 (-n2) command-line flags. Exiting ... ' )
    sys.exit(1)

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
  CompareNetcdfs( netcdf1,netcdf2 )

if __name__ == '__main__':
  main()
