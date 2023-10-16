"""
______________________________________
//////////////        \\\\\\\\\\\\\\\\
////////                     \\\\\\\\\
||||||||  BLOCKTOOLS LIBRARY  ||||||||
\\\\\\\\_____          ______/////////
\\\\\\\\\\\\\\________////////////////

Author: Michele Filippucci on intership at ISAC-CNR (TO)

This library is a set of tools for the analysis of atmospheric blocking in the northern hemisphere.
The index used for atm blocking diagnostic is described in "Davini et al. - 2012 - Bidimensional diagnostics, variability, and 
trends of northern hemisphere blocking". Some differences and features are added: the persistence and area criteria
are applied at the level of tracking. Tracking also allow the user to perform lagrangian analysis.

This library was developed using daily datasets.

The requirements for this library are:
Python              3.8.10
xarray              0.18.2
numpy               1.20.3
scipy               1.6.3
tqdm                4.61.1
matplotlib          3.4.2
Cartopy             0.19.0.post1

______________________________________
//////////////        \\\\\\\\\\\\\\\\
////////                     \\\\\\\\\
||||||||  LIST OF FUNCTIONS:  ||||||||
\\\\\\\\_____          ______/////////
\\\\\\\\\\\\\\________////////////////

_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/

_ _ _ DAV(dataset,fn_out = "",data_return = False,freq_also = False,mer_gradient_filter = False,long_filter = False) _ _ _ 

This function is capable of creating an nc file identical to the input dataset (located in fn_out) with additional attributes:
DAV  and (when freq_also == True) DAV_freq. 
DAV is a matrix with the same shape of the zg matrix limited to 20-70 nord latitudes. It is zero where there is no blocking and it is 1 where 
there is.
As an alternative you can change the flag "data_return" and the function will return a dataset object from the class xarray containing 
the same additional attributes
dataset: input dataset
fn_out: location and name of the new dataset. for example "user/home/dataset.nc"
data_return: when False fn_out is used, otherwise it is returned.
freq_also: if True the frequency of blocking (time mean of DAV) is calculated and stored in the netcdf
mer_gradient_filter: if True a filter described in Davini et al. 2012 is applied. It's purpose is to avoid the contamination of the signal
from geopotential lows at the equator.
long_filter: (EXPERIMENTAL) this filter check if the blocking condition is satisfied by at least three longitudinally neighbour grid boxes.

_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/

_ _ _ TM(dataset,output) _ _ _ 

Tibaldi and Molteni Index This function takes a .nc file containing z500 variable and computes the Tibaldi and Monteni index for the latitude 
60° N. It outputs a .dat file containing the design matrix (features, boolean label) needed for training a neural network.
dataset: input dataset
output: .dat output

_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/

EXPERIMENTAL!

_ _ _ ContourTracking3D(dataset,fn_out = "",var_name = "pIB_boolean",data_return = False) _ _ _ 

This is an experimental tracking function that treats time as a third dimension equivalent to longitude and latitude. This is computationally
more efficient but it treats badly the merging and splitting of blocking events and non-continuous datasets.

_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/

_ _ _ ContourTracking2D(dataset,fn_out = "",var_name = "DAV",data_return = False, pers = 0) _ _ _ 

This is a tracking function. 2D stands for the number of dimension that the method label from scipy takes under consideration.
This function takes a .nc file containing the DAV attribute from the function DAV and creates a new .nc file containing an additional 
attribute called var_name + _tracked which is zero when blocking is not occuring and is n when the nth blocking
event is occuring.

dataset: input dataset
fn_out: output dataset
var_name: name of the boolean matrix in the input dataset
data_return: similar to DAV method data_return. If it is true no output is saved and the dataset is returned in the script.
pers: a persistency filter on pers days is applied.
min_area: minimum area (km^2) required to indentify an area as blocked on a single day.

_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/

_ _ _ StationaryWaves(data_array, collapse_time = True, plev = 50000 #default is 500hPa) _ _ _

This funtion takes a data_array containing the geopotential height at 500hPa and outputs a matrix that is the geopotential height minus
the zonal mean. In this wayh we can look at the Rossby stationary waves.
collapse_time : if True the data_array returned is averaged along time dimension.
plev : pressure level, default is 500hPa

_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/

_ _ _ ZonalMean(data_array,collapse_time=True,plev = 50000,no_plev = True) _ _ _ 

This function takes a data_array and returns its zonal mean.
data_array: can be 4 dimensional (with pressure) or 3 dimensional (time,lon,lat)
collapse_time: if True the data_array returned is averaged along time dimension.
plev: the pressure level if present.
no_plev: the pressure level if you want to select it during the computation

_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/

_ _ _ Bootstrap(datasets=[],phys_quantity = "zg",anomaly = True,quantity_mean = 0,N_sample = 1000,N_composite = 330) _ _ _ 

This function take a list of datasets and a physical quantity contained in the dataset and performs a bootstrap.
The function returns a matrix of 5 percentile values that can be compared with a field obtained from a composite to assert if the 
field is different enough to be considered as a signal.
datasets : list of dataset (for example list of netcdf of a ensemble run)
phys_quantity: physical quantity contained in the dataset
anomaly: boolean variable; True if we want to get the percentile of an anomaly from the mean
quantity_mean :if we want the anomaly we need the mean field. This is not calculated in the function because it would be computationally expensive
and it's often already calculated in the scripts where the function is employed.
N_sample: number of "random composites". In other words the number of extractions that is used to calculated the pdf fr the percentile.
N_composite: number of element of the random composite.

_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/

_ _ _ GeoAreaCoordinates(name="none") _ _ _ 

This function returns the coordinates plon,plat of a region.
name: implemented regions are "Europe","Atlantic_Ocean","Greenland".

_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/

_ _ _ GetIndex(ds,coord="",key="") _ _ _ 

This function find the index of the numpy array correspondent to a ceratin variable value. For example if we want to find the index correspondent
to a 40° latitude then we use coord = "lat" and key = "40".
ds : dataset with desired coordinates
coord : coordinate name
key : coordinate value

_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/

_ _ _ Area(arr,boolarr,lats=[0,90],lons=[-180,180],grid=2.5) _ _ _ 

This method calculate the area of a portion of a matrix in km^2.The portion of the matrix is given as a boolean array which is True 
where the area has to be calculated. lons and lats lists define the boundaries of the original array grid defines the dimension of a grid point.
! note that the array isn't (lon,lat) but (lat,lon) instead

_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/

_ _ _ OrderIndexes(arr) _ _ _ 

This function is necessary for the tracking algorithm to work. Its purpose is to take a tracked matrix that is zero where there is no blocking
and that is = label where there is a labeled blocking event and re-order and re-assign the label so that they are unique and increasingly ordered.
arr : tracked matrix

_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/

_ _ _ CenterofMass(tuple,label,grid=2.5) _ _ _ 

This function is based on the homonym scipy function and returns a list of xs (longitudes) and ys (latitudes) that are the center of mass
coordinates for each time-step of the dataset.
tuple: tracked matrix
label: label of the event
grid: grid dimension needed to convert the index in latitudes and longitudes.

_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/
"""




import numpy as np
import xarray as xr
from scipy.ndimage.measurements import label
from scipy.ndimage.measurements import center_of_mass
from tqdm import tqdm
import matplotlib.pyplot as plt
import cartopy.util as cutil
import sys
import math
# generate random integer values
from random import seed
from random import randint

np.set_printoptions(precision=2,threshold=np.inf)

def StationaryWaves(data_array,
              collapse_time = True,
              plev = 50000 #default is 500hPa
              ):
  lon = data_array["lon"]
  lat = data_array["lat"]
  time = data_array["time"]
  #print(len(data_array.values.shape))
  if len(data_array.values.shape) == 4:
    tuple = data_array.loc[:,plev,:,:].values
  if len(data_array.values.shape) == 3:
    tuple = data_array.values
  if len(data_array.values.shape) != 3 and len(data_array.values.shape) != 4 :
    print("Error : unknown data format")
    return 0
  zonal_mean = np.mean(tuple,axis = (0,2))
  plt.plot(zonal_mean)
  zonal = tuple.copy()
  #for i in range(tuple.shape[0]):
  for j in range(tuple.shape[1]):
    zonal[:,j,:] = tuple[:,j,:]-zonal_mean[j]
  if collapse_time:
    zonal_anomaly = np.mean(zonal,axis=0)
    #print(zonal_mean.shape)
    dataarray_za = xr.DataArray(data=zonal_anomaly,dims=["lat","lon"],coords = dict(lat=lat,lon=lon))
  if not collapse_time:
    dataarray_za = xr.DataArray(data=zonal,dims=["time","lat","lon"],coords = dict(time =time,lat=lat,lon=lon))
  dataarray_za.name = data_array.name
  return dataarray_za

def ZonalMean(data_array,
             collapse_time=True,
             plev = 50000,
             no_plev = True
             ):
  lon = data_array["lon"]
  lat = data_array["lat"]
  if no_plev == False:
    plev = data_array["plev"]
  time = data_array["time"]
  #print(len(data_array.values.shape))
  if no_plev :
    if len(data_array.values.shape) == 4:
      tuple = data_array.loc[:,plev,:,:].values
    if len(data_array.values.shape) == 3:
      tuple = data_array.values
    if len(data_array.values.shape) != 3 and len(data_array.values.shape) != 4 :
      print("Error : unknown data format")
      return 0
    zonal_mean = np.mean(tuple,axis = (0,2))
    dataarray_zm = xr.DataArray(data=zonal_mean,dims=["lat"],coords = dict(lat=lat))
  if not no_plev:
    tuple = data_array.values
    zonal_mean = np.mean(tuple,axis = (0,3))
    dataarray_zm = xr.DataArray(data=zonal_mean,dims=["lev","lat"],coords = dict(lev=plev,lat=lat))
  dataarray_zm.name = data_array.name
  return dataarray_zm
  

def Bootstrap(datasets=[],phys_quantity = "zg",quantity_mean = 0,anomaly = True,N_sample = 1000,N_composite = 330):
  means = []
  for j in range(N_sample):
    seed(j)
    matrix = []
    for i in range(N_composite):
      rand_ensemble = randint(0,len(datasets)-1)
      ds = datasets[rand_ensemble]
      rand_time = randint(0,len(ds["time"].values)-1)
      if len(ds[phys_quantity].values.shape) == 4:
        matrix.append(ds[phys_quantity].values[rand_time,0,:,:])
      if len(ds[phys_quantity].values.shape) == 3:
        matrix.append(ds[phys_quantity].values[rand_time,:,:])
    if(anomaly == True):
      means.append((sum(matrix)/len(matrix))-quantity_mean)
    if(anomaly == False):
      means.append((sum(matrix)/len(matrix)))    
  percentile = np.percentile(means,[2.5,97.5],axis=0)
  return percentile    
          


def GeoAreaCoordinates(name="none"):
  plon = 0
  plat = 0
  if name == "Europe":
    plon = 10
    plat = 55
  if name == "Atlantic_Ocean":
    plon = -40
    plat = 35
  if name == "Greenland":
    plon = -45
    plat = 65
  if name == "none" or (plon,plat) == (0,0):
    print("Error: no valid name recognized")
    
  return plon,plat
  
  

def GetIndex(ds,coord="",key=""):
  index = 0
  for x in ds[coord].values:
    if str(x) == key:
      break
    index += 1
  return index

def Area(arr,boolarr,lats=[0,90],lons=[-180,180],grid=2.5):
  #This method calculate the area of a portion of a matrix in km^2. 
  #The portion of the matrix is given as a boolean array which is True 
  #where the area has to be calculated.
  #lons and lats lists define the boundaries of the original array
  #grid defines the dimension of a grid point.
  #! note that the array isn't (lon,lat) but (lat,lon) instead
  area = 0
  circearth = 40075
  for ilon in range(len(arr[0,:])):
    for jlat in range(len(arr[:,0])):
      if boolarr[jlat,ilon] == True:
        area += np.cos((lats[0] + jlat*grid)*(math.pi/360))*(2.5*circearth/360)**2
        
  return area
  

def OrderIndexes(arr):
  boolarr = arr > 0
  newarr = arr[boolarr]
  newarr = np.unique(np.sort(newarr))
  newval = range(1,len(newarr)+1)
  for i in tqdm(range(0,len(newarr))):
    arr[arr == newarr[i]] = newval[i]
    #np.where(arr == newarr[i], newval[i],arr)
  return arr

"""
This function returns the coordinates in space of the center of mass
of a blocking event. The object returned is a list of the coordinates
at each time step.
"""

def CenterofMass(tuple,label,grid=2.5):
  time = np.shape(tuple)[0]
  x = [] # lats
  y = [] # lons
  bool = tuple==label
  for t in range(time):
    if True in bool[t,:,:] and not True in bool[t,:,0]:
      cm = center_of_mass(bool[t,:,:])
      cm = np.array([cm[0]*grid,cm[1]*grid-180])
      x.append(cm[0])
      y.append(cm[1])
    if True in bool[t,:,0]: #this is the case of an event on the boundary
      #first we double the matrix and center it on the boundary
      shp = np.shape(bool[t,:,:])
      updt = np.zeros(shp*np.array((1,2)))
      updt[:,:shp[1]]=bool[t,:,:]
      updt[:,shp[1]:]=bool[t,:,:]
      bool[t,:,:] = updt[:,int(shp[1]/2):int(shp[1]*3/2)]
      cm = center_of_mass(bool[t,:,:])
      #we then have to translate the matrix 180 eastward
      if cm[1] < 180/grid:
        cm = np.array([cm[0]*grid,cm[1]*grid-180+180])
        #print("minor: " + str(cm))
      else:
        if cm[1] >= 180/grid:
          cm = np.array([cm[0]*grid,cm[1]*grid-180-180])
        #print("major: " + str(cm))
      x.append(cm[0])
      y.append(cm[1])
      #this method is exact for even long-shape. When shape is odd there is
      #an error of 1.25 degrees.
  return x,y

"""
This function is capable of creating an nc file identical (located in fn_out)
to the BlockTools.dataset with additional attributes:
pIB_boolean and (when freq_also == True) pIB_frequencies
As an alternative you can change the flag "data_return" and the function
will return a dataset object from the class xarray containing the same
additional attributes
"""  
def DAV(dataset,fn_out = "",\
        data_return = False,\
        freq_also = False,\
        mer_gradient_filter = False,\
        long_filter = False):
  
  print("__Starting a DAV process__")
  print("input: zg500 , freq_also = " + str(freq_also) + ", mer_gradient_filter = " + str(mer_gradient_filter) )

  if fn_out=="" and data_return==False:
    string = "Specify the kind of output you want"
    print(string)
    return 0
  #checking if dataset is right
  try:
    zg = dataset["zg"]
    """
    if len(zg.shape)>3:
      try:
        dataset = dataset.mean["lev"]
      except:
        dataset = dataset.mean["plev"]
      zg = dataset["zg"]
      print("zg dimension reduced")
    """
    string = "data successfully received"
  except:
    string = "zg variable was not found.\n\ Hint: check the content of your dataset."
    print(string)
    return 0

  #define XArray using the costructor for an appropriate output
  times = zg.coords["time"]
  lon = zg.coords["lon"]
  lat = zg.coords["lat"]

  #____CHECK GEOP HEIGHT____
  #ERA5 dataset uses geopotential, not height
  if zg.values[0,0,0] > 10000:
      zg = zg/9.80665

  #.values gives tuples
  #compute GHGS/GHGN
  GHGS = (+ zg.loc[:,30.0:75.0,:].values - zg.loc[:,15.0:60.0,:].values)/15.0
  GHGN = (- zg.loc[:,30.0:75.0,:].values + zg.loc[:,45.0:90.0,:].values)/15.0

  #look for grid points where the conditions for pIB are satisfied
  #using where function from xarray
  #first term is GHGN condition, second term is GHGS condition. Tuples are multiplied
  #(no matrix product)
  if mer_gradient_filter == False:
    TuplepIB = xr.where(GHGN < -10.0, 1.0, 0.0) * xr.where(GHGS > 0., 1.0 , 0.0)
  #filter for meridional gradient
  if mer_gradient_filter == True:
    GHGS2 = (+ zg.loc[:,15:60,:].values - zg.loc[:,0:45,:].values)/15.0
    TuplepIB = xr.where(GHGS2 < -5,1.0,0.0)*xr.where(GHGN < -10.0, 1.0, 0.0)*\
               xr.where(GHGS > 0., 1.0 , 0.0)
  #check = TuplepIB
  #15 degrees continuous longitude filter
  if long_filter == True:
    temp1 = TuplepIB
    temp2 = temp1
    shift = 3
    for i in range(len(lon)):
      #shift 3 to the right
      if i < shift:
        temp1[:,:,i] = TuplepIB[:,:,len(lon)-shift+i]
      else:
        temp1[:,:,i] = TuplepIB[:,:,i-shift]
      #shift 3 to the left
      if i < len(lon)-shift:
        temp2[:,:,i] = TuplepIB[:,:,i+shift]
      else:
        temp2[:,:,i] = TuplepIB[:,:,i-len(lon)+shift]
    TuplepIB = temp1*temp2

  #define XArray using the costructor for an appropriate output
  DAV = xr.DataArray(data=np.zeros(zg.shape), 
                     dims=["time","lat","lon"],
                     coords = dict(time=dataset["time"],lat=dataset["lat"],lon=dataset["lon"]))
  #print(pIB)
  DAV.loc[:,30:75,:] = TuplepIB[:,:,:]  
  #print(pIB)
  	
  dataset = dataset.assign(DAV=DAV)
  
  if freq_also == True:
    DAV_freq = sum(DAV)*100/DAV.values.shape[0]
    dataset = dataset.assign(DAV_freq = DAV_freq)
    #dataset["DAV_freq"].name = "DAV_freq"
  if data_return == False:
    print("saving file in: " + fn_out)
    try:
      print(dataset)
      dataset.to_netcdf(fn_out)
    except:
      print("something went wrong")
  if data_return == True:
    return dataset
  else:
    return 0

"""
Tibaldi and Molteni Index
This function takes a .nc file containing z500 variable and computes the Tibaldi
and Monteni index for the latitude 60° N.
It outputs a .dat file containing the design matrix (features, boolean label) needed
for training a neural network.
"""

def TM(dataset,
       output):
  #checking if dataset is right
  try:
    zg = dataset["zg"]
    string = "data successfully received"
  except:
    string = "zg variable was not found.\n\ Hint: use read() to load data."
    print(string)
    return 0
  #____CHECK GEOP HEIGHT____
  #ERA5 dataset uses geopotential, not height
  if zg.values[0,0,0] > 10000:
      zg = zg.values/9.80665
  #.values gives tuples
  print(dataset["lat"])
  N = GetIndex(dataset,"lat","75.0") #north
  C = GetIndex(dataset,"lat","60.0") #center
  S = GetIndex(dataset,"lat","45.0") #south
  print(N,C,S)
  file = open(output, "a")
  for i in range(len(dataset["time"])):
    for j in range(len(dataset["lon"])):
      string = ""
      for k in range(12):
        string += str(zg[i,0,S+k,j]) + " "
      GHGS = (+ zg[i,C,k] - zg[i,S,k])/15.0
      GHGN = (- zg[i,C,k] + zg[i,N,k])/15.0
      flag = int(GHGN < -10.0) * int(GHGS > 0.)
      string += str(flag)
#        print(string)
      file.write(string + "\n")
  file.close()
  return 0

"""
Contour Tracking
This function takes a .nc file containing the pIB_boolean attribute from the function
boolean_pIB and creates a new .nc file containing an additional attribute called 
pIB_tracked which is zero when blocking is not occuring and is n when the nth blocking
event is occuring.
3D version uses label method on a 3D matrix (time,longitude,latitude)
2D version uses label method in a foor loop (on t) on a 2D matrix (longitude,latitude)
"""
def ContourTracking3D(dataset,fn_out = "",var_name = "pIB_boolean",data_return = False):
  if fn_out=="" and data_return==False:
    string = "Specify the kind of output you want"
    print(string)
    return 0
  try:
    pIB_boolean = dataset[var_name]
  except:
    print("Error Code 1: dataset not valid. The variable " + var_name + " cannot be found")
    return 1
  arr = pIB_boolean.values[:,0,:,:] 

  #filtra circa 9 punti griglia

  #label method from scipy.ndimage.measurements is used
  #structure = np.ones((3,3,3))
  structure = [[[0,0,0],[0,1,0],[0,0,0]],\
               [[0,1,0],[1,1,1],[0,1,0]],\
               [[0,0,0],[0,1,0],[0,0,0]]] #this matrix defines what is defined as neighbour
  #neighbour points are labeled with the same sequential number
  arr,ncomponents=label(arr,structure)
  #applying some filters
  for t in np.arange(0,len(dataset.time.values-1)):
    bool = arr[t,:,:] > 0
    list = np.unique(arr[t,bool])
    for l in list:
      boolarr = arr[t,:,:] == l
      n = np.count_nonzero(boolarr)
      #filtering cluster dimension
      if n < 9:
        arr[t,:,:] = xr.where(boolarr, 0,arr[t,:,:])

  arr = OrderIndexes(arr)
  #initialize coords for new .nc
  times = pIB_boolean.coords["time"].values
  #plev = pIB_boolean.coords["plev"].values
  lon = pIB_boolean.coords["lon"].values
  lat = pIB_boolean.coords["lat"].values
  #initialize dataset object for the new .nc
  pIB_tracked = xr.DataArray(0,coords=[times,lat,lon],dims = ['time','lat','lon'])
  pIB_tracked[:,:,:] = 0
  pIB_tracked[:,0,:,:] = arr
  #assign dataset to dataset which is now updated
  dataset = dataset.assign(pIB_tracked = pIB_tracked)

  #output data
  if data_return == False:
    dataset.to_netcdf(fn_out)
  if data_return == True:
    return dataset
  else:
    return 0

def ContourTracking2D(dataset,fn_out = "",var_name = "DAV",data_return = False, pers = 0,min_area = 500000,max_length=28):
  if fn_out=="" and data_return==False:
    string = "Specify the kind of output you want"
    print(string)
    return 0
  try:
    pIB_boolean = dataset[var_name]
  except:
    print("Error Code 1: dataset not valid. The variable " + var_name + " cannot be found")
    return 1
  
  print("__Starting a Tracking process__")
  print("input: " + var_name + ", pers = " + str(pers))

  #loop over time
  max = 0
  lastlat = len(dataset.lat.values) -1
  lastlon = len(dataset.lon.values) -1
  times = dataset.time.values
  print("connected component analysis")
  for t in tqdm(np.arange(0,len(times)-1)):
    if t > 0:
      tmp = np.amax(arr[t-1,:,:])
      if max < tmp: #update maximum value in matrix
        max = tmp
    arr = pIB_boolean.values[:,:,:]
    #label method from scipy.ndimage.measurements is used
    structure = [[0,1,0],\
                 [1,1,1],\
                 [0,1,0]] #this matrix defines what is defined as neighbour

    #neighbour points are labeled with the same sequential number
    arr[t,:,:],ncomponents=label(arr[t,:,:],structure)
    arr[t,:,:] = xr.where(arr[t,:,:] > 0, arr[t,:,:] + max , arr[t,:,:])

    #making it periodic in longitude
    for j in range(0,lastlat):
      if arr[t,j,lastlon] > 0 and arr[t,j,0] > 0:
        arr[t,:,:] = xr.where(arr[t,:,:] == arr[t,j,lastlon], arr[t,j,0], arr[t,:,:])

    #applying some filters
    bool = arr[t,:,:] > 0
    comp = np.unique(arr[t,bool])
    for l in comp:
      boolarr = arr[t,:,:] == l
      length = 0
      for i in range(boolarr.shape[1]):
        if np.any(boolarr[:,i]):
          length += 1
      area = Area(arr[t,:,:],boolarr)
      #filtering cluster dimension
      if area < min_area or length > max_length: 
        #the minimum area value has been chosen looking at present litterature
        #the maxmimum length has been set to 60 degrees, hence 24 cells with 2.5 deg resolution.
        arr[t,:,:] = xr.where(boolarr, 0,arr[t,:,:])
    """
    TRACKING IN TIME
    """
    if t > 0:
      lbl = 1
      bool1 = arr[t-1,:,:] > 0
      bool2 = arr[t,:,:] > 0
      comp1 = np.unique(arr[t-1,bool1])
      comp2 = np.unique(arr[t,bool2])
      for l1 in comp1:
        boolarr1 = arr[t-1,:,:] == l1
        for l2 in comp2:
          #first we use a filter for avoiding lagrangian tracking between distant days
          diff = times[t]-times[t-1]
          diff = int(diff)
          diff = diff/(1e9*60*60*24) #conversion from ms to days
          if diff > 1:
            break
          #then we find link between clusters
          boolarr2 = arr[t,:,:] == l2
          boolarr = boolarr1*boolarr2
          n = np.count_nonzero(boolarr)
          n_ex = np.count_nonzero(boolarr1)
          n_new = np.count_nonzero(boolarr2)
          if n > n_ex/2 or n > n_new/2: #50% overlap
            #new label which is always different
            arr[t,:,:] = xr.where(boolarr2,l1,arr[t,:,:])
        if diff > 1:
          break

  arr[:,:,:] = OrderIndexes(arr[:,:,:])

  """
  PERSISTENCE MODULE
  """
  if pers > 0:
    counter = 0
    safe = []
    print("persistence analysis")
    for t in tqdm(np.arange(0,len(dataset.time.values))):
      bool1 = arr[t,:,:] > 0
      try:
        bool2 = arr[t+pers,:,:] > 0
      except:
        arr[t:,:,:] = 0 #not possible to check so out of output
        print("exited with " + str(len(dataset.time.values)-t-pers) + " elements remaining")
        print(str(counter) + " blocking events where found")
        break
      comp1 = np.unique(arr[t,bool1]) #labels at day t
      comp2 = np.unique(arr[t+pers,bool2]) #labels at day t+pers
      for l1 in comp1:
        if not l1 in safe:
          if not l1 in comp2: #if pers days after there is no l1 the event is deleted
            arr[:,:,:] = xr.where(arr[:,:,:]==l1,0,arr[:,:,:])
          else:
            safe.append(l1) #if pers days after there is l1 the event is saved
            counter += 1
  print("ordering indexes")
  arr = OrderIndexes(arr)
  print("number of labels: " + str(np.amax(arr)))

  #initialize coords for new .nc
  times = pIB_boolean.coords["time"].values
  lon = pIB_boolean.coords["lon"].values
  lat = pIB_boolean.coords["lat"].values

  #initialize dataset object for the new .nc
  DAV_tracked = xr.DataArray(0,coords=[times,lat,lon],dims = pIB_boolean.dims)
  DAV_tracked[:,:,:] = arr

  #assign new data_array to dataset
  dataset["DAV_tracked"] = DAV_tracked
  
  """
  Update DAV and DAV_freq (if present) after area and persistence filter are applied.
  """
  dataset[var_name] = DAV_tracked > 0
  dataset[var_name + "_freq"] = dataset[var_name].mean(dim="time")

  #output data
  if data_return == False:
    print("saving netcdf in: " + fn_out)
    dataset.to_netcdf(fn_out)
  if data_return == True:
    return dataset
  else:
    return 0
