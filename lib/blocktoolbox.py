"""
______________________________________
//////////////        \\\\\\\\\\\\\\\\
////////                     \\\\\\\\\
||||||||  BLOCKTRACK LIBRARY  ||||||||
\\\\\\\\_____          ______/////////
\\\\\\\\\\\\\\________////////////////

Author: Michele Filippucci, UniTN - IUSS Pavia
With the help and avice of: Paolo Davini, CNR-Isac

This library is a set of tools for the analysis of atmospheric blocking in the northern hemisphere.
The index used for atm blocking diagnostic is described in "Davini et al. - 2012 - Bidimensional diagnostics, variability, and trends of northern hemisphere blocking". Some differences and features are added: the persistence and area criteriamare applied at the level of tracking. Tracking also allow the user to perform lagrangian analysis.

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
DAV is a matrix with the same shape of the zg matrix limited to 20-70 nord latitudes. It is False (or "0") where there is no blocking 
and it is True (or "1") where there is.
As an alternative it is possible to change the flag "data_return" and the function will return a dataset object from the class xarray containing 
the same additional attributes
dataset: input dataset
fn_out: location and name of the new dataset. for example "user/home/dataset.nc"
data_return: when False fn_out is used, otherwise it is returned.
freq_also: if True the frequency of blocking (time mean of DAV) is calculated and stored in the netcdf
mer_gradient_filter: if True a filter for avoiding the detection of equatorial cut of lows is applied. The filter and its functioning is
described in Davini et al. 2012.

_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/

_ _ _ TM(dataset,output) _ _ _ 

Tibaldi and Molteni Index This function takes a .nc file containing z500 variable and computes the Tibaldi and Monteni index for the latitude 
60° N. It outputs a .dat file containing the design matrix (features, boolean label).
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

np.set_printoptions(precision=2,threshold=np.inf)

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
        area += np.cos(np.deg2rad(lats[0] + jlat*grid))*(grid*circearth/360)**2
        
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
        ):
  
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
    string = "zg variable was not found."
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
        #the feature matrix is composed of 12 geopotential values 15 deg north and south
        #of the ispected grid point.
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


def ContourTracking2D(dataset,fn_out = "",var_name = "DAV",data_return = False,
                      char_dic=True,save_track = True, fn_dic= ''):
  if fn_out=="" and data_return==False:
    string = "Specify the kind of output you want"
    print(string)
    return 0
  try:
    pIB_boolean = dataset[var_name]
    #initialize the array for performing the tracking
    arr = pIB_boolean.values
  except:
    print("Error Code 1: dataset not valid. The variable " + var_name + " cannot be found")
    return 1
  
  #initialiazing a blocking characteristics dictionary
  dic = []

  print("__Starting a Tracking process__")
  print("input: " + var_name)

  #loop over time
  #initialize some varaibles
  max = 0
  times = dataset.time.values
  #store the lat and lon dimensions for later computations.
  #more robust than arr.shape, as the order of coordinates may vary
  lastlat = len(dataset.lat.values) -1
  lastlon = len(dataset.lon.values) -1
  print("connected component analysis")
  for t in tqdm(np.arange(0,len(times)-1)):
    if t > 0:
      tmp = np.amax(arr[t-1,:,:])
      if max < tmp: #update maximum value in matrix
        max = tmp
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
        arr[t,:,:] = np.where(arr[t,:,:] == arr[t,j,lastlon], arr[t,j,0], arr[t,:,:])

    """
    TRACKING IN TIME
    """
    if t > 0:
      diff = times[t]-times[t-1]
      lbl = 1
      bool1 = arr[t-1,:,:] > 0
      bool2 = arr[t,:,:] > 0
      comp1 = np.unique(arr[t-1,bool1])
      comp2 = np.unique(arr[t,bool2])
      for l1 in comp1:
        boolarr1 = arr[t-1,:,:] == l1
        for l2 in comp2:
          #first we use a filter for avoiding lagrangian tracking between distant days
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

  print('rearranging indexes')
  arr[:,:,:] = OrderIndexes(arr[:,:,:])
  #create a dictionary where the number of items is equal to the number of labels
  for l in np.unique(arr):
    dic.append({})

  if char_dic == True:
    #compute blocking events characteristics
    """
    FILTERS MODULE
    """
    #initialize dictionary
    for l in np.unique(arr).astype(int):
      dic[l]['persistence'] = 0
      dic[l]['avg_area'] = 0
      dic[l]['avg_aspect_ratio'] = 0
      dic[l]['distance_traveled'] = 0
      dic[l]['track'] = []
      dic[l]['date'] = ''
    print('calculating characteristics')
    print('persistence, track, date:')
    #PERSISTENCE MODULE
    past_events = []
    for t in tqdm(np.arange(0,len(times)-1)):
      bool = arr[t,:,:] > 0
      today_events = np.unique(arr[t,bool]).astype(int) #labels at day t
      for l in today_events:
        if l not in past_events:
          past_events.append(l)
          if save_track == True:
            dic[l]['track'] = CenterofMass(arr[t:,:,:],l,grid=2.5) #center of mass traj
            xs,ys = dic[l]['track']
            dist = 0
            for i in range(len(xs)-1):
              dist += ((xs[i+1]-xs[i])**2 + (ys[i+1]-ys[i])**2)**0.5
            dic[l]['distance_traveled'] = dist
            dic[l]['date'] = times[t]
        dic[l]['persistence'] += 1

    print('area and longitudinal extent')
    #AREA and LONGITUDINAL EXTENT MODULE
    for t in tqdm(np.arange(0,len(times)-1)):
      bool = arr[t,:,:] > 0
      today_events = np.unique(arr[t,bool]).astype(int)
      for l in today_events:
        boolarr = arr[t,:,:] == l
        area = Area(arr[t,:,:],boolarr)
        lon_ext = 0
        for i in range(boolarr.shape[1]):
          if np.any(boolarr[:,i]):
            lon_ext += 1
        lat_ext = 0
        for i in range(boolarr.shape[0]):
          if np.any(boolarr[i,:]):
            lat_ext += 1
        #updating aspect_ratio
        lat_pos = sum(dic[l]['track'][0])/len(dic[l]['track'][0])
        dic[l]['avg_aspect_ratio'] += ((lon_ext*np.cos(np.deg2rad(lat_pos)))/lat_ext)/dic[l]['persistence']#the longitudinal dimension of a single grid diminuish with the latitude
        #updating area
        dic[l]['avg_area'] += area/dic[l]['persistence'] 

    #save dictionary
    np.save(fn_dic,dic)

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
  Update DAV_freq (if present) after area and persistence filter are applied.
  """
  dataset[var_name + "_freq"] = xr.where(dataset['DAV_tracked']>0,1,0).mean(dim="time")

  #output data
  if data_return == False:
    print("saving netcdf in: " + fn_out)
    dataset.to_netcdf(fn_out)
  if data_return == True:
    return dataset
  else:
    return 0


def FilterEvents(ds,fn_dic='',fn_out = "",fn_dic_out="",var_name = "DAV_tracked",data_return = False,
                  pers_min = 5,pers_max = 25,min_area = 500000,max_area=4e6,max_ar=5):

  print("__Starting a Filtering process__")
  print("input: " + var_name)

  try:
    pIB_boolean = ds[var_name]
    #initialize the array for performing the tracking
    arr = pIB_boolean.values
  except:
    print("Error Code 1: dataset not valid. The variable " + var_name + " cannot be found")
    return 1

  #import dictionary
  dic = np.load(fn_dic,allow_pickle=True)

  print('applying filters')
  to_pop = []
  to_retain = []
  l = 0
  for event in tqdm(dic):
    if event['persistence'] < pers_min or event['persistence'] > pers_max or event['avg_area'] < min_area or event['avg_area'] > max_area or event['avg_aspect_ratio'] > max_ar:
      to_pop.append(l)
      #I should find a better way to do this, as it isn't efficient at all. For example I could use the time information stored in the dictionary
      arr = np.where(arr==l,0,arr)
    else:
      to_retain.append(l)
    l+=1 #didn't use enumerate to use the progress bar.


  #didn't find another way to update the dictionary. It is a bit strange
  dic_filtered = []
  for l,event in enumerate(dic):
    if l in to_retain:
      dic_filtered.append(event)

  print('rearranging indexes')
  arr = OrderIndexes(arr)

  #save dictionary
  if fn_dic != '':
    np.save(fn_dic,dic_filtered)

  print("number of labels: " + str(np.amax(arr)))

  #initialize coords for new .nc
  times = ds["time"].values
  lon = ds["lon"].values
  lat = ds["lat"].values

  #initialize dataset object for the new .nc
  DAV_tracked = xr.DataArray(0,coords=[times,lat,lon],dims = pIB_boolean.dims)
  DAV_tracked[:,:,:] = arr

  #assign new data_array to dataset
  ds["DAV_tracked"] = DAV_tracked
  
  """
  Update DAV_freq (if present) after area and persistence filter are applied.
  """
  ds["DAV_freq"] = xr.where(ds['DAV_tracked']>0,1,0).mean(dim="time")

  #output data
  if data_return == False:
    print("saving netcdf in: " + fn_out)
    ds.to_netcdf(fn_out)
  if data_return == True:
    return ds
  else:
    return 0



