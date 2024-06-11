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

def Area(boolarr,lat_lim=[0,90],lon_lim=[-180,180],grid=2.5):
  #This method calculate the area of a portion of a matrix in km^2.
  #This function can also compute the latitudinal and longitudinal extent of the blocking event. 
  #The portion of the matrix is given as a boolean array which is True 
  #where the area has to be calculated.
  #lon_lim and lat_lim lists define the boundaries of the original array
  #grid defines the dimension of a grid point.
  #! note that the array isn't (lon,lat) but (lat,lon) instead
  area = 0
  
  #quantities needed to compute the longitudinal and latitudinal extent as well
  lats = []
  lons = []
  
  lon_len = len(boolarr[0,:])
  lat_len = len(boolarr[:,0])

  circearth = 40075
  for ilon in range(lon_len):
    for jlat in range(lat_len):
      if boolarr[jlat,ilon] == True:
        #extent
        lats.append(lat_lim[0] + jlat*grid)
        lons.append(lon_lim[0] + ilon*grid)
        #area        
        area += np.cos(np.deg2rad(lat_lim[0] + jlat*grid))*(grid*circearth/360)**2
  #finding the extent
  #here the latitudinal position of the event is identified as the mean of the latitudes list, coherently with the center of mass method
  lon_ext = (np.amax(lons) - np.amin(lons))*np.cos(np.deg2rad(np.mean(lats))) + grid
  lat_ext = np.amax(lats) - np.amin(lats) + grid
        
  return area,lon_ext,lat_ext
  

def OrderIndexes(arr):
  examined_idxs = []
  for t in tqdm(range(arr.shape[0])):
    today_idxs = np.unique(arr[t,:,:])
    #remove '0' from today indexes
    bool = today_idxs > 0
    today_idxs = today_idxs[bool]
    for idx in today_idxs:
      if idx in examined_idxs:
        arr[t,:,:][arr[t,:,:]==idx] = examined_idxs.index(idx)+1
      else:
        arr[t,:,:][arr[t,:,:]==idx] = len(examined_idxs)+1
        examined_idxs.append(idx)

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
to the BlockTools.dataset with additional attributes
"""  
def DAV(dataset,
        mer_gradient_filter = False
        ):
  
  print("__Starting a DAV process__")
  print("input: zg500 , mer_gradient_filter = " + str(mer_gradient_filter) )
  
  #checking if dataset is right
  try:
    zg = dataset["zg"]
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
  
  DAV_freq = sum(DAV)*100/DAV.values.shape[0]
  dataset = dataset.assign(DAV_freq = DAV_freq)
  return dataset
  
"""
Geopotential Height Anomaly (GHA) index.
This is a simple implementation of a Geopotential Height Anomaly index for istantaneous blocking detection. This computes the daily geopotential height anomaly for each grid point of a given dataset by comparing the grid value at each day with its mean value over a 90 days time window centered on the same day. Moreover, the algorithm computes the standard deviation of the geopotential height anomaly over the same window. Once these quantities are calculated, a grid point identied as blocked if the anomaly exceeds the standard deviation moltiplied by a multiplicative threshold that is given to the index as a input variable.
"""

def GHA(dataset,
        multiplicative_threshold = 1.26, 
        eulerian_persistence=0,
        ):
  print("__Starting a GHA process__")
  print("input: zg500 , multiplicative threshold = " + str(multiplicative_threshold) + ' , eulerian_persistence = ' + str(eulerian_persistence))
  bound_up=90
  bound_down=30
  #checking if dataset is right
  
  try:
    zg = dataset["zg"]
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
  zg_reduced=zg.loc[:,bound_down:bound_up,:].values
  gha_reduced=np.zeros(zg_reduced.shape)
  for t in tqdm(range(45,len(times)-45)):
    if eulerian_persistence==0:
      zg_I = zg_reduced[t,:,:]
      zg_tmp = zg_reduced[t-45:t+45,:,:]
      zg_mean = np.mean(zg_tmp,axis=0)
      zg_std = np.std(zg_tmp,axis=0)
      gha_reduced[t,:,:] = np.where(zg_I-zg_mean>multiplicative_threshold*zg_std,1,0)
    if eulerian_persistence>0:
      zg_I = zg_reduced[t:t+eulerian_persistence,:,:]
      zg_tmp = zg_reduced[t-45:t+45,:,:]
      zg_mean = np.mean(zg_tmp,axis=0)
      zg_std = np.std(zg_tmp,axis=0)
      gha_reduced[t,:,:] = np.prod(np.where(zg_I-np.repeat(np.expand_dims(zg_mean,0),eulerian_persistence,axis=0)>multiplicative_threshold*zg_std,1,0), axis=0)
      #we also want to considered as blocked the last days of blocking, hence:
      gha_reduced[t,:,:] = np.logical_or(gha_reduced[t,:,:],(zg_I[0,:,:]-zg_mean>multiplicative_threshold*zg_std)*gha_reduced[t-1,:,:])
    gha = xr.DataArray(data=np.zeros(zg.shape), 
                     dims=["time","lat","lon"],
                     coords = dict(time=dataset["time"],lat=dataset["lat"],lon=dataset["lon"]))
  gha.loc[:,bound_down:bound_up,:] = gha_reduced
  dataset = dataset.assign(GHA=gha)
  
  gha_freq = sum(gha)*100/gha.values.shape[0]
  dataset = dataset.assign(GHA_freq = gha_freq)
  return dataset
  
    

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

def ContourTracking2D(dataset,var_name = "DAV",save_track = True, fn_dic= ''):
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
  for t in tqdm(range(len(times))):
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

  #create a dictionary where the number of items is equal to the number of labels + 1
  #add an empty first element to create correspondence between idx and label
  dic.append({})
  #find the idxs
  idxs = np.unique(arr)
  bool = idxs > 0
  idxs = idxs[bool].astype(int)
  #initialize the remaining dictionaries
  for l in idxs:
    dic.append({})
    #compute blocking events characteristics
  """
  CHARACTERISTICS MODULE
  """
  #initialize dictionary
  for l in idxs:
    dic[l]['persistence'] = 0
    dic[l]['avg_area'] = 0
    dic[l]['avg_aspect_ratio'] = 0
    dic[l]['distance_traveled'] = 0
    dic[l]['track'] = []
    dic[l]['date'] = ''
    dic[l]['time'] = 0
  print('calculating characteristics')
  print('persistence, track, date:')
  #PERSISTENCE MODULE
  past_events = []
  len_time = len(times)
  for t in tqdm(range(len_time)):
    bool = arr[t,:,:] > 0
    today_events = np.unique(arr[t,bool]).astype(int) #labels at day t
    for l in today_events:
      dic[l]['persistence'] += 1
      if l not in past_events:
        past_events.append(l)
        if save_track == True:
          #100 is a safe value for making the algorithm a little more efficient, as there is no event longer than 100 days.
          dic[l]['track'] = CenterofMass(arr[t:min([t+100,len_time]),:,:],l,grid=2.5) #center of mass traj
          ys,xs = dic[l]['track']
          dist = 0
          for i in range(len(xs)-1):
            lon2km_coeff = np.cos(np.deg2rad(np.mean([ys[i+1],ys[i]])))*111.320
            lat2km_coeff = 110.574
            if xs[i+1]*xs[i] > 0: #same sign  
              dist += (((xs[i+1]-xs[i])*lon2km_coeff)**2 + ((ys[i+1]-ys[i])*lat2km_coeff)**2)**0.5
            if xs[i+1]*xs[i] <= 0 and abs(xs[i])>100 and xs[i] > 0: #different sign->boundary of the domain
              dist += (((xs[i+1]-xs[i]+360)*lon2km_coeff)**2 + ((ys[i+1]-ys[i])*lat2km_coeff)**2)**0.5
            if xs[i+1]*xs[i] <= 0 and abs(xs[i])>100 and xs[i] <= 0: #different sign->boundary of the domain
              dist += (((xs[i+1]-xs[i]-360)*lon2km_coeff)**2 + ((ys[i+1]-ys[i])*lat2km_coeff)**2)**0.5
          dic[l]['distance_traveled'] = dist
          dic[l]['avg_dist_traveled'] = dist/len(xs) 
          dic[l]['date'] = times[t]
          dic[l]['time'] = t

  print('area and longitudinal extent')
  #AREA and LONGITUDINAL EXTENT MODULE
  for t in tqdm(range(len_time)):
    bool = arr[t,:,:] > 0
    today_events = np.unique(arr[t,bool]).astype(int)       #remove '0' from the list
    for l in today_events:
      boolarr = arr[t,:,:] == l
      area,lon_ext,lat_ext = Area(boolarr)
      #updating aspect_ratio
      dic[l]['avg_aspect_ratio'] += (lon_ext/lat_ext)/dic[l]['persistence']
      #updating area
      dic[l]['avg_area'] += area/dic[l]['persistence']  

  print("number of labels: " + str(np.amax(arr)))

  #initialize coords for new .nc
  times = pIB_boolean.coords["time"].values
  lon = pIB_boolean.coords["lon"].values
  lat = pIB_boolean.coords["lat"].values

  #initialize dataset object for the new .nc
  DAV_tracked = xr.DataArray(0,coords=[times,lat,lon],dims = pIB_boolean.dims)
  DAV_tracked[:,:,:] = arr

  #assign new data_array to dataset
  dataset[var_name+"_tracked"] = DAV_tracked
  dataset[var_name] = xr.where(DAV_tracked>0,1,0)
  
  """
  Update DAV_freq (if present) after area and persistence filter are applied.
  """
  dataset[var_name + "_freq"] = xr.where(dataset[var_name+'_tracked']>0,1,0).mean(dim="time")*100

  #output data
  return dataset,dic



def FilterEvents(ds,dic,var_name = "DAV",
                 pers_min = 5,min_area = 500000,max_avg_dist=1000):

  print("__Starting a Filtering process__")
  print("input: " + var_name + "_tracked")

  try:
    pIB_boolean = ds[var_name + "_tracked"]
    #initialize the array for performing the tracking
    arr = pIB_boolean.values
  except:
    print("Error Code 1: dataset not valid. The variable " + var_name + "_tracked" + " cannot be found")
    return 1

  #import dictionary
  if dic == 0:
      dic = np.load(fn_dic,allow_pickle=True)

  print('applying filters')
  to_retain = []
  l = 1
  for event in tqdm(dic[1:]): #skip the first empty element
    if event['persistence'] < pers_min or event['avg_area'] < min_area or event['avg_dist_traveled'] > max_avg_dist:
      ti = event['time']
      tf = event['time'] + event['persistence']
      arr[ti:tf,:,:] = np.where(arr[ti:tf,:,:]==l,0,arr[ti:tf,:,:])
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

  print("number of labels: " + str(np.amax(arr)))

  #initialize coords for new .nc
  times = ds["time"].values
  lon = ds["lon"].values
  lat = ds["lat"].values

  #initialize dataset object for the new .nc
  DAV_tracked = xr.DataArray(0,coords=[times,lat,lon],dims = pIB_boolean.dims)
  DAV_tracked[:,:,:] = arr

  #assign new data_array to dataset
  ds[var_name+"_tracked"] = DAV_tracked
  
  """
  Update DAV_freq (if present) after area and persistence filter are applied.
  """
  ds[var_name + "_freq"] = xr.where(ds[var_name+"_tracked"]>0,1,0).mean(dim="time")*100

  #output data
  return ds,dic_filtered




