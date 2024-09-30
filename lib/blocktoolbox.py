"""
______________________________________
//////////////        \\\\\\\\\\\\\\\\
////////                     \\\\\\\\\
||||||||      BLOCKTRACK      ||||||||
\\\\\\\\_____          ______/////////
\\\\\\\\\\\\\\________////////////////

Author: Michele Filippucci, UniTN - IUSS Pavia
With the help and advice of: Paolo Davini, CNR-Isac

This code is a set of tools for the Lagrangian analysis of atmospheric blocking in the Northern Hemisphere.
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

'''
_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/

Area: this method calculate the area of a portion of a matrix in km^2.The portion of the matrix is given as a boolean array which is True 
where the area has to be calculated. lons and lats lists define the boundaries of the original array grid defines the dimension of a grid point.
! note that the array is (lat,lon)

Inputs:
boolarr: Portion of the matrix, Boolean Numpy Array
lats: Latitudinal boundaries of the boolarr, List
lons: Longitudinal boundaries of the boolarr, List
grid: grid dimension needed to convert the index in latitudes and longitudes, Float

Returns:
Area: The area of the given region expressed in km^2, Float
lon_ext: The longitudinal extent of the given region expressed in km^2, Float
lat_ext: The latitudinal extent of the given region expressed in km^2, Float

'''

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
  
'''
_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/

OrderIndexes: This function is necessary for the tracking algorithm to work. Its purpose is to take a tracked matrix that is zero where there is no blocking
and that is = label where there is a labeled blocking event and re-order and re-assign the label so that they are unique and increasingly ordered.

Inputs:
arr : tracked matrix, Numpy Array

Returns:
arr: reordered tracked matrix, Numpy Array

'''

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

'''
_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/

CenterofMass: This function is based on the homonym scipy function and returns a list of xs (longitudes) and ys (latitudes) that are the center of mass
coordinates for each time-step of the dataset.
NOTE: this method doesn't consider the sphericity of Earth, so the center of mass results slightly displaced north from its actual position.

Inputs:
tuple: tracked matrix, Numpy Array
label: label of the event, Integer
grid: grid dimension needed to convert the index in latitudes and longitudes., Float

Returns:
x: list of latitudinal positions of the center of mass, List
y: list of longitudinal positions of the center of mass, List
'''


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

'''
_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/

DAV: This function computes the blocked grid points in a gridded dataset following the geopotential height gradient reversal index
described in Davini et al. (2012). The blocked grid points are marked with a '1' boolean value in a matrix with the same dimension
as the geopotential height that is '0' everywhere else. The matrix is stored in a copy dataset of the input dataset, which is then
given as an output.

Inputs:
dataset: input dataset that must contain the daily geopotential height at 500hPa in the area [-180,180]°lon [0,90]°lat. The rank of
the input data must be 2, Xarray Dataset
mer_gradient_filter: flag that determines whether the meridional gradient filter is applied, as described in Davini et al. 2012, 
Boolean

Returns:
dataset: the input dataset + the matrix defining the blocked grid points

'''

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
  
'''
_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/

GHA: This function computes the blocked grid points in a gridded dataset following the geopotential height anomaly index
described in Woollings et al. 2018. The blocked grid points are marked with a '1' boolean value in a matrix with the same dimension
as the geopotential height that is '0' everywhere else. The matrix is stored in a copy dataset of the input dataset, which is then
given as an output.

Inputs:
dataset: input dataset that must contain the daily geopotential height at 500hPa in the area [-180,180]°lon [0,90]°lat. The rank of
the input data must be 2, Xarray Dataset
multiplicative_threshold: a number that determines how large the anomaly should be in terms of sigmas 
(anomaly > multiplicative_threshold*sigma), Float

Returns:
dataset: the input dataset + the matrix defining the blocked grid points
'''

def GHA(dataset,
        multiplicative_threshold = 1.26, 
        ):
  print("__Starting a GHA process__")
  print("input: zg500 , multiplicative threshold = " + str(multiplicative_threshold))
  bound_up=80
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
  #define zg and gha in the index domain.
  #converting everything to numpy array to optimize computational expense.
  zg_reduced=zg.loc[:,bound_down:bound_up,:].values
  gha_reduced=np.zeros(zg_reduced.shape)
  #define the domain of the reference distribution.
  zg_ref=zg.loc[:,45:80,:].values
  zg = zg.values
  
  for t in tqdm(range(45,len(times)-45)):
    #compute the anomalies with respect to a three months period for the index domain
    anomalies = zg_reduced[t,:,:] - np.mean(zg_reduced[t-45:t+45,:,:],axis=0)
    #compute the reference anomaly distribution for the three months period
    anomalies_ref = zg_ref[t-45:t+45,:,:] - np.repeat(np.expand_dims(np.mean(zg_ref[t-45:t+45,:,:],axis=0),axis=0),90,axis=0)
    ref_mean = np.mean(anomalies_ref,axis=(0,1,2))
    ref_std = np.std(anomalies_ref,axis=(0,1,2))
    #compare the anomalies at time t with the reference distribution.
    gha_reduced[t,:,:] = np.where(anomalies-ref_mean>multiplicative_threshold*ref_std,1,0)
    
  gha = xr.DataArray(data=np.zeros(zg.shape), 
                     dims=["time","lat","lon"],
                     coords = dict(time=dataset["time"],lat=dataset["lat"],lon=dataset["lon"]))
  gha.loc[:,bound_down:bound_up,:] = gha_reduced
  dataset = dataset.assign(GHA=gha)
  
  gha_freq = sum(gha)*100/gha.values.shape[0]
  dataset = dataset.assign(GHA_freq = gha_freq)
  return dataset
  
'''
_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/

LWAA: This function computes the blocked grid points in a gridded dataset following an index similar to GHA but taking the 
Local Wave Activity as an input variable. The blocked grid points are marked with a '1' boolean value in a matrix with the same dimension
as the geopotential height that is '0' everywhere else. The matrix is stored in a copy dataset of the input dataset, which is then
given as an output.

Inputs:
dataset: input dataset that must contain the daily Local Wave Activity in the area [-180,180]°lon [0,90]°lat. The rank of
the input data must be 2, Xarray Dataset
multiplicative_threshold: a number that determines how large the anomaly should be in terms of sigmas 
(anomaly > multiplicative_threshold*sigma), Float

Returns:
dataset: the input dataset + the matrix defining the blocked grid points
'''

def LWAA(dataset,
        multiplicative_threshold = 1.26, 
        ):
  print("__Starting a LWA anomaly process__")
  print("input: zg500 , multiplicative threshold = " + str(multiplicative_threshold))
  bound_up=80
  bound_down=30
  #checking if dataset is right
  
  try:
  #add sign - to keep computatoins similar to zg
    lwa = - dataset["lwa"]
    string = "data successfully received"
  except:
    string = "lwa variable was not found.\n\ Hint: check the content of your dataset."
    print(string)
    return 0

  #define XArray using the costructor for an appropriate output
  times = lwa.coords["time"]
  lon = lwa.coords["lon"]
  lat = lwa.coords["lat"]

  #define lwa and LWAA in the index domain.
  #converting everything to numpy array to optimize computational expense.
  lwa_reduced=lwa.loc[:,bound_down:bound_up,:].values
  lwaa_reduced=np.zeros(lwa_reduced.shape)
  #define the domain of the reference distribution.
  lwa_ref=lwa.loc[:,45:80,:].values
  lwa = lwa.values
  
  for t in tqdm(range(45,len(times)-45)):
    #compute the anomalies with respect to a three months period for the index domain
    anomalies = lwa_reduced[t,:,:] - np.mean(lwa_reduced[t-45:t+45,:,:],axis=0)
    #compute the reference anomaly distribution for the three months period
    anomalies_ref = lwa_ref[t-45:t+45,:,:] - np.repeat(np.expand_dims(np.mean(lwa_ref[t-45:t+45,:,:],axis=0),axis=0),90,axis=0)
    ref_mean = np.mean(anomalies_ref,axis=(0,1,2))
    ref_std = np.std(anomalies_ref,axis=(0,1,2))
    #compare the anomalies at time t with the reference distribution.
    lwaa_reduced[t,:,:] = np.where(anomalies-ref_mean>multiplicative_threshold*ref_std,1,0)
    
  lwaa = xr.DataArray(data=np.zeros(lwa.shape), 
                     dims=["time","lat","lon"],
                     coords = dict(time=dataset["time"],lat=dataset["lat"],lon=dataset["lon"]))
  lwaa.loc[:,bound_down:bound_up,:] = lwaa_reduced
  dataset = dataset.assign(LWAA=lwaa)
  
  lwaa_freq = sum(lwaa)*100/lwaa.values.shape[0]
  dataset = dataset.assign(LWAA_freq = lwaa_freq)
  return dataset

'''
_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/

ContourTracking2D: This function performs the tracking of the blocking events. It takes as input the output dataset of any of 
the functions DAV, GHA or LWAA. It then returns a dataset similar to the input dataset but with an additional variable that
has the same dimensions as the geopotential height field and identifies each blocked grid cell with a label that is unique for
each blocking event. In addition to this, this function returns a dictionary object where all the characteristics of the blocking
events are stored.

Inputs:
dataset: input dataset that must be the output of any of the functions DAV, GHA or LWAA, Xarray Dataset
var_name: the name of the variable contained inside of the dataset ('DAV', 'GHA' or 'LWAA'), String
geop_name: the name of the geopotential height contained inside of the dataset. This is needed for computing the intensity, String
overlap: the overlap criteria sets the portion of area that two consecutive blocking days must share to be considered as the same
block. It must be in the range [0,1], Float
save_track: flag to set whether the center of mass has to be computed.
save_intensity: flag to set whehter the intensity has to be computed.

Returns:
dataset: the input dataset + the tracked events matrix
dic: a list of dictionaries containing the features of the tracked events. The element 0 of the list is null and the index of the
other elements correspond to the label of the blocking event. The structure of one element of the list is the following:
dic[event_label]:
{'persistence': the number of days the blocking event lasts. Int,
'avg_area': the average area of the blocking event during its life cycle. Float [km^2], 
'avg_aspect_ratio': the average aspect ratio (lon_ext/lat_ext) of the blocking event during its life-cycle. Float, 
'distance_traveled': the total distance traveled during the event life-cycle. Float [km], 
'track': the track of the center of mass of the blocking event. (xs (List),ys(List)), 
'date': the date of the first day of the blocking event. np.DataTime64 (first blocking day),
'intensity': the average magnitude of the geop height anomaly in the blocked area during the blocking life-cycle. Float [m], 
'avg_dist_traveled': distance_traveled/persistence. Float [km]}

'''


def ContourTracking2D(dataset,var_name = "DAV",geop_name='zg',overlap=0.5,save_track = True,save_intesity=True):
  try:
    pIB_boolean = dataset[var_name]
    zg = dataset[geop_name].values
    zg_clim = zg.mean(axis=0)
    #initialize the array for performing the tracking
    arr = pIB_boolean.values
  except:
    print("Error Code 1: dataset not valid. The variable " + var_name + " or " + geop_name + " cannot be found")
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
      try:
        diff = (times[t]-times[t-1])
        diff = int(diff)/(1e9*60*60*24)
      except:
        diff = (times[t]-times[t-1]).days
      bool1 = arr[t-1,:,:] > 0
      bool2 = arr[t,:,:] > 0
      comp1 = np.unique(arr[t-1,bool1])
      comp2 = np.unique(arr[t,bool2])
      for l1 in comp1:
        boolarr1 = arr[t-1,:,:] == l1
        for l2 in comp2:
          #first we use a filter for avoiding lagrangian tracking between distant days
          #diff = diff.days
          #diff = diff/(1e9*60*60*24) #conversion from ms to days
          if diff > 1:
            break
          #then we find link between clusters
          boolarr2 = arr[t,:,:] == l2
          boolarr = boolarr1*boolarr2
          n = np.count_nonzero(boolarr)
          n_ex = np.count_nonzero(boolarr1)
          n_new = np.count_nonzero(boolarr2)
          if n > n_ex*overlap:# or n > n_new*overlap: #overlap criterium
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
    dic[l]['area'] = []
    dic[l]['aspect_ratio'] = []
    dic[l]['distance_traveled'] = []
    dic[l]['track'] = []
    dic[l]['date'] = ''
    dic[l]['intensity'] = []
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
          dist = []
          for i in range(len(xs)-1):
            lon2km_coeff = np.cos(np.deg2rad(np.mean([ys[i+1],ys[i]])))*111.320
            lat2km_coeff = 110.574
            if xs[i+1]*xs[i] > 0: #same sign  
              dist.append((((xs[i+1]-xs[i])*lon2km_coeff)**2 + ((ys[i+1]-ys[i])*lat2km_coeff)**2)**0.5)
            if xs[i+1]*xs[i] < 0 and abs(xs[i])<=100: #different sign near 0°lon
              dist.append((((xs[i+1]-xs[i])*lon2km_coeff)**2 + ((ys[i+1]-ys[i])*lat2km_coeff)**2)**0.5)
            if xs[i+1]*xs[i] <= 0 and abs(xs[i])>100 and xs[i] > 0: #different sign->boundary of the domain
              dist.append((((xs[i+1]-xs[i]+360)*lon2km_coeff)**2 + ((ys[i+1]-ys[i])*lat2km_coeff)**2)**0.5)
            if xs[i+1]*xs[i] <= 0 and abs(xs[i])>100 and xs[i] <= 0: #different sign->boundary of the domain
              dist.append((((xs[i+1]-xs[i]-360)*lon2km_coeff)**2 + ((ys[i+1]-ys[i])*lat2km_coeff)**2)**0.5)
          dic[l]['distance_traveled'] = dist
          dic[l]['date'] = times[t]
          dic[l]['time'] = t

  print('area, intensity and longitudinal extent')
  #AREA and LONGITUDINAL EXTENT MODULE
  for t in tqdm(range(len_time)):
    bool = arr[t,:,:] > 0
    today_events = np.unique(arr[t,bool]).astype(int)       #remove '0' from the list
    for l in today_events:
      boolarr = arr[t,:,:] == l
      area,lon_ext,lat_ext = Area(boolarr)
      if save_intesity==True:
        #the intensity is computed through the average anomaly associated to the blocked grid cells.
        an_tmp = abs(np.reshape(zg[t,:,:]-zg_clim,boolarr.shape))
        dic[l]['intensity'].append(np.sum(an_tmp[boolarr])/(len(an_tmp[boolarr])))
      #updating aspect_ratio
      dic[l]['aspect_ratio'].append((lon_ext/lat_ext))
      #updating area
      dic[l]['area'].append(area)
      
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

'''
_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/

FilterEvents: This function applies a set of filters to the blocking events detected through ContourTracking2D. It takes as an input
both the dataset and the dictionaries produced by the tracking algorithm and it updates them following a persistence, area and
distance_traveled criteria.

Inputs:
dataset: input dataset that must be the output ContourTracking2D, Xarray Dataset
var_name: the name of the variable contained inside of the dataset ('DAV', 'GHA' or 'LWAA'), String
pers_min: the minimum persistence that a blocking event must have, Int
min_avg_area: the minimum area of a blocking event in km^2, Float
max_avg_dist: the maximum per day distance traveled by a blocking event in km, Float

Returns:
dataset: the fitered input dataset
dic: the filtered input dictionary
'''

def FilterEvents(ds,dic,var_name = "DAV",
                 pers_min = 5,min_avg_area = 500000,max_avg_dist=1000):

  print("__Starting a Filtering process__")
  print("input: " + var_name + "_tracked")

  try:
    pIB_boolean = ds[var_name + "_tracked"]
    #initialize the array for performing the tracking
    arr = pIB_boolean.values
  except:
    print("Error Code 1: dataset not valid. The variable " + var_name + "_tracked" + " cannot be found")
    return 1

  print('applying filters')
  to_retain = []
  l = 1
  for event in tqdm(dic[1:]): #skip the first empty element
    try:
      if event['persistence'] < pers_min or sum(event['area'])/len(event['area']) < min_avg_area or sum(event['distance_traveled'])/len(event['distance_traveled']) > max_avg_dist:
        ti = event['time']
        tf = event['time'] + event['persistence']
        arr[ti:tf,:,:] = np.where(arr[ti:tf,:,:]==l,0,arr[ti:tf,:,:])
  #      if event['persistence'] > 4:
  #        print(l)
  #        print(event)
  #        try:
  #          print(sum(event['distance_traveled'])/len(event['distance_traveled']))
  #          print(sum(event['area'])/len(event['area']))
  #        except:
  #          print('1 day event')
      else:
        to_retain.append(l)
      l+=1 #didn't use enumerate to use the progress bar.
    except:
      print(l)
  #   print(event)

  #didn't find another way to update the dictionary. It is a bit strange
  dic_filtered = []
  dic_filtered.append({})
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




