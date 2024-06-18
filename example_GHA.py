#my own class
import sys
sys.path.append("/home/michele/github/blocktrack/lib") #path to your '/blocktrack/lib' folder
import blocktoolbox as BT

#other libraries needed
import xarray as xr
import numpy as np
import time
start_time = time.time()

print("Starting operation")

data_dir='/archive/michele/data/ERA5/'
fn = data_dir + "ERA5_northem_2.5x2.5_zg_daily_1960-1965.nc"
work_dir='/home/michele/DATA/blocktrack_debug/'
fn_out_unfilt = work_dir + "ERA5_northem_2.5x2.5_zg_DAV_tracked_unfilt_daily_1960-1965.nc"
fn_out = work_dir + "ERA5_northem_2.5x2.5_zg_DAV_tracked_daily_1960-1965.nc"
fn_dic_unfilt = work_dir + "dictionaries/ERA5_northem_2.5x2.5_zg_DAV_tracked_unfilt_daily_1960-1965_DIC.npy"
fn_dic = work_dir + "dictionaries/ERA5_northem_2.5x2.5_zg_DAV_tracked_daily_1960-1965_DIC.npy"

ds = xr.load_dataset(fn)
print("Data correctly read")
ds = BT.GHA(ds,multiplicative_threshold = 1.26)
print(np.amax(ds['GHA_freq']))
print("GHA function correctly executed")
ds,dic_unfilt = BT.ContourTracking2D(ds,var_name='GHA')
#saving temporary data
ds.to_netcdf(fn_out_unfilt)
np.save(fn_dic_unfilt,dic_unfilt)
ds,dic_filt = BT.FilterEvents(ds,dic_unfilt,var_name='GHA',pers_min = 5,min_area = 5e5,max_avg_dist=1000)
print("ContourTracking2D function correclty executed")
#saving data
np.save(fn_dic,dic_filt)
ds.to_netcdf(fn_out, encoding={'time': {'dtype': 'i4'}})
print("Data created")
  
print("--- %s seconds ---" % (time.time() - start_time))