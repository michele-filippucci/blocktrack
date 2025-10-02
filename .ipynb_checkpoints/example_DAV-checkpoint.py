#my own class
import sys
sys.path.append("/home/michele/prog/github/blocktrack/lib") #path to your '/blocktrack/lib' folder
import blocktoolbox as bt

#other libraries needed
import xarray as xr
import time
start_time = time.time()

print("Starting operation")

data_dir='/archive/michele/data/ERA5/' #data directory
fn = data_dir + "ERA5_northem_2.5x2.5_zg_daily_1959-2021.nc"
work_dir='/home/michele/DATA/blocktrack_debug/' #output directory
fn_out_unfilt = work_dir + "ERA5_northem_2.5x2.5_zg_DAV_tracked_unfilt_daily_1940-2023.nc"
fn_out = work_dir + "ERA5_northem_2.5x2.5_zg_DAV_tracked_daily_1940-2023.nc"
#dictionaries for saving the tracking data. Here we create a separated directory to avoid confusion.
fn_dic_unfilt = work_dir + "dictionaries/ERA5_northem_2.5x2.5_zg_DAV_tracked_unfilt_daily_1940-2023_DIC.npy" 
fn_dic = work_dir + "dictionaries/ERA5_northem_2.5x2.5_zg_DAV_tracked_daily_1940-2023_DIC.npy"

ds = xr.open_dataset(fn)
print("Data correctly read")
ds = bt.DAV(ds,mer_gradient_filter = True)
print(np.amax(ds['DAV_freq']))
print("DAV function correctly executed")
ds,dic_unfilt = bt.ContourTracking2D(ds,grid=2.5)
#saving temporary data
ds.to_netcdf(fn_out_unfilt)
np.save(fn_dic_unfilt,dic_unfilt)
ds,dic_filt = bt.FilterEvents(ds,dic_unfilt,var_name = 'DAV',
                              pers_min = 5,min_avg_area = 500000)
print("ContourTracking2D function correclty executed")
#saving data
np.save(fn_dic,dic_filt)
ds.to_netcdf(fn_out)
print("Data created")
  
print("--- %s seconds ---" % (time.time() - start_time))