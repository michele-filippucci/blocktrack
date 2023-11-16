#my own class
import sys
sys.path.append("/home/michele/prog/github/blocktrack/lib") #path to your '/blocktrack/lib' folder
import blocktoolbox as bt

#other libraries needed
import xarray as xr
import time
start_time = time.time()

print("Starting operation")

fn = "/home/michele/DATA/ERA5/ERA5_northem_2.5x2.5_zg_daily_1959-2021_mod.nc" #path to your .nc input file
#Remember: it must contain a "zg" variable at 500 hPa pressure level, 2.5°x2.5° horizontal resolution and
#daily time-stepping. The computation is performed for the Northern Hemisphere only.
fn_out = "/home/michele/DATA/ERA5/ERA5_northem_2.5x2.5_DAV_daily_1959-2021.nc" #path to your .nc output file

print("Importing data from: " + fn)
ds = xr.load_dataset(fn)
print("Dataset correctly loaded")
ds = bt.DAV(ds,mer_gradient_filter = True,data_return=True)
print("DAV function correctly executed")
#ds = bt.ContourTracking2D(ds,pers=4,min_area = 500000,max_length=28,data_return=True)
#print("ContourTracking2D function correclty executed")
ds.to_netcdf(fn_out)
print("Data created and exported in: " + fn_out)
  
print("--- %s seconds ---" % (time.time() - start_time))
