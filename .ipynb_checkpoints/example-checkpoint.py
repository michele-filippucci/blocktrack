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
ds = ds = bt.DAV(ds,mer_gradient_filter = True,data_return=True)
print("DAV function correctly executed")
ds,dic = bt.ContourTracking2D(ds,var_name='DAV',data_return=True,char_dic=True,exp_dic=False)
#a dictionary object is created which contains the features of the tracked events.
print("ContourTracking2D function correclty executed")
ds,dic_filt = bt.FilterEvents(ds,dic,data_return=True,var_name='GHA',fn_dic_out=string_dic_out,
              pers_min = 5,pers_max = 25,min_area = 2e6,max_area=10e6,max_ar=100,max_avg_dist=10000)
#the filtering is performed based on the features stored in the dictionary. The new set of events is then
#stored is a new dictionary (dic_filt)
ds.to_netcdf(fn_out)
print("Data created and exported in: " + fn_out)
  
print("--- %s seconds ---" % (time.time() - start_time))
