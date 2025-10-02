# blocktrack
A set of functions for the detection of atmospheric blocking events in a gridded dataset and the Lagrangian tracking of atmospheric blocking trajectories.

If you intend to use this set of functions for your analysis do not hesitate and email me at: michele.filippucci@unitn.it
The documentation for the algorithm may not be updated or exaustive and you may need guidance.

You do not need to install anything to use the library. You just have to include the directory of the cloned github repository in your code as done in the example provided.

libraries needed for the algorithm ro work:
xarray
numpy
scipy
tqdm
math
sys

We tested different versions of these libraries and the algorithm should work with the latest releases. If you encounter problems please contact me.

The dataset processed by the functions is a NETCDF gridded dataset of daily geopotential height at 500 hPa. The dataset must refer to the Northern Hemisphere and the cordinates should be in the range [-180,180] degrees of longitude and [0,90] degrees of latitude. The geopotential height must be defined in meter. The grid size can be set by the user while intializing the functions. Not every resolution is compatible. I suggest using 2.5,1.5 or 1 degrees of spatial resolution. The temporal resolution must be 1 day. 

A series of indices based on the 500hPa geopotential height are implemented: 
- DAV: the index introduced by Davini et al 2012.
- GHA: a simple implementation of an anomaly based index similar to the one adopted in Woollings et al. 2018.
- MIX: similar to Peings et al. 2021. It mixes the previous two indices getting the good and the bad of both.
We also implemented two additional indices that use as input variable the local wave activity computed through the falwa package by Huang and Nakaumura (https://github.com/csyhuang/hn2016_falwa). In this case the input variable is the density weighted vertical integral of local wave activity (again, a single level dataset). 
- LWAA: the index works exactly like the GHA index but uses the local wave activity (vertically integrated) as input variable.
- LWAA_2: this index follows the definition of blocking index adopted in Barpanda and Nakamura 2025.

For further details see the comments in blocktoolbox.py

The output has the same dimensions and characteristics of the input dataset, but a series of additional variables are included:
- A boolean matrix with the same dimensions of the geopotential height field named as the adopted blocking index (DAV, GHA, LWAA ...) that is 1 for blocked grid points and 0 for unblocked grid points.
- A similar matrix named as the adopted blocking index plus the desinence [...]_tracked that is similar to the previous but instead of 1 it contains a integer label that is different for each blocking event.
- A matrix named named as the adopted blocking index plus the desinence [...]_freq that is the frequency of blocking in the Northern Hemisphere.

Moreover, the algorithm outputs a list of dictionaries as a numpy object.
- the i-th element of the list is a dictionary containing the information of the i-th blocking event.

You can save a import the datasets and dictionaries using numpy save and load and xarrays routines.

By providing both the dictionaries and the dataset to the filter_events function it is possible to apply filters to the blocking dataset. The output of the filtering function is a new dataset where the events not meeting the filter criteria are deleted and the events are re-labeled from 1 to n where n is the new total number of blocking events. A new dictionary is also created. All the variables of the input datasets are updated to take into account the filtering (e.g. DAV_tracked, DAV, DAV_freq).

The algorithm is rather fast, taking around 15 minutes to track the blocking events for the whole ERA5 reanalysis (80 years) on a regular laptop (Macbook with A2 processor). The filtering of the events is also rather fast, taking just a few seconds. If you have any suggestion for improvement you are more than welcome to contact me.

The reference article for citing the algorithm is:
Filippucci, M., Bordoni, S., & Davini, P. (2024). Impact of stochastic physics on the representation of atmospheric blocking in EC-Earth3. Weather and Climate Dynamics, 5(4), 1207-1222.
The zenodo reference for the algorithm is:
https://doi.org/10.5281/zenodo.13837897

Please, cite me if you decide to use my python algorithm

