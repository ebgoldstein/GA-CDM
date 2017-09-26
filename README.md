# Genetic Algorithm — CDM tuning routine

This bunch of code is used to tune the Coastal Dune Model — specifically CDM-L, avaialble here:  https://github.com/ebgoldstein/Coastal-Dune-Model-Lateral, which is a variant of CDM v2.0 (https://github.com/csdms-contrib/Coastal-Dune-Model). Nick Cohn wrote the first version of the 'cdm_params.m' file in the CDM folder. The calibration routine uses the genetic algorithm routine supplied in Matlab through the Global Optimization toolbox (https://www.mathworks.com/products/global-optimization.html). 

All Genetic Algorithms require a fitness function. The fitness funciton used here is a displacement based error metric presented by Bosboom and Reniers (https://www.adv-geosci.net/39/37/2014/adgeo-39-37-2014.html), which utilizes the algorithm of Kroon and Slump (2009) available on the Mathworks page (https://www.mathworks.com/matlabcentral/fileexchange/21451-multimodality-non-rigid-demon-algorithm-image-registration)

To make this algorithm work, there are 2 additional components needed:
  1. compile CDM in the CDM folder
  2. Download the Kroon and Slump algorithm and place the files in a folder entitled 'demon_reg' (within this directory).
  
The Scirpt to operate the code is 'CDMGA.m', and the actual function to be optimized (which runs CDM and computes the error) is 'CDMfunction.m'
  
The 3 sets of topo and vegetation data (April, September and November, all 2016) are rescaled to 1m x 1m domains from 3 point cloud files derived from Kite-based structure from motion photogrammetry at Fort Fisher State Recreation Area, NC,USA. 

And the data files are here:
  1. April 2016: https://doi.org/10.6084/m9.figshare.3370648.v2
  2. September 2016: https://doi.org/10.6084/m9.figshare.3856248.v2
  3. November 2016: https://doi.org/10.6084/m9.figshare.4746613.v1
 
The technique capture technique is described here: https://dx.doi.org/10.7287/peerj.preprints.1444v1 
 
If you donwload the files and process them (using the included script SfMdense2topogrid.m) then you need to rename the dense cloud files as  FF1dense.txt, FF2dense.txt, FF3dense.txt, respectively.
