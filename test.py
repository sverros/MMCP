from mpi4py import MPI
import numpy as np
import time
import sys
import argparse
from mapio.shake import ShakeGrid
from mapio.multiple import MultiGrid
from readstation import readStation
from setup import initialize
from loop import main
from realizations import realizations

"""
Parallel code for computing the spatial correlation for a ShakeMap,
adding to a ShakeMap grid, and computing multiple realizations
VARIABLES:
    voi - variable of interest, i.e. PGA
    num_realization- integer for desired number of realizations
    method- string for distance measure
    input data- grid.xml, uncertainty.xml, and stationlist.xml
        stored in Inputs directory. fault.txt and event.xml should also
        be stored in this directory.
mpi4py is used for parallelization
File may be run using:
mpiexec -n # python test.py path imt distance_measure N
where # is the desired number of processors, path is the path to the input folder, 
imt is the intensity measure, distance_measure is the appropriate ShakeMap distance measure, 
and N is the number of realizations
"""

start_time = time.time()

# Start MPI
comm = MPI.COMM_WORLD
size = comm.Get_size()
my_rank = comm.Get_rank()

# Get variables
parser = argparse.ArgumentParser(description='Run Method of Successive Conditional Simulations to compute realizations of a spatially correlated random field')
parser.add_argument('direc', metavar='path', type=str, 
                   help='path to the inputs directory')
parser.add_argument('voi', metavar='imt', type=str, 
                   help='intensity measure, i.e., pga')
parser.add_argument('method', metavar='distance_measure', type=str, 
                   help='distance measure, i.e., rjb, rrup')
parser.add_argument('num_realizations', metavar='N', type = int, 
                   help='number of realizations')
args = parser.parse_args()
direc = args.direc
num_realizations = args.num_realizations
voi = args.voi
method = args.method

# Constant variables
radius = [15]
# For now the correlation model is constant, as openquake only include this correlation model. 
corr_model = 'JB2009'
vscorr = True

# Get shakemap for desired variable, PGA, uncertainty grid and stationdata
shakegrid = ShakeGrid.load(direc + 'grid.xml', adjust='res')
shakemap = shakegrid.getLayer(voi)
unc_grid = ShakeGrid.load(direc + 'uncertainty.xml', adjust = 'res')
unc_INTRA = unc_grid.getLayer('gmpe_intra_std%s' % voi)
unc_INTER = unc_grid.getLayer('gmpe_inter_std%s' % voi)
stationlist = direc+'stationlist.xml'
stationdata = readStation(stationlist)

# Initialize the grid
# In this step we use the ShakeMap outputs to determine the grid points, grid spacing, site collections, 
# station data, and other initial values
variables = initialize(shakegrid, unc_INTRA, unc_INTER, stationdata, direc, voi, method)
if my_rank == 0:
    print(variables['K'], 'stations', variables['M']*variables['N'], 'data points')
initialization_time = time.time() - start_time
if my_rank == 0:
    print('Initialization time', initialization_time)

# Compute the grid, mu, and sigma arrays
# In this step, we use the correlation model to compute the covariance matrices for each point on the 
# ShakeMap grid. This computation is done in parallel
out = main(variables, radius, voi, corr_model, vscorr)
main_time = time.time()-  start_time - initialization_time
if my_rank == 0:
    print('Main time', main_time)

# Compute realizations of the random field
# After computing the covariance matrices for each point, we can compute realizations of the random
# fields. If multiple cores are used, each core needs the data for every point on the grid.
if num_realizations == 1:
    # Master will compute this single realization
    if my_rank == 0:
        data = realizations(1, 1, radius, variables, 
                            out['grid_arr'], out['mu_arr'], 
                            out['sigma_arr'], out['list_sizes_grid'], out['list_sizes_mu'],
                            shakegrid, voi, comm, direc, method)
else:
    # Master broadcasts the arrays to the other cores
    if my_rank == 0:
        grid_arr = out['grid_arr']
        mu_arr = out['mu_arr']
        sigma_arr = out['sigma_arr']
        list_sizes_grid = out['list_sizes_grid']
        list_sizes_mu = out['list_sizes_mu']
    else:
        grid_arr = None
        mu_arr = None
        sigma_arr = None
        list_sizes_grid = None
        list_sizes_mu = None        

    grid_arr = comm.bcast(grid_arr, root = 0)
    mu_arr = comm.bcast(mu_arr, root = 0)
    sigma_arr = comm.bcast(sigma_arr, root = 0)
    list_sizes_grid = comm.bcast(list_sizes_grid, root = 0)
    list_sizes_mu = comm.bcast(list_sizes_mu, root = 0)

    my_reals = np.arange(my_rank, num_realizations, size) 

    # Each core does a set of realizations
    data = realizations(num_realizations, my_reals, radius, variables,
                        grid_arr, mu_arr, sigma_arr, list_sizes_grid, list_sizes_mu,
                        shakegrid, voi, comm, direc, method)

realization_time = time.time() - start_time - initialization_time - main_time
if my_rank == 0:
    print('Realization time', realization_time)

# Output is stored in data_files folder
