"""
Parallel code for computing the spatial correlation for a ShakeMap,
adding to a ShakeMap grid, and computing multiple realizations
VARIABLES:
    voi - variable of interest, i.e. PGA
    r - radius of influence
    num_realization- integer for desired number of realizations
    corr_model- JB2009 or GA2010
    vs_corr- Vs30 correlated bool, see JB2009
    plot_on- boolean plotting value
    compute_loss- boolean for loss computation
    multiple_maps- integer, for shakemaps that are situated exactly
    ontop of one another
    method- string for distance measure
    input data- grid.xml, uncertainty.xml, and stationlist.xml
        stored in Inputs directory
        
mpi4py is used for parallelization
File may be run using:
mpiexec -n # python test.py
where # is the desired number of processors
"""
from mpi4py import MPI
import numpy as np
import time
import sys
from mapio.shake import ShakeGrid
from mapio.multiple import MultiGrid
from neicio.readstation import readStation
#from neicio.shake import ShakeGrid
from setup import initialize
from loop import main
from realizations import realizations

# Start MPI
start_time = time.time()
comm = MPI.COMM_WORLD
size = comm.Get_size()
my_rank = comm.Get_rank()

# Variables
voi = 'pga'
radius = [45]
num_realizations = 10
corr_model = 'JB2009'
vscorr = True
plot_on = True
compute_loss = True
multiple_maps = 4
method = 'rjb'
direc = '/Users/sverros/Documents/Modules/MM_CP/input/'

# Get shakemap for desired variable, PGA, uncertainty grid and stationdata
shakegrid = ShakeGrid.load(direc + 'grid.xml', adjust='res')
shakemap = shakegrid.getLayer(voi)
unc_grid = ShakeGrid.load(direc + 'uncertainty.xml', adjust = 'res')
unc_INTRA = unc_grid.getLayer('gmpe_intra_std%s' % voi)
unc_INTER = unc_grid.getLayer('gmpe_inter_std%s' % voi)
### UPDATE!
stationlist = direc+'stationlist.xml'
stationdata = readStation(stationlist)

# Initialize the grid
if my_rank == 0:
    print 'Calling initialize'

variables = initialize(shakegrid, unc_INTRA, unc_INTER, stationdata, direc, voi, method)
if my_rank == 0:
    print variables['K'], 'stations', variables['M']*variables['N'], 'data points'

initialization_time = time.time() - start_time
if my_rank == 0:
    print 'Initialization time', initialization_time

# Compute the grid, mu, and sigma arrays
if my_rank == 0:
    print 'Calling main'
out = main(variables, radius, voi, corr_model, vscorr)

main_time = time.time()-  start_time - initialization_time
if my_rank == 0:
    print 'main time', main_time

if num_realizations == 1:
    # Master will compute this single realization
    if my_rank == 0:
        print 'Computing realizations'
        data = realizations(1, 1, radius, variables, 
                            out['grid_arr'], out['mu_arr'], 
                            out['sigma_arr'], out['list_sizes_grid'], out['list_sizes_mu'],
                            compute_loss, shakemap, voi, comm, multiple_maps, plot_on, direc, method)
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
                        compute_loss, shakemap, voi, comm, multiple_maps, plot_on, direc, method)

realization_time = time.time() - start_time - initialization_time - main_time
if my_rank == 0:
    print 'realization time', realization_time
