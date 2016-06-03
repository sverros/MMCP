from mpi4py import MPI
import numpy as np
import time
import sys
import os
import argparse
from mapio.shake import ShakeGrid
from mapio.multiple import MultiGrid
from readstation import readStation
from init_grid import initialize
from loop import main
from realizations import realizations

def run_method(direc, voi, num_realizations, radius, corr_model, vscorr, output_dir):
    """
    Parallel code for computing the spatial correlation for a ShakeMap,
    adding to a ShakeMap grid, and computing multiple realizations.
    File may be run using:
    mpiexec -n # python test.py path imt distance_measure N
    where # is the desired number of processors. Required command line parameters are listed below:
    :param path:
        string, path to directory containing grid, stationlist, uncertainty, event xmls and fault.txt
    :param imt:
        string, intensity measures to use, i.e., 'pga pgv psa03'
    :param distance_measure:
        string, appropriate distance measure for the ShakeMap, i.e., 'rjb'
    :param N:
        integer, number of realizations to compute
    :param output_dir:
        path to directory where output is stored
    """
    start_time = time.time()
    
    # Start MPI
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    my_rank = comm.Get_rank()
    
    # Get shakemap, uncertainty grid, and stationdata
    shakegrid = ShakeGrid.load(os.path.join(direc,'grid.xml'), adjust='res')
    unc_grid = ShakeGrid.load(os.path.join(direc,'uncertainty.xml'), adjust = 'res')
    stationlist = os.path.join(direc,'stationlist.xml')
    stationdata = readStation(stationlist)
    
    # Initialize the grid
    # In this step we use the ShakeMap outputs to determine the grid points, grid spacing, site collections, 
    # station data, and other initial values
    variables = initialize(shakegrid, unc_grid, stationdata, direc, voi)
    if my_rank == 0:
        print(variables['K'], 'stations', variables['M']*variables['N'], 'data points')
    initialization_time = time.time() - start_time
    if my_rank == 0:
        print('Initialization time', initialization_time)
    sys.stdout.flush()

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
    for ii in range(0, np.size(voi)):
        if num_realizations == 1:
            # Master will compute this single realization
            if my_rank == 0:
                data = realizations(1, 1, radius, variables, 
                                    out['grid_arr'], out['mu_arr'][voi[ii]], 
                                    out['sigma_arr'][voi[ii]], out['list_sizes_grid'], out['list_sizes_mu'],
                                    shakegrid, voi[ii], comm, direc, method, output_dir)
        else:
            # Master broadcasts the arrays to the other cores
            if my_rank == 0:
                grid_arr = out['grid_arr']
                mu_arr = out['mu_arr'][voi[ii]]
                sigma_arr = out['sigma_arr'][voi[ii]]
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
                                shakegrid, voi[ii], comm, direc, output_dir)

    realization_time = time.time() - start_time - initialization_time - main_time
    if my_rank == 0:
        print('Realization time', realization_time)

    # Output is stored in data_files folder
