import numpy as np
import copy
import os
import h5py
from mapio.shake import ShakeGrid
from openquake.hazardlib.geo.geodetic import distance
from readstation import readStation
from setup import initialize


def realizations(total_real, my_reals, radius, variables, grid_arr, mu_arr, sigma_arr, list_sizes_grid, list_sizes_mu,
                 shakegrid, voi, comm, dir, method, output_dir):    
    '''
    Function realizations uses output from the main function in loop.py to compute realizations of the spatially variable random field.
    :param total_real:
        integer, total number of realizations assigned to each core
    :param my_reals:
        numpy array, which realizations each core is computing
    :param radius:
        float, radius of influence
    :param variable:
        dict, output from initialize function in setup.py
    :param grid_arr:
        numpy array of all grid array values, note that these are all combined into one large array
        these are indices that each grid point depends on
    :param mu_arr:
        numpy array of all mu arrays, note that these are all combined into one large array
        Sig12.T*Sig11inv
    :param sigma_arr:
        numpy array of R values
    :param list_sizes_grid
        numpy array, the number of elements of grid_arr belonging to each grid point
    :param list_sizes_mu:
        numpy array, the number of elements of mu)arr belonging to each grid point
    :param shakegrid:
        shakegrid object
    :param voi:
        string, intensity measure
    :param comm:
        mpi communicator
    :param dir:
        string, path to inputs folder
    :param method:
        string, method for computing distances
    Outputs are saved to a file. If multiple grid.xml files are used, the epsilon matrices will be saved to file. Otherwise realizations
    of the spatially variable ShakeMap will be saved.
    '''
    num_realizations = np.size(my_reals)
    if num_realizations == 0:
        return
    my_rank = comm.Get_rank()
    size = comm.Get_size()

    # Determine if multiple grid.xml files are to be used.
    multiple_maps = 0
    isd = True
    while isd == True:
        isd = os.path.isdir(os.path.join(dir,'%i'%(multiple_maps+1)))
        if isd == True:
            multiple_maps += 1
            
    # Set data file names. If multiple realizations are used, store epsilon
    if multiple_maps > 0:
        write_correlation = True
        filename = os.path.join(output_dir,'Epsilon_%s_%i.hdf5'%(voi,my_rank))
    else:
        write_correlation = False
        filename = os.path.join(output_dir,'SVSM_%s_%i.hdf5'%(voi, my_rank))

    shakemap = shakegrid.getLayer(voi)
    N = variables['N']
    M = variables['M']

    event_attr = shakegrid.getEventDict()
    grid_attr =  shakegrid.getGeoDict()

    # Set up dictionaries to store data
    uncertaintydata, data, data_new, sm_dict = {},{}, {}, {}
    uncertaintydata['map0'] = variables['uncertaintydata']
    stationlist = os.path.join(dir,'stationlist.xml')
    stationdata = readStation(stationlist)
    data['map0'] = variables['data']
    sm_dict['map0'] = shakemap

    # If there are multiple maps, store other maps data
    for i in range(1, multiple_maps+1):
        folder = '%i/'%i
        sm_grid = ShakeGrid.load(os.path.join(dir,folder,'grid.xml'), adjust = 'res')
        sm_dict['map%i'%i] = sm_grid.getLayer(voi)
        event_attr = sm_grid.getEventDict()
        unc_grid = ShakeGrid.load(os.path.join(dir,folder,'uncertainty.xml'), adjust = 'res')
        stationlist = os.path.join(dir,folder,'stationlist.xml')
        stationdata = readStation(stationlist)
        voi_list = []
        voi_list.append(voi)
        variables = initialize(sm_grid, unc_grid, stationdata, dir, voi_list, method)
        uncertaintydata["map{0}".format(i)] = variables['uncertaintydata']
        data["map{0}".format(i)] = variables['data']

    list_size_mu = np.reshape(list_sizes_mu, [M*N,1])
    list_size_grid = np.reshape(list_sizes_grid, [M*N,1])
    sigma_arr = np.reshape(sigma_arr, [M*N,1])

    # Set header information for file. Change if neccesary
    f = h5py.File(filename, 'w')
    f.attrs['Conventions'] = 'COARDS, CF-1.5'
    f.attrs['title'] = 'filename'
    f.attrs['history'] = 'Created with python MultiHazardGrid.save(%s)' % filename
    f.attrs['GMT_version'] = 'NA'

    xvar = np.linspace(grid_attr.xmin,grid_attr.xmax,grid_attr.nx)
    yvar = np.linspace(grid_attr.ymin,grid_attr.ymax,grid_attr.ny)
    x = f.create_dataset('x',data=xvar,compression='gzip',shape=xvar.shape,dtype=str(xvar.dtype))
    x.attrs['CLASS'] = 'DIMENSION_SCALE'
    x.attrs['NAME'] = 'x'
    x.attrs['_Netcdf4Dimid'] = 0 #no idea what this is
    x.attrs['long_name'] = 'x'
    x.attrs['actual_range'] = np.array((xvar[0],xvar[-1]))
    
    y = f.create_dataset('y',data=yvar,compression='gzip',shape=yvar.shape,dtype=str(yvar.dtype))
    y.attrs['CLASS'] = 'DIMENSION_SCALE'
    y.attrs['NAME'] = 'y'
    y.attrs['_Netcdf4Dimid'] = 1 #no idea what this is
    y.attrs['long_name'] = 'y'
    y.attrs['actual_range'] = np.array((yvar[0],yvar[-1]))

    # Compute realizations of the field, COR
    for j in range(0, num_realizations):
        X = np.zeros([M*N,1])
        for i in range(0,M*N):
            st_g = np.sum(list_size_grid[0:i])
            st_m = np.sum(list_size_mu[0:i])
            end_g = st_g + list_size_grid[i]
            end_m = st_m + list_size_mu[i]
            rand_arr = np.random.randn()
            nzeros = list_size_mu[i] - list_size_grid[i]
            x = np.append(np.zeros(nzeros), X[np.array(grid_arr[st_g:end_g], dtype = 'i')])            
            mu = np.dot(mu_arr[st_m:end_m], x)
            X[i] = mu + rand_arr * sigma_arr[i]

        COR = np.reshape(X, [M,N])
        layerkey = 'realization_%i'%j

        # Write data to file
        if write_correlation == True:
            dset = f.create_dataset(layerkey,data=COR,compression='gzip')
            dset.attrs['long_name'] = layerkey

        else:
            for i in range(0, multiple_maps+1):
                xx = 'map%i'%i
                X = np.multiply(COR, uncertaintydata[xx][voi])
                DATA_NEW = data[xx][voi]*np.exp(X)
                dset = f.create_dataset(layerkey,data=DATA_NEW,compression='gzip')
                dset.attrs['long_name'] = layerkey

        if np.mod(j+1, 25) == 0:
            print('Done with', j+1, 'of', num_realizations, 'iterations.')

    f.close()
    
