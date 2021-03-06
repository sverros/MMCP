from mpi4py import MPI
import sys
import numpy as np
import time
from openquake.hazardlib.correlation import JB2009CorrelationModel, BaseCorrelationModel
from openquake.hazardlib.geo.geodetic import geodetic_distance
from openquake.hazardlib.imt import from_string

def main(var, r, voi, cor_model, vs_corr):
    '''
    Uses grid information and correlation model to compute covariance matrices. This function encapsualtes the
    successive conditional simulation method. 
    :param var:
        dict, output from setup.initialize function
    :param r:
        float, radius of influence
    :param voi:
        list of strings, variable of interest, i.e., ['pga', 'pgv']
    :param cor_model:
        string, correlation model, currently 'JB2009'
    :param vs_corr:
        boolean, associated with cor_model
    :returns:
        dict, the following keys
        grid_arr- array of all grid array values, note that these are all combined into one large array
            these are indices that each grid point depends on
        mu_arr- array of all mu arrays, note that these are all combined into one large array
            Sig12.T*Sig11inv
        sigma_arr- array of R values
        list_sizes_grid- The number of elements of grid_arr belonging to each grid point 
        list_sizes_mu- The number of elements of mu)arr belonging to each grid point

    EXAMPLE: If your grid array for the first 5 entries were as follows:
           [0], [0,1], [0, 1, 2], [0, 1, 2, 3], [1, 2, 3, 4]
           your grid_arr would be [0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 1, 2, 3, 4]
           and list_sizes_grid would be [1, 2, 3, 4, 4]
    '''
    
    start_time = time.time()

    comm = MPI.COMM_WORLD
    my_rank = comm.Get_rank()
    size = comm.Get_size()

    M = var['M']
    N = var['N']
    K = var['K']

    if cor_model == 'JB2009':
        CM = JB2009CorrelationModel(vs_corr)

    # Compute which rows to compute
    my_rows = dist_rows(M, my_rank, size)

    # Initialize vectors for storing data
    grid_arr = [None]*(np.size(my_rows)*N)
    list_sizes_grid = np.zeros([np.size(my_rows)*N,1], dtype = 'int')
    list_sizes_mu = np.zeros([np.size(my_rows)*N,1], dtype = 'int')

    mu_arr = {}
    sigma_arr = {}

    for i in range(0, np.size(voi)):
        mu_arr[voi[i]] = [None]*(np.size(my_rows)*N)
        sigma_arr[voi[i]] = np.zeros([np.size(my_rows)*N,1])

    # Get spcing of horozontal and vertical points
    ld  = set_up_grid_dist(M,N,var['site_collection_SM'])

    set_up_time = time.time() - start_time

    base = {}
    other = {}

    for it in range(0, np.size(my_rows)):
        i = my_rows[it]

        # Find the number of points in radius horozontally and vertically for each row
        vhva = calc_vert_hor(i, r, ld['l'], ld['d'])
        
        # Calculate the full distance matrix for each row
        dist = calc_full_dist(i, vhva['vert'], vhva['hor'], N, var['site_collection_SM'])
        first_time_per_row = 1

        for j in range(0,N):
            num = i*N+j
            
            # Find the reduced distance matrix 
            dist_calc = reduce_distance(
                j, vhva['vert'], vhva['hor'], vhva['added_vert'], N, 
                dist['distance_matrix'], dist['grid_indices'])

            # Include stations in distance matrix and find the 
            # indices of the points within the radius
            out = inc_stations(j, i, N, K, r, var['site_collection_SM'], 
                               var['site_collection_station'], dist_calc['dist_mat'], 
                               dist_calc['inc_ind'])

            if np.size(dist_calc['inc_indices']) == 1:
                # Correlation value is not dependent on anything, result is random
                grid_arr [it*N+j]= np.zeros(0)
                for ii in range(0, np.size(voi)):
                    mu_arr   [voi[ii]][it*N+j]= np.zeros(0)
                    sigma_arr[voi[ii]][it*N+j]= 1
                list_sizes_grid[it*N+j] = 1
                list_sizes_mu[it*N+j] = 1
            else:
                # Check if reduced distance matrix is full distance matrix
                if ((vhva['vert'] == 1 and dist_calc['num_indices'] == vhva['hor']+1)or\
                        (vhva['vert'] != 1 and dist_calc['num_indices'] == 2*vhva['hor']+1))and\
                        (np.size(out['inc_sta_indices']) == 0):
                    # If this is the first full distance matrix per row, calculate base case
                    if first_time_per_row == 1:
                        for ii in range(0, np.size(voi)):
                            base[voi[ii]] = calculate_corr(out['dist_mat'], voi[ii], CM)
                        first_time_per_row = 0

                    # Use the base case matrices
                    grid_arr [it*N+j] = np.array(dist_calc['inc_ind'][0:-1])
                    for ii in range(0, np.size(voi)):
                        mu_arr   [voi[ii]][it*N+j] = np.array(base[voi[ii]]['Sig12'].T*base[voi[ii]]['Sig11inv'])
                        sigma_arr[voi[ii]][it*N+j] = base[voi[ii]]['R']
                    list_sizes_grid[it*N+j] = int(np.size(grid_arr [it*N+j]))
                    list_sizes_mu[it*N+j] = int(np.size(mu_arr[voi[0]][it*N+j]))

                else:
                    # Need to compute the correlation matrix for this individual case
                    for ii in range(0, np.size(voi)):
                        other[voi[ii]] = calculate_corr(out['dist_mat'], voi[ii], CM)

                    grid_arr [it*N+j] = np.array(dist_calc['inc_ind'][0:-1])
                    for ii in range(0, np.size(voi)):
                        mu_arr   [voi[ii]][it*N+j] = np.array(other[voi[ii]]['Sig12'].T*other[voi[ii]]['Sig11inv'])
                        sigma_arr[voi[ii]][it*N+j] = other[voi[ii]]['R']
                    list_sizes_grid[it*N+j] = int(np.size(grid_arr [it*N+j]))
                    list_sizes_mu[it*N+j] = int(np.size(mu_arr[voi[0]][it*N+j]))

            if np.mod(i*N+j,5000) == 0:
                print('Finishing step:', i*N+j)


    # Reshape grid and mu_arr into vector of values. Store number of entries in each in new vector, list_sizes
    grid_vec = np.zeros([np.sum(list_sizes_grid)])
    mu_vec = {}

    for ii in range(0, np.size(voi)):
        mu_vec[voi[ii]] = np.zeros([np.sum(list_sizes_mu)])
    
    index1 = 0
    index2 = 0
    for it in range(0, np.size(my_rows)):
        i = my_rows[it]
        for j in range(0,N):
            if np.size(grid_arr[it*N+j]) == 0:
                grid_vec[index1:index1+list_sizes_grid[it*N+j]] = -1
            else:    
                grid_vec[index1:index1+list_sizes_grid[it*N+j]] = np.reshape(grid_arr[it*N+j], [list_sizes_grid[it*N+j]])
            index1 += list_sizes_grid[it*N+j]
            if np.size(mu_arr[voi[0]][it*N+j]) == 0:
                for ii in range(0, np.size(voi)):
                    mu_vec[voi[ii]][index2:index2+list_sizes_mu[it*N+j]] = -1
            else:
                for ii in range(0, np.size(voi)):
                    mu_vec[voi[ii]][index2:index2+list_sizes_mu[it*N+j]] = np.reshape(mu_arr[voi[ii]][it*N+j], [list_sizes_mu[it*N+j]])
            index2 += list_sizes_mu[it*N+j]
    loop_time = time.time() - set_up_time - start_time

    if my_rank == 0:
        print('Loop time', loop_time)

    # Send gathered information to master, order sending must match order of recieving.
    for rank in range(1,size):                                                                                                                                 
       if my_rank == rank:
           comm.Send([list_sizes_grid, MPI.INT], dest = 0, tag = my_rank*11)
           comm.Send([list_sizes_mu, MPI.INT], dest = 0, tag = my_rank*111)
           for ii in range(0, np.size(voi)):
               comm.Send([sigma_arr[voi[ii]], MPI.DOUBLE], dest = 0, tag = my_rank*100*(ii+100))
           comm.Send([grid_vec, MPI.INT], dest = 0, tag = my_rank*12)
           for ii in range(0, np.size(voi)):
               comm.Send([mu_vec[voi[ii]], MPI.DOUBLE], dest = 0, tag = my_rank*13*(ii+100))
   # Master gathers all information
    if my_rank == 0:
        full_list_sizes_grid = np.zeros([M,N])
        full_list_sizes_mu = np.zeros([M,N])
        temp_list_sizes_grid = np.zeros(size)
        temp_list_sizes_mu = np.zeros(size)
        full_sigma_arr = {}
        for ii in range(0, np.size(voi)):
            full_sigma_arr[voi[ii]] = np.zeros([M,N])
        ii = 0
        # Fill in masters known entries
        for i in np.arange(0, M, size):
            for j in range(0, N):
                full_list_sizes_grid[i,j] = list_sizes_grid[ii]
                full_list_sizes_mu[i,j] = list_sizes_mu[ii]
                ii += 1
        # Recieve all other cores size information
        temp_sigma = {}
        for rank in range(1, size):
            num_rows = dist_rows(M, rank, size)
            temp_list_grid = np.zeros([np.size(num_rows)*N,1], dtype = 'int')
            comm.Recv([temp_list_grid, MPI.DOUBLE], source = rank, tag = rank*11)
            temp_list_mu = np.zeros([np.size(num_rows)*N,1], dtype = 'int')
            comm.Recv([temp_list_mu, MPI.DOUBLE], source = rank, tag = rank*111)
            for kk in range(0, np.size(voi)):
                temp_sigma[voi[kk]] = np.zeros([np.size(num_rows)*N,1])
                comm.Recv([temp_sigma[voi[kk]], MPI.DOUBLE], source = rank, tag = rank*100*(kk+100))

            temp_list_sizes_grid[rank] = np.sum(temp_list_grid)
            temp_list_sizes_mu[rank] = np.sum(temp_list_mu)
            
            ii = 0
            for i in np.arange(rank, M, size):
                for j in range(0, N):
                    full_list_sizes_grid[i,j] = temp_list_grid[ii]
                    full_list_sizes_mu[i,j] = temp_list_mu[ii]
                    for kk in range(0, np.size(voi)):
                        full_sigma_arr[voi[kk]][i,j] = temp_sigma[voi[kk]][ii]
                    ii += 1


        full_grid_arr = np.zeros([np.sum(full_list_sizes_grid)])
        full_mu_arr = {}
        for kk in range(0, np.size(voi)):
            full_mu_arr[voi[kk]] = np.zeros([np.sum(full_list_sizes_mu)])
        
        # Fill in grid entries
        ii_g = 0
        ii_m = 0
        for i in np.arange(0, M, size):
            for j in range(0, N):
                if i == 0:
                    st_g = np.sum(full_list_sizes_grid[i,0:j])
                    st_m = np.sum(full_list_sizes_mu[i,0:j])
                else:
                    st_g = np.sum(full_list_sizes_grid[0:i,:]) + np.sum(full_list_sizes_grid[i,0:j])
                    st_m = np.sum(full_list_sizes_mu[0:i,:]) + np.sum(full_list_sizes_mu[i,0:j])
                end_g = st_g + full_list_sizes_grid[i, j]
                end_m = st_m + full_list_sizes_mu[i, j]
                full_grid_arr[st_g:end_g] = grid_vec[ii_g:ii_g+full_list_sizes_grid[i,j]]
                for kk in range(0, np.size(voi)):
                    full_mu_arr[voi[kk]][st_m:end_m] = mu_vec[voi[kk]][ii_m:ii_m+full_list_sizes_mu[i,j]]
                ii_g += full_list_sizes_grid[i,j]
                ii_m += full_list_sizes_mu[i,j]
        
        temp_mu_arr = {}
        for rank in range(1, size):
            num_rows = dist_rows(M, rank, size)
            temp_grid_arr = np.zeros(temp_list_sizes_grid[rank])
            comm.Recv([temp_grid_arr, MPI.INT], source = rank, tag = rank*12)
            for kk in range(0, np.size(voi)):
                temp_mu_arr[voi[kk]] = np.zeros(temp_list_sizes_mu[rank])
                comm.Recv([temp_mu_arr[voi[kk]], MPI.DOUBLE], source = rank, tag = rank*13*(kk+100))
        
            ii_g = 0
            ii_m = 0
            for i in np.arange(rank, M, size):
                for j in range(0, N):
                    st_g = np.sum(full_list_sizes_grid[0:i,:]) + np.sum(full_list_sizes_grid[i,0:j])
                    end_g = st_g + full_list_sizes_grid[i, j]
                    st_m = np.sum(full_list_sizes_mu[0:i,:]) + np.sum(full_list_sizes_mu[i,0:j])
                    end_m = st_m + full_list_sizes_mu[i, j]
                    full_grid_arr[st_g:end_g] = temp_grid_arr[ii_g:ii_g+full_list_sizes_grid[i,j]]
                    for kk in range(0, np.size(voi)):
                        full_mu_arr[voi[kk]][st_m:end_m] = temp_mu_arr[voi[kk]][ii_m:ii_m+full_list_sizes_mu[i,j]]
                    ii_g += full_list_sizes_grid[i,j]
                    ii_m += full_list_sizes_mu[i,j]
                    
    sys.stdout.flush()
    comm.Barrier()
    send_time = time.time() - loop_time - set_up_time - start_time
    if my_rank == 0:
        print('Send time', send_time)
        return {'grid_arr':full_grid_arr, 'mu_arr':full_mu_arr, 'sigma_arr':full_sigma_arr, 
                'list_sizes_grid':full_list_sizes_grid, 'list_sizes_mu':full_list_sizes_mu}

def dist_rows(M, my_rank, size):
    """
    Distributes the rows as evenly as possible for the number of cores
    :param M:
        integer, number of rows
    :param my_rank:
        integer, id number of core
    :param size:
        integer, number of processes
    :returns:
        numpy array, array of rows for core to compute
    """
    my_rows = np.arange(my_rank,M, size)

    return my_rows

def calculate_corr(dist_mat, voi, CM):
    """
    Calculates correlation model for distance matrix and voi
    :param dist_mat:
        numpy matrix, reduced distance matrix
    :param voi:
        string, variable of interest
    :param JB_cor_model:
        correlation model object, see correlation.py in oq-hazardlib
    :returns:
        dict, with the following keys
        Sig12, Sig11inv- partitions of correlation matrix
        R - Sqrt of sigma
    """
    if voi.startswith('psa'):
        voi_str = voi[1:3]
        voi_num = voi[3:]
        voi = voi_str+'('+voi_num[0]+'.'+voi_num[1]+')'
    correlation_model = CM._get_correlation_model(dist_mat, from_string(voi.upper()))
    
    Sig11 = np.mat(correlation_model[0:-1, 0:-1])
    Sig12 = np.mat(correlation_model[0:-1, -1]).T
    Sig22 = np.mat(correlation_model[-1,-1])

    Sig11inv = np.mat(np.linalg.pinv(Sig11))
    sigma = Sig22 - (Sig12.T*Sig11inv*Sig12)

    R = np.sqrt(sigma)
    
    return {'Sig12':Sig12, 'Sig11inv':Sig11inv, 'R':R}

def inc_stations(j, i, N, K, r, site_collection_SM, 
                 site_collection_station, dist_mat, inc_ind):
    """
    If there are stations included within the radius for a point, 
    this function will add those stations to the 
    distance matrix and determine the array of points included in the radius, x
    :param i,j:
        integer, row and column 
    :param N:
        integer, number of points in row
    :param K:
        integer, number of stations
    :param r:
        float, radius 
    :param site_collection_SM:
        site collection object, for shakemap grid points
    :param site_collection_station:
        site collection object, for station data
    :param dist_mat:
        numpy array, reduced distance matrix
    :param X:-
        numpy array, array of previously calculated correlation values
    :param inc_ind:-
        numpy array, indices of included points
    :param inc_indices:
         integer, total number of points in the top most row of distance matrix
    :returns:
        dict, with following keys
        dist_mat- reduced distance matrix, modified to include stations
        inc_sta_indices- array of indices of stations included in the radius
    """
    
    # Compute the distances for all stations to the grid point we're looking at
    dist_sta_sit = np.array(
        geodetic_distance(site_collection_SM.lons[j+i*N], site_collection_SM.lats[j+i*N],
                          site_collection_station.lons[0:K], site_collection_station.lats[0:K]))
    
    # Find which of those stations are in the radius we are considering
    inc_sta_indices = np.where(dist_sta_sit < r)[0]
    if np.size(inc_sta_indices) != 0:
        sta_to_sta_dist = np.zeros([np.size(inc_sta_indices), np.size(inc_sta_indices)])
        sta_to_grd_dist = np.zeros([np.size(inc_sta_indices), np.size(inc_ind)])
        
        iinds = np.array(inc_ind).T[0]
        # Calculate distance between each included station and all included grid points,
        # then calculate the distance from each included station to every other included station
        for eta in range(0, np.size(inc_sta_indices)):
            sta_to_grd_dist[eta, :] = geodetic_distance(
                site_collection_station.lons[inc_sta_indices[eta]], 
                site_collection_station.lats[inc_sta_indices[eta]],
                site_collection_SM.lons[iinds], site_collection_SM.lats[iinds])
            sta_to_sta_dist[eta, eta+1:] = geodetic_distance(
                site_collection_station.lons[inc_sta_indices[eta]], 
                site_collection_station.lats[inc_sta_indices[eta]],
                site_collection_station.lons[inc_sta_indices[eta+1:]], 
                site_collection_station.lats[inc_sta_indices[eta+1:]])
            
        sta_to_sta_dist = sta_to_sta_dist + sta_to_sta_dist.T
        station_distance_matrix = np.concatenate((sta_to_sta_dist, sta_to_grd_dist), axis=1)

        # Concatenate the station distance matrix with the modified distance matrix, dist_mat
        dist_mat = np.concatenate((station_distance_matrix[:, np.size(inc_sta_indices):], dist_mat), axis=0)
        dist_mat = np.concatenate((station_distance_matrix.T, dist_mat), axis=1)

    return {'dist_mat':dist_mat, 'inc_sta_indices':inc_sta_indices}

def reduce_distance(j, vert, hor, added_vert, N, distance_matrix, grid_indices):
    """
    Find which columns/rows in the distance matrix to keep
    :param j:
        integer, column
    :param vert:
        integer, number of rows included in the radius
    :param hor:
        integer, number of columns included in radius
    :param added_vert:
        integer, number of rows in between first row and first included row 
    :param N:
        integer, number of points in row
    :param distance_matrix:
        numpy array, full distance matrix
    :param grid_indices:
        numpy array, indices included in full distance matrix
    :returns:
        dict, with following keys
        dist_mat- reduced distance matrix
        inc_indices- number of points in top most row of dist_mat
        inc_ind - indices of points included in dist_mat
        num_indices- number of points in top most row of distance_matrix
    """
    inc_ind = [None]*(vert*(2*hor+1))
    n_grid_indices = 0
    num_indices = 0
    
    if j < hor:
        # On left end of grid
        inc_indices = range(0,np.size(grid_indices))
        if vert == 1:
            num_indices = (j+1)
            inc_indices = np.where(np.mod(inc_indices,hor+1) >=hor+1 - num_indices)
        else:
            num_indices = (j+hor+1)
            inc_indices = np.where(np.mod(inc_indices,2*hor+1) >=2*hor+1 - num_indices)
            
        for k in range(0,vert):
            if k == vert-1:
                for eta in range(0,j+1):
                    inc_ind[n_grid_indices] = [eta+N*k+ N*added_vert]
                    n_grid_indices += 1
            else:
                for eta in range(0,j+hor+1):
                    inc_ind[n_grid_indices] = [eta+N*k+N*added_vert]
                    n_grid_indices += 1
        del inc_ind[n_grid_indices:]
            
    elif j > N-hor-1:
        # On right end of grid
        inc_indices = range(0,np.size(grid_indices))
        if vert == 1:
            num_indices = (1+hor)
            inc_indices = np.where(np.mod(inc_indices,hor+1) >=hor+1 - num_indices)
        else:
            num_indices = (N - j+hor)
            inc_indices = np.where(np.mod(inc_indices,2*hor+1) <=num_indices-1)
        for k in range(0,vert):
            if k == vert-1:
                for eta in range(j-hor,j+1):
                    inc_ind[n_grid_indices] = [eta+N*k+N*added_vert]
                    n_grid_indices += 1
            else:
                for eta in range(j-hor,N):
                    inc_ind[n_grid_indices] = [eta+N*k+N*added_vert]
                    n_grid_indices += 1
        del inc_ind[n_grid_indices:]
            
    else:
        # In the middle of the grid, all points included
        inc_indices = range(0,np.size(grid_indices))
        if vert == 1:
            num_indices = (1+hor)
            inc_indices = np.where(np.mod(inc_indices,hor+1) >=hor+1 - num_indices)
        else:
            num_indices = (2*hor+1)
            inc_indices = np.where(np.mod(inc_indices,2*hor+1) >=2*hor+1 - num_indices)
        for k in range(0,vert):
            if k == vert-1:
                for eta in range(j-hor,j+1):
                    inc_ind[n_grid_indices] = [eta+N*k+N*added_vert]
                    n_grid_indices += 1
            else:
                for eta in range(j-hor,j+hor+1):
                    inc_ind[n_grid_indices] = [eta+N*k+N*added_vert]
                    n_grid_indices += 1
        del inc_ind[n_grid_indices:]

    inc_indices = np.array(inc_indices).flatten()
        
    # dist_mat: the distance matrix modified for the current point
    dist_mat = distance_matrix[:,inc_indices]
    dist_mat = dist_mat[inc_indices, :]

    return {'dist_mat':dist_mat, 'inc_indices':inc_indices, 'inc_ind':inc_ind, 'num_indices':num_indices}

def calc_full_dist(row, vert, hor, N, site_collection_SM):
    """
    Calculates full distance matrix. Called once per row.
    INPUTS:
    :param vert:
        integer, number of included rows
    :param hor:
        integer, number of columns within radius 
    :param N:
        integer, number of points in row
    :param site_collection_SM:
        site collection object, for ShakeMap data
    :returns:
        dict, with following keys
        grid_indices- indices of points included in distance matrix
        distance_matrix- full distance matrix
    """

    # gathers indices for full distance matrix for each row
    grid_indices = [None]*(vert*(2*hor+1))
    n_grid_indices = 0
    for k in range(row-vert+1, row+1):
        if k == row:
            for j in range(0,hor+1):
                grid_indices[n_grid_indices] = j + N*k
                n_grid_indices += 1
        else:
            for j in range(0,2*hor+1):
                grid_indices[n_grid_indices] = j + N*k
                n_grid_indices += 1 
    del grid_indices[n_grid_indices:]

    distance_matrix = np.zeros([np.size(grid_indices), np.size(grid_indices)])

    # Create full distance matrix for row
    for k in range(0, np.size(grid_indices)):
        distance_matrix[k, k:] = geodetic_distance(
            site_collection_SM.lons[grid_indices[k ]], site_collection_SM.lats[grid_indices[k]],
            site_collection_SM.lons[grid_indices[k:]], site_collection_SM.lats[grid_indices[k:]]).flatten()
    
    distance_matrix = distance_matrix + distance_matrix.T
    
    return {'grid_indices':grid_indices, 'distance_matrix':distance_matrix}

def calc_vert_hor(i, r, l, d):
    """
    Calculates the number of vertical of horozontal points in the full distance matrix
    :param i:
        integer, current row
    :param r:
        float, radius
    :param l:
        numpy array, vector if distances between points vertically
    :param d:
        numpy array, vector if distances between points horozontally
    :returns:
        dict, with following keys
        vert- number of rows in radius above current point
        hor- number of columns in radius to the left of current point
        added_vert- number of rows in between first row and the first row included in radius
    """

    hor = int(np.floor(r/d[i]))

    if i*l[1] <= r:
        vert = i+1
        added_vert = 0
    else:
        vert = int(np.floor(r/l[1])) + 1
        added_vert = i - vert + 1

    return {'vert':vert, 'hor':hor, 'added_vert':added_vert}

def set_up_grid_dist(M,N, site_collection_SM):
    """
    Calculates the vertical and horozontal spacing between points for each row
    INPUTS:
    :param M,N:
        integer, number of points in grid vertically and horozontally
    :param site_collection_SM:
        site collection object, for ShakeMap data
    :returns:
        dict, with following keys
        l- vector of distances between points vertically
        d- vector of distance between points horozontally
    """
    l = np.zeros([M-1])
    d = np.zeros([M])
    
    # Calculate vertical and horozonal spacing between points for each row
    l[:] = geodetic_distance(
        site_collection_SM.lons[range(N,N*M,N)],   site_collection_SM.lats[range(N,N*M,N)],
        site_collection_SM.lons[range(0,N*M-N,N)], site_collection_SM.lats[range(0,N*M-N,N)])
    d[:] = geodetic_distance(
        site_collection_SM.lons[range(0, M*N, N)], site_collection_SM.lats[range(0, M*N, N)],
        site_collection_SM.lons[range(1, M*N, N)], site_collection_SM.lats[range(1, M*N, N)])
    return {'l':l, 'd':d}
