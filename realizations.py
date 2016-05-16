import numpy as np
import copy
from model_losses import model_loss, initialize_loss, model_econ, MMI_conv2
from neicio.shake import ShakeGrid
from openquake.hazardlib.geo.geodetic import distance
from neicio.readstation import readStation
from setup import initialize
from plotting import plot, plot_fats, plot_econ, plot_intensity
import time

def realizations(total_real, my_reals, radius, variables, grid_arr, mu_arr, sigma_arr, list_sizes_grid, list_sizes_mu,
                 compute_loss, shakemap, voi, comm, multiple_maps, plot_on, dir, method):

    num_realizations = np.size(my_reals)
    if num_realizations == 0:
        return

    my_rank = comm.Get_rank()
    size = comm.Get_size()

    if multiple_maps > 0:
        plot_fault = True
    else:
        plot_fault = False

    N = variables['N']
    M = variables['M']

    attr = shakemap.getAttributes()
    event_attr = attr['event']
    grid_attr =  attr['grid_specification']

    # Set up dictionaries to store data
    uncertaintydata, data, data_new, sm_dict = {},{}, {}, {}
    uncertaintydata['map0'] = variables['uncertaintydata']
    stationlist = dir+'stationlist.xml'
    stationdata = readStation(stationlist)
    data['map0'] = variables['data']
    sm_dict['map0'] = shakemap

    if compute_loss == True:
        init_dict, ccodes, fats, bestest = {},{},{},{}
        init_dict['map0'] = init_loss_outputs(sm_dict['map0'], voi, my_rank, variables['dist_to_fault'], event_attr['magnitude'])
        ccodes ['map0'] = [None]*num_realizations
        fats   ['map0'] = [None]*num_realizations
        bestest['map0'] = [None]*num_realizations
    
    for i in range(1, multiple_maps+1):
        folder = '%i/'%i
        sm_dict['map%i'%i] = ShakeGrid(dir+folder+'grid.xml', variable = '%s' % voi)
        attr = sm_dict['map%i'%i].getAttributes()
        event_attr = attr['event']
        unc_INTRA = ShakeGrid(dir+folder+'uncertainty.xml', variable = 'GMPE_INTRA_STD%s' % voi)
        unc_INTER = ShakeGrid(dir+folder+'uncertainty.xml', variable = 'GMPE_INTER_STD%s' % voi)
        stationlist = dir+folder+'stationlist.xml'
        stationdata = readStation(stationlist)
        variables = initialize(sm_dict['map%i'%i], unc_INTRA, unc_INTER, stationdata, dir, method)
        uncertaintydata["map{0}".format(i)] = variables['uncertaintydata']
        data["map{0}".format(i)] = variables['data']

        if compute_loss == True:
            init_dict['map%i'%i] = init_loss_outputs(sm_dict['map%i'%i], voi, my_rank, variables['dist_to_fault'], event_attr['magnitude'])
            ccodes ['map%i'%i] = [None]*num_realizations
            fats   ['map%i'%i] = [None]*num_realizations
            bestest['map%i'%i] = [None]*num_realizations


    g = open('../random_files/rand%i.txt' % (my_rank+1), 'r')
    
    list_size_mu = np.reshape(list_sizes_mu, [M*N,1])
    list_size_grid = np.reshape(list_sizes_grid, [M*N,1])
    sigma_arr = np.reshape(sigma_arr, [M*N,1])

    for j in range(0, num_realizations):
        X = np.zeros([M*N,1])
        for i in range(0,M*N):
            st_g = np.sum(list_size_grid[0:i])
            st_m = np.sum(list_size_mu[0:i])
            end_g = st_g + list_size_grid[i]
            end_m = st_m + list_size_mu[i]
            rand_arr = float(g.readline())
            #nzeros = np.size(mu_arr[i]) - np.size(grid_arr[i])
            #x = np.append(np.zeros(nzeros), X[np.array(grid_arr[i], dtype = 'i')])
            nzeros = list_size_mu[i] - list_size_grid[i]
            x = np.append(np.zeros(nzeros), X[np.array(grid_arr[st_g:end_g], dtype = 'i')])            
            mu = np.dot(mu_arr[st_m:end_m], x)
            X[i] = mu + rand_arr * sigma_arr[i]

        COR = np.reshape(X, [M,N])
        
        for i in range(0, multiple_maps+1):
            xx = 'map%i'%i
            X = np.multiply(COR, uncertaintydata[xx])
            DATA_NEW = data[xx]*np.exp(X)
            if j == num_realizations-1:
                data_new["map{0}".format(i)] = DATA_NEW
            if compute_loss == True:
                expobj = copy.deepcopy(init_dict[xx]['expobj'])
                mmi_data = MMI_conv2(DATA_NEW, voi, variables['dist_to_fault'], event_attr['magnitude'])

                fat_loss = model_loss(mmi_data, sm_dict[xx], expobj)
                expobj = copy.deepcopy(init_dict[xx]['expobj'])
                econ_est = model_econ(mmi_data, sm_dict[xx], expobj)
                ccodes_temp, fats_temp = [], []

                ccodes[xx][j] = fat_loss['ccodes']
                fats[xx][j] = np.sum(fat_loss['fats'])
                bestest[xx][j] = econ_est['BestEstimate']

        if np.mod(j+1, 25) == 0:
            print "Done with", j+1, "of", num_realizations, "iterations."

    g.close()

    
    g = open('fats_econs_my_rank%i_%i.txt'%(my_rank, num_realizations), 'w')
    for j in range(0, num_realizations):
        g.write(str(fats['map0'][j]) + '\n')
        g.write(str(bestest['map0'][j]) + '\n')

    g.close()


    if plot_on == True:
        if compute_loss == True:
            loss_data = {}
            for i in range(0, multiple_maps+1):
                xx = 'map%i'%i
                data_dict = {'ccodes':ccodes[xx], 'fats':fats[xx], 'bestest':bestest[xx]}
                loss_data[xx] = gather_results(total_real, my_reals, my_rank, data_dict, size, comm)

            if my_rank == 0:
                econ_xy_lim = np.zeros(4)
                fats_xy_lim = np.zeros(4)
                fats_xy_lim[0] = np.min(loss_data['map0']['fats'])
                econ_xy_lim[0] = np.min(loss_data['map0']['bestest'])

                for i in range(0, multiple_maps+1):
                    xx = 'map%i'%i
                # compute x,y ranges for histogram plots
                    if np.min(loss_data[xx]['fats']) < fats_xy_lim[0]:
                        fats_xy_lim[0] = np.min(loss_data[xx]['fats'])
                    if init_dict[xx]['loss_dict']['lkm20'] < fats_xy_lim[0]:
                        fats_xy_lim[0] = init_dict[xx]['loss_dict']['lkm20']-2
                    if np.max(loss_data[xx]['fats']) > fats_xy_lim[1]:
                        fats_xy_lim[1] = np.max(loss_data[xx]['fats'])
                    if init_dict[xx]['loss_dict']['lkp20'] > fats_xy_lim[1]:
                        fats_xy_lim[1] = init_dict[xx]['loss_dict']['lkp20']+2
                    if np.min(loss_data[xx]['bestest']) < econ_xy_lim[0]:
                        econ_xy_lim[0] = np.min(loss_data[xx]['bestest'])
                    if init_dict[xx]['loss_dict']['ekbm20'] < econ_xy_lim[0]:
                        econ_xy_lim[0] = init_dict[xx]['loss_dict']['ekbm20']
                    if np.max(loss_data[xx]['bestest']) > econ_xy_lim[1]:
                        econ_xy_lim[1] = np.max(loss_data[xx]['bestest'])
                    if init_dict[xx]['loss_dict']['ekbp20'] > econ_xy_lim[1]:
                        econ_xy_lim[1] = init_dict[xx]['loss_dict']['ekbp20']
                econ_xy_lim[3] = total_real/3
                fats_xy_lim[3] = total_real/3
                fats_int = (fats_xy_lim[1]-fats_xy_lim[0])/30
                econ_int = (econ_xy_lim[1]-econ_xy_lim[0])/30
                econ_xy_lim[0] -= econ_int
                econ_xy_lim[1] += econ_int

                for i in range(0, multiple_maps+1):
                    xx = 'map%i'%i
                    plot_fats(fats_xy_lim, fats_int, loss_data[xx]['fats'],
                              init_dict[xx]['loss_dict'], i)

                    plot_econ(econ_xy_lim, econ_int, loss_data[xx]['bestest'],
                              init_dict[xx]['loss_dict'], i)

        for i in range(0, multiple_maps+1):
            if my_rank == 0:
                plotting_data = {'cor':COR, 'data':sm_dict['map%i'%i], 'data_new':data_new['map%i'%i]}
                plot(plotting_data, variables, voi, sm_dict['map%i'%i], stationdata, i, plot_fault)     
                MMI_sm = MMI_conv2(sm_dict['map%i'%i].griddata, voi, variables['dist_to_fault'], event_attr['magnitude'])
                MMI_dn = MMI_conv2(data_new['map%i'%i], voi, variables['dist_to_fault'], event_attr['magnitude'])
                plot_intensity(MMI_sm, MMI_dn, sm_dict['map%i'%i], i, plot_fault)
                

def gather_results(total_real, my_reals, my_rank, data, size, comm):

    if my_rank == 0:
        ccodes = [None]*total_real
        fats = [None]*total_real
        bestest = [None]*total_real
        for xx in range(0, np.size(my_reals)):
            ccodes[my_reals[xx]]  = data['ccodes'][xx]
            fats[my_reals[xx]]    = data['fats'][xx]
            bestest[my_reals[xx]] = data['bestest'][xx]
        for xx in range(1, size):
            for realz in np.arange(xx, total_real, size):
                ccodes[realz] = comm.recv(source = xx, tag = 11*xx*realz)
                fats[realz]   = comm.recv(source = xx, tag = 13*xx*realz)
                bestest[realz]   = comm.recv(source = xx, tag = 15*xx*realz)
        # For now, just summing all fatalities and economic losses
#        for xx in range(0, np.shape(ccodes)[0]):
#            fats[xx] = sum(fats[xx])
            
        print 'Mean/std. dev. fatalities:', np.mean(fats), np.std(fats), '\n'
        print 'Economic Losses, average (dollars):'
        print 'Best estimate:', np.mean(bestest)

        loss_data= {'ccodes':ccodes, 'fats':fats, 'bestest':bestest}
    else:
        for xx in range(0, np.size(my_reals)):
            comm.send(data['ccodes'][xx], dest = 0, tag = 11*my_rank*my_reals[xx])
            comm.send(data['fats'][xx], dest = 0, tag = 13*my_rank*my_reals[xx])
            comm.send(data['bestest'][xx], dest = 0, tag = 15*my_rank*my_reals[xx])
        loss_data = None

    return loss_data

def init_loss_outputs(shakemap, voi, my_rank, dist, mag):
    expobj, outgrid = initialize_loss(shakemap, voi, dist, mag)
    if my_rank == 0:
        smdata = shakemap.griddata.copy()
        data = MMI_conv2(smdata, voi,  dist, mag)
        loss = model_loss(data, shakemap, expobj)
        lk = sum([int(loss['fats'][i]) for i in range(0, np.size(loss['fats']))])
        print lk
        loss = model_loss(data*1.2, shakemap, expobj)
        lkp20 = sum([int(loss['fats'][i]) for i in range(0, np.size(loss['fats']))])
        loss = model_loss(data*0.8, shakemap, expobj)
        lkm20 = sum([int(loss['fats'][i]) for i in range(0, np.size(loss['fats']))])
        loss = model_loss(data*1.1, shakemap, expobj)
        lkp10 = sum([int(loss['fats'][i]) for i in range(0, np.size(loss['fats']))])
        loss = model_loss(data*0.9, shakemap, expobj)
        lkm10 = sum([int(loss['fats'][i]) for i in range(0, np.size(loss['fats']))])
        loss = model_econ(data, shakemap, expobj)
        ekb = loss['BestEstimate']
        loss = model_econ(data*1.2, shakemap, expobj)
        ekbp20 = loss['BestEstimate']
        loss = model_econ(data*0.8, shakemap, expobj)
        ekbm20 = loss['BestEstimate']
        loss = model_econ(data*1.1, shakemap, expobj)
        ekbp10 = loss['BestEstimate']
        loss = model_econ(data*0.9, shakemap, expobj)
        ekbm10 = loss['BestEstimate']
        
        loss_dict = {'lk':lk, 'lkp20':lkp20,'lkm20':lkm20,'lkp10':lkp10,'lkm10':lkm10,
                     'ekb':ekb,'ekbp20':ekbp20, 'ekbm20':ekbm20, 'ekbp10':ekbp10,'ekbm10':ekbm10}
 
    else:
        loss_dict = {}
    expobj, outgrid = initialize_loss(shakemap, voi, dist, mag)

    return {'loss_dict':loss_dict, 'expobj':expobj}
