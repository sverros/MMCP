import sys
import numpy as np
import matplotlib.pyplot as plt
import cartopy
from matplotlib import cm
import matplotlib.colors as clors
from matplotlib.colors import LinearSegmentedColormap

WATER_COLOR = [.47,.60,.81]

def plot_econ(lim, int, data, ld, number):
    xmin, xmax, ymin, ymax = lim
    bins = np.arange(xmin,xmax,int)
    fig = plt.figure(figsize = (10,10))
    plt.title('Histogram of estimated economic losses (%i)\nR = 45 km, mean = %1.2e, stddev = %1.2e (billions)' % ( number, (np.mean(data)/10e9), (np.std(data))/10e9), y = 1.02)
    plt.xlabel('estimated economic losses')
    plt.ylabel('number of realizations')
    n, bins, patches = plt.hist(data, bins, normed=0, histtype='bar')                                                                               
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.setp(patches, 'facecolor', 'g', 'alpha', 0.75)
    plt.axvline(np.median(data), lw = 4,label = 'median of SCS loss')
    tag = 'ekb'
    plt.axvline(ld[tag],          color = 'r', lw = 4,label = 'median of SM loss, %1.2e'%(ld[tag]/10e9))
    plt.axvline(ld[tag + 'p20'],  color = 'k', linestyle = '--', lw = 4, label = 'SM +%%20 gm Loss, %1.2e'%(ld[tag+'p20']/10e9))
    plt.axvline(ld[tag + 'm20'],  color = 'k', linestyle = '--', lw = 4, label = 'SM -%%20 gm Loss, %1.2e'%(ld[tag+'m20']/10e9))
    plt.axvline(ld[tag + 'p10'],  color = 'k', linestyle = ':', lw = 4, label = 'SM +%%10 gm Loss, %1.2e'%(ld[tag+'p10']/10e9))
    plt.axvline(ld[tag + 'm10'],  color = 'k', linestyle = ':', lw = 4, label = 'SM -%%10 gm Loss, %1.2e'%(ld[tag+'m10']/10e9))
    plt.legend(framealpha=1)
    plt.savefig('figures/economic_losses%i.png'%(number))
    plt.close(fig)
    return
 
def plot_fats(lim, int, data, ld, number):
    xmin, xmax, ymin, ymax = lim
    bins = np.arange(xmin,xmax,int)
    fig = plt.figure(figsize = (10,10))
    plt.title('Histogram of estimated fatalities (%i)\nR = 45 km, mean = %f, stddev = %f' % ( number, np.mean(data), np.std(data)), y = 1.02)
    plt.xlabel('estimated fatalities')
    plt.ylabel('number of realizations')
    n, bins, patches = plt.hist(data, bins, normed=0, histtype='bar')                                                                               
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.setp(patches, 'facecolor', 'g', 'alpha', 0.75)
    plt.axvline(np.median(data), lw = 4,label = 'median of SCS loss')
    tag = 'lk'
    plt.axvline(ld[tag],          color = 'r', lw = 4,label = 'median of SM loss, %i'%ld[tag])
    plt.axvline(ld[tag + 'p20'],  color = 'k', linestyle = '--', lw = 4, label = 'SM +%%20 gm Loss, %i'%ld[tag+'p20'])
    plt.axvline(ld[tag + 'm20'],  color = 'k', linestyle = '--', lw = 4, label = 'SM -%%20 gm Loss, %i'%ld[tag+'m20'])
    plt.axvline(ld[tag + 'p10'],  color = 'k', linestyle = ':', lw = 4, label = 'SM +%%10 gm Loss, %i'%ld[tag+'p10'])
    plt.axvline(ld[tag + 'm10'],  color = 'k', linestyle = ':', lw = 4, label = 'SM -%%10 gm Loss, %i'%ld[tag+'m10'])
    plt.legend(framealpha=1)
    plt.savefig('figures/fatalities%i.png'%(number))
    plt.close(fig)
 
    return                    

def plot(out, variables, voi, shakemap, stationdata, number, plot_fault):
    if voi == 'PGA':
        units = 'pctg'
    elif voi == 'PGV':
        units = 'cms'


    attributes = shakemap.getAttributes()
    intensity = stationdata['name']
    SM = []
    IN = []
    for value in enumerate(intensity):
        if ((value == 'UNCERTAINTY')or(value == 'DYFI')or(value == 'MMI')or(value == 'CIIM')):
            IN.append(value[0])
        else:
            SM.append(value[0])

    sm_station_lons = [stationdata['lon'][j] for j in SM]
    sm_station_lats = [stationdata['lat'][j] for j in SM]
    in_station_lats = [stationdata['lat'][j] for j in IN]
    in_station_lons = [stationdata['lon'][j] for j in IN]            
    palette = cm.jet

    fig = plt.figure(figsize=(10,10))
    proj = cartopy.crs.PlateCarree()

    ax = plt.axes(projection=proj)
    cartopy.feature.COASTLINE.scale = '50m'
    cartopy.feature.LAND.scale = '50m'
    cartopy.feature.OCEAN.scale = '50m'
    ax.add_feature(cartopy.feature.OCEAN,facecolor=WATER_COLOR)
    ax.add_feature(cartopy.feature.COASTLINE)
    ax.add_feature(cartopy.feature.BORDERS, linestyle=':',zorder=10)
    ax.gridlines(crs=proj, draw_labels=True, linestyle='-')
    ax.set_extent(shakemap.getRange())
    map = ax.imshow(out['cor'],extent=shakemap.getRange(), origin='upper',cmap=palette)
    plt.plot(sm_station_lons, sm_station_lats, 'g>', markersize = 6)
    plt.plot(in_station_lons, in_station_lats, 'r^', markersize = 6)
    if plot_fault == True:
        points = []
        f = open('input/fault.txt', 'r')
        for line in f:
            if line[0] != '#':
                line_parts = line.split()
                points.append(float(line_parts[0]))
                points.append(float(line_parts[1]))
        
        f.close()
        points = np.asarray(points)
        points = np.reshape(points, [np.size(points)/2, 2])

        plt.plot(points[:,1], points[:,0], color = 'k', linewidth = 2, transform=proj)
        
        if number != 0:
            f2 = open('input/fault_points.txt', 'r')
            for i, line in enumerate(f2):    
                if i == number-1:
                    plot_point = line.split()
            f2.close()
            print float(plot_point[1]), float(plot_point[0])
            plt.plot(float(plot_point[1]), float(plot_point[0]), marker='D', color = 'k', transform=proj)

    locstr = attributes['event']['event_description']
    mag = attributes['event']['magnitude']
    datestr = attributes['event']['event_timestamp'].strftime('%b %d, %Y %H:%M:%S')
    th = plt.title('Correlation Matrix (%i)\n %s - M%.1f\n %s, units epsilon' % (number, locstr,mag,voi), y = 1.08)
    ch=plt.colorbar(map, shrink=0.6, pad = 0.1)
    plt.savefig('figures/epsilon%i.png'%number)
    plt.close(fig)

    fig = plt.figure(figsize = (10,10))
    proj = cartopy.crs.PlateCarree()
    ax = plt.axes(projection=proj)
    cartopy.feature.COASTLINE.scale = '50m'
    cartopy.feature.LAND.scale = '50m'
    cartopy.feature.OCEAN.scale = '50m'
    ax.add_feature(cartopy.feature.OCEAN,facecolor=WATER_COLOR)
    ax.add_feature(cartopy.feature.COASTLINE)
    ax.add_feature(cartopy.feature.BORDERS, linestyle=':',zorder=10)
    ax.gridlines(crs=proj, draw_labels=True, linestyle='-')
    ax.set_extent(shakemap.getRange())
    map = ax.imshow(np.log(shakemap.griddata),extent=shakemap.getRange(), origin='upper',cmap=palette, vmin = 0, vmax = 5)

    if plot_fault == True:
        plt.plot(points[:,1], points[:,0], color = 'k', linewidth = 2, transform=proj)
        if number != 0:
            plt.plot(float(plot_point[1]), float(plot_point[0]), marker='D', color = 'k', transform=proj)


    plt.plot(sm_station_lons, sm_station_lats, 'g>', markersize = 6)
    plt.plot(in_station_lons, in_station_lats, 'r^', markersize = 6)
    locstr = attributes['event']['event_description']
    mag = attributes['event']['magnitude']
    datestr = attributes['event']['event_timestamp'].strftime('%b %d, %Y %H:%M:%S')
    th = plt.title('ShakeMap (%i)\n %s - M%.1f\n%s, units ln(%s)' % (number, locstr,mag, voi,units), y = 1.08)
    ch=plt.colorbar(map, shrink=0.6, pad = 0.1)
    plt.savefig('figures/shakemap%i.png'%number)
    plt.close(fig)

    fig = plt.figure(figsize = (10,10))
    proj = cartopy.crs.PlateCarree()
    ax = plt.axes(projection=proj)
    cartopy.feature.COASTLINE.scale = '50m'
    cartopy.feature.LAND.scale = '50m'
    cartopy.feature.OCEAN.scale = '50m'
    ax.add_feature(cartopy.feature.OCEAN,facecolor=WATER_COLOR)
    ax.add_feature(cartopy.feature.COASTLINE)
    ax.add_feature(cartopy.feature.BORDERS, linestyle=':',zorder=10)
    ax.gridlines(crs=proj, draw_labels=True, linestyle='-')
    ax.set_extent(shakemap.getRange())
    map = ax.imshow(np.log(out['data_new']),extent=shakemap.getRange(), origin='upper',cmap=palette, vmin = 0, vmax = 5)
    if plot_fault == True:
        plt.plot(points[:,1], points[:,0], color = 'k', linewidth = 2, transform=proj)
        if number != 0:
            plt.plot(float(plot_point[1]), float(plot_point[0]), marker='D', color = 'k', transform=proj)

    plt.plot(sm_station_lons, sm_station_lats, 'g>', markersize = 6)
    plt.plot(in_station_lons, in_station_lats, 'r^', markersize = 6)
    locstr = attributes['event']['event_description']
    mag = attributes['event']['magnitude']
    datestr = attributes['event']['event_timestamp'].strftime('%b %d, %Y %H:%M:%S')
    th = plt.title('ShakeMap with added variability (%i)\n %s - M%.1f\n%s, units ln(%s)' % (number, locstr,mag, voi, units), y = 1.08)
    ch=plt.colorbar(map, shrink=0.6, pad = 0.1)
    plt.savefig('figures/scm%i.png'%number)
    plt.close(fig)


    return


def plot_intensity(output1, output2, shakemap, number, plot_fault):

    attributes = shakemap.getAttributes()
    palette = make_intensity_cmap()
    bounds = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    fig = plt.figure(figsize=(10,10))
    proj = cartopy.crs.PlateCarree()
    ax = plt.axes(projection=proj)
    cartopy.feature.COASTLINE.scale = '50m'
    cartopy.feature.LAND.scale = '50m'
    cartopy.feature.OCEAN.scale = '50m'
    ax.add_feature(cartopy.feature.OCEAN,facecolor=WATER_COLOR)
    ax.add_feature(cartopy.feature.COASTLINE)
    ax.add_feature(cartopy.feature.BORDERS, linestyle=':',zorder=10)
    ax.gridlines(crs=proj, draw_labels=True, linestyle='-')
    ax.set_extent(shakemap.getRange())
    map = ax.imshow(output1,extent=shakemap.getRange(), origin='upper',cmap=palette, vmin = 0, vmax = 10)
    if plot_fault == True:
        points = []
        f = open('input/fault.txt', 'r')
        for line in f:
            if line[0] != '#':
                line_parts = line.split()
                points.append(float(line_parts[0]))
                points.append(float(line_parts[1]))

        f.close()
        points = np.asarray(points)
        points = np.reshape(points, [np.size(points)/2, 2])

        plt.plot(points[:,1], points[:,0], color = 'k', linewidth = 2, transform=proj)

        if number != 0:
            f2 = open('input/fault_points.txt', 'r')
            for i, line in enumerate(f2):
                if i == number-1:
                    plot_point = line.split()
            f2.close()
            print float(plot_point[1]), float(plot_point[0])
            plt.plot(float(plot_point[1]), float(plot_point[0]), marker='D', color = 'k', transform=proj)

    locstr = attributes['event']['event_description']
    mag = attributes['event']['magnitude']
    datestr = attributes['event']['event_timestamp'].strftime('%b %d, %Y %H:%M:%S')
    th = plt.title('ShakeMap (%i)\n %s - M%.1f\n, units MI' % (number, locstr,mag), y = 1.08)
    ch=plt.colorbar(map, shrink=0.6, pad = 0.1, cmap = palette, ticks=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    plt.savefig('figures/sm_intensity%i.png'%number)
    plt.close(fig)

    fig = plt.figure(figsize=(10,10))
    proj = cartopy.crs.PlateCarree()
    ax = plt.axes(projection=proj)
    cartopy.feature.COASTLINE.scale = '50m'
    cartopy.feature.LAND.scale = '50m'
    cartopy.feature.OCEAN.scale = '50m'
    ax.add_feature(cartopy.feature.OCEAN,facecolor=WATER_COLOR)
    ax.add_feature(cartopy.feature.COASTLINE)
    ax.add_feature(cartopy.feature.BORDERS, linestyle=':',zorder=10)
    ax.gridlines(crs=proj, draw_labels=True, linestyle='-')
    ax.set_extent(shakemap.getRange())
    map = ax.imshow(output2,extent=shakemap.getRange(), origin='upper',cmap=palette, vmin = 0, vmax = 10)
    if plot_fault == True:
        plt.plot(points[:,1], points[:,0], color = 'k', linewidth = 2, transform=proj)
        if number != 0:
            plt.plot(float(plot_point[1]), float(plot_point[0]), marker='D', color = 'k', transform=proj)
    locstr = attributes['event']['event_description']
    mag = attributes['event']['magnitude']
    datestr = attributes['event']['event_timestamp'].strftime('%b %d, %Y %H:%M:%S')
    th = plt.title('ShakeMap with Variability (%i)\n %s - M%.1f\n, units MI' % (number, locstr,mag), y = 1.08)
    ch=plt.colorbar(map, shrink=0.6, pad = 0.1, cmap = palette, ticks = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    plt.savefig('figures/sm_wd_intensity%i.png'%number)
    plt.close(fig)

def make_intensity_cmap():

    cdict = { 'red' : ((0,   255.0/255.0, 255.0/255.0),
                       (0.1, 255.0/255.0, 255.0/255.0),
                       (0.2, 191.0/255.0, 191.0/255.0),
                       (0.3, 160.0/255.0, 160.0/255.0),
                       (0.4, 128.0/255.0, 128.0/255.0),
                       (0.5, 122.0/255.0, 122.0/255.0),
                       (0.6, 255.0/255.0, 255.0/255.0),
                       (0.7, 255.0/255.0, 255.0/255.0),
                       (0.8, 255.0/255.0, 255.0/255.0),
                       (0.9, 255.0/255.0, 255.0/255.0),
                       (1.0, 200.0/255.0, 200.0/255.0)),
              'green':((0,   255.0/255.0, 255.0/255.0),
                       (0.1, 255.0/255.0, 255.0/255.0),
                       (0.2, 204.0/255.0, 204.0/255.0),
                       (0.3, 230.0/255.0, 230.0/255.0),
                       (0.4, 255.0/255.0, 255.0/255.0),
                       (0.5, 255.0/255.0, 255.0/255.0),
                       (0.6, 255.0/255.0, 255.0/255.0),
                       (0.7, 200.0/255.0, 200.0/255.0),
                       (0.8, 145.0/255.0, 145.0/255.0),
                       (0.9,   0,   0),
                       (1.0,   0,   0)),
              'blue': ((0,   255.0/255.0, 255.0/255.0),
                       (0.1, 255.0/255.0, 255.0/255.0),
                       (0.2, 255.0/255.0, 255.0/255.0),
                       (0.3, 255.0/255.0, 255.0/255.0),
                       (0.4, 255.0/255.0, 255.0/255.0),
                       (0.5, 147.0/255.0, 147.0/255.0),
                       (0.6, 0, 0),
                       (0.7, 0, 0),
                       (0.8, 0, 0),
                       (0.9, 0, 0),
                       (1.0,0, 0))}

    my_cmap = LinearSegmentedColormap('my_cmap', cdict)
    return my_cmap
