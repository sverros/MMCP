import numpy as np
import math
import os
from openquake.hazardlib.site import Site, SiteCollection
from openquake.hazardlib.geo import Point
from openquake.hazardlib.geo.geodetic import distance
from openquake.hazardlib.geo.mesh import Mesh
from xml.dom import minidom
from shakemap.grind.fault import Fault
from shakemap.grind.distance import get_distance
from shakemap.grind.source import Source, read_event_file


def initialize(shakegrid, unc_grid, stationdata, dir, voi, method, dm=1, dn=1):
    """
    Set up grid spacing, site collections, and other data values.
    :param shakegrid:
        shakegrid object, shakegrid associated with grid.xml. See mapio.shake
    :param unc_grid:
        shakegrid object, shakegrid associated with uncertainty.xml
    :param stationdata:
        dict, data from stationlist.xml. See readstation.py
    :param dir:
        string, path to inputs folder
    :param method:
        string, distance measure method, i.e. 'rjb'
    :param dn dm:
        integer, specifies vertical and horozontal sampling interval.
    :returns:
        dict with the following keys,
        M,N - integer, number of points vertically and horozontally
    	K- integer, number of stations
    	uncertaintydata- dict, array of uncertainty data at each point and for each intensity measure
    	site_collection_SM- site collection object, site collection for ShakeMap data
    	site_collection_station- site collection object, site collection for station data
    	location_lat/lon_g- numpy array, lat and lons for grid data
    	data- dict, ShakeMap data at each point on grid and for each intensity measure
    	site_collection_both- site collection for both SM and stations
    	dist_to_fault- numpy array, distance from each point to fault
    """
    grid_attr = shakegrid.getGeoDict()
    event_attr= shakegrid.getEventDict()

    # Determine the size of the grid                                                                
    m = grid_attr.ny
    n = grid_attr.nx
    M = int(math.floor(m/dm))
    N = int(math.floor(n/dn))

    # Determine number of stations
    K = np.size(stationdata['lat'])

    uncertainty = {}
    DATA = {}

    # Compute the bias and use to determine uncertainty
    bias = ComputeBias(os.path.join(dir,'uncertainty.xml'))
    for i in range(0, np.size(voi)):
        unc_INTRA = unc_grid.getLayer('gmpe_intra_std%s' % voi[i])
        unc_INTER = unc_grid.getLayer('gmpe_inter_std%s' % voi[i])
        if bias == True:
            uncertainty1 = unc_INTRA.getData()
        else:
            rand = np.random.randn()
            uncertainty1 = np.sqrt(np.power(rand*unc_INTER.griddata,2) + np.power(unc_INTRA.griddata,2))

        shakemap = shakegrid.getLayer(voi[i])
        uncertainty[voi[i]] = np.zeros([M,N])
        DATA[voi[i]] = np.zeros([M,N])

        # Allocate matrices for storing data and locations
        DATA1 = shakemap.getData()
        for k in range(0, M):
            for j in range(0, N):
                uncertainty[voi[i]][k, j] = uncertainty1[dm*k,dn*j]
                DATA[voi[i]][k, j] = DATA1[dm*k,dn*j]

    site_SM = []
    location_lat_g = np.zeros([M,N])
    location_lon_g = np.zeros([M,N])

    for i in range(0, M):
        for j in range(0, N):
            location_lat_g[i/dm, j/dn], location_lon_g[i/dm, j/dn] = grid_attr.getLatLon(i,j)

    # Puts griddata into Site class, then turns the sites into a SiteCollection
    for i in range(0,M):
        for j in range(0,N):
            site_SM.append(Site(location = Point(location_lon_g[i,j], location_lat_g[i,j]), 
                                vs30 = 760, vs30measured = True, z1pt0 = 100, z2pt5 = 1))
    site_collection_SM = SiteCollection(site_SM)

    # Store lat lon points
    location_lat_s = np.empty([K])
    location_lon_s = np.empty([K])
    site_station = []

    # Puts stationdata into Site class, then turns the sites into a SiteCollection
    for i in range(0,K):
        location_lat_s[i] = stationdata['lat'][i]
        location_lon_s[i] = stationdata['lon'][i]

        site = Site(location=Point(location_lon_s[i], location_lat_s[i]),
                    vs30 = 760, vs30measured = True, z1pt0 = 100, z2pt5 = 1)
        site_station.append(site)

    site_collection_station = SiteCollection(site_station)
    site_both = site_station + site_SM
    site_collection_both = SiteCollection(site_both)

    fault = Fault.readFaultFile(os.path.join(dir,'fault.txt'))
    g = open(os.path.join(dir,'event.xml'), 'r')
    event = read_event_file(g)
    g.close()
    source = Source(event, fault)
    
    dist = get_distance(method, np.reshape(location_lat_g, [M*N,1]), np.reshape(location_lon_g, [M*N,1]), 
                                 np.zeros([M*N,1]), source)
    dist_to_fault = np.reshape(dist[method], [M,N])

    return {'M': M, 'N': N, 'K': K, 'uncertaintydata':uncertainty, \
                'site_collection_SM':site_collection_SM, \
                'site_collection_station':site_collection_station, \
                'location_lat_g': location_lat_g, 'location_lon_g':location_lon_g, \
                'data':DATA, \
                'site_collection_both':site_collection_both, \
                'dist_to_fault':dist_to_fault}

def ComputeBias(uncertainty):
    '''
    Determined whether bias is computed in ShakeMap
    :param uncertainty:
        string, path to uncertainty.xml
    :returns:
        boolean
    '''
    root = minidom.parse(uncertainty)
    evu = root.getElementsByTagName('event_specific_uncertainty')
    for instance in evu:
        value = instance.getAttribute('value')
        if float(value) == -1.0:
            return False
        else:
            return True

            
