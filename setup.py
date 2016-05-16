import numpy as np
import math
from openquake.hazardlib.site import Site, SiteCollection
from openquake.hazardlib.geo import Point
from openquake.hazardlib.geo.geodetic import distance
from xml.dom import minidom
import sys
sys.path.insert(0, '/Users/sverros/Documents/Modules/shakemap/shakemap')
from shakelib.fault import Fault
from shakelib.distance import getDistance
from openquake.hazardlib.geo.mesh import Mesh

def initialize(shakemap, unc_INTRA, unc_INTER, stationdata, dir, method, dm=1, dn=1):
    """
    Set up grid spacing, site collections, and other data values
    
    INPUTS: 
    shakemap- shake grid of grid.xml
    uncertainty- shake grid of uncertainty.xml
    stationdata- data from stationlist.xml
    vscorr - boolean, determines if Vs30 are correlation. See
    JB2009
    dm, dn- vertical and horozontal spacing, respectively. Default value is 1
    OUT:
    Returns a dictionary with the following keys:    
    M,N - number of points vertically and horozontally
    K- number of stations
    uncertaintydata- array of uncertainty data at each point
    site_collection_SM- site collection for ShakeMap data
    site_collection_station- site collection for station data
    location_lat/lon_g- lat and lons for grid data, for plotting
    data- ShakeMap data at each point
    itensity- factor for non-native data
    site_collection_both- site collection for both SM and stations
    """
    attributes = shakemap.getAttributes()
    grid_attr = attributes['grid_specification']
    event_attr= attributes['event']

    # Determine the size of the grid                                                                
    m = grid_attr['nlat']
    n = grid_attr['nlon']

    # M,N number of vertical, horizontal points considered
    M = int(math.floor(m/dm))
    N = int(math.floor(n/dn))

    # Determine number of stations
    K = np.size(stationdata['lat'])

    bias = ComputeBias(dir+'uncertainty.xml')

    if bias == True:
        uncertainty = unc_INTRA.griddata
    else:
        rand = np.random.randn()
#        if abs(rand) > 3:
#            rand = np.sign(rand)*(6 - np.sign(rand)*rand)
        uncertainty = np.sqrt(np.power(rand*unc_INTER.griddata,2) + np.power(unc_INTRA.griddata,2))

    # Allocate matrices for storing data and locations
    DATA = shakemap.griddata
    location_lat_g = np.transpose(np.ones([N,M])*np.linspace(grid_attr['lat_max'], grid_attr['lat_min'], M))
    location_lon_g = np.ones([M,N])*np.linspace(grid_attr['lon_min'], grid_attr['lon_max'], N)
    site_SM = []

    

    if dm != 1:
        uncertainty1 = np.zeros([M,N])
        DATA1 = np.zeros([M,N])
        for i in range(0, M):
            for j in range(0, N):
                uncertainty1[i, j] = uncertainty[dm*i,dn*j]
                DATA1[i/dm, j] = DATA[dm*i,dn*j]
                #location_lat_g1[i/dm, j/dn] = location_lat_g[i,j]
                #location_lon_g1[i/dm, j/dn] = location_lon_g[i,j]
        uncertainty = uncertainty1
        DATA = DATA1


    # Puts griddata into Site class, then turns the sites into a SiteCollection
    for i in range(0,M):
        for j in range(0,N):
            site_SM.append(Site(location = Point(location_lon_g[i,j], location_lat_g[i,j]), 
                                vs30 = 760, vs30measured = True, z1pt0 = 100, z2pt5 = 1))

    site_collection_SM = SiteCollection(site_SM)

    # Store lat lon points
    location_lat_s = np.empty([K])
    location_lon_s = np.empty([K])
    intensity = np.empty([K])
    site_station = []

    # Puts stationdata into Site class, then turns the sites into a SiteCollection
    for i in range(0,K):
        location_lat_s[i] = stationdata['lat'][i]
        location_lon_s[i] = stationdata['lon'][i]
        if stationdata['name'][i] == 'DERIVED':
            intensity[i] = 1
        else:
            intensity[i] = 0

        site = Site(location=Point(location_lon_s[i], location_lat_s[i]),
                    vs30 = 760, vs30measured = True, z1pt0 = 100, z2pt5 = 1)
        site_station.append(site)

    site_collection_station = SiteCollection(site_station)
    site_both = site_station + site_SM
    site_collection_both = SiteCollection(site_both)


    #mesh = Mesh(location_lon_g, location_lat_g, np.zeros([M,N]))
    #fault = Fault.readFaultFile(dir+'fault.txt')
    #quads = fault.getQuadrilaterals()
    #
    #dist_to_fault = np.reshape(getDistance(method, mesh, quads), [M,N])

    dist_to_fault = distance(event_attr['lat'], event_attr['lon'], event_attr['depth'], location_lat_g, location_lon_g, np.zeros([M,N]))

    return {'M': M, 'N': N, 'K': K, 'uncertaintydata':uncertainty, \
                'site_collection_SM':site_collection_SM, \
                'site_collection_station':site_collection_station, \
                'location_lat_g': location_lat_g, 'location_lon_g':location_lon_g, \
                'data':DATA, 'intensity':intensity, \
                'site_collection_both':site_collection_both, \
                'dist_to_fault':dist_to_fault}

def ComputeBias(uncertainty):
    root = minidom.parse(uncertainty)
    evu = root.getElementsByTagName('event_specific_uncertainty')
    for instance in evu:
        value = instance.getAttribute('value')
        if float(value) == -1.0:
            return False
        else:
            return True

            
