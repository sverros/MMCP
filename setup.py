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

def initialize(shakegrid, unc_INTRA, unc_INTER, stationdata, dir, voi, method, dm=1, dn=1):
    """
    Set up grid spacing, site collections, and other data values
    
    INPUTS: 
        shakemap- shake grid of grid.xml
        unc_INTRA- shake grid of intraevent uncertainty
        unc_INTER- shake grid og interevent uncertainty
        stationdata- data from stationlist.xml
        dir- directory of Inputs folder
        method- distance measure
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

    shakemap = shakegrid.getLayer(voi)

    #attributes = shakemap.getHeaderData()
    grid_attr = shakegrid.getGeoDict()
    event_attr= shakegrid.getEventDict()

    # Determine the size of the grid                                                                
    m = grid_attr.ny
    n = grid_attr.nx
    M = int(math.floor(m/dm))
    N = int(math.floor(n/dn))

    # Determine number of stations
    K = np.size(stationdata['lat'])

    # Compute the bias and use to determine uncertainty
    bias = ComputeBias(dir+'uncertainty.xml')
    if bias == True:
        uncertainty1 = unc_INTRA.getData()
    else:
        rand = np.random.randn()
        uncertainty1 = np.sqrt(np.power(rand*unc_INTER.griddata,2) + np.power(unc_INTRA.griddata,2))

    # Allocate matrices for storing data and locations
    DATA = shakemap.getData()
    site_SM = []
    location_lat_g = np.zeros([M,N])
    location_lon_g = np.zeros([M,N])
    uncertainty = np.zeros([M,N])
    DATA = np.zeros([M,N])

    for i in range(0, M):
        for j in range(0, N):
            uncertainty[i, j] = uncertainty1[dm*i,dn*j]
            DATA[i/dm, j] = DATA[dm*i,dn*j]
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

    mesh = Mesh(location_lon_g, location_lat_g, np.zeros([M,N]))
    fault = Fault.readFaultFile(dir+'fault.txt')
    quads = fault.getQuadrilaterals()
    
    dist_to_fault = np.reshape(getDistance(method, mesh, quads), [M,N])

    #dist_to_fault = distance(event_attr['lat'], event_attr['lon'], event_attr['depth'], location_lat_g, location_lon_g, np.zeros([M,N]))

    return {'M': M, 'N': N, 'K': K, 'uncertaintydata':uncertainty, \
                'site_collection_SM':site_collection_SM, \
                'site_collection_station':site_collection_station, \
                'location_lat_g': location_lat_g, 'location_lon_g':location_lon_g, \
                'data':DATA, 'intensity':intensity, \
                'site_collection_both':site_collection_both, \
                'dist_to_fault':dist_to_fault}

def ComputeBias(uncertainty):
    ''' ComputeBias
    INPUTS:
    uncertainty- file name of uncertainty.xml
    OUTPUTS: Boolean
    '''
    root = minidom.parse(uncertainty)
    evu = root.getElementsByTagName('event_specific_uncertainty')
    for instance in evu:
        value = instance.getAttribute('value')
        if float(value) == -1.0:
            return False
        else:
            return True

            
