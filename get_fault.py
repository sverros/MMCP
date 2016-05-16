import numpy as np
import math
from scipy.special import sph_harm, lpmv
import sys
sys.path.insert(0, '/Users/sverros/Documents/Modules/shakemap/shakemap')
from shakelib.fault import Fault
from shakelib.distance import getDistance
from openquake.hazardlib.geo.mesh import Mesh
from neicio.readstation import readStation
from neicio.shake import ShakeGrid

voi = 'PGA'
direc = '/Users/sverros/Documents/Modules/MM_CP/input/'
shakemap = ShakeGrid(direc+'grid.xml', variable = '%s' % voi)
unc_INTRA = ShakeGrid(direc+'uncertainty.xml', variable = 'GMPE_INTRA_STD%s' % voi)
unc_INTER = ShakeGrid(direc+'uncertainty.xml', variable = 'GMPE_INTER_STD%s' % voi)
stationlist = direc+'stationlist.xml'
stationdata = readStation(stationlist)
attr = shakemap.getAttributes()
event_attr = attr['event']
grid_attr =  attr['grid_specification']
M = grid_attr['nlat']
N = grid_attr['nlon']
lats = np.transpose(np.ones([N,M])*np.linspace(grid_attr['lat_max'], grid_attr['lat_min'], M))
lons = np.ones([M,N])*np.linspace(grid_attr['lon_min'], grid_attr['lon_max'], N)

mesh = Mesh(lons, lats, np.zeros([M,N]))
method = 'rjb'
fault = Fault.readFaultFile(direc+'fault.txt')
quads = fault.getQuadrilaterals()

R_jb = getDistance(method, mesh, quads)

#g = open(direc+'fault.txt', 'r')
#points = []
#usd = 0.0
#lsd = 0.0
#
#for i,line in enumerate(g):
#    if i != 0:
#        lines = line.split()
#        points.append(Point(float(lines[1]), float(lines[0]), 0.0))
#        if i == 1:
#            usd = float(lines[2])
#            lsd = float(lines[2])
#        else:
#            if  usd > float(lines[2]):
#                usd = float(lines[2])
#            if  lsd < float(lines[2]):
#                lsd = float(lines[2])
#g.close()
#
#line = Line(points)
#fault = SimpleFaultSurface.from_fault_data(line, usd, lsd, 15, 1.0)
#
#
#R_jb = fault.get_joyner_boore_distance(mesh)


print R_jb





