from mpi4py import MPI
import argparse
import os
import sys
from main import run_method
import h5py
import matplotlib.pyplot as plt

# Start MPI                                                                                                                                                 
comm = MPI.COMM_WORLD
size = comm.Get_size()
my_rank = comm.Get_rank()

# Get variables
parser = argparse.ArgumentParser(description='Run Method of Successive Conditional Simulations to compute realizations of a spatially correlated random field')
parser.add_argument('direc', metavar='path', type=str,
                   help='path to the inputs directory or directory containing grid, stationlist, uncertainty, and event xml files and fault.txt')
parser.add_argument('voi', metavar='imt', type=str,
                   help='string, intensity measures, i.e., pga pgv')
parser.add_argument('num_realizations', metavar='N', type = int,
                   help='integer, number of realizations')
parser.add_argument('-r', dest='rad', type=int,
                    default = 45,
                    help='optional argument, radius, defaults to 45 km')
parser.add_argument('-out', dest='output_directory', type=str,
                    default = 'data_files/',
                    help='optional argument, path to output directory, defaults to data_files folder')

args = parser.parse_args()
direc = args.direc
num_realizations = args.num_realizations
voi = args.voi
voi = voi.split()
output_dir = args.output_directory
radius = [args.rad]

# Constant variables
# For now the correlation model is constant, as openquake only include this correlation model.
corr_model = 'JB2009'
vscorr = True

run_method(direc, voi, num_realizations, radius, corr_model, vscorr, output_dir)

sys.stdout.flush()

if my_rank == 0:
    
    f1 = h5py.File(os.path.join(output_dir,'Epsilon_pgv_0.hdf5'), 'r')
    
    
    dset1 = f1['realization_0']
    data1 = dset1[:]
    
    fig = plt.figure(figsize=(10,10))
    map = plt.imshow(data1)
    plt.savefig(os.path.join(output_dir,'visualization_check_pgv.png'))
    plt.close(fig)
    f1.close()
    
    f2 = h5py.File(os.path.join(output_dir,'Epsilon_pga_0.hdf5'), 'r')
    dset2 = f2['realization_0']
    data2 = dset2[:]
    
    fig = plt.figure(figsize=(10,10))
    map = plt.imshow(data2)
    plt.savefig(os.path.join(output_dir, 'visualization_check_pga.png'))
    plt.close(fig)
    
    f2.close()
