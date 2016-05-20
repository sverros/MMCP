from mpi4py import MPI
import argparse
from main import run_method

# Get variables
parser = argparse.ArgumentParser(description='Run Method of Successive Conditional Simulations to compute realizations of a spatially correlated random field')
parser.add_argument('direc', metavar='path', type=str, 
                   help='path to the inputs directory or directory containing grid, stationlist, uncertainty, and event xml files and fault.txt')
parser.add_argument('voi', metavar='imt', type=str, 
                   help='string, intensity measures, i.e., pga pgv')
parser.add_argument('method', metavar='distance_measure', type=str, 
                   help='string, distance measure, i.e., rjb, rrup')
parser.add_argument('num_realizations', metavar='N', type = int, 
                   help='integer, number of realizations')
parser.add_argument('-out', dest='output_directory', type=str,
                    default = 'data_files/', 
                    help='path to output directory, defaults to data_files folder')
args = parser.parse_args()
direc = args.direc
num_realizations = args.num_realizations
voi = args.voi
voi = voi.split()
method = args.method
output_dir = args.output_directory

# Constant variables
radius = [45]
# For now the correlation model is constant, as openquake only include this correlation model. 
corr_model = 'JB2009'
vscorr = True

run_method(direc, voi, method, num_realizations, radius, corr_model, vscorr, output_dir)

 
