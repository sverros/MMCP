# MMCP

Dependencies:

numpy:  http://www.numpy.org/

mpi4py: http://pythonhosted.org/mpi4py/
	Install with conda:
	conda install --channel mpi4py mpich mpi4py 

mapio:  http://github.com/usgs/MapIO
	Install by using pip with git:
       	pip install git+git://github.com/usgs/MapIO.git

shakemap: https://github.com/usgs/shakemap
	Install by using pip with git:
	pip install git+git://github.com/usgs/shakemap.git

openquake: http://www.globalquakemodel.org/openquake/about/
	Install by using pip with git:
       	pip install git+git://github.com/gem/oq-hazardlib.git

geos (required for openquake):
     Install with:
     conda install geos

h5py: http://docs.h5py.org/en/latest/quick.html#quick
      Install with:
      conda install h5py


Installation:

Install with pip with git:
        pip install git+git://github.com/sverros/MMCP.git


Running Code:

Inputs:

xml files (grid.xml, stationlist.xml, uncertainty.xml, event.xml) and fault file (fault.txt) should be in the input folder. 
If multiple ShakeMaps are to be used (multiple sets of grid, stationlist, and uncertainty files on the EXACT SAME grid), store one set
of files in the input folder, and store each subsequent set (grid, stationlist, uncertainty) in a subfolder named '1', '2', etc.

The variables that need to be specified when running this code are as follows:
    
path: The path to the input directory

imt: The intensity measure of interest. Should be a string, i.e., 'pga'

distance_measure: The distance measure for the site-to-fault distance, i.e., 'rjb'

N: The number of realizations to be computed.

This code is written with MPI and must be executed with the mpiexec command. The number of cores is specified at run-time.

Execute with: mpiexec -n xx python test.py path imt distance_measure N
where xx is the number of cores.

Outputs:

The outputs are stored in hdf5 files in the data_files folder. If a single grid.xml is used, then the ShakeMap with added variability will be written to file (SVSM_*.hdf5), otherwise the random fields will be written to file (Epsilon_*.hdf5). The * is used to denote the rank of the core that wrote the output.

Visualization of the output:

To plot the data, the visualizer.py script may be run. The script only plots one field at a time, specified using the file and realization key in the script. Dependencies include h5py and matplotlob, both of which can be installed using conda. The field is outputted as a png saved as visualization_check.png.


 
