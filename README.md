# MMCP

Dependencies:

mpi4py: http://pythonhosted.org/mpi4py/
	Install with conda:
	conda install --channel mpi4py mpich mpi4py 
numpy:  http://www.numpy.org/
mapio:  http://github.com/usgs/MapIO
	Install by using pip with git:
       	pip install git+git://github.com/usgs/MapIO.git
neicio (to be removed): https://github.com/usgs/neicio
       Install by using pip with git:
       pip install git+git://github.com/usgs/neicio.git
shakemap.shakelib: ??


openquake: http://www.globalquakemodel.org/openquake/about/
	Install by using pip with git:
       	pip install git+git://github.com/gem/oq-hazardlib.git



Installation:

Install with pip with git:
        pip install git+git://github.com/sverros/MMCP.git


Running Code:

xml files (grid.xml, stationlist.xml, uncertainty.xml) and fault file (fault.txt) should be in folder input

Main code is in test.py file, modify variables here.

Execute with: mpiexec -n xx python test.py
where xx is the number of cores.
