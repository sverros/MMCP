import h5py
import matplotlib.pyplot as plt

f = h5py.File('data_files/SVSM_0.hdf5', 'r')

dset = f['realization_0']
data = dset[:]

fig = plt.figure(figsize=(10,10))
map = plt.imshow(data)
plt.savefig('visualization_check.png')
plt.close(fig)
