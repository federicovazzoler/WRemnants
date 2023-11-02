import narf
import h5py
import pdb

file = h5py.File("wremnants-data/data/angularCoefficients/w_z_coeffs_absY_scetlib_dyturboCorr.hdf5", "r")
res = narf.ioutils.pickle_load_h5py(file["results"])
print(res["Z"])