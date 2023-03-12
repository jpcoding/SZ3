import numpy as np
from pathlib import Path
from pysz import SZ
import sys

# prepare your data in numpy array format
HOME = str(Path.home())
data = np.fromfile(HOME + '/data/hurricane_100x500x500/Uf48.bin.f32', dtype=np.float32)
data = np.reshape(data, (100, 500, 500))

# init SZ (both SZ2 and SZ3 are supported)
# Please change the path to the SZ dynamic library file in your system
lib_extention = {
    "darwin": "libSZ3c.dylib",
    "windows": "SZ3c.dll",
}.get(sys.platform, "libSZ3c.so")

sz = SZ("../../install/lib/{}".format(lib_extention))

# compress, both input and output data are numpy array
data_cmpr, cmpr_ratio = sz.compress(data, 0, 1e-6, 0, 0,config_file="/Users/pjiao/git/sz3-orig/tools/sz3/sz3.config")
print("compression ratio = {:5G}".format(cmpr_ratio))

# decompress, both input and output data are numpy array
data_dec = sz.decompress(data_cmpr, data.shape, data.dtype)

# verify
sz.verify(data, data_dec)
