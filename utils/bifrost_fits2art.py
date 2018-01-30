# This is how a referece file look like:

# dens                     Dataset {1, 25, 15, 269}
# dx                       Dataset {15}
# dy                       Dataset {25}
# temperature              Dataset {1, 25, 15, 269}
# xne                      Dataset {1, 25, 15, 269}
# z                        Dataset {1, 25, 15, 269}

# Imports

import h5py
import numpy as np
from astropy.io import fits
from astropy import units as u 

from IPython import embed

b_dens        = "BIFROST_en024048_hion2_lgr_385.fits"
b_temperature = "BIFROST_en024048_hion2_lgtg_385.fits"
b_xne         = "BIFROST_en024048_hion2_lgne_385.fits"
b_Pressure    = "BIFROST_en024048_hion2_lgp_385.fits"

b_dens = fits.open(b_dens)
b_temperature = fits.open(b_temperature)
b_xne = fits.open(b_xne)
b_Pressure = fits.open(b_Pressure)

dens =  np.einsum('ijk->jki',np.power(10,b_dens[0].data)) * (u.kg / (u.m**3))

dx = [ b_dens[0].header['DX'] * x for x in range(dens.shape[0]) ] * u.Mm
dy = [ b_dens[0].header['DY'] * y for y in range(dens.shape[1]) ] * u.Mm

temperature =  np.einsum('ijk->jki',np.power(10,b_temperature[0].data)) * u.K # in K
xne =  np.einsum('ijk->jki',np.power(10,b_xne[0].data)) / u.m**3 
pressure =  np.einsum('ijk->jki',np.power(10,b_Pressure[0].data)) * (u.N / u.m**2)

dz_vec = b_dens[1].data[:]
dz = np.broadcast_to(dz_vec, (dx.shape[0], dy.shape[0], dz_vec.shape[0])) * u.Mm

ny, nx, nz = dens.shape

embed()

# we might set cuts

# ny = 20
# nx = 10

# write h5 file

b = h5py.File("bifrost.h5","w")

h5_dens        = b.create_dataset("dens", (1,ny,nx,nz), dtype='float32')
h5_dx          = b.create_dataset("dx", (nx,), dtype='float32')
h5_dy          = b.create_dataset("dy", (ny,), dtype='float32')
h5_temperature = b.create_dataset("temperature", (1,ny,nx,nz), dtype='float32')
h5_xne         = b.create_dataset("xne", (1,ny,nx,nz), dtype='float32')
h5_z           = b.create_dataset("z", (1,ny,nx,nz), dtype='float32')
h5_Pgas        = b.create_dataset("Pgas", (1,ny,nx,nz), dtype='float32')

# fill with values

h5_dens[0,:,:,:]        = dens[0:ny,0:nx,::-1].cgs.value
h5_dx[:]                = dx[0:nx].cgs.value
h5_dy[:]                = dy[0:ny].cgs.value
h5_temperature[0,:,:,:] = temperature[0:ny,0:nx,::-1].cgs.value
h5_xne[0,:,:,:]         = xne[0:ny,0:nx,::-1].cgs.value
h5_z[0,:,:,:]           = dz[0:ny,0:nx,::-1].cgs.value
h5_Pgas[0,:,:,:]        = pressure[0:ny,0:nx,::-1].cgs.value

# fill with empty attributes

b.attrs["units"] = 0


b.flush()
b.close()
