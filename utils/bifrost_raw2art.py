#!/usr/bin/env python

# This is how a referece file look like:

# dens                     Dataset {1, 25, 15, 269}
# dx                       Dataset {15}
# dy                       Dataset {25}
# temperature              Dataset {1, 25, 15, 269}
# xne                      Dataset {1, 25, 15, 269}
# z                        Dataset {1, 25, 15, 269}

# Imports

import getopt, sys
import h5py
import numpy as np

from astropy.io import fits
from astropy import units as u
from helita.sim import bifrost

help_message     = '''

    ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    :: SolarALMA: Bifrost raw to ART intput                   ::
    ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

     -i / --input  - root_name
     -n / --nsnap  - snap number
     -h / --help   - print this help msg
    
    ............................................................

'''

def main(argv=None):

	if argv is None:
		argv = sys.argv

	opts, args = getopt.getopt(argv[1:], "hi:n:", ["help","input=","nsnap="]) # double dots means that opt require input

	# option processing

	print(help_message)

	for option, value in opts:
		if option in ("-h", "--help"):
			print(help_message)
			return 0
		if option in ("-i", "--input"):
			root_name = value
			print('[!] root name : ', root_name)
		if option in ("-n", "--nsnap"):
			snap_no = int(value)
			print('[!] snap_no   : ', snap_no)

	snap = bifrost.BifrostData(root_name,snap_no,verbose=False)

	snap_has_hion = (snap.params['do_hion'] == 1)

	if snap_has_hion:
		print('[!] Will read xne and temperature from non-LTE solver')

	b_dens = snap.r * snap.params['u_r'] # now in * (u.g/u.cm**2)

	try:
		b_Pressure = snap.p * snap.params['u_p'] # now in * (u.dyn/u.cm**2)
	except:
		print('[!] Pressure will be read from eostable')
		eostab = bifrost.Rhoeetab()
		b_Pressure = eostab.tab_interp(b_dens, (snap.e * snap.params['u_e'])/(snap.r * snap.params['u_r']), out='pg') 

	if snap_has_hion:
		b_temperature = snap.hiontg
	else:
		try:
			b_temperature = snap.tg
		except:
			print('[!] Temperature will be read from eostable')
			eostab = bifrost.Rhoeetab()
			b_temperature = eostab.tab_interp(b_dens, (snap.e * snap.params['u_e'])/(snap.r * snap.params['u_r']), out='tg')

	if snap_has_hion:
		b_xne = snap.hionne # in * (1.0/u.cm**3)
	else:
		try:
			b_xne = snap.ne
		except:
			print('[!] Electron density number will be read from eostable')
			eostab = bifrost.Rhoeetab()
			b_xne = eostab.tab_interp(b_dens, (snap.e * snap.params['u_e'])/(snap.r * snap.params['u_r']), out='ne') 

	dens =  np.einsum('ijk->jik',b_dens)[:,:,::-1] * (u.g/u.cm**3) 

	dx = [ snap.dx * x for x in range(dens.shape[0]) ] * u.Mm
	dy = [ snap.dy * y for y in range(dens.shape[1]) ] * u.Mm

	temperature =  np.einsum('ijk->jik',b_temperature)[:,:,::-1] * u.K # in K
	xne =  np.einsum('ijk->jik',b_xne)[:,:,::-1] / u.cm**3 
	pressure =  np.einsum('ijk->jik',b_Pressure)[:,:,::-1] * (u.dyn/u.cm**2)

	dz_vec = - snap.z[::-1]
	dz = np.broadcast_to(dz_vec, (dx.shape[0], dy.shape[0], dz_vec.shape[0])) * u.Mm

	ny, nx, nz = dens.shape

	# we might set cuts

	# ny = 20
	# nx = 10

	# write h5 file

	final_file = "%s_%03d.h5" % (root_name,snap_no)

	b = h5py.File( final_file ,"w")

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

	print('[!] File %s saved!' % final_file)

if __name__ == '__main__':
	sys.exit(main())