#!/usr/bin/env python

# bifrost_to_art.py --- Script to make ART input from Bifrost raw snapshots
# M. Szydlarski [miko@astro.uio.no] - RoCS' 2019

# This is how a referece file look like:

# dens                     Dataset {1, 25, 15, 269}
# dx                       Dataset {15}
# dy                       Dataset {25}
# temperature              Dataset {1, 25, 15, 269}
# xne                      Dataset {1, 25, 15, 269}
# z                        Dataset {1, 25, 15, 269}

import getopt
import sys
import re
import os
import glob
import warnings
from os.path import relpath
import numpy as np
import getopt
import sys
from scipy.interpolate import interp1d


if (sys.version_info < (3, 0, 0)):
    sys.stderr.write("You need python 3.0 or later to run this script\n")
    sys.stderr.write('Your Python version: %d.%d.%d %s %d\n' %
                     sys.version_info)
    exit(1)
    pass

try:
    import h5py
    import astropy.units as u
except ImportError as error:
    print("")
    print("[!] " + error.__class__.__name__ + ": " + error.msg)
    print("")
    print("[!] No astropy or h5py ?")
    print("")

try:
    from helita.sim import bifrost
except ImportError as error:
    print("")
    print("[!] " + error.__class__.__name__ + ": " + error.msg)
    print("")
    print("[!] Looks like you don't have helita package")
    print("[!] Get it from here: https://github.com/ITA-Solar/helita.git")
    print("")

warnings.filterwarnings("ignore")

# --- Global variables -------------
global verbose, with_bfield

prog_name = os.path.basename(sys.argv[0])
prog_ver = 0.0
verbose = True
with_bfield = False  # include b-field
with_interpolation = False  # interpolate bifrost cubes to regular grid
zs_Mm = 14.0  # start is in corona
ze_Mm = -1.0  # end is in corona
high_res_nz = 512  # number of grid points in 'z'

help_message = '''
   :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
   :: Script to generate ART/*.h5 input from Bifrost raw snapshot
   :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    You can setup some parameters from command line:

    -i / --input    - "base_name" point to the snap file
                      e.g., -i 'cb24bih'
    -s / --snaps    - select snaps range [start:end] (seperate by ":")
                      e.g., 100:200
    -h / --help     - Print this help msg

    Examples
    --------
    This makes *.h5 from snaphots "refine_xxx.aux"
    where xxx is in range [200:382] and save it as 'refine_xxx.h5'.

    >>> ./bifrost_make_art.py -i refine -s 200:382
   ...........................................................................
'''


def get_indx_from_axis_val(axis_vec, val):
    return (np.abs(axis_vec-val)).argmin()


def main(argv=None):

    global help_message
    global config_f_only

    config_f_only = False

    if argv is None:
        argv = sys.argv

        try:
            opts, args = getopt.getopt(
                argv[1:], "hi:s:c", ["help", "input=", "snaps=", "cfonly"])
        except getopt.GetoptError:
            print(help_message)
            sys.exit(2)

        # option processing
        for option, value in opts:
            if option in ("-h", "--help"):
                print(help_message)
                return 0
            if option in ("-c", "--cfonly"):
                config_f_only = True
            if option in ("-i", "--input"):
                base_name = value
            if option in ("-s", "--snaps"):
                srange = value
                if (':' in srange):
                    snap_range = range(
                        int(srange.split(':')[0]), int(srange.split(':')[1]))
                else:
                    snap_range = [int(srange)]
            pass

    # Read object
    bobj = bifrost.BifrostData(base_name, snap=snap_range[0], verbose=False)

    # check if hion
    has_hion = (bobj.params['do_hion'][0] == 1)

    # Read geometry of the cube (once is enough)
    x = bobj.x
    y = bobj.y
    z = - bobj.z  # makes positive values == above photospere

    # At least for ALMA we don't need whole box so we can cut corona & bottom
    # TODO: make it as argument

    zs = get_indx_from_axis_val(z, zs_Mm)
    ze = get_indx_from_axis_val(z, ze_Mm)

    if with_interpolation:
        hrz = np.linspace(-0.5, 11.0, high_res_nz)[::-1]

    print(
        'Takes only z range from: %d [%f Mm] --> %d [%f Mm]' % (zs, z[zs], ze, z[ze]))

    for snap in snap_range:

        name_template = 'bifrost_%s_%04d_%s.%s'

        out_file = name_template % (base_name, snap, 'model', 'h5')
        cfg_file = name_template % (base_name, snap, 'art', 'cfg')
        int_file = name_template % (base_name, snap, 'int', 'h5')

        if config_f_only:

            # Create config file for ART

            print("[Write *.cfg file for ART]")

            with open(cfg_file, 'w') as f:
                f.write('input_model = %s\n' % out_file)
                f.write('output_profiles = %s\n' % int_file)
                f.write('\n')
                f.write('rt_solver = 0\n')
                f.write('verbose = 1\n')
                f.write('mu = 1.0\n')
                f.write('\n')
                f.write('#can be witt (default) or pisk (slower)\n')
                f.write('#eos = witt\n')
                f.write('\n')
                f.write('# Cut the model from above until T = temperature_cut\n')
                f.write('temperature_cut = 500000.0\n')
                f.write('logfile = null\n')
                f.write('get_contribution_funtion = 0\n')
                f.write('\n')
                # Photosphere
                f.write('region =     5000.000000, 1.0, 1, 1.0, none, none\n')
                # Band 3
                # lambda = 3.258614 [mm]
                f.write('region = 32586136.739130, 1.0, 1, 1.0, none, none\n')
                # lambda = 3.189281 [mm]
                f.write('region = 31892814.680851, 1.0, 1, 1.0, none, none\n')
                # lambda = 3.122838 [mm]
                f.write('region = 31228381.041667, 1.0, 1, 1.0, none, none\n')
                # lambda = 2.882620 [mm]
                f.write('region = 28826197.884615, 1.0, 1, 1.0, none, none\n')
                # lambda = 2.828231 [mm]
                f.write('region = 28282307.358491, 1.0, 1, 1.0, none, none\n')
                # lambda = 2.775856 [mm]
                f.write('region = 27758560.925926, 1.0, 1, 1.0, none, none\n')
                # Band 6
                # lambda = 1.309137 [mm]
                f.write('region = 13091373.711790, 1.0, 1, 1.0, none, none\n')
                # lambda = 1.297803 [mm]
                f.write('region = 12978028.484848, 1.0, 1, 1.0, none, none\n')
                # lambda = 1.286663 [mm]
                f.write('region = 12866629.098712, 1.0, 1, 1.0, none, none\n')
                # lambda = 1.223643 [mm]
                f.write('region = 12236426.857143, 1.0, 1, 1.0, none, none\n')
                # lambda = 1.213735 [mm]
                f.write('region = 12137346.477733, 1.0, 1, 1.0, none, none\n')
                # lambda = 1.203986 [mm]
                f.write('region = 12039857.751004, 1.0, 1, 1.0, none, none\n')
                # Band 7
                # lambda = 0.885388 [mm]
                f.write('region =  8853882.398110, 1.0, 1, 1.0, none, none\n')
                # lambda = 0.880189 [mm]
                f.write('region =  8801892.483852, 1.0, 1, 1.0, none, none\n')
                # lambda = 0.875051 [mm]
                f.write('region =  8750509.573847, 1.0, 1, 1.0, none, none\n')
                # lambda = 0.855084 [mm]
                f.write('region =  8550840.216771, 1.0, 1, 1.0, none, none\n')
                # lambda = 0.850234 [mm]
                f.write('region =  8502338.570618, 1.0, 1, 1.0, none, none\n')
                # lambda = 0.845438 [mm]
                f.write('region =  8454384.038353, 1.0, 1, 1.0, none, none\n')

            print("[Done] for: %s" % cfg_file)

        else:

            bobj.set_snap(snap)
            r = bobj.get_var('r')[:, :, zs:ze] * bobj.params['u_r'][0]

            if with_interpolation:
                fint = interp1d(z[zs:ze],  r, axis=-1, kind='cubic')
                r = fint(hrz)

            if has_hion:
                tg = bobj.get_var('hiontg')[:, :, zs:ze]
            else:
                try:
                    tg = bobj.get_var('tg')[:, :, zs:ze]
                except:
                    print('[!] Temperature will be read from eostable')
                    eostab = bobj.Rhoeetab()
                    ienergy = (bobj.get_var('e')[:, :, zs:ze]
                               * bobj.params['u_e'][0]) / r
                    tg = eostab.tab_interp(rho, ienergy, out='tg')

            if with_interpolation:
                fint = interp1d(z[zs:ze],  tg, axis=-1, kind='cubic')
                tg = fint(hrz)

            try:
                p = bobj.get_var('p')[:, :, zs:ze] * bobj.params['u_p'][0]
            except:
                print('[!] Pressure will be read from eostable')
                eostab = bobj.Rhoeetab()
                ienergy = (bobj.get_var('e')[:, :, zs:ze]
                           * bobj.params['u_e'][0]) / r
                p = eostab.tab_interp(
                    rho, ienergy, out='p') * bobj.params['u_p'][0]

            if with_interpolation:
                fint = interp1d(z[zs:ze],  p, axis=-1, kind='cubic')
                p = fint(hrz)

            if has_hion:
                nel = bobj.get_var('hionne')[:, :, zs:ze]
            else:
            	 nel = (bobj.get_electron_density().value)[:, :, zs:ze]

            if with_interpolation:
                fint = interp1d(z[zs:ze],  nel, axis=-1, kind='cubic')
                nel = fint(hrz)

            nx, ny, nz = tg.shape

            if with_bfield:
                # move field to cell center
                bx = bifrost.cstagger.xup(bobj.get_var('bx'))[
                    :, :, zs:ze] * bobj.params['u_b'][0]
                by = bifrost.cstagger.yup(bobj.get_var('by'))[
                    :, :, zs:ze] * bobj.params['u_b'][0]
                bz = bifrost.cstagger.zup(bobj.get_var('bz'))[
                    :, :, zs:ze] * bobj.params['u_b'][0]

            vx = np.asarray(range(nx)) * bobj.params['dx'][0] * u.Mm
            vy = np.asarray(range(ny)) * bobj.params['dy'][0] * u.Mm
            if with_interpolation:
                vz = np.broadcast_to(hrz, (ny, nx, nz)) * u.Mm
            else:
                vz = np.broadcast_to(z[zs:ze], (ny, nx, nz)) * u.Mm

            if os.path.exists(out_file):
                os.remove(out_file)

            with h5py.File(out_file, 'w') as f:
                h5_dens = f.create_dataset(
                    "dens", (1, ny, nx, nz), dtype='float32')
                h5_dx = f.create_dataset("dx", (nx,), dtype='float32')
                h5_dy = f.create_dataset("dy", (ny,), dtype='float32')
                h5_temperature = f.create_dataset(
                    "temperature", (1, ny, nx, nz), dtype='float32')
                h5_xne = f.create_dataset(
                    "xne", (1, ny, nx, nz), dtype='float32')
                h5_z = f.create_dataset("z", (1, ny, nx, nz), dtype='float32')
                h5_Pgass = f.create_dataset(
                    "Pgas", (1, ny, nx, nz), dtype='float32')

                if with_bfield:
                    h5_bx = f.create_dataset(
                        "bx", (1, ny, nx, nz), dtype='float32')
                    h5_by = f.create_dataset(
                        "by", (1, ny, nx, nz), dtype='float32')
                    h5_bz = f.create_dataset(
                        "bz", (1, ny, nx, nz), dtype='float32')

                # fill with values

                h5_dens[0, :, :, :] = np.einsum('ijk->jik', r)
                h5_dx[:] = vx.cgs.value
                h5_dy[:] = vy.cgs.value
                h5_temperature[0, :, :, :] = np.einsum('ijk->jik', tg)
                h5_xne[0, :, :, :] = np.einsum('ijk->jik', nel)
                h5_z[0, :, :, :] = vz.cgs.value
                h5_Pgass[0, :, :, :] = np.einsum('ijk->jik', p)

                if with_bfield:
                    h5_bx[0, :, :, :] = np.einsum('ijk->jik', bx)
                    h5_by[0, :, :, :] = np.einsum('ijk->jik', by)
                    h5_bz[0, :, :, :] = np.einsum('ijk->jik', bz)

                # fill with empty attributes

                f.attrs["units"] = 0
                f.flush()
                f.close()
                print("[DONE] for: %s" % out_file)
            print("[CHECKING] final *.h5 file")
            # Check output
            chc = h5py.File(out_file, "r")
            for i in chc.items():
                print(i)

        print("[ALL][DONE] for snap: %d" % snap)


if __name__ == '__main__':
    sys.exit(main())
