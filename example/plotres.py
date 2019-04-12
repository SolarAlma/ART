import matplotlib.gridspec as gridspec
import numpy  as np

plt.close("all")
plt.ion()

# Open synthetic file

f = h5py.File('synthetic.h5', 'r')

# Define plot stuff

a = np.float32(f["Stokes_I"][:])
tit = ['1 mm', '3.1 mm']
f = plt.figure(figsize=(6,5))
ax = []
extent=(0,a.shape[2]*0.048, 0, a.shape[1]*0.048)
gs = gridspec.GridSpec(1,2)
ax.append(f.add_subplot(gs[0, 0]))
ax.append(f.add_subplot(gs[0, 1]))
gs.update( wspace=0.02, hspace=0.0,left=0.06, right=0.99, top=0.95, bottom=0.1)



# Make plot

ax[0].imshow(a[0,:,:,0], interpolation='nearest', aspect=1, extent=extent)
ax[1].imshow(a[0,:,:,1], interpolation='nearest', aspect=1, extent=extent)

# Change axes

ax[0].set_ylabel('y [Mm]')
ax[0].set_xlabel('x [Mm]')
ax[1].set_xlabel('x [Mm]')
ax[0].set_aspect(1.0)
ax[1].set_aspect(1.0)
ax[1].set_yticklabels([])
ax[0].axes.minorticks_on()
ax[1].axes.minorticks_on()
ax[0].set_title(tit[0])
ax[1].set_title(tit[1])


# Save figure

f.savefig('fig_test.pdf', format='pdf', bbox='standard', compression=1)

f.show()
