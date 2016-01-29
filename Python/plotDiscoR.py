import sys
import h5py as h5
import numpy as np
import matplotlib.pyplot as plt

def loadCheckpoint(filename):

    f = h5.File(filename, "r")

    piph = f['Data']['Cells'][:,-1][...]
    prim = f['Data']['Cells'][:,:-1][...]
    index = f['Grid']['Index'][...]
    idPhi0 = f['Grid']['Id_phi0'][...]
    nphi = f['Grid']['Np'][...]
    t = f['Grid']['T'][0]
    riph = f['Grid']['r_jph'][...]
    ziph = f['Grid']['z_kph'][...]

    r = np.zeros(piph.shape)
    R = 0.5*(riph[1:] + riph[:-1])
    for i in xrange(index.shape[0]):
        for k in xrange(index.shape[1]):
            ind0 = index[i,k]
            ind1 = ind0 + nphi[i,k]
            r[ind0:ind1] = R[i]

    return t, r, prim

def plotCheckpoint(file):
    
    print("Loading {0:s}...".format(file))

    t, r, prim = loadCheckpoint(file)

    print("   Plotting...")
    nq = prim.shape[1]

    fig, ax = plt.subplots(2,3,figsize=(14,9))
    plotAx(ax[0,0], r, prim[:,0], "linear", "linear", r"$r$", r"$\rho$", 'k+')
    plotAx(ax[0,1], r, prim[:,1], "linear", "linear", r"$r$", r"$P$", 'k+')
    plotAx(ax[1,0], r, prim[:,2], "linear", "linear", r"$r$", r"$u_r$", 'k+')
    plotAx(ax[1,1], r, prim[:,3], "linear", "linear", r"$r$", r"$u_\phi$",'k+')
    plotAx(ax[1,2], r, prim[:,4], "linear", "linear", r"$r$", r"$u_z$", 'k+')
    if nq > 5:
        plotAx(ax[0,2], r, prim[:,5], "linear", "linear", r"$r$", r"$q$", 'k+')

    title = "DISCO t = {0:.3g}".format(t)

    #fig.suptitle(title, fontsize=18)

    plt.tight_layout()

    name = file.split('.')[0].split('_')[-1]
    plotname = "plot_r_{0:s}.png".format(name)
    
    print("   Saving {0:s}...".format(plotname))
    fig.savefig(plotname)

    plt.close(fig)

def plotAx(ax, x, y, xscale, yscale, xlabel, ylabel, *args, **kwargs):
    ax.plot(x, y, *args, **kwargs)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)


if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("Makes plots of Disco prims as a function of r.")
        print("usage: python plotDiscoR.py <checkpoint.h5 ...>")
        sys.exit()

    files = sys.argv[1:]
    for f in files:
        plotCheckpoint(f)
