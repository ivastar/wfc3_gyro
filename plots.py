import matplotlib
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
import matplotlib.gridspec as gridspec
import pylab
import unicorn
import astropy
from astropy.table import Table as table

def mag_radius():
    
    from matplotlib.ticker import ScalarFormatter
    
    gabe = table.read('../REF/cosmos-wide_v0.3_drz_bkg_sci.cat', format='ascii.sextractor')
    orig = table.read('/3DHST/Photometry/Release/v4.0/COSMOS/Detection/cosmos_3dhst.v4.0.F160W_orig.cat', format='ascii.sextractor')
    cat = table.read('test2_drz_sci.cat', format='ascii.sextractor')
    
    top = 0.075
    bottom = 0.1
    left = 0.15
    fig = unicorn.catalogs.plot_init(xs=10,aspect=0.5, left=left, right=0.1, bottom=bottom, top=top, NO_GUI=False)
    fig.subplots_adjust(wspace=0.15)
    fig.subplots_adjust(hspace=0.1)
    fs = 10
    matplotlib.rc('xtick',labelsize=fs)
    matplotlib.rc('ytick',labelsize=fs)
        
    gs1 = gridspec.GridSpec(1,2)
    
    ax1 = fig.add_subplot(gs1[0,0], ylim=[1, 20], xlim=[14,26])
    ax1.plot(orig['MAG_AUTO'],orig['FLUX_RADIUS'], '.', color='0.7',markersize=0.5)
    ax1.plot(cat['MAG_AUTO'],cat['FLUX_RADIUS']*0.1/0.06, 'o', color='black', markersize=1, alpha=0.5)
    ax1.set_xlabel('MAG_AUTO')
    ax1.set_ylabel('FLUX_RADIUS')
    ax1.set_yscale('log')
    ax1.set_title('ALIGNED READS')
    ax1.yaxis.set_major_formatter(ScalarFormatter())
    
    ax2 = fig.add_subplot(gs1[0,1], ylim=[1, 20], xlim=[14,26])
    ax2.plot(orig['MAG_AUTO'],orig['FLUX_RADIUS'], '.', color='0.7',markersize=0.5)
    ax2.plot(gabe['MAG_AUTO'],gabe['FLUX_RADIUS']*0.1/0.06, 'o', color='red', markersize=1., alpha=0.5)
    ax2.set_xlabel('MAG_AUTO')
    ax2.set_ylabel('FLUX_RADIUS')
    ax2.set_yscale('log')
    ax2.set_title('PIPELINE FLTS (GABE)')
    ax2.yaxis.set_major_formatter(ScalarFormatter())
    
    plt.show(block=False)