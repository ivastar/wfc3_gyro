import matplotlib
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
import matplotlib.gridspec as gridspec
import pylab
import unicorn
import astropy
from astropy.table import Table as table
from astropy.io import fits
import numpy as np
import glob, os

def mag_radius():
    
    from matplotlib.ticker import ScalarFormatter
    import my_python.mk_region_file
    
    gabe = table.read('../REF/cosmos-wide_v0.3_drz_bkg_sci.cat', format='ascii.sextractor')
    orig = table.read('/3DHST/Photometry/Release/v4.0/COSMOS/Detection/cosmos_3dhst.v4.0.F160W_orig.cat', format='ascii.sextractor')
    cat = table.read('test7_drz_sci.cat', format='ascii.sextractor')
    
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
    
    ax1 = fig.add_subplot(gs1[0,0], ylim=[0.8, 20], xlim=[14,26])
    ax1.plot(orig['MAG_AUTO'],orig['FLUX_RADIUS']*0.06/0.1, '.', color='0.7',markersize=0.5)
    ax1.plot(cat['MAG_AUTO'],cat['FLUX_RADIUS'], 'o', color='black', markersize=1, alpha=0.5)
    ax1.set_xlabel('MAG_AUTO F160W')
    ax1.set_ylabel('FLUX_RADIUS')
    ax1.set_yscale('log')
    cr = (cat['FLUX_RADIUS']*0.1/0.06 < 2.)
    stars = (cat['MAG_AUTO'] > 15.) & (cat['MAG_AUTO'] < 22.) & (cat['FLUX_APER_5']/cat['FLUX_APER'] > 1.1) & (cat['FLUX_APER_5']/cat['FLUX_APER'] < 1.2)
    #stars = (cat['MAG_AUTO'] > 15.) & (cat['MAG_AUTO'] < 23.) & (cat['FLUX_RADIUS']*0.1/0.06 < 5.15 - 0.115*cat['MAG_AUTO']) & (cat['FLUX_RADIUS']*0.1/0.06 > 2.1)
    
    #ax1.plot(cat['MAG_AUTO'][cr],cat['FLUX_RADIUS'][cr], 'o', color='blue', markersize=2, alpha=1.0)
    ax1.plot(cat['MAG_AUTO'][stars],cat['FLUX_RADIUS'][stars], 'o', color='red', markersize=2, alpha=1.0, markeredgecolor='red')
    print 'STARS: mean: {} / median{}'.format(np.mean(cat['FWHM_IMAGE'][stars]), np.median(cat['FWHM_IMAGE'][stars]))
    
    #yy_select = -0.115*np.array([14.,23.]) +5.15
    #ax1.plot(np.array([14.,23.]), yy_select,'--',color='red')
    #ax1.plot([23.,23.], np.array([2.1, yy_select[-1]]), '--',color='red')
    #ax1.plot([14.,23.], [2.1,2.1],'--',color='red')
    
    ax1.yaxis.set_major_formatter(ScalarFormatter())
    
    ax2 = fig.add_subplot(gs1[0,1], ylim=[0.5,4], xlim=[14,26])
    #ax2.plot(orig['MAG_AUTO'],orig['FLUX_RADIUS'], '.', color='0.7',markersize=0.5)
    #ax2.plot(gabe['MAG_AUTO'],gabe['FLUX_RADIUS']*0.1/0.06, 'o', color='red', markersize=1., alpha=0.5)
    ax2.plot(orig['MAG_AUTO'], orig['FLUX_APER_5']/orig['FLUX_APER'], '.', color='0.7',markersize=0.5)
    ax2.plot(cat['MAG_AUTO'], cat['FLUX_APER_5']/cat['FLUX_APER'], 'o', color='black', markersize=1., alpha=0.5)
    ax2.plot(cat['MAG_AUTO'][stars], cat['FLUX_APER_5'][stars]/cat['FLUX_APER'][stars], 'o', color='red', markersize=2., alpha=1.0, markeredgecolor='red')
    ax2.plot(cat['MAG_AUTO'][cr], cat['FLUX_APER_5'][cr]/cat['FLUX_APER'][cr], 'o', color='blue', markersize=2., alpha=1.0)
    ax2.set_xlabel('MAG_AUTO F160W')
    ax2.set_ylabel('Flux (2.0\")/Flux (0.5\")')
    #ax2.set_title('PIPELINE FLTS (GABE)')
    ax2.yaxis.set_major_formatter(ScalarFormatter())
    
    my_python.mk_region_file.mk_region_file_from_lists(cat['X_WORLD'][stars],cat['Y_WORLD'][stars],outfile = 'stars', printids='no', color='cyan')
    my_python.mk_region_file.mk_region_file_from_lists(cat['X_WORLD'][cr],cat['Y_WORLD'][cr],outfile = 'cr', printids='no', color='yellow')
    
    plt.show(block=False)
    
    fig.savefig('mag_radius.pdf', dpi=200, transparent=False)
    fig.savefig('mag_radius.png', dpi=200, transparent=False)
    
def noise_distribution(master_root='icxe15010'):
        
    import threedhst
    import stwcs
    import scipy
    import scipy.optimize
    import unicorn
    import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator
    import matplotlib.gridspec as gridspec
    import drizzlepac
    from drizzlepac import astrodrizzle
        
    master_asn = threedhst.utils.ASNFile('{}_asn.fits'.format(master_root))
    
    #print 'Read files...'
    ref = fits.open('{}_drz_sci.fits'.format(master_root))
    ref_wcs = stwcs.wcsutil.HSTWCS(ref, ext=0)

    seg = fits.open('{}_drz_seg.fits'.format(master_root))    
    seg_data = np.cast[np.float32](seg[0].data)
          
    yi, xi = np.indices((1014,1014))
    
    # read candels
    candels = table.read('candels_noise.txt', format='ascii')
   
    # Make plot
    fs=8.
    matplotlib.rc('xtick',labelsize=fs)
    matplotlib.rc('ytick',labelsize=fs)
    fig = unicorn.catalogs.plot_init(xs=10., aspect=0.5, fontsize=fs, left=0.1, right=0.1, top=0.1, bottom=0.1)
    gs = gridspec.GridSpec(2,4,top=0.95, bottom=0.1)
    fig.subplots_adjust(wspace=0.05)
    fig.subplots_adjust(hspace=0.1)
    
    fp = open('sigma_table.txt','a')

    #### Loop through the pointings in this orbit
    for j, root in enumerate(master_asn.exposures):
        
        asn = threedhst.utils.ASNFile('{}_asn.fits'.format(root))
        ax = fig.add_subplot(gs[j])
        
        if j == 0:
            flt_color='blue'
        else:
            flt_color='red'
          
        #### Loop through FLTs
        for exp in asn.exposures:
        
            NSAMP = len(asn.exposures)
        
            flt = fits.open('{}_flt.fits'.format(exp))
            flt_wcs = stwcs.wcsutil.HSTWCS(flt, ext=1)
         
            if flt[1].header['BG_SUB'] ==  'No':
                flt[1].data -= flt[1].header['MDRIZSKY']          
        
            if exp == asn.exposures[0]:
                print 'Segmentation image: {}_blot.fits'.format(exp)
                blotted_seg = astrodrizzle.ablot.do_blot(seg_data, ref_wcs, flt_wcs, 1, coeffs=True, interp='nearest', sinscl=1.0, stepsize=10, wcsmap=None)         
            
            mask = (blotted_seg == 0) & (flt['DQ'].data == 0) & (xi > 10) & (yi > 10) & (xi < 1004) & (yi < 1004)
        
            n, bin_edge, patch  = ax.hist(flt[1].data[mask], bins=300, range=[-3,3], color='0.5', alpha=0.5, histtype='step', normed=True)
        
            if exp == asn.exposures[0]:
                centers = (bin_edge[:-1] + (bin_edge[1:] - bin_edge[:-1])/2)
                flt_coeff, flt_var = scipy.optimize.curve_fit(gauss, centers, n, p0=[np.max(n),0.,1.])
            
        drz = fits.open('{}_drz_sci.fits'.format(root))
        drz_wcs = stwcs.wcsutil.HSTWCS(drz, ext=0)
        drz_sh = np.shape(drz[0].data)
        dyi, dxi = np.indices(drz_sh)
        blotted_seg = astrodrizzle.ablot.do_blot(seg_data, ref_wcs, drz_wcs, 1, coeffs=False, interp='nearest', sinscl=1.0, stepsize=10, wcsmap=None)
        mask = (blotted_seg == 0) & (dxi > 10) & (dyi > 10) & (dxi < drz_sh[1]-10) & (dyi < drz_sh[0]-10)
        n, bin_edge, patch = ax.hist(drz[0].data[mask], bins=100, range=[-1.,1.], color='black', alpha=1.0, histtype='step', normed=True)
        centers = (bin_edge[:-1] + (bin_edge[1:] - bin_edge[:-1])/2)
        coeff, var = scipy.optimize.curve_fit(gauss, centers, n, p0=[np.max(n),0.,1.])

        flt = fits.open('{}_flt.fits'.format(root))
        flt_wcs = stwcs.wcsutil.HSTWCS(flt, ext=1)
        blotted_seg = astrodrizzle.ablot.do_blot(seg_data, ref_wcs, flt_wcs, 1, coeffs=True, interp='nearest', sinscl=1.0, stepsize=10, wcsmap=None)
        mask = (blotted_seg == 0) & (flt['DQ'].data == 0) & (xi > 10) & (yi > 10) & (xi < 1004) & (yi < 1004)
        n, bin_edge, patch  = ax.hist(flt[1].data[mask], bins=100, range=[-1.,1.], color=flt_color, alpha=1.0, histtype='step', normed=True)
        centers = (bin_edge[:-1] + (bin_edge[1:] - bin_edge[:-1])/2)
        flt_orig_coeff, flt_orig_var = scipy.optimize.curve_fit(gauss, centers, n, p0=[np.max(n),0.,1.])

        if root == master_asn.exposures[0]:
            flag = 0
        else:
            flag = 1

        fp.write('{}\t{:7.4f}\t{:7.4f}\t{:7.2f}\t{:7.4f}\t{}\n'.format(root, np.abs(flt_orig_coeff[2]), np.abs(coeff[2]), flt[1].header['SAMPTIME'], flt[1].header['MDRIZSKY'],flag))

        ### plot candels
        #ax.bar(candels['centers'], 4.5*candels['n']/np.max(candels['n']), width=0.001, align='center', edgecolor='orange', color='orange', alpha=0.3)
    
        ax.set_ylim([0, 5.0])
        ax.set_xlim([-3.0, 3.0])
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        if j == 0 or j==4:
            ax.set_ylabel('Normalized N', fontsize=fs)
        else:
            ax.set_yticklabels([])
        if j > 3:    
            ax.set_xlabel('e$^-$/s', fontsize=fs)
        else:
            ax.set_xticklabels([])
        ax.set_title(root, fontsize=fs)
        ax.text(-2.5, 4.5,'$\sigma_{reads}$ = %4.3f' % (flt_coeff[2]), color='0.5', fontsize=fs)
        ax.text(-2.5, 4.,'$\sigma_{FLT}$ = %4.3f' % (flt_orig_coeff[2]), color=flt_color, fontsize=fs)    
        ax.text(-2.5, 3.5,'$\sigma_{DRZ}$ = %4.3f' % (coeff[2]), color='black', fontsize=fs)
        #ax.text(-2.5, 3.0,'$\sigma_{CAND}$ = %4.3f' % (4.44509564e-03), color='orange', fontsize=fs)
    
    plt.show(block=False)
    
    fig.savefig('{}_noise.png'.format(master_root), dpi=200, transparent=False)
    
def candels_mosaic_noise():
    
    import scipy
    
    mos = fits.open('/3DHST/Photometry/Release/v4.0/COSMOS/HST_Images/cosmos_3dhst.v4.0.F160W_orig_sci.fits')
    wht = fits.open('/3DHST/Photometry/Release/v4.0/COSMOS/HST_Images/cosmos_3dhst.v4.0.F160W_orig_wht.fits')
    seg = fits.open('/3DHST/Photometry/Release/v4.0/COSMOS/Detection/cosmos_3dhst.v4.0.F160W_seg.fits.gz')
    
    mask = (seg[0].data == 0) & (wht[0].data != 0)
    n, bin_edge = np.histogram(mos[0].data[mask], bins=100, range=[-0.05,0.05], normed=True)
    centers = (bin_edge[:-1] + (bin_edge[1:] - bin_edge[:-1])/2)
    coeff, var = scipy.optimize.curve_fit(gauss, centers, n, p0=[np.max(n),0.,0.01])
    print 'Gaussian coeffs: {}'.format(coeff)
    
    t = table([centers,n], names =('centers','n'))
    t.write('candels_noise.txt', format='ascii')
    
def overall_offsets():
    
    import glob
    import matplotlib.colors as colors
    import matplotlib.cm as cmx
    
    files = glob.glob('shifts_icxe*010.txt')
    
    fig = unicorn.catalogs.plot_init(xs=10.,aspect=0.24, left=0.12, right=0.05, bottom=0.15, top=0.15, NO_GUI=False)    
    gs = gridspec.GridSpec(1,4,top=0.9, bottom=0.15)
    fig.subplots_adjust(wspace=0.17)
    fig.subplots_adjust(hspace=0.1)
    fs = 8

    jet = cm = plt.get_cmap('jet_r')
    cNorm = colors.Normalize(vmin=0, vmax=9)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
    
    for j, file in enumerate(files):
        print file, '{}'.format(file.split('_')[1].split('.')[0])
        ax = fig.add_subplot(gs[j])
        if file.startswith('shifts_icxe15'):
            ax.set_xlim([-5.,70.])
            ax.set_ylim([-5.,20.])
        if file.startswith('shifts_icxe16'):
            ax.set_xlim([-5.,50.])
            ax.set_ylim([-10.,5.])
        if file.startswith('shifts_icxe17'):
            ax.set_xlim([-50.,5.])
            ax.set_ylim([-10.,5.])
        if file.startswith('shifts_icxe18'):
            ax.set_xlim([-40.,5.])
            ax.set_ylim([-5.,40.])
        ax.set_xlabel('$\Delta$ x [pix]', fontsize=fs, labelpad=0.1)
        if j == 0:
            ax.set_ylabel('$\Delta$ y [pix]', fontsize=fs)
        ax.set_title('{}'.format(file.split('_')[1].split('.')[0]), fontsize=fs)
            
        cc = 0.
        with open(file) as f:
            for line in f:
                if not line.startswith('#'):
                    data = line.split()
                    cc += 1.
                    color = scalarMap.to_rgba(cc)
                    #ax.arrow(0.,0., float(data[1]), float(data[2]), head_width=0., head_length=0., fc=color, ec = color)
                    ax.annotate("",xy=(float(data[1]), float(data[2])), xytext=(0,0), 
                        arrowprops=dict(arrowstyle='->', color=color))
        
    plt.show(block=False)
    fig.savefig('overall_offsets.png', dpi=200, transparent=False)


def overall_offsets_vs_distance():
    
    import glob
    import matplotlib.colors as colors
    import matplotlib.cm as cmx
    
    files = glob.glob('shifts_icxe*010.txt')
    
    fig = unicorn.catalogs.plot_init(xs=10.,aspect=0.47, left=0.2, right=0.05, bottom=0.15, top=0.15, NO_GUI=False)    
    gs = gridspec.GridSpec(1,2,top=0.9, bottom=0.15)
    fig.subplots_adjust(wspace=0.2)
    fig.subplots_adjust(hspace=0.2)
    fs = 10
    matplotlib.rc('xtick',labelsize=fs)
    matplotlib.rc('ytick',labelsize=fs)

    jet = cm = plt.get_cmap('jet_r')
    cNorm = colors.Normalize(vmin=0, vmax=8)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)


    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])
    lines = ['-','--','-.',':']
    labels = ['COSMOS-15', 'COSMOS-16', 'COSMOS-17','COSMOS-18']
    
    for j, file in enumerate(files):
        print file
        root = '{}'.format(file.split('_')[1].split('.')[0])
        ax1.set_xlim([-0.5,9.])
        ax1.set_ylim([-0.005,0.15])
        ax1.set_xlabel('Commanded Offset Relative to First Position\n[arcmin]', fontsize=fs, labelpad=0.1)
        ax1.set_ylabel('Offset from Commanded Position\n[arcmin]', fontsize=fs)
        root = '{}'.format(file.split('_')[1].split('.')[0])
            
        data = table.read(file, format='ascii',  names=('file','x','y','rot','scale','x_rms','y_rms'))    
        small_off  = np.sqrt(data['x']**2 + data['y']**2)*0.12/60.
        #colors = (scalarMap.to_rgba(cc) for cc in range(len(small_off)))
        big_off = [0.0]
        drift_rate = [0.0]
        
        origin = fits.open(data['file'][0])
        x_orig, y_orig = origin[1].header['CRVAL1'], origin[1].header['CRVAL2'] 
        
        for k, file in enumerate(data['file'][1:]):
            flt = fits.open(file)
            x, y = flt[1].header['CRVAL1'], flt[1].header['CRVAL2']
            big_off.append(60.*np.sqrt((x_orig-x)**2 + (y_orig-y)**2))
            
            if k+1 < 4:
                t_exp = 255
            else:
                t_exp = 277
            drift = table.read(file.split('_')[0]+'_shifts.txt', format='ascii', 
                names=('file','x','y','rot','scale','x_rms','y_rms'))
            drift_rate.append(25.*np.sqrt((drift['x'][0]-drift['x'][-1])**2 + (drift['y'][0]-drift['y'][-1])**2)/t_exp)
            
        print len(small_off), len(big_off)
        ax1.plot(big_off, small_off,linestyle=lines[j], color='0.5', label=labels[j], zorder=0)
        ax1.scatter(big_off, small_off, c=range(len(small_off)), cmap='jet', s=35., edgecolors='black', alpha=0.7)
        ax1.legend(loc='upper left', frameon=False, labelspacing=0.8, fontsize=9)        
        
        ax2.set_xlim([-0.5, 9.])
        ax2.set_ylim([-0.01,0.5])
        ax2.set_xlabel('Commanded Offset Relative to First Position\n[arcmin]', fontsize=fs, labelpad=0.1)
        ax2.set_ylabel('Drift Rate During Observations\n[pix per 25 seconds]', fontsize=fs)
        
        ax2.plot(big_off, drift_rate,linestyle=lines[j], color='0.5', zorder=0)
        ax2.scatter(big_off, drift_rate, c=range(len(small_off)), cmap='jet', s=35., edgecolors='black', alpha=0.7)
        
        
        
    plt.show(block=False)
    fig.savefig('overall_offsets_vs_distance.png', dpi=200, transparent=False)


def overall_offsets_vs_time():
    
    import glob
    import matplotlib.colors as colors
    import matplotlib.cm as cmx
    
    files = glob.glob('shifts_icxe*010.txt')
    
    fig = unicorn.catalogs.plot_init(xs=10.,aspect=0.47, left=0.2, right=0.05, bottom=0.15, top=0.15, NO_GUI=False)    
    gs = gridspec.GridSpec(1,2,top=0.9, bottom=0.15)
    fig.subplots_adjust(wspace=0.2)
    fig.subplots_adjust(hspace=0.2)
    fs = 10
    matplotlib.rc('xtick',labelsize=fs)
    matplotlib.rc('ytick',labelsize=fs)

    jet = cm = plt.get_cmap('jet_r')
    cNorm = colors.Normalize(vmin=0, vmax=8)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)


    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])

    time = np.array([0, 368, 713, 1060, 1405, 1777, 2147, 2519]) + 333.
    lines = ['-','--','-.',':']
    labels = ['COSMOS-15', 'COSMOS-16', 'COSMOS-17','COSMOS-18']
    
    for j, file in enumerate(files):
        root = '{}'.format(file.split('_')[1].split('.')[0])
        ax1.set_xlim([-5., 3100.])
        ax1.set_ylim([-0.005,0.15*60.])
        ax1.set_xlabel('Time Since Beginning of Orbit\n[seconds]', fontsize=fs)
        ax1.set_ylabel('Offset from Commanded Position\n[arcsec]', fontsize=fs)
        root = '{}'.format(file.split('_')[1].split('.')[0])
            
        data = table.read(file, format='ascii',  names=('file','x','y','rot','scale','x_rms','y_rms'))    
        small_off  = np.sqrt(data['x']**2 + data['y']**2)*0.12
        #colors = (scalarMap.to_rgba(cc) for cc in range(len(small_off)))
        drift_rate = [0.0]
                
        for k, file in enumerate(data['file'][1:]):
            if k+1 < 4:
                t_exp = 255
            else:
                t_exp = 277
            drift = table.read(file.split('_')[0]+'_shifts.txt', format='ascii', 
                names=('file','x','y','rot','scale','x_rms','y_rms'))
            drift_rate.append(25.*np.sqrt((drift['x'][0]-drift['x'][-1])**2 + (drift['y'][0]-drift['y'][-1])**2)/t_exp)
            
        ax1.plot(time, small_off, linestyle=lines[j], color='0.5', label=labels[j], zorder=0)
        ax1.scatter(time, small_off, c=range(len(small_off)), cmap='jet', s=35., edgecolors='black', alpha=0.7)
        ax1.legend(loc='upper left', frameon=False, labelspacing=0.8, fontsize=9)
        
        
        ax2.set_xlim([-5., 3100.])
        ax2.set_ylim([-0.01,0.5])
        ax2.set_xlabel('Time Since Beginning of Orbit\n[seconds]', fontsize=fs)
        ax2.set_ylabel('Drift Rate During Pointing\n[pix per 25 seconds]', fontsize=fs)
        
        ax2.plot(time, drift_rate,linestyle=lines[j], color='0.5', zorder=0)
        ax2.scatter(time, drift_rate, c=range(len(small_off)), cmap='jet', s=35., edgecolors='black', alpha=0.7)
        
    plt.show(block=False)
    fig.savefig('overall_offsets_vs_time.png', dpi=200, transparent=False)
    fig.savefig('overall_offsets_vs_time.pdf', dpi=200, transparent=False)

        
def gyro_drift():
    
    import glob
    import matplotlib.colors as colors
    import matplotlib.cm as cmx
    import threedhst
    
    files = glob.glob('icxe*010_asn.fits')
    
    fig = unicorn.catalogs.plot_init(xs=10.,aspect=0.24, left=0.12, right=0.05, bottom=0.15, top=0.15, NO_GUI=False)    
    gs = gridspec.GridSpec(1,4,top=0.9, bottom=0.15)
    fig.subplots_adjust(wspace=0.17)
    fig.subplots_adjust(hspace=0.1)
    fs = 8

    jet = cm = plt.get_cmap('jet_r')
    cNorm = colors.Normalize(vmin=0, vmax=9)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
    
    for j, file in enumerate(files):
        print file
        ax = fig.add_subplot(gs[j])
        ax.set_xlabel('$\Delta$ x [pix]', fontsize=fs, labelpad=0.1)
        if j == 0:
            ax.set_ylabel('$\Delta$ y [pix]', fontsize=fs)
        ax.set_title('{}'.format(file.split('_')[0]), fontsize=fs)
        ax.set_xlim([-3., 3.])
        ax.set_ylim([-3., 3.])
        
        asn = threedhst.utils.ASNFile(file)
        cc = 1.
        for exp in asn.exposures[1:]:
            
            data = table.read(exp+'_shifts.txt', format='ascii', names=('file','x','y','rot','scale','x_rms','y_rms'))
            cc += 1.
            color = scalarMap.to_rgba(cc)
            ax.plot(data['x'], data['y'], '-', color=color, alpha=0.5, linewidth=1.5)
            ax.plot(data['x'], data['y'], 'o', color=color, markersize=3., alpha=0.5, markeredgecolor=None)
            
    plt.show(block=False)
    fig.savefig('gyro_drift.png', dpi=200, transparent=False)                    

def footprints_plot(root='icxe15010'):
    
    import unicorn.survey_paper as sup
    import matplotlib.colors as colors
    import matplotlib.cm as cmx
    
    if root == 'icxe15010':
        aspect = 1.75
        xlim = [150.265, 150.157]
        ylim = [2.45, 2.64]
        xticklab = [r'$10^\mathrm{h}01^\mathrm{m}00^\mathrm{s}$', r'$10^\mathrm{h}00^\mathrm{m}45^\mathrm{s}$']
        xtickv = [sup.degrees(10,01,00, hours=True),sup.degrees(10,00,45, hours=True)]
        yticklab = [r'$+02^\circ30^\prime00^{\prime\prime}$',r'$+02^\circ35^\prime00^{\prime\prime}$']
        ytickv = [sup.degrees(2, 30, 00, hours=False),sup.degrees(2, 35, 00, hours=False)]
        label = 'COSMOS-15'
        factor=10.

    if root == 'icxe16010':
        aspect=0.9
        xlim = [150.265, 150.1]
        ylim = [2.607, 2.74]
        xticklab = [r'$10^\mathrm{h}01^\mathrm{m}00^\mathrm{s}$', r'$10^\mathrm{h}00^\mathrm{m}45^\mathrm{s}$',r'$10^\mathrm{h}00^\mathrm{m}30^\mathrm{s}$']
        xtickv = [sup.degrees(10,01,00, hours=True),sup.degrees(10,00,45, hours=True),sup.degrees(10,00,30, hours=True)]
        yticklab = [r'$+02^\circ38^\prime00^{\prime\prime}$',r'$+02^\circ40^\prime00^{\prime\prime}$', r'$+02^\circ42^\prime00^{\prime\prime}$', r'$+02^\circ44^\prime00^{\prime\prime}$']
        ytickv = [sup.degrees(2, 38, 00, hours=False),sup.degrees(2, 40, 00, hours=False),sup.degrees(2, 42, 00, hours=False),sup.degrees(2, 44, 00, hours=False)]
        label='COSMOS-16'
        factor=20.
    
    if root == 'icxe17010':
        aspect=1.4
        xlim = [150.2, 150.06]
        ylim = [2.52, 2.72]
        xticklab = [r'$10^\mathrm{h}00^\mathrm{m}45^\mathrm{s}$', r'$10^\mathrm{h}00^\mathrm{m}30^\mathrm{s}$',r'$10^\mathrm{h}00^\mathrm{m}15^\mathrm{s}$']
        xtickv = [sup.degrees(10,00,45, hours=True),sup.degrees(10,00,30, hours=True),sup.degrees(10,00,15, hours=True)]
        yticklab = [r'$+02^\circ35^\prime00^{\prime\prime}$',r'$+02^\circ40^\prime00^{\prime\prime}$']
        ytickv = [sup.degrees(2, 35, 00, hours=False),sup.degrees(2, 40, 00, hours=False)]
        label='COSMOS-17'
        factor=240.

    if root == 'icxe18010':
        aspect=1.577
        xlim = [150.14, 150.01]
        ylim = [2.53, 2.735]
        xticklab = [r'$10^\mathrm{h}00^\mathrm{m}30^\mathrm{s}$', r'$10^\mathrm{h}00^\mathrm{m}20^\mathrm{s}$',r'$10^\mathrm{h}00^\mathrm{m}10^\mathrm{s}$']
        xtickv = [sup.degrees(10,00,30, hours=True),sup.degrees(10,00,20, hours=True),sup.degrees(10,00,10, hours=True)]
        yticklab = [r'$+02^\circ35^\prime00^{\prime\prime}$',r'$+02^\circ40^\prime00^{\prime\prime}$']
        ytickv = [sup.degrees(2, 35, 00, hours=False),sup.degrees(2, 40, 00, hours=False)]
        label='COSMOS-18'
        factor=240.
    
    
    
    fig = unicorn.catalogs.plot_init(square=True, xs=5., aspect=aspect, 
        fontsize=8, left=0.18, right=0.02, bottom=0.10, top=0.10)
    ax = fig.add_subplot(111)
    jet = cm = plt.get_cmap('jet')
    cNorm = colors.Normalize(vmin=0, vmax=9)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)    
    
    reg_file = root+'_asn.reg'
    
    poly = []
    with open(reg_file) as f:
        for line in f:
            if not line.startswith('fk5'):
                region = line.split('#')[0]
                poly.append(sup.polysplit(region=region, get_shapely=True))

    shifts = table.read('shifts_{}.txt'.format(root), format='ascii', 
        names=('file','x','y','rot','scale','x_rms','y_rms'))
        
    cc = 0
    xcen_all = []
    ycen_all = []
    for j,(pp, x_off, y_off, file) in enumerate(zip(poly, shifts['x'], shifts['y'], shifts['file'])):
        cc += 1.
        color = scalarMap.to_rgba(cc)
        x, y = pp.exterior.xy
        flt = fits.open(file)
        xcen = flt[1].header['CRVAL1O']
        ycen = flt[1].header['CRVAL2O']
        x_off = (flt[1].header['CRVAL1B']-flt[1].header['CRVAL1O'])*20.
        y_off = (flt[1].header['CRVAL2B']-flt[1].header['CRVAL2O'])*20.
        print file, xcen, xcen+x_off, ycen, ycen+y_off
        #xcen = (np.mean(x[:-1]))
        #ycen = (np.mean(y[:-1]))
        xcen_all.append(xcen)
        ycen_all.append(ycen)
        ax.plot(x,y,'-', color=color)
        #ax.annotate("",xy=(xcen+(x_off*0.12)/factor, ycen+(y_off*0.12)/factor), xytext=(xcen, ycen), 
        #    arrowprops=dict(arrowstyle='->', color=color))
        #ax.plot([xcen, xcen+x_off], [ycen, ycen+y_off], '-')
        ax.annotate("",xy=(xcen+x_off, ycen+y_off), xytext=(xcen, ycen), 
            arrowprops=dict(arrowstyle='->', color=color))

    ax.plot(xcen_all, ycen_all, '+:', markersize=10., color='0.5', alpha=0.5) 
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)     
    ax.set_xticklabels(xticklab)
    xtick = ax.set_xticks(xtickv)
    ax.set_yticklabels(yticklab)
    ytick = ax.set_yticks(ytickv)
    ax.set_title(label)       
    plt.show(block=False)
    
    fig.savefig('footprint_{}.png'.format(label.lower()), dpi=200, transparent=False)          
    
def footprints_all():
    
    import unicorn.survey_paper as sup
    import matplotlib.colors as colors
    import matplotlib.cm as cmx
    
    fig = unicorn.catalogs.plot_init(square=True, xs=10., aspect=1.1, 
        fontsize=8, left=0.18, right=0.02, bottom=0.10, top=0.10)
    ax = fig.add_subplot(111)
    jet = cm = plt.get_cmap('jet')
    cNorm = colors.Normalize(vmin=0, vmax=9)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)    
    
    xlim = [150.270, 149.99]
    ylim = [2.45, 2.751]
    xticklab = [r'$10^\mathrm{h}01^\mathrm{m}00^\mathrm{s}$', r'$10^\mathrm{h}00^\mathrm{m}45^\mathrm{s}$',r'$10^\mathrm{h}00^\mathrm{m}30^\mathrm{s}$',r'$10^\mathrm{h}00^\mathrm{m}15^\mathrm{s}$',r'$10^\mathrm{h}00^\mathrm{m}00^\mathrm{s}$']
    xtickv = [sup.degrees(10,01,00, hours=True),sup.degrees(10,00,45, hours=True),sup.degrees(10,00,30, hours=True),sup.degrees(10,00,15, hours=True),sup.degrees(10,00,00, hours=True)]
    yticklab = [r'$+02^\circ30^\prime00^{\prime\prime}$',r'$+02^\circ35^\prime00^{\prime\prime}$', r'$+02^\circ40^\prime00^{\prime\prime}$',r'$+02^\circ45^\prime00^{\prime\prime}$']
    ytickv = [sup.degrees(2, 30, 00, hours=False),sup.degrees(2, 35, 00, hours=False),sup.degrees(2, 40, 00, hours=False),sup.degrees(2, 45, 00, hours=False)]
    factor=20.
    
    roots = ['icxe15010', 'icxe16010', 'icxe17010','icxe18010']
    labels = ['COSMOS-15', 'COSMOS-16','COSMOS-17','COSMOS-18']
    lines = ['-','--','-.',':']
    
    for root, label, linestyle in zip(roots, labels, lines):
        
        reg_file = root+'_asn.reg'
    
        poly = []
        with open(reg_file) as f:
            for line in f:
                if not line.startswith('fk5'):
                    region = line.split('#')[0]
                    poly.append(sup.polysplit(region=region, get_shapely=True))

        shifts = table.read('shifts_{}.txt'.format(root), format='ascii', 
            names=('file','x','y','rot','scale','x_rms','y_rms'))
        
        cc = 0
        xcen_all = []
        ycen_all = []
        for j,(pp, file) in enumerate(zip(poly, shifts['file'])):
            cc += 1.
            color = scalarMap.to_rgba(cc)
            x, y = pp.exterior.xy
            flt = fits.open(file)
            xcen = flt[1].header['CRVAL1O']
            ycen = flt[1].header['CRVAL2O']
            x_off = (flt[1].header['CRVAL1B']-flt[1].header['CRVAL1O'])*factor
            y_off = (flt[1].header['CRVAL2B']-flt[1].header['CRVAL2O'])*factor
            xcen_all.append(xcen)
            ycen_all.append(ycen)
            ax.plot(x,y,'-', color=color)
            ax.annotate("",xy=(xcen+x_off, ycen+y_off), xytext=(xcen, ycen), 
                arrowprops=dict(arrowstyle='->', color=color))
            ax.text(xcen, ycen+0.005, file.split('_')[0], fontsize=9, va='top', ha='center')
        
        ax.plot(xcen_all, ycen_all, linestyle=linestyle, markersize=10., 
            color='0.5', alpha=0.7, label=label) 
        ax.plot(xcen_all, ycen_all, '+', markersize=10., color='0.5', alpha=0.7) 
        ax.text(xcen_all[0], ycen_all[0]-0.005, label, fontsize=10, va='top', ha='center')

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)     
    ax.set_xticklabels(xticklab)
    xtick = ax.set_xticks(xtickv)
    ax.set_yticklabels(yticklab)
    ytick = ax.set_yticks(ytickv)
    ax.legend(loc='lower right', frameon=False, labelspacing=0.8, fontsize=9, handlelength=10, borderpad=5.)
    #ax.set_title(label)       
    plt.show(block=False)
        
    fig.savefig('footprints_all.png'.format(label.lower()), dpi=200, transparent=False)          
    
def make_hist_plot_bw():
    
    import astropy.io.ascii as ascii
    from matplotlib.ticker import MultipleLocator, FormatStrFormatter,AutoMinorLocator, ScalarFormatter

    minorLocator   = MultipleLocator(10)
    
    matplotlib.rc('xtick',labelsize=10)
    matplotlib.rc('ytick',labelsize=10)
    
    data = ascii.read('/3DHST/Spectra/Work/HOPR/GYRO/gyro_shifts.txt')
        
    
    fig = unicorn.catalogs.plot_init(square=True, xs=6., aspect=0.9, fontsize=10., left=0.15, right=0.05, top=0.1, bottom=0.15)
    ax = fig.add_subplot(111)

    ax.fill_between([0.,0.5], 0, 30, color='0.85', alpha=0.9)
    nn, bins, patches = ax.hist(data['rate']*25., bins=60, range=[0.,6.], color='black', histtype='step', linestyle='dotted')    

    test_val = data['rate'][data['root'] == 'ib2u22prq']
        
    ax.arrow(test_val*25., 17.1, 0., -1., head_length=0.5, head_width=0.1, fc='red', ec='red')
    ax.text(0., 18.0, 'Fig. 1 Example', ha='left', va='center', fontsize=13, fontweight='heavy')
    
    ax.arrow(0., 21., 0.5, 0., head_length=0.05, head_width=0.30, fc='black', ec='black', length_includes_head=True)
    ax.arrow(0.5, 21., -0.5, 0., head_length=0.05, head_width=0.30, fc='black', ec='black', length_includes_head=True)
    ax.text(0.0, 22.5, 'Expected Drift', ha='left', va='center', fontsize=15, fontweight='heavy')
    ax.text(1.1, 5, 'gyro\nproblem', fontsize=12.5, fontweight=0.1, multialignment='center')
    ax.text(4.3, 5, 'failed\nguide star\nacquisition', fontsize=12.5, fontweight=0.1, multialignment='center')
    
    shift_files = glob.glob('shifts_icxe*010.txt')
    drift_rate = []
    for shift_file in shift_files:
        data = table.read(shift_file, format='ascii',  names=('file','x','y','rot','scale','x_rms','y_rms'))    
                
        for k, file in enumerate(data['file'][1:]):
            if k+1 < 4:
                t_exp = 255
            else:
                t_exp = 277
            drift = table.read(file.split('_')[0]+'_shifts.txt', format='ascii', 
                names=('file','x','y','rot','scale','x_rms','y_rms'))
            drift_rate.append(25.*np.sqrt((drift['x'][0]-drift['x'][-1])**2 + (drift['y'][0]-drift['y'][-1])**2)/t_exp)
    
    nn_sh, bins_sh, patches_sh = ax.hist(drift_rate, bins=5, range=[0.,0.5], color='red', histtype='step', lw=3) 
    print "Mean drift: {}\nMedian drift: {}".format(np.mean(drift_rate), np.median(drift_rate))   

    ax.set_xlim([-0.15, 2.7])
    ax.set_ylim([0., 25])
    ax.set_xlabel('Drift in WFC3/IR Pixels per 25 Seconds', fontsize=12)
    ax.set_ylabel('Number of Images')
    ax.tick_params(axis='both', which='major',labelsize=12)
    minorLocator   = AutoMinorLocator()
    ax.xaxis.set_minor_locator(minorLocator)
    minorLocator   = AutoMinorLocator(5)
    ax.yaxis.set_minor_locator(minorLocator)
    ax.tick_params(which='minor', length=2)   
    plt.show(block=False)
    
    plt.savefig('drift_hist_bw.png',dpi=100,transparent=False)
    plt.savefig('drift_hist_bw.pdf',dpi=100,transparent=False)

def xy_shifts():
    import scipy
    
    # x gaussian: [ 3.41423523  1.29066909  6.73498474] or uniform [-70, 70np.]
    # y gaussian: [ 5.27724105 -0.78413847 -4.80223471]
    
    tab = table.read('allshifts.txt', format='ascii')

    nn, bins, patches = plt.hist(tab['col2'], bins=50, range=[-70,70], color='red', alpha=0.5)
    centers = 0.5*(bins[1:]+bins[:-1])
    coeff_g, var_g = scipy.optimize.curve_fit(gauss, centers , nn, p0=[np.max(nn),0.,100.])
    fit_g = gauss(centers, *coeff_g)
    plt.plot(centers, fit_g, color='red')
    print coeff_g, np.min()
    
    nn1, bins1, patches = plt.hist(tab['col3'], bins=bins, range=[-70,70], color='blue', alpha=0.5)
    centers1 = 0.5*(bins1[1:]+bins1[:-1]) 
    coeff_g, var_g = scipy.optimize.curve_fit(gauss, centers1 , nn1, p0=[np.max(nn1),0.,100.])
    fit_g = gauss(centers, *coeff_g)
    plt.plot(centers1, fit_g, color='blue')
    print coeff_g
    
    plt.show(block=False)

def gauss(x, *p):
    
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))
    

def area_mag():
    
    import unicorn
    import unicorn.survey_paper as sup
    from my_python.phot_paper import determine_outline as determine_outline
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MultipleLocator, FormatStrFormatter,AutoMinorLocator
    
    fig = unicorn.catalogs.plot_init(square=True, xs=5, aspect=1.0, fontsize=10, left=0.15, right=0.02, top=0.03, bottom=0.12)
    ax = fig.add_subplot(111)
       
    WFC3_IR_AREA = 4.6 # arcmin2
    imaging = {'UDF/XDF':(236.1e3/3000, 1*WFC3_IR_AREA, 150, 0.3), 'HUDF09-1': (13, 1*WFC3_IR_AREA,350,0.25), 'HUDF09-2':(19, 1*WFC3_IR_AREA, 350,0.25), 'CANDELS-Deep': (12, 0.04*3600, 350,0.07), 'CANDELS-Wide': (2./3., 0.2*3600, 0, 0.3), 'CLASH': (2, 25*WFC3_IR_AREA, 0,0.3), 'FF clusters': (24, 6*WFC3_IR_AREA, 150, 0.3)} 
    
    
    for key in imaging.keys():
        mag = 2.5*np.log10(np.sqrt(imaging[key][0]/0.125))+25.5
        ax.plot(imaging[key][1],mag, 's', markersize=7, color='black')
        ax.text(imaging[key][1]+imaging[key][2], mag+imaging[key][3], key, ha='center', va='top', fontsize=9)
        print key, mag
        if key == 'HUDF09-1':
            ax.plot([13+50, 13+130], [mag, mag+0.17], color='black', lw=1.)
        if key == 'HUDF09-2':
            ax.plot([19+50, 19+130], [mag, mag+0.17], color='black', lw=1.)
        if key == 'FF clusters':
            ax.plot()

    
    ax.plot([2000.], [25.0], 's', markersize=12, color='red')
    ax.text(2000., 25.3, 'GO-14114', multialignment='center', color='red',fontsize=10, ha='center', va='top')
    ax.text(1200, 28.5, 'HST/WFC3 F160W\n imaging', fontsize =12, multialignment='center')
    
    ax.set_xlim([-100, 2490.])
    ax.set_ylim([24.65, 29.5])
    ax.set_xlabel('Area [arcmin$^2$]', fontsize=12)
    ax.set_ylabel('5$\sigma$ Point Source Depth [AB mag]', fontsize=12)
    minorLocator   = AutoMinorLocator()
    ax.xaxis.set_minor_locator(minorLocator)
    minorLocator   = AutoMinorLocator(5)
    ax.yaxis.set_minor_locator(minorLocator)
    yticklab = ['25.0', '26.0', '27.0', '28.0', '29.0']
    ytickv = [25.0, 26.0, 27.0, 28.0, 29.0]
    ax.set_yticklabels(yticklab)
    ytick = ax.set_yticks(ytickv)
    
    plt.show(block=False)
    
    fig.savefig('area_mag.png', dpi=200, transparent=False)
    fig.savefig('area_mag.pdf', dpi=200, transparent=False)

def mosaic_demo():
    
    mosaic = fits.open('test7_drz_sci.fits')
    ### open COSMOS-15 mosaics
    cos15 = fits.open('icxe15010_test_drz_sci.fits')
    ### open CANDELS mosaic
    candels = fits.open('cosmos_3dhst_cutout.fits')  

    fig = unicorn.catalogs.plot_init(square=True, xs=10., aspect=0.835, 
        fontsize=8, left=0.02, right=0.02, bottom=0.02, top=0.02)

    gs1 = gridspec.GridSpec(1,1)
    gs2 = gridspec.GridSpec(2,1, left=0.65, right=0.98, top=0.8, bottom=0.4)
    gs3 = gridspec.GridSpec(1,2, left=0.40, right=0.88, top = 0.3, bottom=0.05)
    fig.subplots_adjust(wspace=0.05)
    fig.subplots_adjust(hspace=0.05)

    ax1 = fig.add_subplot(gs1[0,0], xlim=[0,12200], ylim=[0,10200])
    im1 = ax1.imshow(mosaic[0].data, cmap = pylab.cm.Greys, vmin=-0.01, 
        vmax=0.5, interpolation='None')
    ax1.plot([2150,3050,3050,2150,2150],[3250,3250,3850,3850,3250],'-',color='black')
    ax1.plot([1900,2125,2125,1900,1900],[5250,5250,5400,5400,5250],'-',color='black')
    ax1.tick_params(axis='both',which='both',top='off',bottom='off', right='off', left='off')
    ax1.axis('off')
    ax1.set_xticklabels([])
    ax1.set_yticklabels([])
     
    ### plot mosaic with boxes at the correct positions
    
    ### plot candels outline
    
    ### plot zooms of the unsmeared and smeared mosaics
    ax2 = fig.add_subplot(gs2[0,0], xlim=[0,225], ylim=[0,150])
    im2 = ax2.imshow(mosaic[0].data[5250:5400,1900:2125], cmap = pylab.cm.Greys, vmin=-0.01, 
        vmax=0.4, interpolation='bicubic')
    ax2.tick_params(axis='both',which='both',top='off',bottom='off', right='off', left='off')
    ax2.text(20,15,'COSMOS-DASH Mosaic', fontsize=10)
    ax2.set_xticklabels([])
    ax2.set_yticklabels([])

    ax3 = fig.add_subplot(gs2[1,0], xlim=[0,225], ylim=[0,150])
    im3 = ax3.imshow(cos15[0].data[5250:5400,1900:2125], cmap = pylab.cm.Greys, vmin=-0.01, 
        vmax=0.4, interpolation='bicubic')
    ax3.text(20,15,'Uncorrected (smeared) image', fontsize=10)    
    ax3.tick_params(axis='both',which='both',top='off',bottom='off', right='off', left='off')
    ax3.set_xticklabels([])
    ax3.set_yticklabels([])
    
    ### plot zooms of our mosaic and CANDELS
    ax4 = fig.add_subplot(gs3[0,0], xlim=[0,900], ylim=[0,600])
    im4 = ax4.imshow(mosaic[0].data[3250:3850,2150:3050], cmap = pylab.cm.Greys, vmin=-0.01, 
        vmax=0.5, interpolation='bicubic')
    ax4.text(50,50,'COSMOS-DASH Mosaic', fontsize=10)
    ax4.tick_params(axis='both',which='both',top='off',bottom='off', right='off', left='off')
    ax4.set_xticklabels([])
    ax4.set_yticklabels([])

    ax5 = fig.add_subplot(gs3[0,1], xlim=[0,1500], ylim=[0,1000])
    im5 = ax5.imshow(candels[0].data[880:1880,610:2110], cmap = pylab.cm.Greys, vmin=-0.01, 
        vmax=0.1, interpolation='bicubic')
    ax5.text(50,83,'CANDELS Mosaic', fontsize=10)    
    ax5.tick_params(axis='both',which='both',top='off',bottom='off', right='off', left='off')
    ax5.set_xticklabels([])
    ax5.set_yticklabels([])

    plt.show(block=False)
    
    fig.savefig('mosaic_demo.png', dpi=200, transparent=False)
    fig.savefig('mosaic_demo.pdf', dpi=200, transparent=False)
    

def plot_galfit():
    """
    Makes the plot comparing the GALFIT derived parameters from the new mosaics to CANDELS.
    """    

    import threedhst
    from my_python.pyspherematch import spherematch
    from matplotlib.ticker import ScalarFormatter, MultipleLocator, FormatStrFormatter
    
    fs=9
    matplotlib.rc('xtick',labelsize=fs)
    matplotlib.rc('ytick',labelsize=fs)
    
    mag_limit=21.
    
    ### read in CANDELS size measurements
    old = table.read('F160W_galfit_v4.0.cat.FITS')
    big_cat = table.read('/3DHST/Photometry/Release/v4.1/COSMOS/Catalog/cosmos_3dhst.v4.1.cat.FITS')
    old = old[(old['mag'] <= mag_limit) & (old['f'] < 2) & (big_cat['use_phot'] == 1)]
    
    ### read in new data from Arjen: H_v2.cat
    print 'Reading H_v4.cat'
    new = table.read('H_v4.cat', format='ascii')
    new = new[(new['mag'] <= mag_limit) & (new['f'] < 2)]
    
    ### cross-match coordinates
    idx_n, idx_o, d = spherematch(new['RA'], new['DEC'], old['ra'], old['dec'],tol = 0.5/3600.)
    print 'Total matches: {}'.format(len(idx_n))
    
    ### make plot
    fig = unicorn.catalogs.plot_init(square=True, xs=8., aspect=1.0, 
        fontsize=8, left=0.15, right=0.1, bottom=0.15, top=0.05)

    gs = gridspec.GridSpec(2,2)
    fig.subplots_adjust(wspace=0.20)
    fig.subplots_adjust(hspace=0.20)

    ax1 = fig.add_subplot(gs[0,0], xlim=[16,22], ylim=[16,22])
    ax1.plot([0,100],[0,100],':', color='gray', alpha=0.8)
    ax1.plot(new['mag'][idx_n],old['mag'][idx_o],'o', color='0.5', alpha=0.5, markersize=4.)
    ax1.set_xlabel('COSMOS-DASH F160W Magnitude', fontsize=fs)
    ax1.set_ylabel('CANDELS F160W Magnitude', fontsize=fs)
    diff = new['mag'][idx_n] - old['mag'][idx_o]
    print 'Magnitude: \nMEDIAN: {}\nNMAD: {}\n\n'.format(np.median(diff), threedhst.utils.nmad(diff))
    
    ax2 = fig.add_subplot(gs[0,1], xlim=[0.1,3.5], ylim=[0.1,3.5])
    ax2.plot([0.01,5],[0.01,5],':', color='gray', alpha=0.8)
    ax2.plot(new['re'][idx_n],old['re'][idx_o],'o', color='0.5', alpha=0.5, markersize=4.)
    ax2.set_xlabel('COSMOS-DASH $R_{eff}$ [arcsec]', fontsize=fs)
    ax2.set_ylabel('CANDELS $R_{eff}$ [arcsec]', fontsize=fs, labelpad=0)

    ax2.set_xscale('log')
    ax2.xaxis.set_major_formatter(ScalarFormatter())
    ax2.xaxis.set_major_locator(MultipleLocator(1))
    ax2.xaxis.set_major_formatter(FormatStrFormatter('%d'))
    ax2.set_yscale('log')
    ax2.yaxis.set_major_formatter(ScalarFormatter())
    ax2.yaxis.set_major_locator(MultipleLocator(1))
    ax2.yaxis.set_major_formatter(FormatStrFormatter('%d'))
    
    diff = np.log10(new['re'][idx_n]) - np.log10(old['re'][idx_o])
    print 'log(Reff): \nMEDIAN: {}\nNMAD: {}\n\n'.format(np.median(diff), threedhst.utils.nmad(diff))
    print np.median(new['re'][idx_n]/old['re'][idx_o])

    ax3 = fig.add_subplot(gs[1,0], xlim=[0.6, 8.9], ylim=[0.6,8.9])
    ax3.plot([-1,100],[-1,100],':', color='gray', alpha=0.8)
    ax3.plot(new['n'][idx_n],old['n'][idx_o],'o', color='0.5', alpha=0.5, markersize=4.)
    ax3.set_xlabel("COSMOS-DASH Sersic Index", fontsize=fs)
    ax3.set_ylabel("CANDELS Sersic Index", fontsize=fs, labelpad=0)
    ax3.set_xscale('log')
    ax3.xaxis.set_major_formatter(ScalarFormatter())
    ax3.xaxis.set_major_locator(MultipleLocator(1))
    ax3.xaxis.set_major_formatter(FormatStrFormatter('%d'))
    ax3.set_yscale('log')
    ax3.yaxis.set_major_formatter(ScalarFormatter())
    ax3.yaxis.set_major_locator(MultipleLocator(1))
    ax3.yaxis.set_major_formatter(FormatStrFormatter('%d'))

    diff = (new['n'][idx_n]) - (old['n'][idx_o])
    print 'log(n): \nMEDIAN: {}\nNMAD: {}\n\n'.format(np.median(diff), threedhst.utils.nmad(diff))

    ax4 = fig.add_subplot(gs[1,1], xlim=[0,1.0], ylim=[0,1.0])
    ax4.plot([0,100],[0,100],':', color='gray', alpha=0.8)
    ax4.plot(new['q'][idx_n],old['q'][idx_o],'o', color='0.5', alpha=0.5, markersize=4.)
    ax4.set_xlabel('COSMOS-DASH Axis Ratio', fontsize=fs)
    ax4.set_ylabel('CANDELS Axis Ratio', fontsize=fs)
    diff = new['q'][idx_n] - old['q'][idx_o]
    print 'Axis ratio: \nMEDIAN: {}\nNMAD: {}\n\n'.format(np.median(diff), threedhst.utils.nmad(diff))
    
    plt.show(block=False)
    fig.savefig('galfit_gyro_comp.pdf', dpi=200)
    fig.savefig('galfit_gyro_comp.png', dpi=200)
    
    
def mag_depth():
    
    orig = table.read('/3DHST/Photometry/Release/v4.0/COSMOS/Detection/cosmos_3dhst.v4.0.F160W_orig.cat', format='ascii.sextractor')
    cat = table.read('test7_drz_sci.cat', format='ascii.sextractor')
    
    fig = unicorn.catalogs.plot_init(square=True, xs=6., aspect=0.9, fontsize=10., left=0.15, right=0.15, top=0.1, bottom=0.15)
    ax = fig.add_subplot(111)

    n_orig, bins, patches = ax.hist(orig['MAG_AUTO'], bins=25, 
        range=[16., 30.], color='0.5', alpha=0.5, lw=2, histtype='step')
    n_cat, bins, patches = ax.hist(cat['MAG_AUTO'], bins=bins, 
        range=[16., 30.], color='red', alpha=0.5, histtype='step', lw=2)
        
    ax.set_xlim([17.,28.])
    ax.set_ylim([0, 8e3])
    ax.set_xlabel('F160W Magnitude', fontsize=12)
    ax.set_ylabel('N')
    
    ax2 = ax.twinx()#fig.add_subplot(111, sharex=ax, frameon=False)
    bin_c = bins[:-1] + (bins[1:] - bins[:-1])/2
    width = np.mean(bins[1:] - bins[:-1])
    ax2.plot(bin_c+width/2, (n_cat/142.429)/(n_orig/183.9), color='black', lw=2, alpha=0.8)
    ax2.plot([15.,31.], [1.0, 1.0], linestyle='dashed', color='0.8', alpha=0.5)
    ax2.plot([15.,31.], [0.9, 0.9], linestyle='dashed', color='0.8', alpha=0.5)
    ax2.plot([15.,31.], [0.75,0.75], linestyle='dashed', color='0.8', alpha=0.5)
    ax2.plot([15.,31.], [0.5, 0.5], linestyle='dashed', color='0.8', alpha=0.5)
    ax2.set_ylim([0.,1.3])
    ax2.set_xlim([17.,28.])
    ax2.set_ylabel('Fraction')
    plt.show(block=False)
    
    fig.savefig('mag_depth.pdf', dpi=200)
    fig.savefig('mag_depth.png', dpi=200)    
    
def psf_plot():
    
    psf = fits.open('test7_psf_v2.fits')
    wht = fits.open('test7_wht_v2.fits')
    
    ### make plot
    fig = unicorn.catalogs.plot_init(square=True, xs=6., aspect=0.6, 
        fontsize=8, left=0.05, right=0.05, bottom=0.1, top=0.00)

    gs = gridspec.GridSpec(1,2)
    fig.subplots_adjust(wspace=0.10)
    fig.subplots_adjust(hspace=0.10)
    
    ax1 = fig.add_subplot(gs[0], xlim=[0,68], ylim=[0,68])
    ax1.axis('off')
    ax2 = fig.add_subplot(gs[1], xlim=[0,68], ylim=[0,68])
    ax2.axis('off')
    
    im1 = ax1.imshow(psf[0].data/np.max(psf[0].data), cmap = pylab.cm.Greys_r, vmin=0.0, 
        vmax=0.015, interpolation='None')
    cbar_ax1 = fig.add_axes([0.06, 0.12, 0.40, 0.07])
    cbar1 = plt.colorbar(im1, cax = cbar_ax1,orientation='horizontal', ticks = [0.0, 0.005, 0.010, 0.015])
    cbar1.ax.set_xticklabels([0.0, 0.005, 0.010, 0.015])
    
    im2 = ax2.imshow(wht[0].data/np.max(wht[0].data), cmap = pylab.cm.Greys_r, vmin=0.7, 
        vmax=1., interpolation='None')
    cbar_ax2 = fig.add_axes([0.54, 0.12, 0.40, 0.07])
    cbar2 = plt.colorbar(im2, cax = cbar_ax2,orientation='horizontal', ticks = [0.7,0.8,0.9,1.0])
    cbar2.ax.set_xticklabels([0.7,0.8,0.9,1.0])
        
    plt.show(block=False)
    
    fig.savefig('psf_plot.pdf', dpi=200)
    fig.savefig('psf_plot.png', dpi=200)    
    
def rotation_check():
    
    import unicorn.survey_paper as sup
    from my_python.pyspherematch import spherematch as psh
    import matplotlib.colors as colors
    import matplotlib.cm as cmx
    
    os.chdir('/3DHST/Spectra/Work/14114/PREPARE/')
    
    REF_CAT = '../REF/IPAC_ACS.fits'
    NEW_CAT = 'test7_drz_sci.cat'
    
    ref_cat = table.read(REF_CAT)
    new_cat = table.read(NEW_CAT, format='ascii.sextractor')
    
    idx1, idx2, dd = psh(new_cat['X_WORLD'], new_cat['Y_WORLD'], ref_cat['ra'], ref_cat['dec'], tol=1./3600.)

    fig = unicorn.catalogs.plot_init(square=True, xs=10., aspect=1.1, 
        fontsize=8, left=0.18, right=0.02, bottom=0.10, top=0.10)
    ax = fig.add_subplot(111)
    jet = cm = plt.get_cmap('jet')
    cNorm = colors.Normalize(vmin=0, vmax=9)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)    
    
    xlim = [150.270, 149.99]
    ylim = [2.45, 2.751]
    xticklab = [r'$10^\mathrm{h}01^\mathrm{m}00^\mathrm{s}$', r'$10^\mathrm{h}00^\mathrm{m}45^\mathrm{s}$',r'$10^\mathrm{h}00^\mathrm{m}30^\mathrm{s}$',r'$10^\mathrm{h}00^\mathrm{m}15^\mathrm{s}$',r'$10^\mathrm{h}00^\mathrm{m}00^\mathrm{s}$']
    xtickv = [sup.degrees(10,01,00, hours=True),sup.degrees(10,00,45, hours=True),sup.degrees(10,00,30, hours=True),sup.degrees(10,00,15, hours=True),sup.degrees(10,00,00, hours=True)]
    yticklab = [r'$+02^\circ30^\prime00^{\prime\prime}$',r'$+02^\circ35^\prime00^{\prime\prime}$', r'$+02^\circ40^\prime00^{\prime\prime}$',r'$+02^\circ45^\prime00^{\prime\prime}$']
    ytickv = [sup.degrees(2, 30, 00, hours=False),sup.degrees(2, 35, 00, hours=False),sup.degrees(2, 40, 00, hours=False),sup.degrees(2, 45, 00, hours=False)]
    factor=20.
    
    ax.quiver(new_cat['X_WORLD'][idx1], new_cat['Y_WORLD'][idx1], (new_cat['X_WORLD'][idx1]-ref_cat['ra'][idx2])*3600.*2,
        (new_cat['Y_WORLD'][idx1]-ref_cat['dec'][idx2])*3600.*2,
        width=0.8, linewidth=0.3, units='dots', headwidth=4, minshaft=1, headlength=5, pivot='tail', scale=0.025)
    
    print np.median(new_cat['X_WORLD'][idx1]-ref_cat['ra'][idx2])*3600., np.median(new_cat['Y_WORLD'][idx1]-ref_cat['dec'][idx2])*3600.

    roots = ['icxe15010', 'icxe16010', 'icxe17010','icxe18010']
    labels = ['COSMOS-15', 'COSMOS-16','COSMOS-17','COSMOS-18']
    lines = ['-','--','-.',':']
    
    for root, label, linestyle in zip(roots, labels, lines):
        
        reg_file = root+'_asn.reg'
    
        poly = []
        with open(reg_file) as f:
            for line in f:
                if not line.startswith('fk5'):
                    region = line.split('#')[0]
                    poly.append(sup.polysplit(region=region, get_shapely=True))

        shifts = table.read('shifts_{}.txt'.format(root), format='ascii', 
            names=('file','x','y','rot','scale','x_rms','y_rms'))
        
        cc = 0
        xcen_all = []
        ycen_all = []
        for j,(pp, file) in enumerate(zip(poly, shifts['file'])):
            cc += 1.
            color = scalarMap.to_rgba(cc)
            x, y = pp.exterior.xy
            flt = fits.open(file)
            xcen = flt[1].header['CRVAL1O']
            ycen = flt[1].header['CRVAL2O']
            x_off = (flt[1].header['CRVAL1B']-flt[1].header['CRVAL1O'])*factor
            y_off = (flt[1].header['CRVAL2B']-flt[1].header['CRVAL2O'])*factor
            xcen_all.append(xcen)
            ycen_all.append(ycen)
            ax.plot(x,y,'-', color=color)
            #ax.annotate("",xy=(xcen+x_off, ycen+y_off), xytext=(xcen, ycen), 
            #    arrowprops=dict(arrowstyle='->', color=color))
            #ax.text(xcen, ycen+0.005, file.split('_')[0], fontsize=9, va='top', ha='center')
        
        ax.plot(xcen_all, ycen_all, linestyle=linestyle, markersize=10., 
            color='black', alpha=0.7, label=label) 
        ax.plot(xcen_all, ycen_all, '+', markersize=10., color='0.5', alpha=0.7) 
        #ax.text(xcen_all[0], ycen_all[0]-0.005, label, fontsize=10, va='top', ha='center')


    ax.set_xlim(xlim)
    ax.set_ylim(ylim)     
    ax.set_xticklabels(xticklab)
    xtick = ax.set_xticks(xtickv)
    ax.set_yticklabels(yticklab)
    ytick = ax.set_yticks(ytickv)
    ax.legend(loc='lower right', frameon=False, labelspacing=0.8, fontsize=9, handlelength=10, borderpad=5.)
    plt.show(block=False)
        
    fig.savefig('rotation_check.png'.format(label.lower()), dpi=200, transparent=False)          

def rotation_fits():

    import glob
    import matplotlib.colors as colors
    import matplotlib.cm as cmx
    
    files = glob.glob('PREP_ROT/shifts_icxe*010.txt')
    
    fig = unicorn.catalogs.plot_init(xs=10.,aspect=0.45, left=0.2, right=0.05, bottom=0.15, top=0.15, NO_GUI=False)    
    gs = gridspec.GridSpec(1,2,top=0.9, bottom=0.15)
    fig.subplots_adjust(wspace=0.2)
    fig.subplots_adjust(hspace=0.2)
    fs = 10
    matplotlib.rc('xtick',labelsize=fs)
    matplotlib.rc('ytick',labelsize=fs)

    jet = cm = plt.get_cmap('jet_r')
    cNorm = colors.Normalize(vmin=0, vmax=8)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)


    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])

    time = np.array([0, 368, 713, 1060, 1405, 1777, 2147, 2519]) + 333.
    lines = ['-','--','-.',':']
    labels = ['COSMOS-15', 'COSMOS-16', 'COSMOS-17','COSMOS-18']
    
    for j, file in enumerate(files):
        root = '{}'.format(file.split('_')[1].split('.')[0])
        ax1.set_xlim([-0.5, 13])
        ax1.set_ylim([-0.075,0.075])
        ax1.set_xlabel('Read Number', fontsize=fs)
        ax1.set_ylabel('Rotation [deg]', fontsize=fs)
        root = '{}'.format(file.split('_')[1].split('.')[0])
            
        data = table.read(file, format='ascii',  names=('file','x','y','rot','scale','x_rms','y_rms'))    
        rot_rate = [0.0]
                
        for k, file in enumerate(data['file'][1:]):
            if k+1 < 4:
                t_exp = 255
            else:
                t_exp = 277
            drift = table.read('PREP_ROT/'+file.split('_')[0]+'_shifts.txt', format='ascii', 
                names=('file','x','y','rot','scale','x_rms','y_rms'))
            
            rot = drift['rot']
            tt = (drift['rot'] > 180.)
            rot[tt] = drift['rot'][tt] - 360.
            rot_rate.append(rot[0]-rot[-1])
            ax1.plot(np.arange(len(drift['rot']))+1, drift['rot'], linestyle=lines[j], 
                color='0.5', label=labels[j], zorder=0)
            ax1.scatter(np.arange(len(drift['rot']))+1, drift['rot'], c='0.5', s=10., edgecolors='black', alpha=0.7)
        
        
        ax2.set_xlim([-5., 3100.])
        ax2.set_ylim([-0.075,0.075])
        ax2.set_xlabel('Time Since Beginning of Orbit\n[seconds]', fontsize=fs)
        ax2.set_ylabel('Rotation During Pointing [deg]', fontsize=fs)
        
        ax2.plot(time, np.array(rot_rate),linestyle=lines[j], color='0.5', zorder=0)
        ax2.scatter(time, np.array(rot_rate), c=range(len(rot_rate)), cmap='jet', s=35., edgecolors='black', alpha=0.7)
        print rot_rate
        
    plt.show(block=False)
    fig.savefig('rotation_fits.png', dpi=200, transparent=False)
    #fig.savefig('overall_offsets_vs_time.pdf', dpi=200, transparent=False)

def mag_radius_test():
    
    from matplotlib.ticker import ScalarFormatter
    import my_python.mk_region_file
    
    orig = table.read('test7_drz_sci_1.8.cat', format='ascii.sextractor')
    cat = table.read('test7_drz_sci_0.8.cat', format='ascii.sextractor')
    
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
    
    ax1 = fig.add_subplot(gs1[0,0], ylim=[0.8, 20], xlim=[14,26])
    ax1.plot(orig['MAG_AUTO'],orig['FLUX_RADIUS'], '.', color='0.4',markersize=0.9)
    ax1.plot(cat['MAG_AUTO'],cat['FLUX_RADIUS'], 'o', color='black', markersize=1, alpha=0.5)
    ax1.set_xlabel('MAG_AUTO F160W')
    ax1.set_ylabel('FLUX_RADIUS')
    ax1.set_yscale('log')
    cr = (cat['FLUX_RADIUS']*0.1/0.06 < 2.)
    stars = (cat['MAG_AUTO'] > 15.) & (cat['MAG_AUTO'] < 22.) & (cat['FLUX_APER_5']/cat['FLUX_APER'] > 1.1) & (cat['FLUX_APER_5']/cat['FLUX_APER'] < 1.2)
    
    #ax1.plot(cat['MAG_AUTO'][stars],cat['FLUX_RADIUS'][stars], 'o', color='red', markersize=2, alpha=1.0, markeredgecolor='red')
    print 'STARS: mean: {} / median{}'.format(np.mean(cat['FWHM_IMAGE'][stars]), np.median(cat['FWHM_IMAGE'][stars]))
    
    ax1.yaxis.set_major_formatter(ScalarFormatter())
    
    ax2 = fig.add_subplot(gs1[0,1], ylim=[0.5,4], xlim=[14,26])
    ax2.plot(orig['MAG_AUTO'], orig['FLUX_APER_5']/orig['FLUX_APER'], '.', color='0.4',markersize=0.9)
    ax2.plot(cat['MAG_AUTO'], cat['FLUX_APER_5']/cat['FLUX_APER'], 'o', color='black', markersize=1., alpha=0.5)
    #ax2.plot(cat['MAG_AUTO'][stars], cat['FLUX_APER_5'][stars]/cat['FLUX_APER'][stars], 'o', color='red', markersize=2., alpha=1.0, markeredgecolor='red')
    ax2.plot(cat['MAG_AUTO'][cr], cat['FLUX_APER_5'][cr]/cat['FLUX_APER'][cr], 'o', color='blue', markersize=2., alpha=1.0)
    ax2.set_xlabel('MAG_AUTO F160W')
    ax2.set_ylabel('Flux (2.0\")/Flux (0.5\")')
    ax2.yaxis.set_major_formatter(ScalarFormatter())
    
    my_python.mk_region_file.mk_region_file_from_lists(cat['X_WORLD'][stars],cat['Y_WORLD'][stars],outfile = 'stars', printids='no', color='cyan')
    my_python.mk_region_file.mk_region_file_from_lists(cat['X_WORLD'][cr],cat['Y_WORLD'][cr],outfile = 'cr', printids='no', color='yellow')
    
    plt.show(block=False)
    
    #fig.savefig('mag_radius_test.pdf', dpi=200, transparent=False)
    #fig.savefig('mag_radius_test.png', dpi=200, transparent=False)

def signal_to_noise_test():

    import threedhst.utils

    top = 0.05
    bottom = 0.15
    left = 0.15
    fig = unicorn.catalogs.plot_init(xs=10,aspect=0.5, left=left, right=0.05, bottom=bottom, top=top, NO_GUI=False)
    fig.subplots_adjust(wspace=0.01)
    fig.subplots_adjust(hspace=0.01)
    fs = 9
    matplotlib.rc('xtick',labelsize=fs)
    matplotlib.rc('ytick',labelsize=fs)
    
    etc_mag = [18.,19.,20.,21.,22.,23.,24.,25.]
    etc_sn = np.array([509.,318.,195.8,116.7,65.1,33.1,15.2,6.5])
    etc_sn_high = np.array([493.,295.,166.9,85.9,39.97,17.2,7.09,2.87])
    
    factor = (1.-(0.8/3/2.75)) #0.678
    
    gs = gridspec.GridSpec(2,4)
    
    asn_files = glob.glob('icxe1*010_asn.fits')
    
    for file in asn_files:
        asn = threedhst.utils.ASNFile(file)
        for ii, exp in enumerate(asn.exposures):

            cat = table.read('/3DHST/Spectra/Work/14114/TEST_SEXTR/{}_drz_sci.cat'.format(exp), format='ascii.sextractor')
            
            ax = fig.add_subplot(gs[ii], ylim=[3., 200], xlim=[21,25.5])
            ax.plot(cat['MAG_APER'], cat['FLUX_APER']/(cat['FLUXERR_APER']/factor), '.', color='0.6',markersize=0.9)
    
            stars = (cat['MAG_AUTO'] > 15.) & (cat['MAG_AUTO'] < 22.) & (cat['FLUX_APER_6']/cat['FLUX_APER_1'] > 1.1) & (cat['FLUX_APER_6']/cat['FLUX_APER_1'] < 1.2)
            ax.plot(cat['MAG_APER'][stars], cat['FLUX_APER'][stars]/(cat['FLUXERR_APER'][stars]/factor), '^', color='black', markersize=2, alpha=1.0, markeredgecolor='0.3')
    
            if file == 'icxe18010_asn.fits':
                ax.text(24,70,'Pointing {}'.format(ii+1), fontsize=9, va='top', ha='center')
                #ax.plot(etc_mag,etc_sn_high, 'o', color='#b10026',markersize=2)
                ax.plot(etc_mag,etc_sn_high, '--', color='#b10026',markersize=2)
                ax.plot(etc_mag,etc_sn_high/1.07, '-', color='#b10026',markersize=2)
                #ax.plot(etc_mag,etc_sn, 'o', color='black',markersize=2)
                ax.plot(etc_mag,etc_sn, '--', color='black',markersize=2)
                ax.plot(etc_mag,etc_sn/1.37, '-', color='black',markersize=2)
                ax.plot([10,30],[5.,5.], ':', color='0.5')
                #ax.plot([10,30],[4.,4.], ':', color='0.5')
                ax.plot([25.,25.],[1.,3000.], ':', color='0.5')
                ax.set_yscale('log')
                if (ii == 0) or (ii == 4):
                    ax.set_ylabel('S/N [r = 0.2" aperture]', fontsize=fs)
                else:
                    ax.set_yticklabels([])
                if (ii < 3):
                    ax.set_xticklabels([])
                else:
                    ax.set_xlabel('AB Magnitude\n[r = 0.2" aperture]', fontsize=fs)
                    ax.set_xticklabels([22,23,24,25], fontsize=fs)
                    xtick = ax.set_xticks([22,23,24,25])
                    
    
    plt.show(block=False)
    fig.savefig('signal_to_noise_test.pdf', dpi=200, transparent=False)
    fig.savefig('signal_to_noise_test.png', dpi=200, transparent=False)
    
def sources_of_noise():
    
    top = 0.075
    bottom = 0.1
    left = 0.15
    fig = unicorn.catalogs.plot_init(xs=10,aspect=0.5, left=left, right=0.1, bottom=bottom, top=top, NO_GUI=False)
    fig.subplots_adjust(wspace=0.15)
    fig.subplots_adjust(hspace=0.1)
    fs = 8
    matplotlib.rc('xtick',labelsize=fs)
    matplotlib.rc('ytick',labelsize=fs)
        
    gs1 = gridspec.GridSpec(1,2)
    
    ax1 = fig.add_subplot(gs1[0,0], ylim=[0.0,0.27], xlim=[0.1,2.5])

    exptime=275.
    rdnoise = 15.
    xx = np.linspace(0.2,2.4,10)
    noise_15 = np.sqrt(rdnoise**2+xx*exptime + (xx*exptime*0.002)**2 + (0.17*exptime))/exptime
    rd_15 = np.sqrt(rdnoise**2)/exptime
    bg_15 = np.sqrt(xx*exptime)/exptime

    rdnoise = np.sqrt(11)*15.
    noise_21_11 = np.sqrt(rdnoise**2+xx*exptime + (xx*exptime*0.002)**2 + (0.17*exptime))/exptime
    rd_21_11 = np.sqrt(rdnoise**2)/exptime

    rdnoise = np.sqrt(2)*21.
    noise_21_2 = np.sqrt(rdnoise**2+xx*exptime + (xx*exptime*0.002)**2 + (0.17*exptime))/exptime
    rd_21_2 = np.sqrt(rdnoise**2)/exptime

    rdnoise = np.sqrt(4)*21.
    noise_21_4 = np.sqrt(rdnoise**2+xx*exptime + (xx*exptime*0.002)**2 + (0.17*exptime))/exptime
    rd_21_4 = np.sqrt(rdnoise**2)/exptime


    ax1.plot(xx, noise_21_11, '-',color='#b10026',lw=1, label=r'RDN = 15$\times\sqrt{11}$ e$^{-}$: (11$\times$25 sec exp (w/ ramp fits))')
    ax1.plot([0.,3.], [rd_21_11,rd_21_11], ':',color='#b10026',lw=1)

    ax1.plot(xx, noise_21_4, '-',color='#fc4e2a',lw=1, label=r'RDN = 21$\times\sqrt{4}$ e$^{-}$')
    ax1.plot([0.,3.], [rd_21_4,rd_21_4], ':',color='#fc4e2a',lw=1)

    ax1.plot(xx, noise_21_2, '-',color='#feb24c',lw=1, label=r'RDN = 21$\times\sqrt{2}$ e$^{-}$')
    ax1.plot([0.,3.], [rd_21_2,rd_21_2], ':',color='#feb24c',lw=1)

    ax1.plot(xx, noise_15, '-',color='black',lw=1, label=r'RDN = 15 e$^{-}$: one 275 sec exp (/w ramp fit)')
    ax1.plot([0.,3.], [rd_15,rd_15], ':',color='black',lw=1)
    ax1.plot(xx, bg_15, '--',color='black',lw=1)

    ax1.text(2.0,0.105,'full noise model', rotation=11, va='center', ha='center', fontsize=fs, backgroundcolor='white', alpha=0.7)
    ax1.text(2.0,0.085,'sky noise', rotation=12, va='center', ha='center', fontsize=fs, backgroundcolor='white', alpha=0.7)
    ax1.text(2.0,0.053,'read noise',va='center', ha='center', fontsize=fs, backgroundcolor='white', alpha=0.7)
    #ax1.fill_between([0.6,0.8], 0.,2.,color='0.5', alpha=0.5)
    ax1.legend(loc='upper left', fontsize=fs, frameon=False)
    ax1.set_ylabel('$\sigma [e^-/s]$', fontsize=fs)
    ax1.set_xlabel('Background [e$^-$/s]', fontsize=fs)
    
    ax2 = fig.add_subplot(gs1[0,1], ylim=[0.2,1.25], xlim=[0.1,2.5])
    ax2.plot(xx, 2.5*np.log10(noise_21_11/noise_15), '-',color='#b10026',lw=1)
    ax2.plot(xx, 2.5*np.log10(noise_21_4/noise_15), '-',color='#fc4e2a',lw=1)
    ax2.plot(xx, 2.5*np.log10(noise_21_2/noise_15), '-',color='#feb24c',lw=1)
    ax2.grid(True, linestyle=':',linewidth=0.5)
    ax2.fill_between([0.6,0.8], 0.,2.,color='0.5', alpha=0.5)
    ax2.fill_between([2.1,2.3], 0.,2.,color='0.5', alpha=0.5)
    ax2.text(0.7, 1.1,'nominal\nbackground', rotation=90., va='center', ha='center', fontsize=8)
    ax2.text(2.2, 0.9,'Oct. 2015\nbackground', rotation=90., va='center', ha='center', fontsize=8)
    ax2.set_ylabel('2.5*log$_{10}$($\sigma/\sigma_{15}$) [mag]', fontsize=fs)
    ax2.set_xlabel('Background [e$^-$/s]', fontsize=fs)
    
    
    plt.show(block=False)
    
    fig.savefig('sources_of_noise.pdf', dpi=200, transparent=False)
    fig.savefig('sources_of_noise.png', dpi=200, transparent=False)
    
def plot_empty_apertures():
    
    import glob
    import matplotlib.colors as colors
    import matplotlib.cm as cmx
    import threedhst.utils
    from scipy.interpolate import interp1d
    import scipy

    growth = table.read('COSMOS_growth_curve.dat', format='ascii')
    fp = interp1d(growth['radius_arcsec'],growth['fraction'])

    fig = unicorn.catalogs.plot_init(xs=10.,aspect=0.45, left=0.2, right=0.05, bottom=0.15, top=0.15, NO_GUI=False)    
    gs = gridspec.GridSpec(1,2,top=0.9, bottom=0.15)
    fig.subplots_adjust(wspace=0.2)
    fig.subplots_adjust(hspace=0.2)
    fs = 10
    matplotlib.rc('xtick',labelsize=fs)
    matplotlib.rc('ytick',labelsize=fs)

    jet = cm = plt.get_cmap('jet')
    cNorm = colors.Normalize(vmin=0, vmax=8)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

    ax1 = fig.add_subplot(gs[0], ylim=[0.1,1.], xlim=[0.1,0.5])
    ax2 = fig.add_subplot(gs[1], ylim=[10.,15.], xlim=[0.1,2])

    files = glob.glob('icxe*0_asn.fits')
    for file in files:
        asn = threedhst.utils.ASNFile('{}'.format(file))
        for j, exp in enumerate(asn.exposures):
            data = table.read('{}_drz_sci.apertures.F160W.v4.0.dat'.format(exp), format='ascii')    
            color = scalarMap.to_rgba(j+1)
            ax1.plot(data['radius_pix']*0.12825, data['sigma'], 'o-', color=color, alpha=0.8, mfc=color, mec=color)
            fp1 = interp1d(data['radius_pix']*0.12825, data['sigma'])
            if (file == 'icxe15010_asn.fits') & ((j == 0) | (j == 7)):
                coeff, var = scipy.optimize.curve_fit(sigma, data['radius_pix'], data['sigma'], p0=[0.1, 1])
                ax1.plot(data['radius_pix']*0.12825, 0.175*coeff[0]*data['radius_pix']**2,'--',color='black')
                ax1.plot(data['radius_pix']*0.12825, 0.175*coeff[0]*data['radius_pix'],'--',color='black')

            growth_y = fp(data['radius_pix']*0.12825)
            ax2.plot(data['radius_pix']*0.12825, -2.5*np.log10(5*(data['sigma']**2/fp(data['radius_pix']*0.12825))*15369.)+25. ,'o-', color=color, alpha=0.8, mfc=color, mec=color)
            
            
    ax1.set_ylabel('$\sigma$ [e-/s]')
    ax1.set_xlabel('Radius [arcsec]')
    

    #ax2.plot(growth['radius_arcsec'],1/growth['fraction'],'-')
    #ax2.plot(growth['radius_arcsec'], 1/fp(growth['radius_arcsec']),'o')
            
    guided = glob.glob('../TEST_HOPR_COOPER/REDUCE/ibt355*_drz_sci.apertures.F160W.v4.0.dat')
    dash = glob.glob('../TEST_HOPR_COOPER/REDUCE/ibt305*_drz_sci.apertures.F160W.v4.0.dat')
    #ax1.set_yscale('log')
    
    for file in guided:
        data = table.read(file, format='ascii')
        #ax1.plot(data['radius_pix']*0.128, data['sigma'], 'o-', color='black', alpha=0.8)
        
    for file in dash:
        data = table.read(file, format='ascii')
        #ax1.plot(data['radius_pix']*0.128, data['sigma'], 'o-', color='gray', alpha=0.8)

    

    plt.show(block=False)
    
def gauss(x, *p):
    
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))
    
def sigma(x, *p):
    
    alpha, beta = p
    sigma1 = 0.175
    return sigma1*alpha*x**beta
            