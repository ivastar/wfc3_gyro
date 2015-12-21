import glob
import os
import astropy
from astropy.io import fits
from astropy.table import Table as table
import numpy as np
import threedhst
import drizzlepac
from drizzlepac import astrodrizzle


def prepare_test_exp(RAW_PATH='../RAW/'):
    
    import threedhst.prep_flt_astrodrizzle 
    
    all_asn = glob.glob('*asn.fits')

    for file in all_asn:
        master_root = file.split('_asn')[0]
        asn = threedhst.utils.ASNFile('{}_asn.fits'.format(master_root))
        for root in asn.exposures:
            os.system('rsync -av {}/{}_flt.fits.gz .'.format(RAW_PATH, root))
            os.system('gunzip -f {}_flt.fits.gz'.format(root))
        
        threedhst.prep_flt_astrodrizzle.subtract_flt_background(root=master_root)    
        drizzlepac.astrodrizzle.AstroDrizzle(master_root+'_asn.fits', clean=False, final_pixfrac=1.0, context=False, final_bits=576, resetbits=0, preserve=False, driz_cr_snr='8.0 5.0', driz_cr_scale = '2.5 0.7')
        

def make_lists():
    
    all_asn = glob.glob('*asn.fits')
    asn_252 = []
    asn_277 = []
    for file in all_asn:
        asn =  threedhst.utils.ASNFile(file)
        flt = fits.open(asn.exposures[0]+'_flt.fits')            
        if flt[1].header['SAMPTIME'] < 260:
            asn_252.append(asn.product.lower())
        else:
            asn_277.append(asn.product.lower())
    t1 = table([asn_252,], names=('files',))   
    t1.write('asn_252.txt', format='ascii') 
    t2 = table([asn_277,], names=('files',))
    t2.write('asn_277.txt', format='ascii')

def split_the_reads_and_drizzle():
    
    import wfc3_gyro.prepare
    
    all_asn = glob.glob('*asn.fits')
    for file in all_asn:
        asn = threedhst.utils.ASNFile(file) 
        for exp in asn.exposures:
            wfc3_gyro.prepare.split_IMA(root=exp)
            wfc3_gyro.prepare.make_pointing_asn(root=exp,master_root=asn.product.lower())
            wfc3_gyro.prepare.subtract_background_reads(root=exp,master_root=asn.product.lower())
            drizzlepac.astrodrizzle.AstroDrizzle('{}_asn.fits'.format(exp), output=exp+'_split', 
                clean=True, final_scale=0.12, final_pixfrac=0.8, context=False, 
                resetbits=0, final_bits=576, driz_cr = True, driz_sep_bits=576, 
                preserve=False, driz_cr_snr='8.0 5.0', driz_cr_scale = '2.5 0.7')
        

def sigma_table():
    
    import stwcs
    import scipy
    import scipy.optimize
    
    all_asn = glob.glob('*asn.fits')
    
    ### open file: filename, sigma_flt, sigma_drz, exptime
    fp = open('sigma_table.txt','w')
    fp.write('#exp sigma_flt sigma_flt3 sigma_flt5 sigma_drz sigma_drz3 sigma_drz5 exptime bg flag\n')
    
    for asn_file in all_asn:
        
        if asn_file.endswith('0_asn.fits'):
            
            root = asn_file.split('_asn')[0]
            asn = threedhst.utils.ASNFile('{}'.format(asn_file))

            ref = fits.open('{}_drz_sci.fits'.format(root))
            ref_wcs = stwcs.wcsutil.HSTWCS(ref, ext=0)

            seg = fits.open('{}_drz_seg.fits'.format(root))    
            seg_data = np.cast[np.float32](seg[0].data)
          
            yi, xi = np.indices((1014,1014))
            
            fp3 = np.array([[1,1,1],[1,1,1],[1,1,1]], dtype=int)
            fp5 = np.ones([5,5], dtype=int)
            
            for exp in asn.exposures:
                
                flt = fits.open('{}_flt.fits'.format(exp))
                flt_wcs = stwcs.wcsutil.HSTWCS(flt, ext=1)
                
                ### Original FLT
                blotted_seg = astrodrizzle.ablot.do_blot(seg_data, ref_wcs, flt_wcs, 1, coeffs=True, 
                    interp='nearest', sinscl=1.0, stepsize=10, wcsmap=None) 
                mask = (blotted_seg == 0) & (flt['DQ'].data == 0) & (xi > 10) & (yi > 10) & (xi < 1004) & (yi < 1004)   
                n, bin_edge  = np.histogram(flt[1].data[mask], bins=300, range=[-1,1], normed=True)
                centers = (bin_edge[:-1] + (bin_edge[1:] - bin_edge[:-1])/2)
                flt_coeff, flt_var = scipy.optimize.curve_fit(gauss, centers, n, p0=[np.max(n),0.,1.])
                
                ### Smooth by 3
                #flt3 = scipy.ndimage.median_filter(flt[1].data, size=3, mode='nearest')
                flt3 = scipy.ndimage.filters.generic_filter(flt[1].data, function=np.mean, footprint=fp3)
                n, bin_edge  = np.histogram(flt3[mask], bins=300, range=[-1,1], normed=True)
                centers = (bin_edge[:-1] + (bin_edge[1:] - bin_edge[:-1])/2)
                flt_coeff3, flt_var3 = scipy.optimize.curve_fit(gauss, centers, n, p0=[np.max(n),0.,1.])
                

                ### Smooth by 5
                #flt5 = scipy.ndimage.median_filter(flt[1].data, size=5, mode='nearest')
                flt5 = scipy.ndimage.filters.generic_filter(flt[1].data, function=np.mean, footprint=fp5)
                n, bin_edge  = np.histogram(flt5[mask], bins=300, range=[-1,1], normed=True)
                centers = (bin_edge[:-1] + (bin_edge[1:] - bin_edge[:-1])/2)
                flt_coeff5, flt_var5 = scipy.optimize.curve_fit(gauss, centers, n, p0=[np.max(n),0.,1.])
                
                ### Original DRZ
                drz = fits.open('{}_split_drz_sci.fits'.format(exp))
                drz_wcs = stwcs.wcsutil.HSTWCS(drz, ext=0)
                drz_sh = np.shape(drz[0].data)
                dyi, dxi = np.indices(drz_sh)
                blotted_seg = astrodrizzle.ablot.do_blot(seg_data, ref_wcs, drz_wcs, 1, coeffs=False, 
                    interp='nearest', sinscl=1.0, stepsize=10, wcsmap=None)
                mask = (blotted_seg == 0) & (dxi > 10) & (dyi > 10) & (dxi < drz_sh[1]-10) & (dyi < drz_sh[0]-10)
                n, bin_edge = np.histogram(drz[0].data[mask], bins=300, range=[-1.,1.], normed=True)
                centers = (bin_edge[:-1] + (bin_edge[1:] - bin_edge[:-1])/2)
                drz_coeff, drz_var = scipy.optimize.curve_fit(gauss, centers, n, p0=[np.max(n),0.,1.])
                
                ### Smooth by 3
                #drz3 = scipy.ndimage.median_filter(drz[0].data, size=3, mode='nearest')
                drz3 = scipy.ndimage.filters.generic_filter(drz[0].data, function=np.mean, footprint=fp3)
                n, bin_edge = np.histogram(drz3, bins=300, range=[-1.,1.], normed=True)
                centers = (bin_edge[:-1] + (bin_edge[1:] - bin_edge[:-1])/2)
                drz_coeff3, drz_var3 = scipy.optimize.curve_fit(gauss, centers, n, p0=[np.max(n),0.,1.])
                
                ### Smooth by 5
                #drz5 = scipy.ndimage.median_filter(drz[0].data, size=5, mode='nearest')
                drz5 = scipy.ndimage.filters.generic_filter(drz[0].data, function=np.mean, footprint=fp5)
                n, bin_edge = np.histogram(drz5, bins=300, range=[-1.,1.], normed=True)
                centers = (bin_edge[:-1] + (bin_edge[1:] - bin_edge[:-1])/2)
                drz_coeff5, drz_var5 = scipy.optimize.curve_fit(gauss, centers, n, p0=[np.max(n),0.,1.])

                flag = 2
                
                print '{}\t{:7.4f}\t{:7.4f}\t{:7.4f}\t{:7.4f}\t{:7.4f}\t{:7.4f}\t{:7.2f}\t{:7.4f}\t{}\n'.format(exp, np.abs(flt_coeff[2]), np.abs(flt_coeff3[2]), np.abs(flt_coeff5[2]), np.abs(drz_coeff[2]), np.abs(drz_coeff3[2]), np.abs(drz_coeff5[2]), flt[1].header['SAMPTIME'], flt[0].header['BGCOMP1'],flag)
                fp.write('{}\t{:7.4f}\t{:7.4f}\t{:7.4f}\t{:7.4f}\t{:7.4f}\t{:7.4f}\t{:7.2f}\t{:7.4f}\t{}\n'.format(exp, np.abs(flt_coeff[2]), np.abs(flt_coeff3[2]), np.abs(flt_coeff5[2]), np.abs(drz_coeff[2]), np.abs(drz_coeff3[2]), np.abs(drz_coeff5[2]), flt[1].header['SAMPTIME'], flt[0].header['BGCOMP1'],flag))
                
    fp.close()    
            
def plot_sigma():
    
    import unicorn
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    import scipy
    
    data = table.read('sigma_table.txt',format='ascii')
    
    fig = unicorn.catalogs.plot_init(xs=10,aspect=0.5, left=0.12, right=0.05, bottom=0.1, top=0.1, NO_GUI=False)
    gs = gridspec.GridSpec(1,2,top=0.9, bottom=0.15)
    xx = np.linspace(0., 0.2, 30.)

    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])
    
    mask1 = (data['flag'] == 2 ) & (data['exptime'] < 260)
    ax1.plot(data['sigma_flt'][mask1], data['sigma_drz'][mask1], 'o', color='0.5', alpha=0.5, label='252s')
    coeff, var = scipy.optimize.curve_fit(fit_func, data['sigma_flt'][mask1], data['sigma_drz'][mask1], p0=[0.,0.5])
    ff = fit_func(xx, *coeff)
    ax1.plot(xx, ff, ':', color='0.5', alpha=0.5)
    
    ax2.plot(data['sigma_flt3'][mask1], data['sigma_drz3'][mask1], 'o', color='0.5', alpha=0.5, label='252s')
    ax2.plot(data['sigma_flt5'][mask1], data['sigma_drz5'][mask1], 'o', color='0.5', alpha=0.5, label='252s')

    mask2 = (data['flag'] == 2 ) & (data['exptime'] > 260)
    ax1.plot(data['sigma_flt'][mask2], data['sigma_drz'][mask2], '^', color='0.5', alpha=0.5, label='277s')
    coeff, var = scipy.optimize.curve_fit(fit_func, data['sigma_flt'][mask2], data['sigma_drz'][mask2], p0=[0.,0.5])
    ff = fit_func(xx, *coeff)
    ax1.plot(xx, ff, ':', color='0.5', alpha=0.5)

    ax2.plot(data['sigma_flt3'][mask2], data['sigma_drz3'][mask2], '^', color='0.5', alpha=0.5, label='277s')
    ax2.plot(data['sigma_flt5'][mask2], data['sigma_drz5'][mask2], '^', color='0.5', alpha=0.5, label='277s')

    mask3 = (data['flag'] == 0 )
    ax1.plot(data['sigma_flt'][mask3], data['sigma_drz'][mask3], 'o', color='blue', alpha=0.5, label='COSMOS, 252s, guided')
    ax2.plot(data['sigma_flt3'][mask3], data['sigma_drz3'][mask3], 'o', color='blue', alpha=0.5, label='COSMOS, 252s, guided')
    ax2.plot(data['sigma_flt5'][mask3], data['sigma_drz5'][mask3], 'o', color='blue', alpha=0.5, label='COSMOS, 252s, guided')

    mask4 = (data['flag'] == 1 ) & (data['exptime'] < 260)
    ax1.plot(data['sigma_flt'][mask4], data['sigma_drz'][mask4], 'o', color='red', alpha=0.5, label='COSMOS, 252s')
    ax2.plot(data['sigma_flt3'][mask4], data['sigma_drz3'][mask4], 'o', color='red', alpha=0.5, label='COSMOS, 252s')
    ax2.plot(data['sigma_flt5'][mask4], data['sigma_drz5'][mask4], 'o', color='red', alpha=0.5, label='COSMOS, 252s')

    mask5 = (data['flag'] == 1 ) & (data['exptime'] > 260)
    ax1.plot(data['sigma_flt'][mask5], data['sigma_drz'][mask5], '^', color='red', alpha=0.5, label='COSMOS, 277s')
    ax2.plot(data['sigma_flt3'][mask5], data['sigma_drz3'][mask5], '^', color='red', alpha=0.5, label='COSMOS, 277s')
    ax2.plot(data['sigma_flt5'][mask5], data['sigma_drz5'][mask5], '^', color='red', alpha=0.5, label='COSMOS, 277s')

    ax1.set_xlim([0.05,0.15])
    ax1.set_ylim([0.05,0.15])
    ax1.plot([0,0.2],[0.,0.2],'-', color='black',alpha=0.5)
    ax1.set_xlabel('$\sigma_{FLT}$')
    ax1.set_ylabel('$\sigma_{DRZ}$')
    ax1.legend(loc=2, frameon=False, numpoints=1, fontsize=9)

    ax2.set_xlim([0.0,0.06])
    ax2.set_ylim([0.0,0.06])
    ax2.plot([0,0.2],[0.,0.2],'-', color='black',alpha=0.5)
    ax2.set_xlabel('$\sigma_{FLT}$')
    ax2.set_ylabel('$\sigma_{DRZ}$')
    ax2.text(0.04, 0.035, '3x3 median top hat')
    ax2.text(0.03, 0.025, '5x5 median top hat')

    plt.show(block=False)
    


def plot_flt_bg():
    
    import unicorn
    import matplotlib.pyplot as plt
    fs= 10
    
    data = table.read('sigma_table.txt',format='ascii')
    
    fig = unicorn.catalogs.plot_init(xs=5,aspect=1.0, left=0.12, right=0.05, bottom=0.1, top=0.1, NO_GUI=False)
    xx = np.linspace(0., 0.2, 30.)
    ax = fig.add_subplot(111)

    mask1 = (data['flag'] == 2 ) & (data['exptime'] < 260)
    ax.plot(data['bg'][mask1], data['sigma_flt'][mask1], 'o', color='0.5', alpha=0.5, label='252s')

    mask2 = (data['flag'] == 2 ) & (data['exptime'] > 260)
    ax.plot(data['bg'][mask2], data['sigma_flt'][mask2], '^', color='0.5', alpha=0.5, label='277s')
    
    
    mask4 = (data['flag'] <= 1 ) & (data['exptime'] < 260)
    ax.plot(data['bg'][mask4], data['sigma_flt'][mask4], 'o', color='red', alpha=0.5, label='COSMOS 252s')
    
    mask5 = (data['flag'] <= 1 ) & (data['exptime'] > 260)
    ax.plot(data['bg'][mask5], data['sigma_flt'][mask5], '^', color='red', alpha=0.5, label='COSMOS 277s')
    
    xx = np.linspace(0.2,2.4,10)
    for exptime, rdnoise, mask in zip([252., 277.],[15.25,14.87], [mask4, mask5]):
        ax.plot(xx, np.sqrt(rdnoise**2+xx*exptime + (xx*exptime*0.002)**2 + (0.15*exptime))/exptime, '-',color='0.5',lw=0.5, alpha=0.5)
        yy = data['bg'][mask]
        tmp = np.sqrt(rdnoise**2+yy*exptime + (yy*exptime*0.002)**2 + (0.15*exptime))/exptime
        print  data['sigma_flt'][mask] - tmp
        
    ax.legend(loc=2, frameon=False, numpoints=1, fontsize=9)

    ax.set_xlim([0.,2.5])
    ax.set_ylim([0.038,0.15])
    ax.set_ylabel('$\sigma$', fontsize=fs)
    ax.set_xlabel('Background [e$^-$/s]', fontsize=fs)
    
    plt.show(block=False)

    fig.savefig('flt_noise_vs_background.png', dpi=200, transparent=False)

def plot_tests():
    
    """
    Plots showing the noise distributions in all the test pointings.
    """
    import stwcs
    import scipy
    import scipy.optimize
    import unicorn
    import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator
    import matplotlib.gridspec as gridspec
    import matplotlib.colors as colors
    import matplotlib.cm as cmx
    import drizzlepac
    from drizzlepac import astrodrizzle
    from astropy.table import Table as table
    
        
    # Make plot
    fs=8.
    matplotlib.rc('xtick',labelsize=fs)
    matplotlib.rc('ytick',labelsize=fs)
    fig = unicorn.catalogs.plot_init(xs=6., aspect=1.0, fontsize=fs, left=0.11, right=0.15, top=0.1, bottom=0.1)
    gs = gridspec.GridSpec(2,2,top=0.95, bottom=0.1)
    fig.subplots_adjust(wspace=0.1)
    fig.subplots_adjust(hspace=0.1)
    
    jet = cm = plt.get_cmap('jet')
    cNorm = colors.Normalize(vmin=0.3, vmax=2.5)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
    
    tab_252 = table.read('asn_252.txt', format='ascii')
    tab_277 = table.read('asn_277.txt', format='ascii')
    

    for j,table in enumerate([tab_252,tab_277]):
        
        bg = []
        sigma = []            
        ax = fig.add_subplot(gs[0,j])
        if j == 0:
            exptime = 252.937439
        else:
            exptime = 277.937958
        
        for root in table['files']:
            
            asn = threedhst.utils.ASNFile('{}_asn.fits'.format(root))
            
            if root.startswith('../'): 
                root = root[:-3]+'010'
            
            ref = fits.open('{}_drz_sci.fits'.format(root))
            ref_wcs = stwcs.wcsutil.HSTWCS(ref, ext=0)

            seg = fits.open('{}_drz_seg.fits'.format(root))    
            seg_data = np.cast[np.float32](seg[0].data)
          
            yi, xi = np.indices((1014,1014))
                        
            for exp in asn.exposures:
                
                if root.startswith('../'):
                    if exp == asn.exposures[0]:
                        flt = fits.open('../../14114/PREPARE/{}_flt.fits'.format(exp[:-3]))
                    else:
                        continue
                else:
                    flt = fits.open('{}_flt.fits'.format(exp))
                flt_wcs = stwcs.wcsutil.HSTWCS(flt, ext=1)
                color = scalarMap.to_rgba(flt[0].header['BGCOMP1'])
                blotted_seg = astrodrizzle.ablot.do_blot(seg_data, ref_wcs, flt_wcs, 1, coeffs=True, 
                    interp='nearest', sinscl=1.0, stepsize=10, wcsmap=None)          
                mask = (blotted_seg == 0) & (flt['DQ'].data == 0) & (xi > 10) & (yi > 10) & (xi < 1004) & (yi < 1004)   

                n, bin_edge, patch  = ax.hist(flt[1].data[mask], bins=300, range=[-1,1], color=color,
                     alpha=0.5, histtype='step', normed=True)
                centers = (bin_edge[:-1] + (bin_edge[1:] - bin_edge[:-1])/2)
                flt_coeff, flt_var = scipy.optimize.curve_fit(gauss, centers, n, p0=[np.max(n),0.,0.1])
                sigma.append(flt_coeff[2])
                bg.append(flt[0].header['BGCOMP1'])
                #print flt[1].header['SAMPTIME'], flt[0].header['BGCOMP1'], flt_coeff
                #fit_g = gauss(centers, *flt_coeff)
                #ax.plot(centers, fit_g, color='red', alpha=0.5, linewidth=0.5)
            
        cbar_ax = fig.add_axes([0.88, 0.545, 0.02, 0.41])
        cbar = matplotlib.colorbar.ColorbarBase(cbar_ax,cmap=jet,norm=cNorm,orientation='vertical')
        cbar.set_label('Background [e-/s]', fontsize=fs)       
        ax.set_ylim([0, 6.5])
        ax.set_xlim([-0.5, 0.5])
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        if j == 0:
            ax.set_ylabel('Normalized N', fontsize=fs)
        else:
            ax.set_yticklabels([])
        ax.set_xlabel('e$^-$/s', fontsize=fs)        
    
        ax = fig.add_subplot(gs[1,j])
        ax.scatter(bg, sigma,alpha=0.7, s=10., color='0.3', edgecolors='none', lw=0)
        xx = np.linspace(0.2,2.4,10)
        ax.plot(xx, np.sqrt(15.**2+xx*exptime + (xx*exptime*0.002)**2 + (0.15*exptime))/exptime, '-',color='red',lw=0.5, alpha=0.5)
        #ax.text(1.5,0.06,r'$\sqrt{20. + BG \times\; t}/t$')
        if j == 0:
            fit_coeff, fitvar = scipy.optimize.curve_fit(fit_func1, bg, sigma, p0=[15.])
            fit_curve = fit_func1(xx, *fit_coeff)
        else:
            fit_coeff, fitvar = scipy.optimize.curve_fit(fit_func2, bg, sigma, p0=[15])
            fit_curve = fit_func2(xx, *fit_coeff)
        print fit_coeff
        ax.plot(xx, fit_curve, ':',color='red',lw=0.5, alpha=0.5)
        #ax.text(0.75,0.12,r'$\sqrt{20. + BG \times\; t + %4.2f \times\; t}/t$' % (fit_coeff[0]))
        ax.set_xlim([0.,2.5])
        ax.set_ylim([0.,0.2])
        if j == 0:
            ax.set_ylabel('$\sigma$', fontsize=fs)
        else:
            ax.set_yticklabels([])
        ax.set_xlabel('Background [e$^-$/s]', fontsize=fs)
    plt.show(block=False)
    
    #fig.savefig('noise_test.png'.format(master_root), dpi=200, transparent=False)
        
def gauss(x, *p):
    
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))

def fit_func1(x, *p):
    
    A, = p
    return np.sqrt(A**2+x*252.937439 + (x*252.937439*0.002)**2  + (0.15*252.937439))/252.937439
    

def fit_func2(x, *p):
    
    A, = p
    return np.sqrt(A**2+x*277.937958 + (x*277.937958*0.002)**2 + (0.15*277.937958))/277.937958
    
def fit_func(x, *p):
    
    A, B = p
    return A + x*B