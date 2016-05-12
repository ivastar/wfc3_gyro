import threedhst
import astropy
from astropy.io import fits
from astropy.table import Table as table
import drizzlepac
from drizzlepac import astrodrizzle
import os
import glob
import numpy as np

def prep_brodwin_repeat():
    
    import threedhst
    import threedhst.prep_flt_astrodrizzle as init
    
    os.chdir('/3DHST/Spectra/Work/14114/TEST_HOPR_BRODWIN/PREP/')

    init.prep_direct_grism_pair(direct_asn='BRODWIN-72-260-F160W_asn.fits', grism_asn=None, radec=None,
    	raw_path = '../RAW/', mask_grow=8, scattered_light=False, final_scale = 0.06,
    	skip_direct=False, ACS=False, align_threshold=6)
    
    os.system('rsync -av BRODWIN-72-260-F160W_drz* ../REF/')
    ### make an radec cat from sextractor cat
    os.chdir('/3DHST/Spectra/Work/14114/TEST_HOPR_BRODWIN/REF/')
    wfc3_gyro.prepare.run_sextractor(mosaic='BRODWIN-72-260-F160W_drz_sci.fits',weight='BRODWIN-72-260-F160W_drz_wht.fits')
    threedhst.sex.sexcatRegions('BRODWIN-72-260-F160W_drz_sci.cat','BRODWIN-72-260-F160W_drz_sci.reg', format=1)
    
    os.chdir('/3DHST/Spectra/Work/14114/TEST_HOPR_BRODWIN/PREP/')
    
    init.prep_direct_grism_pair(direct_asn='BRODWIN-72-260-F160W_asn.fits', grism_asn=None, radec='../REF/radec.cat',
    	raw_path = '../RAW/', mask_grow=8, scattered_light=False, final_scale = 0.06,
    	skip_direct=False, ACS=False, align_threshold=6)
        
def prep_brodwin_drifted():
    
    import threedhst.prep_flt_astrodrizzle    
            
    os.system('rsync -av ../RAW/ib2u22*_flt.fits.gz .')
    os.system('gunzip -f ib2u22*_flt.fits.gz')        

    root='BRODWIN-22-051-F160W'
    refimage='../REF/BRODWIN-72-260-F160W_drz_sci.fits'
    REF_CAT='../REF/BRODWIN-72-260-F160W_drz_sci.cat'

    outshifts = 'shifts_{}.txt'.format(root)
    outwcs = 'shifts_{}_wcs.fits'.format(root)

    threedhst.prep_flt_astrodrizzle.subtract_flt_background(root=root)
    
    drizzlepac.tweakreg.TweakReg(root+'_asn.fits', refimage=refimage, updatehdr = True, updatewcs = True, 
        writecat = False, clean = True, verbose = True, runfile = 'tweakreg.log', 
        wcsname = 'TWEAK', headerlet = False, shiftfile = True, outshifts = outshifts, outwcs = outwcs, 
        refcat = REF_CAT, refxcol = 13, refycol = 14, refxyunits = 'degrees', minobj = 5, searchrad = 1000.0, 
        searchunits ='pixels', 
        use2dhist = True, see2dplot = False, separation = 0.5, tolerance = 1.0, xoffset = 0.0, yoffset = 0.0, 
        fitgeometry = 'rscale', residplot = 'No plot', interactive=False, nclip = 3, sigma = 3.0, clobber=True) 
    
    drizzlepac.astrodrizzle.AstroDrizzle(root+'_asn.fits', clean=False, final_pixfrac=1.0, context=False, final_bits=576, resetbits=0, preserve=False, driz_cr_snr='8.0 5.0', driz_cr_scale = '2.5 0.7', wcskey= 'TWEAK')
    

def hopr_brodwin():
    
    import unicorn
    import wfc3_gyro.prepare        
    import threedhst.prep_flt_astrodrizzle as init
        
    os.chdir('/3DHST/Spectra/Work/14114/TEST_HOPR_BRODWIN/REDUCE/')
    unicorn.candels.make_asn_files()        
    RAW_PATH = '../RAW'
    
    prep_brodwin_drifted()
    os.system('rsync -av {}/ib2u72*_flt.fits.gz .'.format(RAW_PATH))
    os.system('gunzip -f ib2u72*_flt.fits.gz')
    
    init.prep_direct_grism_pair(direct_asn='BRODWIN-72-260-F160W_asn.fits', grism_asn=None, radec='../REF/radec.cat',
    	raw_path = '../RAW/', mask_grow=8, scattered_light=False, final_scale = 0.06,
    	skip_direct=False, ACS=False, align_threshold=6)

    for master_root in ['BRODWIN-22-051-F160W','BRODWIN-72-260-F160W']:
        asn = threedhst.utils.ASNFile('{}_asn.fits'.format(master_root))    
        for root in asn.exposures:
            wfc3_gyro.prepare.split_IMA(root=root)
            wfc3_gyro.prepare.make_pointing_asn(root=root, master_root=master_root)
            wfc3_gyro.prepare.subtract_background_reads(root=root, master_root=master_root)
            wfc3_gyro.prepare.fix_cosmic_rays(root=root, master_root=master_root)
            wfc3_gyro.prepare.align_reads(root=root, refimage='../REF/BRODWIN-72-260-F160W_drz_sci.fits', final_scale=0.12,
                master_catalog='../REF/BRODWIN-72-260-F160W_drz_sci.cat', refxcol = 13, refycol = 14)


    for root in ['ib2u72','ib2u22']:
        drizzle_flts(root=root)
    
def drizzle_flts(root = 'ib2u72'):
    
    ### this only makes sense for the repeated good FLT!
    
    flt_files = glob.glob('{}???_flt.fits'.format(root))
    
    for file in flt_files:
        
        flt = fits.open(file, mode='update')
        flt[1].header['MDRIZSKY'] =  flt[0].header['BGCOMP1']
        flt[1].data += flt[1].header['MDRIZSKY']
        flt.flush()
        drizzlepac.astrodrizzle.AstroDrizzle(file, output=file.split('.fits')[0], clean=True, final_pixfrac=0.8, resetbits=0, context=False, preserve=False, skysub = True, skywidth = 0., skystat = '', skylower = None, skyupper = None, skyclip = 0, skylsigma = 0.0, skyusigma = 0.0, skyuser = 'MDRIZSKY', skyfile = '', wcskey = 'TWEAK', driz_separate = False, driz_sep_wcs = False, median = False, blot = False, driz_cr = False, driz_combine = True, final_wht_type = 'IVM', final_kernel = 'square', final_wt_scl = 'exptime', final_fillval = 0,final_bits = 576, final_units = 'cps', final_wcs = True, driz_sep_bits = 576, driz_cr_snr='8.0 5.0', driz_cr_scale = '2.5 0.7')
        flt = fits.open(file, mode='update')
        flt[1].data -= flt[1].header['MDRIZSKY']
        flt.flush()

def noise_measure():
    
    os.chdir('/3DHST/Spectra/Work/14114/TEST_HOPR_BRODWIN/REDUCE')
    root = 'BRODWIN-72-260-F160W'
    for flt in ['ib2u72jxq','ib2u72k0q']:
        determine_noise(root=root,flt=flt)
    root = 'BRODWIN-22-051-F160W'
    for flt in ['ib2u22pnq','ib2u22prq']:
        determine_noise(root=root,flt=flt, drift=True)
    os.chdir('/3DHST/Spectra/Work/14114/TEST_HOPR_COOPER/REDUCE')
    root = 'EGS13034445-05-087-F140W'
    for flt in ['ibt305i1q','ibt305i8q','ibt305igq','ibt305itq']:
        determine_noise(root=root,flt=flt, drift=True)
    root = 'EGS13034445-55-081-F140W'
    for flt in ['ibt355neq','ibt355nlq','ibt355ntq','ibt355odq']:
        determine_noise(root=root,flt=flt)
    
    
def determine_noise(root = 'BRODWIN-72-260-F160W', flt = 'ib2u72k0q', drift=False):
    
    import stwcs
    import scipy
    import scipy.optimize
    import wfc3_gyro.noise
    
    if drift:
        flag=2
    else:
        flag=1
    
    fp = open('sigma_hopr.txt','a')
    fp.write('#file mean sigma bg flag shift\n')
    
    tab = table.read('{}_shifts.txt'.format(flt), format='ascii', 
                names=('file','x','y','rot','scale','x_rms','y_rms'))
    shift = np.sqrt((tab['x'][0]-tab['x'][-1])**2 + (tab['y'][0]-tab['y'][-1])**2)
    
    yi, xi = np.indices((1014,1014))

    ref = fits.open('{}_drz_sci.fits'.format(root))
    ref_wcs = stwcs.wcsutil.HSTWCS(ref, ext=0)

    seg = fits.open('{}_drz_seg.fits'.format(root))    
    seg_data = np.cast[np.float32](seg[0].data)

    ### Measure noise in CALWF3 FLT image
    flt_file = fits.open('{}_flt.fits'.format(flt))
    flt_wcs = stwcs.wcsutil.HSTWCS(flt_file, ext=1)
                
    ### Original FLT
    ### ib2u72k0q_drz_sci.fits ib2u72k0q_flt_drz_sci.fits ib2u72k0q_flt.fits BRODWIN-72-260-F160W_drz_sci.fits
    blotted_seg = astrodrizzle.ablot.do_blot(seg_data, ref_wcs, flt_wcs, 1, coeffs=True, 
        interp='nearest', sinscl=1.0, stepsize=10, wcsmap=None) 
    mask = (blotted_seg == 0) & (flt_file['DQ'].data == 0) & (xi > 10) & (yi > 10) & (xi < 1004) & (yi < 1004)   
    n, bin_edge  = np.histogram(flt_file[1].data[mask], bins=300, range=[-1,1], normed=True)
    centers = (bin_edge[:-1] + (bin_edge[1:] - bin_edge[:-1])/2)
    flt_coeff, flt_var = scipy.optimize.curve_fit(wfc3_gyro.noise.gauss, centers, n, p0=[np.max(n),0.,1.])

    print '{}_flt.fits\t{:7.4f}\t{:7.4f}\t{:7.4f}\t{}\t{:7.4f}\n'.format(flt,flt_coeff[1], np.abs(flt_coeff[2]), flt_file[0].header['BGCOMP1'], flag, shift)
    fp.write('{}_flt.fits\t{:7.4f}\t{:7.4f}\t{:7.4f}\t{}\t{:7.4f}\n'.format(flt,flt_coeff[1], np.abs(flt_coeff[2]), flt_file[0].header['BGCOMP1'], flag, shift))
    
    ### Drizzled FLT
    drz = fits.open('{}_flt_drz_sci.fits'.format(flt))
    drz_wcs = stwcs.wcsutil.HSTWCS(drz, ext=0)
    drz_sh = np.shape(drz[0].data)
    dyi, dxi = np.indices(drz_sh)
    blotted_seg = astrodrizzle.ablot.do_blot(seg_data, ref_wcs, drz_wcs, 1, coeffs=False, 
        interp='nearest', sinscl=1.0, stepsize=10, wcsmap=None)
    mask = (blotted_seg == 0) & (dxi > 10) & (dyi > 10) & (dxi < drz_sh[1]-10) & (dyi < drz_sh[0]-10)
    n, bin_edge = np.histogram(drz[0].data[mask], bins=300, range=[-1.,1.], normed=True)
    centers = (bin_edge[:-1] + (bin_edge[1:] - bin_edge[:-1])/2)
    drz_coeff, drz_var = scipy.optimize.curve_fit(wfc3_gyro.noise.gauss, centers, n, p0=[np.max(n),0.,1.])

    print '{}_flt_drz_sci.fits\t{:7.4f}\t{:7.4f}\t{:7.4f}\t{}\t{:7.4f}\n'.format(flt,drz_coeff[1], np.abs(drz_coeff[2]), flt_file[0].header['BGCOMP1'], flag, shift)
    fp.write('{}_flt_drz_sci.fits\t{:7.4f}\t{:7.4f}\t{:7.4f}\t{}\t{:7.4f}\n'.format(flt,drz_coeff[1], np.abs(drz_coeff[2]), flt_file[0].header['BGCOMP1'], flag, shift))

    ### Reads drizzled together
    drz = fits.open('{}_drz_sci.fits'.format(flt))
    drz_wcs = stwcs.wcsutil.HSTWCS(drz, ext=0)
    drz_sh = np.shape(drz[0].data)
    dyi, dxi = np.indices(drz_sh)
    blotted_seg = astrodrizzle.ablot.do_blot(seg_data, ref_wcs, drz_wcs, 1, coeffs=False, 
        interp='nearest', sinscl=1.0, stepsize=10, wcsmap=None)
    mask = (blotted_seg == 0) & (dxi > 10) & (dyi > 10) & (dxi < drz_sh[1]-10) & (dyi < drz_sh[0]-10)
    n, bin_edge = np.histogram(drz[0].data[mask], bins=300, range=[-1.,1.], normed=True)
    centers = (bin_edge[:-1] + (bin_edge[1:] - bin_edge[:-1])/2)
    drz_reads_coeff, drz_reads_var = scipy.optimize.curve_fit(wfc3_gyro.noise.gauss, centers, n, p0=[np.max(n),0.,1.])

    print '{}_drz_sci.fits\t{:7.4f}\t{:7.4f}\t{:7.4f}\t{}\t{:7.4f}\n'.format(flt,drz_reads_coeff[1], np.abs(drz_reads_coeff[2]), flt_file[0].header['BGCOMP1'], flag, shift)
    fp.write('{}_drz_sci.fits\t{:7.4f}\t{:7.4f}\t{:7.4f}\t{}\t{:7.4f}\n'.format(flt,drz_reads_coeff[1], np.abs(drz_reads_coeff[2]), flt_file[0].header['BGCOMP1'], flag, shift))
                
    fp.close()        
                
    
def prep_cooper_repeat():
    
    import threedhst
    import threedhst.prep_flt_astrodrizzle as init
    
    os.chdir('/3DHST/Spectra/Work/14114/TEST_HOPR_COOPER/PREP/')

    unicorn.candels.make_asn_files()
    
    init.prep_direct_grism_pair(direct_asn='EGS13034445-55-081-F140W_asn.fits', grism_asn=None, radec=None,
    	raw_path = '../RAW/', mask_grow=8, scattered_light=False, final_scale = 0.06,
    	skip_direct=False, ACS=False, align_threshold=6)
    
    os.system('rsync -av EGS13034445-55-081-F140W_drz* ../REF/')
    os.chdir('/3DHST/Spectra/Work/14114/TEST_HOPR_COOPER/REF/')
    wfc3_gyro.prepare.run_sextractor(mosaic='EGS13034445-55-081-F140W_drz_sci.fits', 
        weight='EGS13034445-55-081-F140W_drz_wht.fits')
    threedhst.sex.sexcatRegions('EGS13034445-55-081-F140W_drz_sci.cat','EGS13034445-55-081-F140W_drz_sci.reg', format=1)
    ### make an radec cat from sextractor cat
    
    os.chdir('/3DHST/Spectra/Work/14114/TEST_HOPR_COOPER/PREP/')
    
    init.prep_direct_grism_pair(direct_asn='EGS13034445-55-081-F140W_asn.fits', grism_asn=None, radec='../REF/radec.cat',
    	raw_path = '../RAW/', mask_grow=8, scattered_light=False, final_scale = 0.06,
    	skip_direct=False, ACS=False, align_threshold=6)

def prep_cooper_drifted():
    
    import threedhst.prep_flt_astrodrizzle    
            
    os.system('rsync -av ../RAW/ibt305*_flt.fits.gz .')
    os.system('gunzip -f ibt305*_flt.fits.gz')        

    root='EGS13034445-05-087-F140W'
    refimage='../REF/EGS13034445-55-081-F140W_drz_sci.fits'
    REF_CAT='../REF/EGS13034445-55-081-F140W_drz_sci.cat'

    outshifts = 'shifts_{}.txt'.format(root)
    outwcs = 'shifts_{}_wcs.fits'.format(root)

    threedhst.prep_flt_astrodrizzle.subtract_flt_background(root=root)
    
    drizzlepac.tweakreg.TweakReg(root+'_asn.fits', refimage=refimage, updatehdr = True, updatewcs = True, 
        writecat = False, clean = True, verbose = True, runfile = 'tweakreg.log', 
        wcsname = 'TWEAK', headerlet = False, shiftfile = True, outshifts = outshifts, outwcs = outwcs, 
        refcat = REF_CAT, refxcol = 13, refycol = 14, refxyunits = 'degrees', minobj = 5, searchrad = 1000.0, 
        searchunits ='pixels', 
        use2dhist = True, see2dplot = False, separation = 0.5, tolerance = 1.0, xoffset = 0.0, yoffset = 0.0, 
        fitgeometry = 'rscale', residplot = 'No plot', interactive=False, nclip = 3, sigma = 3.0, clobber=True) 
    
    drizzlepac.astrodrizzle.AstroDrizzle(root+'_asn.fits', clean=False, final_pixfrac=1.0, context=False, final_bits=576, resetbits=0, preserve=False, driz_cr_snr='8.0 5.0', driz_cr_scale = '2.5 0.7', wcskey= 'TWEAK')
    

def hopr_cooper():
    
    import unicorn
    import wfc3_gyro.prepare        
    import threedhst.prep_flt_astrodrizzle as init
        
    os.chdir('/3DHST/Spectra/Work/14114/TEST_HOPR_COOPER/REDUCE/')
    unicorn.candels.make_asn_files()        
    RAW_PATH = '../RAW'
    
    prep_cooper_drifted()
    os.system('rsync -av {}/ibt355*_flt.fits.gz .'.format(RAW_PATH))
    os.system('gunzip -f ibt355*_flt.fits.gz')
    
    init.prep_direct_grism_pair(direct_asn='EGS13034445-55-081-F140W_asn.fits', grism_asn=None, radec='../REF/radec.cat',
    	raw_path = '../RAW/', mask_grow=8, scattered_light=False, final_scale = 0.06,
    	skip_direct=False, ACS=False, align_threshold=6)

    for master_root in ['EGS13034445-05-087-F140W','EGS13034445-55-081-F140W']:
        asn = threedhst.utils.ASNFile('{}_asn.fits'.format(master_root))    
        for root in asn.exposures:
            wfc3_gyro.prepare.split_IMA(root=root)
            wfc3_gyro.prepare.make_pointing_asn(root=root, master_root=master_root)
            wfc3_gyro.prepare.subtract_background_reads(root=root, master_root=master_root)
            wfc3_gyro.prepare.fix_cosmic_rays(root=root, master_root=master_root)
            wfc3_gyro.prepare.align_reads(root=root, refimage='../REF/EGS13034445-55-081-F140W_drz_sci.fits', 
                final_scale=0.12, master_catalog='../REF/EGS13034445-55-081-F140W_drz_sci.cat', 
                refxcol = 13, refycol = 14)

    for root in ['ibt305','ibt355']:
        drizzle_flts(root=root)
    
def plot_hopr_noise():
    
    import unicorn
    import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib.figure import Figure
    import matplotlib.gridspec as gridspec
    
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
    ax1 = fig.add_subplot(gs1[0,0], ylim=[0.0,0.12], xlim=[0.1,2.5])
    ax2 = fig.add_subplot(gs1[0,1], ylim=[0.0,0.35], xlim=[0.1,2.5])

    for ax, exptime in zip([ax1,ax2],[650, 203]):
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


        ax.plot(xx, noise_21_11, '-',color='#b10026',lw=1, label=r'RDN = 15$\times\sqrt{11}$ e$^{-}$')
        ax.plot([0.,3.], [rd_21_11,rd_21_11], ':',color='#b10026',lw=1)

        ax.plot(xx, noise_21_4, '-',color='#fc4e2a',lw=1, label=r'RDN = 21$\times\sqrt{4}$ e$^{-}$')
        ax.plot([0.,3.], [rd_21_4,rd_21_4], ':',color='#fc4e2a',lw=1)

        ax.plot(xx, noise_21_2, '-',color='#feb24c',lw=1, label=r'RDN = 21$\times\sqrt{2}$ e$^{-}$')
        ax.plot([0.,3.], [rd_21_2,rd_21_2], ':',color='#feb24c',lw=1)

        ax.plot(xx, noise_15, '-',color='black',lw=1, label=r'RDN = 15 e$^{-}$')
        ax.plot([0.,3.], [rd_15,rd_15], ':',color='black',lw=1)
        ax.plot(xx, bg_15, '--',color='black',lw=1)

        #ax.text(2.0,0.105,'full noise model', rotation=11, va='center', ha='center', fontsize=fs, backgroundcolor='white', alpha=0.7)
        #ax.text(2.0,0.085,'sky noise', rotation=12, va='center', ha='center', fontsize=fs, backgroundcolor='white', alpha=0.7)
        #ax.text(2.0,0.053,'read noise',va='center', ha='center', fontsize=fs, backgroundcolor='white', alpha=0.7)
        ax.legend(loc='upper left', fontsize=fs, frameon=False)
        ax.set_ylabel('$\sigma [e^-/s]$', fontsize=fs)
        ax.set_xlabel('Background [e$^-$/s]', fontsize=fs)
    
    ### BRODWIN DATA
    ax1.set_title('BRODWIN 600 sec NSAMP=8 SPARS50 F160W')
    b_data = table.read('/3DHST/Spectra/Work/14114/TEST_HOPR_BRODWIN/REDUCE/sigma_hopr.txt', format='ascii')
    ff = (b_data['flag'] == 1)
    ax1.plot(b_data['bg'][ff], b_data['sigma'][ff],'o', color='black')
    ff = (b_data['flag'] == 2)
    ax1.plot(b_data['bg'][ff], b_data['sigma'][ff],'o', mec='red', mfc='None', lw=2)
    for ii in range(len(b_data)):
        if (b_data['flag'][ii] == 2) & ('_drz_sci.fits' in b_data['file'][ii]) & ('flt_drz_sci.fits' not in b_data['file'][ii]):
            ax1.plot(b_data['bg'][ii], b_data['sigma'][ii],'o',color='red')
            ax1.text(b_data['bg'][ii]+0.05, b_data['sigma'][ii], '{:4.1f} pix'.format(b_data['shift'][ii]), ha='left', va='center', fontsize=6)    
            print b_data['file'][ii]

    ### COOPER DATA
    ax2.set_title('COOPER 200 sec NSAMP=11 SPARS25 F140W')
    c_data = table.read('/3DHST/Spectra/Work/14114/TEST_HOPR_COOPER/REDUCE/sigma_hopr.txt', format='ascii')
    ff = (c_data['flag'] == 1)
    ax2.plot(c_data['bg'][ff], c_data['sigma'][ff],'o',color='black')
    ff = (c_data['flag'] == 2)
    ax2.plot(c_data['bg'][ff], c_data['sigma'][ff],'o',mec='red', mfc='None', lw=2)
    for ii in range(len(c_data)):
        print c_data['file'][ii]
        if (c_data['flag'][ii] == 2) & ('_drz_sci.fits' in c_data['file'][ii]) & ('flt_drz_sci.fits' not in c_data['file'][ii]):
            ax2.plot(c_data['bg'][ii], c_data['sigma'][ii],'o',color='red')
            ax2.text(c_data['bg'][ii]+0.05, c_data['sigma'][ii], '{:4.1f} pix'.format(c_data['shift'][ii]), ha='left', va='center', fontsize=6)
            print c_data['file'][ii]
       
    plt.show(block=False)
    
    #fig.savefig('sources_of_noise.pdf', dpi=200, transparent=False)
    #fig.savefig('sources_of_noise.png', dpi=200, transparent=False)
    
def plot_growth_curve():
    
    os.chdir('/3DHST/Spectra/Work/14114/TEST_HOPR_COOPER/REDUCE/')
    
    import unicorn
    import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib.figure import Figure
    import matplotlib.gridspec as gridspec
    
    top = 0.075
    bottom = 0.1
    left = 0.15
    fig = unicorn.catalogs.plot_init(xs=5,aspect=1.0, left=left, right=0.1, bottom=bottom, top=top, NO_GUI=False)
    fig.subplots_adjust(wspace=0.15)
    fig.subplots_adjust(hspace=0.1)
    fs = 8
    matplotlib.rc('xtick',labelsize=fs)
    matplotlib.rc('ytick',labelsize=fs)
        
    gs1 = gridspec.GridSpec(1,1)
    ax1 = fig.add_subplot(gs1[0,0], ylim=[0.1,5.], xlim=[0.1,1.5])
    
    guided = glob.glob('ibt355*_drz_sci.apertures.F160W.v4.0.dat')
    dash = glob.glob('ibt305*_drz_sci.apertures.F160W.v4.0.dat')
    
    for file in guided:
        data = table.read(file, format='ascii')
        ax1.plot(data['radius_pix']*0.128, data['sigma'], 'o', color='black', alpha=0.8, ms=2.0)
        
    for file in dash:
        data = table.read(file, format='ascii')
        ax1.plot(data['radius_pix']*0.128, data['sigma'], 'o', color='red', alpha=0.8, ms=2.0)
    
    #ax1.set_yscale('log')
        
    plt.show(block=False)