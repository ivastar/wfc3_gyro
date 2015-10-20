"""
Routines to analyze the data from the COSMOS program.
"""
import threedhst
import astropy
from astropy.io import fits
from astropy.table import Table as table
import drizzlepac
from drizzlepac import astrodrizzle
import os
import glob
import numpy as np

def split_IMA(root='icxe15wwq', PATH='../RAW/'):
    
    """
    Make FLTs from the individual reads.
    """
    
    FLAT_F140W = 1 ### FLATCORR=COMPLETE for all of these
    
    ima = fits.open(PATH+root+'_ima.fits.gz')
    flt = fits.open(root+'_flt.fits')
    
    NSAMP = ima[0].header['NSAMP']
    sh = ima['SCI',1].shape
    
    cube = np.zeros((NSAMP, sh[0], sh[1]))
    dq = np.zeros((NSAMP, sh[0], sh[1]), dtype=np.int)
        
    time = np.zeros(NSAMP)
    for i in range(NSAMP):
        cube[NSAMP-1-i, :, :] = ima['SCI',i+1].data*ima['TIME',i+1].header['PIXVALUE']/FLAT_F140W
        dq[NSAMP-1-i, :, :] = ima['DQ',i+1].data
        time[NSAMP-1-i] = ima['TIME',i+1].header['PIXVALUE']
    
    diff = np.diff(cube, axis=0)
    dt = np.diff(time)
        
    readnoise_2D = np.zeros((1024,1024))
    readnoise_2D[512: ,0:512] += ima[0].header['READNSEA']
    readnoise_2D[0:512,0:512] += ima[0].header['READNSEB']
    readnoise_2D[0:512, 512:] += ima[0].header['READNSEC']
    readnoise_2D[512: , 512:] += ima[0].header['READNSED']
    readnoise_2D = readnoise_2D**2
    
    for j in range(1, NSAMP-1):
        print '{}_{:02d}_flt.fits'.format(root,j)
        sci = diff[j,:,:]
        exptime = dt[j]
        var = readnoise_2D + sci        
        err = np.sqrt(var)/exptime
        
        flt[0].header['EXPTIME'] = exptime
        flt['SCI'].data = sci[5:-5,5:-5]/exptime
        flt['ERR'].data = err[5:-5,5:-5]
        flt['DQ'].data = dq[j][5:-5,5:-5]
        flt['SAMP'].data = np.zeros((1014,1014)) + 1.
        flt['TIME'].data = np.zeros((1014,1014)) + exptime
        flt[0].header['IMA2FLT'] = (1, 'FLT {} extracted from IMA file'.format(j))

        print 'Writing {}_{:02d}_flt.fits'.format(root,j)
        flt.writeto('{}_{:02d}_flt.fits'.format(root,j), clobber=True)
        
def make_pointing_asn(root='icxe15wwq', master_root='icxe15010'):
    
    master_asn = fits.open('{}_asn.fits'.format(master_root))
        
    files = glob.glob('{}_*_flt.fits'.format(root))
    nrows = len(files)
    
    #### Primary HDU
    hdu = master_asn[0].copy()
    tbhdu = fits.new_table(master_asn[1].columns, nrows=nrows+1, fill=True)
    for i in range(nrows):
        tbhdu.data[i] = (files[i].split('_flt')[0].upper(), 'EXP-DTH', True)
     
    tbhdu.data[i+1] = (root.upper(), 'PROD-DTH', True)
            
    tbhdu.header = master_asn[1].header.copy()
    tbhdu.header.update('ASN_ID',root)
    tbhdu.header.update('ASN_TAB','{}_asn.fits'.format(root))
    
    #### Create HDUList and write it to output file
    out_fits = fits.HDUList([hdu,tbhdu])
    
    if 'EXTEND' not in hdu.header.keys():
        hdu.header.update('EXTEND', True, after='NAXIS')
                
    print 'Writing {}_asn.fits'.format(root)
    out_fits.writeto('{}_asn.fits'.format(root), clobber=True)


def subtract_background_reads(root='icxe15wwq', master_root='icxe15010'):
    
    import threedhst
    import stwcs
    import scipy
    import scipy.optimize
        
    asn = threedhst.utils.ASNFile('{}_asn.fits'.format(root))
    
    #print 'Read files...'
    ref = fits.open('{}_drz_sci.fits'.format(master_root))
    ref_wcs = stwcs.wcsutil.HSTWCS(ref, ext=0)

    seg = fits.open('{}_drz_seg.fits'.format(master_root))    
    seg_data = np.cast[np.float32](seg[0].data)
          
    yi, xi = np.indices((1014,1014))
    NCOMP=6
    bg_components = np.ones((NCOMP,1014,1014))
    bg_components[1,:,:] = (xi-507)/507.
    bg_components[2,:,:] = (yi-507)/507.
    bg_components[3,:,:] = ((xi-507)/507.)**2
    bg_components[4,:,:] = ((yi-507)/507.)**2
    bg_components[5,:,:] = (xi-507)*(yi-507)/507.**2
            
    bg_flat = bg_components.reshape((NCOMP,1014**2))            

    #### Loop through FLTs, blotting reference and segmentation
    models = []
    for exp in asn.exposures:
        flt = fits.open('{}_flt.fits'.format(exp)) #, mode='update')
        flt_wcs = stwcs.wcsutil.HSTWCS(flt, ext=1)
        
        if exp == asn.exposures[0]:
            print 'Segmentation image: {}_blot.fits'.format(exp)
            blotted_seg = astrodrizzle.ablot.do_blot(seg_data, ref_wcs, flt_wcs, 1, coeffs=True, interp='nearest', sinscl=1.0, stepsize=10, wcsmap=None)         
            
        mask = (blotted_seg == 0) & (flt['DQ'].data == 0) & (flt[1].data > -1) & (xi > 10) & (yi > 10) & (xi < 1004) & (yi < 1004)
        mask &= (flt[1].data < 5*np.median(flt[1].data[mask]))
        data_range = np.percentile(flt[1].data[mask], [2.5, 97.5])
        mask &= (flt[1].data >= data_range[0]) & (flt[1].data <= data_range[1])
        data_range = np.percentile(flt[2].data[mask], [0.5, 99.5])
        mask &= (flt[2].data >= data_range[0]) & (flt[2].data <= data_range[1])
        
        ### Least-sq fit for component normalizations
        data = flt[1].data[mask].flatten()
        wht = (1./flt[2].data[mask].flatten())**2
        templates = bg_flat[:, mask.flatten()]
        p0 = np.zeros(NCOMP)
        p0[0] = np.median(data)
        obj_fun = threedhst.grism_sky.obj_lstsq
        popt = scipy.optimize.leastsq(obj_fun, p0, args=(data, templates, wht), full_output=True, ftol=1.49e-8/1000., xtol=1.49e-8/1000.)
        xcoeff = popt[0]
        model = np.dot(xcoeff, bg_flat).reshape((1014,1014))
        models.append(model)
        
        # add header keywords of the fit components
        flt = fits.open('{}_flt.fits'.format(exp), mode='update')
        flt[1].data -= model
        for i in range(NCOMP):
            if 'BGCOMP%d' %(i+1) in flt[0].header:
                flt[0].header['BGCOMP%d' %(i+1)] += xcoeff[i]
            else:
                flt[0].header['BGCOMP%d' %(i+1)] = xcoeff[i]                
        
        flt.flush()
        coeff_str = '  '.join(['%.4f' %c for c in xcoeff])
        print 'Background subtraction, {}_flt.fits:\n\n  {}'.format(exp, coeff_str)

def align_reads(root='icxe15wwq', threshold=3, final_scale=0.06, refimage='../REF/cosmos-wide_ACS.fits', master_catalog='../REF/IPAC_ACS.fits'):
    
    from drizzlepac import tweakreg

    asn = threedhst.utils.ASNFile('{}_asn.fits'.format(root))

    catfile = '{}.catfile'.format(root)
    fp = open(catfile,'w')

    for exp in asn.exposures:
        se = threedhst.sex.SExtractor()
        se.aXeParams()
        se.copyConvFile()
        se.options['CHECKIMAGE_TYPE'] = 'NONE'
        se.options['FILTER']    = 'Y'
        se.options['WEIGHT_IMAGE'] = '{}_flt.fits[1]'.format(exp)
        se.options['WEIGHT_TYPE'] = 'NONE'
        
        #
        se.params['X_IMAGE'] = True; se.params['Y_IMAGE'] = True
        se.params['MAG_AUTO'] = True
        #
        se.options['CATALOG_NAME'] = '{}_flt.cat'.format(exp)
        se.options['DETECT_THRESH'] = '{}'.format(threshold)
        se.options['ANALYSIS_THRESH'] = '{}' .format(threshold)
        se.options['DETECT_MINAREA'] = '10'
        #
        se.sextractImage('{}_flt.fits[0]'.format(exp))
        threedhst.sex.sexcatRegions('{}_flt.cat'.format(exp), '{}_flt.reg'.format(exp), format=1)
        
        line = '{0}_flt.fits {0}_flt.cat\n'.format(exp)
        fp.write(line)
        
    fp.close()
        
    #### Make room for TWEAK wcsname
    for exp in asn.exposures:
        threedhst.prep_flt_astrodrizzle.clean_wcsname(flt='{}_flt.fits'.format(exp), wcsname='TWEAK')    
        
    tweakreg.TweakReg('{}_asn.fits'.format(root), refimage=refimage, updatehdr=True, updatewcs=True, catfile=catfile, xcol=2, ycol=3, xyunits='pixels', refcat=master_catalog, refxcol = 5, refycol = 6, refxyunits='degrees', shiftfile=True, outshifts='{}_shifts.txt'.format(root), outwcs='{}_wcs.fits'.format(root), searchrad=10, tolerance=1., minobj = 5, wcsname='TWEAK', interactive=False, residplot='No plot', see2dplot=True, clean=True, headerlet=True, clobber=True)
    
    drizzlepac.astrodrizzle.AstroDrizzle('{}_asn.fits'.format(root), output=root, clean=True, final_scale=final_scale, 
        final_pixfrac=0.8, context=False, resetbits=4096, final_bits=576, driz_sep_bits=576, 
        preserve=False, driz_cr_snr='8.0 5.0', driz_cr_scale = '2.5 0.7')
        
    
def prep_FLTs(root='icxe15010', REF_IMAGE='../REF/cosmos-wide_ACS.fits', REF_CAT='../REF/IPAC_ACS.fits'):
    
    import threedhst.prep_flt_astrodrizzle    
    
    outshifts = 'shifts_{}.txt'.format(root)
    outwcs = 'shifts_{}_wcs.fits'.format(root)
    
    threedhst.prep_flt_astrodrizzle.subtract_flt_background(root='icxe15010')
    
    drizzlepac.tweakreg.TweakReg(root+'_asn.fits', refimage=REF_IMAGE, exclusions = '', updatewcs = True, 
        writecat = False, clean = True, verbose = True, runfile = 'tweakreg.log', updatehdr = True, 
        wcsname = 'WCSNEW', headerlet = False, shiftfile = True, outshifts = outshifts, outwcs = outwcs, 
        catfile = '', refcat = REF_CAT, refxcol = 5, refycol = 6, refxyunits = 'degrees', rfluxcol = 2, 
        rfluxunits = 'mag', refnbright = 20000, minobj = 5, searchrad = 100.0, searchunits = 'pixels', 
        use2dhist = True, see2dplot = True, separation = 0.5, tolerance = 1.0, xoffset = 0.0, yoffset = 0.0, 
        fitgeometry = 'rscale', residplot = 'both', nclip = 5, sigma = 3.0, clobber=True) 
    
    drizzlepac.astrodrizzle.AstroDrizzle(root+'_asn.fits', clean=True, final_scale=None, final_pixfrac=1.0, context=False, final_bits=576, preserve=False, driz_cr_snr='5.0 4.0', driz_cr_scale = '2.5 0.7')
    
def run_orbit(master_root='icxe15010', RAW_PATH = '../RAW/'):
    
    asn = threedhst.utils.ASNFile('{}_asn.fits'.format(master_root))
    
    for root in asn.exposures:
        os.system('rsync -av {}/{}_flt.fits.gz .'.format(RAW_PATH, root))
        os.system('gunzip -f {}_flt.fits.gz'.format(root))
        
    prep_FLTs(root=master_root)
    
    for root in asn.exposures:
        split_IMA(root=root)
        make_pointing_asn(root=root, master_root=master_root)
        subtract_background_reads(root=root, master_root=master_root)
        align_reads(root=root)
