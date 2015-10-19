"""
Routines for the gyro proposal. Some of these depend on unicorn.
"""

import astropy
import astropy.io.fits as fits
import os
import glob
#import my_python
import pylab
import numpy as np

def find_SPARS(spars='SPARS'):
    
    os.chdir('/3DHST/Spectra/Work/HOPR/GYRO/')
    
    files = glob.glob('*_flt.fits.gz')
    
    spars_roots = []
    
    for file in files:
        tmp = fits.open(file)
        if spars in tmp[0].header['SAMP_SEQ']:
            print file, tmp[0].header['SAMP_SEQ']
            spars_roots.append(file.split('_')[0])
            
    return spars_roots
    

def measure_offset(root = 'iacr12afq', search_r = 0, show=False, write_to_list=False):
    
    from pyraf import iraf
    import threedhst
    from scipy.spatial import cKDTree as KDT
    import numpy as np
    import matplotlib
    import matplotlib.pyplot as plt
    import unicorn
    
    """
    iacr11tiq: ~ 2 arcmin offset, skip that
    iacr13o1q: large move
    """
    
    if root in ['iacr11tiq', 'iacr13o1q','ia2103t8q', 'ia2103taq','ichx02e5q', 'ichx02e7q', 'ichx02eaq', 'ichx02edq', 'ichx02ehq', 'ichx02ekq','ichx52caq', 'ichx52ccq', 'ichx52cfq', 'ichx52chq', 'ichx52clq', 'ichx52coq','ib6308i9q','ib6308iaq','ib6308idq','ib6308ifq']:
        return float('nan')
    
    if os.path.exists(root+'_ima.fits.gz'):
        ima = fits.open(root+'_ima.fits.gz')
    else:
        print 'NO SUCH FILE: {}_ima.fits.gz'.format(root)
        return float('nan')
        
    NSAMP = ima[0].header['NSAMP']

    ### if fewer than 5 samples, this is not going to work 
    if NSAMP < 5:
        return float('nan')
        
    EXPTIME = ima[0].header['EXPTIME']
    if ima[0].header['SAMP_SEQ'] in ['SPARS100','STEP50','STEP100','STEP200']:
        DELTA_TIME = ima['SCI',1].header['SAMPTIME'] - ima['SCI',2].header['SAMPTIME']
    else:
        DELTA_TIME = ima['SCI',1].header['SAMPTIME'] - ima['SCI',NSAMP-3].header['SAMPTIME']
    print NSAMP, EXPTIME
    
    ### set a multiplicative search radius depending on exposure length
    if search_r == 0:
        search_r = 2.5*(DELTA_TIME/25.) # maximum search radius is 2.5 pixels in 25 seconds drift
    print search_r
    """
    if EXPTIME/NSAMP < 30:
        serach_r = 1.
    elif (EXPTIME/NSAMP < 60) & (EXPTIME/NSAMP > 30.):
        search_r = 2.
    else:
        search_r = 4.
    """
    
    ### first difference image
    if ima[0].header['SAMP_SEQ']  in ['SPARS100', 'STEP50','STEP100','STEP200']:
        diff_first = ima['SCI',2].data*ima['SCI',2].header['SAMPTIME'] - ima['SCI',3].data*ima['SCI',3].header['SAMPTIME']
    else:
        diff_first = ima['SCI',NSAMP-3].data*ima['SCI',NSAMP-3].header['SAMPTIME'] - ima['SCI',NSAMP-2].data*ima['SCI',NSAMP-2].header['SAMPTIME']
    first = fits.PrimaryHDU(data=diff_first, header=ima[0].header)
    first.header['EXTVER'] = 0
    first.writeto('tmp_first.fits',clobber=True)
    
    #first_err = fits.PrimaryHDU(data=ima['ERR',NSAMP-3].data*ima['SCI',NSAMP-3].header['SAMPTIME']-ima['ERR',NSAMP-2].data*ima['SCI',NSAMP-2].header['SAMPTIME'], header=ima[0].header)
    #first_err.header['EXTVER'] = 0
    #first_err.writeto('tmp_first_err.fits', clobber=True)
    
    ### last difference image
    diff_last = ima['SCI',1].data*ima['SCI',1].header['SAMPTIME']-ima['SCI',2].data*ima['SCI',2].header['SAMPTIME']
    last = fits.PrimaryHDU(data=diff_last, header=ima[0].header)
    last.header['EXTVER'] = 0
    last.writeto('tmp_last.fits',clobber=True)
    #last_err = fits.PrimaryHDU(data=ima['ERR',1].data*ima['SCI',1].header['SAMPTIME']-ima['ERR',2].data*ima['SCI',2].header['SAMPTIME'], header=ima[0].header)
    #last_err.header['EXTVER'] = 0
    #last_err.writeto('tmp_last_err.fits', clobber=True)

    
    ### Run sextractor on each of the tmp images
    s = threedhst.sex.SExtractor()
    s.aXeParams()
    s.options['DETECT_MINAREA'] = '10'
    if root in ['ib8c7hnyq','ib8c7hnxq','ibt305i8q', 'ibfup2qpq','ib5x5cqfq']: ### very empty pointings 'ib5x5cqbq', 'ib5x5cqfq', 'ib5x5fq4q', 'ib5x5fq6q'
        s.options['DETECT_THRESH']    = '3.0'
    else:
        s.options['DETECT_THRESH']    = '5.0'         
    s.options['ANALYSIS_THRESH']  = '3.0'
    s.options['CHECKIMAGE_TYPE'] = 'NONE'
    s.options['WEIGHT_TYPE']     = 'NONE'
    s.options['FILTER']    = 'Y'
    

    s.options['WEIGHT_IMAGE']    = 'tmp_last_err.fits[0]'
    s.sextractImage('tmp_last.fits')
    os.system('mv test.cat tmp_last.cat')
    threedhst.sex.sexcatRegions('tmp_last.cat', 'tmp_last.reg', format=1)
    
    s.options['WEIGHT_IMAGE']    = 'tmp_first_err.fits[0]'
    s.sextractImage('tmp_first.fits')
    os.system('mv test.cat tmp_first.cat')
    threedhst.sex.sexcatRegions('tmp_first.cat', 'tmp_first.reg', format=1)

    s_last = astropy.io.ascii.read('tmp_last.cat', format='sextractor') 
    mask_last = (s_last['X_IMAGE'] > 20) & (s_last['Y_IMAGE'] > 20) & (s_last['X_IMAGE'] < 1000) & (s_last['Y_IMAGE'] < 1000)
    coords_last = np.empty((np.sum(mask_last), 2))
    coords_last[:, 0] = np.array(s_last['X_IMAGE'][mask_last]).ravel()
    coords_last[:, 1] = np.array(s_last['Y_IMAGE'][mask_last]).ravel()
    kdt_last = KDT(coords_last)     
     
    s_first = astropy.io.ascii.read('tmp_first.cat', format='sextractor') 
    mask_first = (s_first['X_IMAGE'] > 50) & (s_first['Y_IMAGE'] > 50) & (s_first['X_IMAGE'] < 950) & (s_first['Y_IMAGE'] < 950)
    coords_first = np.empty((np.sum(mask_first), 2))

    coords_first[:, 0] = np.array(s_first['X_IMAGE'][mask_first]).ravel()
    coords_first[:, 1] = np.array(s_first['Y_IMAGE'][mask_first]).ravel()
    kdt_first = KDT(coords_first)     
    
    dists = []
    f_x = []
    f_y = []
    l_x = []
    l_y = []
    print np.sum(mask_first), np.sum(mask_last)
    
    for i, [x,y] in enumerate(coords_last):
        dist, idx = kdt_first.query([x,y], k=1, distance_upper_bound=search_r)
        if not np.isinf(dist):
            dists.append(dist)
            f_x.append(kdt_first.data[idx][0])
            f_y.append(kdt_first.data[idx][1])
            l_x.append(x)
            l_y.append(y)
   
    dx = np.median(np.array(l_x)-np.array(f_x))
    dy = np.median(np.array(l_y)-np.array(f_y))
    dd = np.sqrt(dx**2 + dy**2)    
            
    top = 0.055
    bottom = 0.09
    left = 0.1
    fig = unicorn.catalogs.plot_init(xs=10,aspect=0.35, left=left, right=0.05, bottom=bottom, top=top, NO_GUI=False)
    fig.subplots_adjust(wspace=0.1)
    fig.subplots_adjust(hspace=0.1)
    fs = 8
    matplotlib.rc('xtick',labelsize=fs)
    matplotlib.rc('ytick',labelsize=fs)
    
    a1 = plt.subplot(1, 3, 1)
    a1.imshow(diff_first, vmin=0, vmax=100., origin='lower', cmap=pylab.cm.Greys_r)
    a1.scatter(f_x, f_y, facecolors='none', edgecolors='#fc8d59', s=15)
    if ima[0].header['SAMP_SEQ'] in ['SPARS100','STEP50','STEP100','STEP200']:
        a1.text(0, 0, 't = {:5.1f} [sec]'.format(ima['SCI',2].header['SAMPTIME']), ha='left', va='bottom', 
            backgroundcolor='1.0', fontsize=fs)
    else:
        a1.text(0, 0, 't = {:5.1f} [sec]'.format(ima['SCI',NSAMP-3].header['SAMPTIME']), ha='left', va='bottom', 
            backgroundcolor='1.0', fontsize=fs)
    a1.set_xlim([-20,1050])
    a1.set_ylim([-20,1050])

    a2 = plt.subplot(1, 3, 2)
    a2.imshow(diff_last, vmin=0, vmax=100., origin='lower', cmap=pylab.cm.Greys_r)
    a2.scatter(l_x, l_y, facecolors='none', edgecolors='#d7301f', s=15)
    a2.set_title('{}: {}'.format(root, ima[0].header['SAMP_SEQ']))
    a2.text(0, 0, 't= {:5.1f} [sec]'.format(ima['SCI',1].header['SAMPTIME']), ha='left', va='bottom', 
        backgroundcolor='1.0', fontsize=fs)
    a2.set_xlim([-20,1050])
    a2.set_ylim([-20,1050])

    a3 = plt.subplot(1, 3, 3)
    a3.plot(f_x, f_y, 'o', color='#fc8d59', markeredgecolor='#fc8d59', ms=3)
    a3.plot(l_x, l_y, 'o', color='#d7301f', markeredgecolor='#d7301f', ms=3)
    a3.quiver(f_x, f_y, np.array(l_x)-np.array(f_x), np.array(l_y)-np.array(f_y), width=0.8, linewidth=0.3, units='dots', headwidth=3, minshaft=1, headlength=5)
    a3.set_xlim([-20,1050])
    a3.set_ylim([-20,1050])
    a3.text(0,100,'$\Delta$ t = {:8.2f} [sec]'.format(DELTA_TIME), ha='left', va='bottom', 
        backgroundcolor='1.0', fontsize=fs)
    a3.text(0,50,'$\Delta$ d = {:8.2f} [pix]'.format(np.median(dists)), ha='left', va='bottom', 
        backgroundcolor='1.0', fontsize=fs)
    a3.text(0,0,'v = {:8.4} pix per 100 sec'.format(100.*np.median(dists)/DELTA_TIME), ha='left', 
        va='bottom', backgroundcolor='1.0', fontsize=fs)
    
    
    print root, dd, dx, dy, DELTA_TIME, (dd*100./DELTA_TIME)
    
    if show:
        plt.show(block=False)
    fig.savefig('{}_gyro_shifts_{}.png'.format(root,ima[0].header['SAMP_SEQ']),dpi=100,transparent=False)

    if write_to_list:
        if not os.path.exists('gyro_shifts.txt'):
            fp = open('gyro_shifts.txt','w')
            fp.write('# root SAMP distance dx dy d_time rate dd\n')
        else:
            fp = open('gyro_shifts.txt','a')
        
        fp.write('{}\t{}\t{:10.3f}\t{:10.3f}\t{:10.3f}\t{:10.3f}\t{:10.5f}\n'.format(root, ima[0].header['SAMP_SEQ'], dd, dx, dy, DELTA_TIME, (np.median(dists)/DELTA_TIME), dd/DELTA_TIME))
        fp.close()
     
    return np.median(dists)
     
def measure_all(delete=True):
    
    os.chdir('/3DHST/Spectra/Work/HOPR/GYRO/')
    
    if delete:
        os.system('rm *gyro_shifts*png')
        os.system('rm gyro_shifts.txt')
    
    for seq in ['SPARS25', 'SPARS50', 'SPARS100']:
        roots = find_SPARS(spars=seq)
        
        for root in roots:
            dd = measure_offset(root=root, show=False)
            if (not np.isnan(dd)) & (dd > 0.5):
                measure_offset(root=root, search_r=dd*1.2, show = False, write_to_list=True)
        #
    measure_offset(root='ib2u22prq', search_r=10., write_to_list=True)
    measure_offset(root='ib2u22pnq', search_r=18., write_to_list=True)
    
    for root in ['ib5x5cqbq', 'ib5x5cqfq', 'ib5x5fq4q', 'ib5x5fq6q']:
        measure_offset(root=root, search_r=30., write_to_list=True)
    
def measure_step(delete=True):
    
    if delete:
        os.system('rm *gyro_shifts*STEP*png')
    
    roots = find_SPARS(spars='STEP50')
        
    for root in roots:
        measure_offset(root=root, show=False, write_to_list=True)
           
            
def make_hist_plot():
    
    import matplotlib
    import astropy.io.ascii as ascii
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MultipleLocator, FormatStrFormatter,AutoMinorLocator
    
    matplotlib.rc('xtick',labelsize=10)
    matplotlib.rc('ytick',labelsize=10)
    
    #os.chdir('/3DHST/Spectra/Work/HOPR/GYRO/')

    data = ascii.read('gyro_shifts.txt')
    
    n_25, bins = np.histogram(data['rate'][data['SAMP'] == 'SPARS25']*25., bins=60, range=[0.,6.])
    n_50, bins = np.histogram(data['rate'][data['SAMP'] == 'SPARS50']*25., bins=60, range=[0.,6.])
    n_100, bins = np.histogram(data['rate'][data['SAMP'] == 'SPARS100']*25., bins=60, range=[0.,6.])
    step_50, bins = np.histogram(data['rate'][data['SAMP'] == 'STEP50']*25., bins=60, range=[0.,6.])
    
    
    centers = 0.5*(bins[1:]+bins[:-1])
    width = (bins[1:]-bins[:-1])
    
    fig = plt.figure(1, figsize=(6,5))
    ax = fig.add_subplot(111)

    ax.fill_between([0.,0.5], 0, 30, color='0.9', alpha=0.9)
    ax.bar(centers, n_100, width, align='center', edgecolor='white', color='#b30000', label='SPARS100 (N={})'.format(np.sum(n_100)))
    ax.bar(centers, n_50, width, bottom=(n_100), align='center', edgecolor='white', color='#e34a33', label='SPARS50 (N={})'.format(np.sum(n_50)))
    ax.bar(centers, n_25, width, bottom=(n_100+n_50), align='center', edgecolor='white', color='#fc8d59', label='SPARS25 (N={})'.format(np.sum(n_25)))
    ax.bar(centers, step_50, width, bottom=(n_100+n_50+n_25), align='center', edgecolor='white', color='#fdcc8a', label='STEP50 (N={})'.format(np.sum(step_50)))
    
    test_val = data['rate'][data['root'] == 'ib2u22prq']
        
    ax.arrow(test_val*25., 17.1, 0., -1., head_length=0.5, head_width=0.06, fc='red', ec='red')
    ax.text(0., 18.0, 'Fig. 1 Example', ha='left', va='center', fontsize=9)
    
    ax.arrow(0., 21., 0.5, 0., head_length=0.05, head_width=0.30, fc='black', ec='black', length_includes_head=True)
    ax.arrow(0.5, 21., -0.5, 0., head_length=0.05, head_width=0.30, fc='black', ec='black', length_includes_head=True)
    ax.text(0.0, 22.5, 'Expected Drift', ha='left', va='center')

    ax.set_xlim([-0.15, 3])
    ax.set_ylim([0., 25])
    ax.set_xlabel('Drift in Pixels per 25 Seconds')
    ax.set_ylabel('N')
    minorLocator   = AutoMinorLocator()
    ax.xaxis.set_minor_locator(minorLocator)
    ax.tick_params(which='minor', length=2)    
    plt.show(block=False)

    #ax.legend(loc=1, frameon=False, prop={'size':9}, borderaxespad=1., labelspacing=1.1)
    
    plt.show(block=False)
    
    plt.savefig('drift_hist.png',dpi=100,transparent=False)
    
def make_hist_plot_bw():
    
    import matplotlib
    import astropy.io.ascii as ascii
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MultipleLocator, FormatStrFormatter,AutoMinorLocator
    from matplotlib import rcParams

    minorLocator   = MultipleLocator(10)
    
    matplotlib.rc('xtick',labelsize=10)
    matplotlib.rc('ytick',labelsize=10)
    rcParams['figure.figsize'] = (6,5)
    
    #os.chdir('/3DHST/Spectra/Work/HOPR/GYRO/')

    data = ascii.read('gyro_shifts.txt')
        
    #centers = 0.5*(bins[1:]+bins[:-1])
    #width = (bins[1:]-bins[:-1])
    
    plt.clf()
    fig = plt.figure(1)
    ax = fig.add_subplot(111)

    ax.fill_between([0.,0.5], 0, 30, color='0.85', alpha=0.9)
    nn, bins, patches = ax.hist(data['rate']*25., bins=60, range=[0.,6.], color='black', histtype='step')    

    #ax.bar(centers, nn, width, align='center', edgecolor='black', color='none')
    
    test_val = data['rate'][data['root'] == 'ib2u22prq']
        
    ax.arrow(test_val*25., 17.1, 0., -1., head_length=0.5, head_width=0.1, fc='red', ec='red')
    ax.text(0., 18.0, 'Fig. 1 Example', ha='left', va='center', fontsize=13, fontweight='heavy')
    
    ax.arrow(0., 21., 0.5, 0., head_length=0.05, head_width=0.30, fc='black', ec='black', length_includes_head=True)
    ax.arrow(0.5, 21., -0.5, 0., head_length=0.05, head_width=0.30, fc='black', ec='black', length_includes_head=True)
    ax.text(0.0, 22.5, 'Expected Drift', ha='left', va='center', fontsize=15, fontweight='heavy')
    #ax.text(1.1, 5, 'gyro\nproblem', fontsize=12.5, fontweight=0.1, multialignment='center')
    #ax.text(4.3, 5, 'failed\nguide star\nacquisition', fontsize=12.5, fontweight=0.1, multialignment='center')

    ax.set_xlim([-0.15, 5.5])
    ax.set_ylim([0., 25])
    ax.set_xlabel('Drift in Pixels per 25 Seconds', fontsize=12)
    ax.set_ylabel('N')
    ax.tick_params(axis='both', which='major',labelsize=12)
    minorLocator   = AutoMinorLocator()
    ax.xaxis.set_minor_locator(minorLocator)
    minorLocator   = AutoMinorLocator(5)
    ax.yaxis.set_minor_locator(minorLocator)
    ax.tick_params(which='minor', length=2)    
    
    new = np.array([0.1026, 0.158, 0.3768, 0.8161, 0.8854, 0.8106, 1.492, 0.4064, 0.5438, 0.6395, 0.7944, 0.5094, 0.6058, 1.356, 0.2935,0.5412,0.3815, 0.8022, 0.4679, 0.4038, 0.8373, 0.1917, 0.6154, 0.7328, 1.215, 1.004, 0.5233, 0.6105])
    new = new/4.
    mm, bins_mm, patches = ax.hist(new, bins=bins, range=[0.,6.], color='red', histtype='step')
    
    plt.show(block=True)
    
    plt.savefig('drift_hist_bw.png',dpi=100,transparent=False)
    plt.savefig('drift_hist_bw.pdf',dpi=100,transparent=False)
    
    
def pointing_strategy():
    
    """ 
    Make a figure showing the pointing strategy
    """
    
    import unicorn
    import unicorn.survey_paper as sup
    from my_python.phot_paper import determine_outline as determine_outline
    import matplotlib.pyplot as plt
    
    plt.rcParams['lines.linewidth'] = 0.3
    wfc3_color='#d7301f'
    acs_color='#2c7bb6'
    ultravista_color = '0.5'
    candels_color= '#fdae61'
    
    width=6.
    fontsize=8
    left=0.15
    right=0.05
    top=0.1
    bottom=0.1
    
    yticklab = None
    
    REGIONS_PATH = '/3DHST/Spectra/Work/HOPR/UVISTA/'
    CANDELS_PATH = '/Users/ivastar/3DHST/Regions/CANDELS/'
    
            
    ### Limits of  mosaic
    x0, x1 = 150.82, 149.28 
    y0, y1 = 1.51, 2.93
    #xticklab = [r'$45^\mathrm{m}00^\mathrm{s}$',r'$10^\mathrm{h}00^\mathrm{m}30^\mathrm{s}$', r'$15^\mathrm{m}00^\mathrm{s}$']
    #xtickv = [sup.degrees(10,00,45, hours=True),sup.degrees(10,00,30, hours=True), sup.degrees(10,00,15, hours=True)]
    #yticklab = [r'$10^\prime00^{\prime\prime}$', r'$15^\prime00^{\prime\prime}$', r'$+2^\circ20^\prime00^{\prime\prime}$', r'$25^\prime00^{\prime\prime}$', r'$30^\prime00^{\prime\prime}$', r'$35^\prime00^{\prime\prime}$']
    #ytickv = [sup.degrees(02, 10, 00, hours=False),sup.degrees(02, 15, 00, hours=False), sup.degrees(02, 20, 00, hours=False), sup.degrees(02, 25, 00, hours=False), sup.degrees(02, 30, 00, hours=False), sup.degrees(02, 35, 00, hours=False)]
    xticklab = [r'$10^\mathrm{h}02^\mathrm{m}$',r'$10^\mathrm{h}00^\mathrm{m}$',r'$9^\mathrm{h}58^\mathrm{m}$']
    xtickv = [sup.degrees(10,02,00, hours=True),sup.degrees(10,00,00, hours=True),sup.degrees(9,58,00, hours=True)]
    
    yticklab = [r'$45^\prime$', r'$+2^\circ00^\prime$', r'$15^\prime$', r'$30^\prime$',r'$45^\prime$']
    ytickv = [sup.degrees(01, 45, 00, hours=False),sup.degrees(02, 00, 00, hours=False),sup.degrees(2, 15, 00, hours=False),sup.degrees(02, 30, 00, hours=False),sup.degrees(02, 45, 00, hours=False)]
    
    #### Make square for given plot dimensions
    dx = np.abs(x1-x0)*np.cos(np.mean([y0,y1])/360*2*np.pi)
    dy = (y1-y0)
            
    fig = unicorn.catalogs.plot_init(square=True, xs=width, aspect=dy/dx, fontsize=fontsize, left=left, right=right, top=top, bottom=bottom)
    ax = fig.add_subplot(111)
        
    ### Outer outline of ULTRAVISTA
    file_uv_box = REGIONS_PATH+'uvista_box.reg'
    x,y = determine_outline([file_uv_box])
    #ax.plot(x,y, '-', alpha=0.9, color=ultravista_color, linewidth=3)
    ax.fill(x,y, alpha=0.15, color=ultravista_color)

    ### Deep stripes of ULTRAVISTA, filles
    polys = []
    fp = open(REGIONS_PATH+'uvista_deep.reg')
    lines = fp.readlines()
    fp.close()
    for line in lines[1:]:
        polys.append(sup.polysplit(line, get_shapely=True))
        sup.polys = polys
        un = polys[0]
        for pp in polys[1:]:
            un = un.union(pp)
        if un.geometryType() is 'MultiPolygon':
            for sub_poly in un.geoms:
                x,y = sub_poly.exterior.xy
        else:        
            x,y = un.exterior.xy
        ax.plot(x,y,'--', alpha=1.0, color=ultravista_color, linewidth=3.5)
        ax.fill(x,y, alpha=0.3, color=ultravista_color)
            
    ### CANDELS
    polys = []
    fp = open(REGIONS_PATH+'COSMOS-ALL-F160W.reg')
    lines = fp.readlines()
    fp.close()
    for line in lines[1:]:
        polys.append(sup.polysplit(line, get_shapely=True)) 
        wfcx, wfcy = sup.polysplit(line)
        fi = ax.plot(wfcx, wfcy, color=candels_color, linewidth=1)
        #fi = ax.plot(wfcx, wfcy,':', alpha=0.9, color=wfc3_color)
            
    ### ACS
    file_acs = REGIONS_PATH+'ACS_Cosmos_poly.reg'
    x,y = determine_outline([file_acs])
    ax.plot(x,y, alpha=0.9, color=acs_color, linewidth=1.5)
    
    ### Planned WFC3
    wfc3_polys = []
    fp = open(REGIONS_PATH+'wide_2.reg')
    lines = fp.readlines()
    fp.close()
    for line in lines[1:]:
        if not line.startswith('polygon'):
            continue
        else:
            wfc3_polys.append(sup.polysplit(line, get_shapely=True)) 
            wfcx, wfcy = sup.polysplit(line)
            fi = ax.plot(wfcx, wfcy, color=wfc3_color, alpha=0.8, linewidth=1.)    
    
    ax.set_xlabel(r'$\alpha$', fontsize=11, labelpad=0)
    ax.set_ylabel(r'$\delta$', fontsize=11, labelpad=0)


    labels = ['Proposed F160W', 'CANDELS F160W','ACS F814W', 'ULTRAVISTA Deep', 'ULTRAVISTA Wide']
    ypos = [0.28, 0.24, 0.20,0.16, 0.12]
    colors = [wfc3_color, candels_color, acs_color, '0.3', '0.8']
    for label, pos, color in zip(labels,ypos,colors):
        ax.text(0.07, pos, label,
            horizontalalignment='left',
            verticalalignment='top',
            transform = ax.transAxes,fontsize=9, color=color)

    
    fsi = '20'
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    ax.set_xticklabels(xticklab)
    xtick = ax.set_xticks(xtickv)
    ax.set_yticklabels(yticklab)
    ytick = ax.set_yticks(ytickv)
    
    #### Get field area from full WFC3 polygons
    #print wfc3_polys
    un_wfc3 = wfc3_polys[0]
    for pp in wfc3_polys[1:]:
        un_wfc3 = un_wfc3.union(pp)
    #
    wfc3_union= []
    
    if un_wfc3.geometryType() is 'MultiPolygon':
        total_area = 0
        xavg, yavg, wht = 0, 0, 0
        for sub_poly in un_wfc3.geoms:
            area_i = sub_poly.area*np.cos(y0/360.*2*np.pi)
            total_area += area_i
            x,y = sub_poly.exterior.xy
            wfc3_union.append(sub_poly)
            xavg += np.mean(x)*area_i**2
            yavg += np.mean(y)*area_i**2
            wht += area_i**2
        xavg, yavg = xavg/wht, yavg/wht
    else:        
        total_area = un_wfc3.area*np.cos(y0/360.*2*np.pi)
        x,y = un_wfc3.exterior.xy
        wfc3_union.append(un_wfc3)
        xavg, yavg = np.mean(x), np.mean(y)
    
    
    wfc3_area = 0.
    for wun in wfc3_union:
        wfc3_area += wun.area*np.cos(y0/360.*2*np.pi)*3600.

    print '== WFC3 ==\n: %.1f ' %(wfc3_area)
    
    
    plt.show(block=False)
    fig.savefig('pointing_strategy.png', dpi=100, transparent=False)
    fig.savefig('pointing_strategy.pdf', dpi=100, transparent=False)
    
    
def area_depth():
    
    import unicorn
    import unicorn.survey_paper as sup
    from my_python.phot_paper import determine_outline as determine_outline
    import matplotlib.pyplot as plt
    
    fig = unicorn.catalogs.plot_init(square=True, xs=5, aspect=1., fontsize=10, left=0.1, right=0.1, top=0.1, bottom=0.1)
    ax = fig.add_subplot(111)
       
    WFC3_IR_AREA = 4.6 # arcmin2
    imaging = {'UDF/XDF':(236.1e3/3000, 1*WFC3_IR_AREA, 0, 3.9), 'HUDF09-1': (13, 1*WFC3_IR_AREA,0,3.9), 'HUDF09-2':(19, 1*WFC3_IR_AREA, 0, 3.9), 'CANDELS-Deep': (12, 0.04*3600, 0,-3.9), 'CANDELS-Wide': (2, 0.2*3600, 0,3.9), 'CLASH': (2, 25*WFC3_IR_AREA, 0,-3.9), 'FF clusters': (24, 6*WFC3_IR_AREA, 0,3.9)} 
    
    
    for key in imaging.keys():
        ax.plot(imaging[key][1],imaging[key][0], 's', markersize=7, color='black')
        ax.text(imaging[key][1]+imaging[key][2], imaging[key][0]+imaging[key][3], key, ha='center', va='top', fontsize=8)
    
    ax.plot([2000.], [0.125], 's', markersize=7, color='red')
    ax.text(2000., 0.225, 'This\n proposal', multialignment='center', color='red',fontsize=8, ha='center', va='top')
    
    ax.set_xlim([-400, 2490.])
    ax.set_ylim([0.08, 500.])
    ax.set_yscale('log')
    ax.set_xlabel('Area [arcmin$^2$]', fontsize=12)
    ax.set_ylabel('N', fontsize=12)
    
    
    ### Log plot
    #ax.set_xlim([0.8, 4000.])
    #ax.set_ylim([0.08, 500.])
    #ax.set_yscale('log')
    #ax.set_xscale('log')

    plt.show(block=False)
            
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
        if key == 'HUDF09-1':
            ax.plot([13+50, 13+130], [mag, mag+0.17], color='black', lw=1.)
        if key == 'HUDF09-2':
            ax.plot([19+50, 19+130], [mag, mag+0.17], color='black', lw=1.)
        if key == 'FF clusters':
            ax.plot()

    
    ax.plot([2000.], [25.3], 's', markersize=12, color='red')
    ax.text(2000., 25.7, 'THIS PROPOSAL', multialignment='center', color='red',fontsize=10, ha='center', va='top')
    ax.text(1200, 28.5, 'HST/WFC3 F160W\n imaging', fontsize =12, multialignment='center')
    
    ax.set_xlim([-100, 2490.])
    ax.set_ylim([24.95, 29.5])
    #ax.set_yscale('log')
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
    
    
    ### Log plot
    #ax.set_xlim([0.8, 4000.])
    #ax.set_ylim([0.08, 500.])
    #ax.set_yscale('log')
    #ax.set_xscale('log')

    plt.show(block=False)
    
    fig.savefig('area_mag.png', dpi=200, transparent=False)
    fig.savefig('area_mag.pdf', dpi=200, transparent=False)
    