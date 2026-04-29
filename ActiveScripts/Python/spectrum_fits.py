"""
Models to fit the (x-ray) spectrum, finde the line position and
its associated resolution (fwhm).
"""
import numpy as np
import uproot as up
import os
import os.path
import time 

from copy import copy

import HErmes as he
import hepbasestack as hep
import dashi as d
import pylab as p

d.visual()

logger = hep.logger.get_logger(10)

from .utils import get_stripname

def model_with_shoulder_with_linear_cutoff(xs,\
                                           mu_p, fwhm_p, amp_p,\
                                           mu_s, fwhm_s, amp_s):
    """
    A full fledged model combining two gaussians - one for the Compton 
    shoulder and another one for the line peak. The shoulder has a linear
    cutoff in the reaion of the peak, ranging from mu_p - fwhm_p to mu_p. 
    for xs > mu_p the contribution of the shoulder is set to zero.

    Args:
        xs (array)       : input data
        mu_p (float)     : line peak
        fwhm_p (float)   : line resolution
        amp_p (float)    : number of events at the line peak
        mu_s (float)     : peak of the (Compton) shoulder
        fwhm_s (float)   : width of the shoulder
        amp_s (float)    : number of events at the peak of (Compton) 
                           shoulder

    Returns:
        ndarray          : model prediction
    """
    peak_ys = he.fitting.functions.fwhm_gauss(xs, mu_p, fwhm_p, amp_p) 
    shoulder_ys = he.fitting.functions.fwhm_gauss(xs, mu_s, fwhm_s, amp_s)
   
    linear_ys = np.zeros(len(peak_ys)) 
    # construct the cutoff-part
    if (np.isfinite(mu_p) and np.isfinite(fwhm_p)):
        start_x = xs[xs <= (mu_p - (abs(fwhm_p)/2))][-1]
        start_y = shoulder_ys[xs <= (mu_p - (abs(fwhm_p)/2))][-1]
        end_x   = mu_p
        assert start_x < end_x
        #end_y   = shoulder_ys[xs < mu_p][-1]
        end_y   = 0
        #print (f'linear interpolation from {start_x} to {end_x}. Yrange {start_y} to {end_y}')
        slope = (abs(start_y))/(abs(end_x) - abs(start_x))
        constant = slope*abs(end_x)

        #print (constant)
        for i, x in enumerate(xs):
            if (x >= start_x and x <= end_x):
                linear_ys[i] = -1*slope*x + constant
        #linear_interp = lambda xs : -1*slope*xs + np.ones(len(xs))*abs(start_y)
        #linear_ys[np.logical_and(xs >= start_x, xs<= end_x)] = linear_interp(xs[np.logical_and(xs >= start_x, xs<=end_x)])
        #print (linear_ys)
        #print (start_x, end_x, start_y)
        shoulder_ys[xs >= mu_p] = 0
    return peak_ys + shoulder_ys - linear_ys

########################################################################

def peak_only_model_with_linear_subtraction(xs, mu_p, fwhm_p, amp_p):
    """
    A model with a single peak, however, we subtract a linear
    contribution coming from a possible shoulder
    
    Args:
        xs (array)       : input data
        mu_p (float)     : line peak
        fwhm_p (float)   : line resolution
        amp_p (float)    : number of events at the line peak

    """
    peak_ys = he.fitting.functions.fwhm_gauss(xs, mu_p, fwhm_p, amp_p) 
    
    # construct the cutoff-part
    if (np.isfinite(mu_p) and np.isfinite(fwhm_p)):
        start_x = xs[xs < (mu_p - abs(fwhm_p)/2)][-1]
        start_y = peak_ys[xs < (mu_p - abs(fwhm_p)/2)][-1]
        end_x   = xs[xs < mu_p][-1]
        end_y   = peak_ys[xs < mu_p][-1]
        slope = (end_y - start_y)/(end_x - start_x)
        linear_interp = lambda xs : slope*xs + np.ones(len(xs))*start_y
        peak_ys[np.logical_and(xs >= start_x, xs<= end_x)] = linear_interp(xs[np.logical_and(xs >= start_x, xs<=end_x)])
        peak_ys[xs > end_x] = 0
    return peak_ys

##########################################################################

def create_model(energy, bins, startparams,\
            use_cutoff_for_shoulder=False,\
            limits = ((None, None),\
                      (None, None),\
                      (None, None),\
                      (None, None),\
                      (None, None),\
                      (None, None)),
            constrain_parameters=False,\
            mode='distribution'):
    """
    The fitting routine. Gaussian for peak + gaussian for Compton shoulder. As a default, energy 
    is all the energy values, and this function will create a distribution with binning as specified 
    in bins. This behaviour can be changed by modifying the mode switch.

    Args:
        energy (ndarray)               : energy as obtained by the digitizer
        bins (ndarray)                 : bins to create the distribution
        startparams (ndarray)          : a set of 6 startparams [mu_x, fwhm_x, amp_x, mu_x, fwhm_x, amp_x]

    Keyword Args:
        use_cutoff_for_shoulder (bool) : use a gaussian model with a cutoff
                                         for the fitting of the shoulder
        limits (tuple)                 : min/max limits for the fitparameters. Are
                                         passed through to 
                                         HErmes.fitting.Model.fit_to_data
        mode (str)                     : default is *distribution* which means that a histogram
                                         will be created out of *energy* with the binning *bins*.
                                         Alternatively, this can be set to *direct*, then *energy*
                                         needs to be bin content and bins need to be the corresponding
                                         bincenters. In this case there will no histogramming by 
                                         this function.
    """
    
    if use_cutoff_for_shoulder:
        # use a model with 2 gaussians and
        # a linear cutoff part at the 'end'
        # of the shoulder
        peakmod = he.fitting.Model(model_with_shoulder_with_linear_cutoff)
    else:
        # the x-ray peak
        peakmod = he.fitting.Model(he.fitting.functions.fwhm_gauss)
        # fit the compton shoulder with another gaussian to increase 
        # fit quality
        shouldermod = he.fitting.Model(he.fitting.functions.fwhm_gauss)
        peakmod += shouldermod
    
    if mode == 'distribution':
        peakmod.add_data(energy, bins=bins, create_distribution=True)
    elif mode == 'direct':
        assert len(energy) == len(bins), "In this mode, 'energy' and 'bins' must be the same length"
        peakmod.add_data(energy, xs=bins)
    else:
        raise ValueError(f'Mode {mode} not understood. Has to be one of "distribution" or "direct"')
    peakmod.startparams = startparams
    # we still constrain the parameters, 
    # condition is that xray peak 
    # and shoulder peak should not interfere
    peakmod.fit_to_data(silent=True,\
                        errors=(1,1,1,1,1,1),\
                        use_minuit=True,\
                        limits = limits)
    
    return peakmod

########################################################################

def first_guess(bincenters, bincontent, roi_right_edge=700):
    """
    Get a reasonable set of startparameters for a certain
    energy. This will look for the peak of the line and 
    then try to set the shoulde peak as far away as possible.

    Args:
        bincenters (ndarray)   : xs
        bincontent (ndarray)   : ys

    Keyword Args:
        roi_right_edge (int)   : the right edge of the roi (region-of-interest)

    """
    assert len(bincontent) == len(bincenters), f"Bincenters and bincontent must have the same length! {len(bincenters)} {len(bincontent)}"

    # we know that the spectrum is the full range of the digitizer, 
    # so first look at a region of interest.
    # typically, the peak will be around 500, 
    # so lets take 0 - 700 (default), then the peak 
    # should be savely in the right half of that part of the spectrum
    bincenters = bincenters[:roi_right_edge]
    bincontent = bincontent[:roi_right_edge]

    # assume we have a fair number of bins, so let's try 
    # to get the peak value for the line
    right_half = int(len(bincenters)/2)
    left_10_percent = int(right_half/5)
    logger.info(f'Found spectrum of len {len(bincenters)}. Dividing in half left and right of {right_half}')
    right_half_spectrum = bincontent[right_half:]
    peakval = np.where(right_half_spectrum == max(right_half_spectrum))[0][0]
    logger.info(f'Found an estimate for the line peak {peakval}')
    peakval += right_half 
    # estimate the fwhm
    peak_amp = max(right_half_spectrum)
    logger.info(f'Found an estimate for the peak amplitude {peak_amp}')
    mod = he.fitting.Model(he.fitting.fwhm_gauss)
    mod.startparams = (peakval, 5, peak_amp)
    mod.add_data(bincontent, xs=bincenters, create_distribution=False)
    mod.fit_to_data(silent=True,\
                    errors=(1,1,1),\
                    limits=((peakval - 0.01*peakval, peakval + 0.01*peakval),
                            (1, 100),
                            (peak_amp - 0.01*peak_amp, peak_amp + 0.01*peak_amp)))
    peak_fwhm = mod.best_fit_params[1]

    # then just ensure that the shoulder peak is far enough from the
    # line peak, maybe just start from the half?
    shoulder_peak = right_half
    # and finally assume the peak hight is 0.1 of the xray peak
    shoulder_amp = peak_amp/10
    
    # try to prefit it
    mod = he.fitting.Model(he.fitting.fwhm_gauss)
    mod.startparams = (shoulder_peak, 10, shoulder_amp) # 10 is arbitrary
    mod.add_data(bincontent, xs=bincenters, create_distribution=False)
    if peak_amp <=1:
        logger.warning("Estimate for the  peak amplitutde is doubtful...")
        peak_amp = 100 
    mod.fit_to_data(silent=True,\
                    errors=(1,1,1,1,1,1),\
                    limits=((left_10_percent, peakval - peak_fwhm),\
                            (1, 100),\
                            (1, peak_amp/2)))
    shoulder_fwhm = mod.best_fit_params[1]

    #first_guess_params = np.array([500,20,500, 300, 40, 50])
    first_guess_params = np.array([peakval, peak_fwhm, peak_amp,\
                                   shoulder_peak, shoulder_fwhm, shoulder_amp])
    limits = ((peakval-peak_fwhm, peakval+peak_fwhm),\
              (0.5*peak_fwhm, 2*peak_fwhm),\
              (1, peak_amp*1.1),\
              (left_10_percent, peakval-2*peak_fwhm),\
              (1, shoulder_fwhm*2),
              (1, peak_amp/2))
    for i,k in enumerate(limits):
        logger.debug(f"Setting limit on parameter {i} of {k}")
    return first_guess_params, limits    


########################################################################

def fit_file(infilename = '143_4000.root',\
             file_type = '.root',\
             channel = 0,\
             ptime = -1,\
             savename_prefix = "",\
             detid = -1,\
             plot_dir = '.',\
             energy = None,\
             #bins = np.linspace(150, 550, 200),\
             bins = np.linspace(1,2**14, 2**14) - 0.5,\
             peakposition = 88.03,\
             debug=False):
    """
    One shot function to fit a datafile with 
    energy values obtained by the Caen N6725 digitizer
    with two gaussians. One for the signal peak,
    another one for the Compton shoulder.

    Keyword Args:
        infilename (str)     : the file to fit. Must have energy  
        file_type (str)      : A specialized string to identify the file type
                               Can be either '.root' (default),'.txt' or None
                               In case of '.txt' one line per digitizer 
                               bin is expected (2**14 lines) with one
                               value each. If None is given, then the parameter energy
                               has to be different than None. In that case,
                               there is no need to read out a file 
        channel (int)        : digitizer channel (0-7)
        energy  (ndarray)    : in case there is no file, but we 
                               have the data already, it can be passed via
                               the energies parameter. If there is a file, 
                               this has to be 'None'
        plot_dir (int)       : directory to save the plots in 
        savename_prefix (str): a prefix to be added to the png file the plot is saved as
        bins (ndarray)       : bins for the histogram
        ptime (float)        : peaking time (just used for the plot title)
        detid (int)          : detector id (just uded for the plot title)
        peakposition (float) : the true value of the x-ray peak in keV
                               this is used to recalibrate the x-axis
        debug        (bool)  : Use to reduce the loglevel, otherwise this 
                               function is as silent as possible.
                               Internally sets the loglevel of HErmes to debug
    """
    if not debug:
        he.set_loglevel(100) 
        logger.setLevel(100)

    # set startparameters
    startparams = np.array([500,20,500, 300, 40, 50])
    metainfo = "" # this currently holds the type of measurement
                  # if it is from a file -> CAEN shaping
                  # mode direct -> mc2 
                  # else -> waveform
    if file_type == '.root':
        mode = 'distribution'
        f = up.open(infilename)
        trigger = f.get('ch' + str(channel)).get('trigger').array()
        energy  = f.get('ch' + str(channel)).get('energy').array()    
        # prebin the histogram to estimate start params
        h = d.factory.hist1d(energy, bins)
        startparams, limits = first_guess(h.bincenters, h.bincontent)
        peakmod = create_model(energy, bins, startparams, limits=limits)
        metainfo = 'TrpzShap'
    elif file_type == '.txt':
        mode='direct'
        energy = []
        with open(infilename) as f:
            for line in f.readlines():
                binenergy = 0
                if len(line.split()) == 2:
                    # we have 2 column data. second line is the bincontent
                    # FIXM: for some reason, in the mc2 data, the first
                    # two bins are a float. This can be bascially ignored
                    # since it is out of our region of interest
                    binenergy = int(float(line.split()[1]))
                else:
                    # one column data, just energy
                    binenergy = int(line.split()[0])
                energy.append(binenergy) 
        metainfo = 'TrpzShapMC2'

        energy = np.array(energy)
        if len(energy) != 2**14:
            raise ValueError(f"Bins missing! Only found {len(energy)} bins in the file!")
        bins = np.linspace(1,2**14, 2**14)
        # make them the bincenter instead of the edge
        bins = bins - 0.50
        startparams, limits = first_guess(bins, energy)
        peakmod = create_model(energy, bins, startparams, mode=mode, limits=limits)

    elif file_type == None:
        if energy is None:
            raise ValueError("Need to give energy or filename and filetype!")
        # direct mode, no file. in that energy and bins must be given directly
        mode = 'distribution'
        # prebin the histogram to estimate start params
        h = d.factory.hist1d(energy, bins)
        logger.info(f'Caluclating start parameters')
        startparams, limits = first_guess(h.bincenters, h.bincontent)
        peakmod = create_model(energy, bins, startparams, mode=mode, limits=limits)
        metainfo = 'GaussShap' 
    else:
        raise ValueError(f'Can not process file type {file_type}')

    logger.info(f'Got chi2/ndf of {peakmod.chi2_ndf}')
    logger.info(f'Found peak at {peakmod.best_fit_params[0]}')
    logger.info(f'Best fit parameters {peakmod.best_fit_params}')

    peak = peakmod.best_fit_params[0]

    # recalibrate the x-axis
    logger.info(f'Calibrating x-axis...{peak} -> {peakposition}')
    conversion = peakposition/peak

    # only convert the energies if they 
    # are truly energies.
    # in the 'direct' approach for mc2 data
    # we 'hijack' these to be the bincontents
    if mode =='direct':
        energyxs = energy
    else:
        energyxs = conversion*energy
    

    # use the previously obtained fit params as start params
    # remember the amplitdue is not converted.
    startparams = [conversion*startparams[0], conversion*startparams[1],\
                   startparams[2], conversion*startparams[3],\
                   conversion*startparams[4], startparams[5]]

    limits = [conversion*np.array(limits[0]),\
              conversion*np.array(limits[1]),\
              np.array(limits[2]),\
              conversion*np.array(limits[3]),\
              conversion*np.array(limits[4]),\
              np.array(limits[5])]
    limits = tuple(limits)
    bins = conversion*bins

    startparams = np.array([abs(s) for s in startparams])
    # get a better fit with a more refined model
    logger.info(f'Refitting...with startparams {startparams}')
    peakmod = create_model(energyxs, bins, startparams,\
                      limits = limits,\
                      use_cutoff_for_shoulder= False,\
                      mode=mode)
    logger.info(f'Got chi2/ndf of {peakmod.chi2_ndf}')
    logger.info(f'Found peak at {peakmod.best_fit_params[0]} with a FWHM of {peakmod.best_fit_params[1]}')

    # create a nice figure
    #fig = p.figure(figsize=SINGLE_FIG, dpi=DPI)
    fig = p.figure()
    #fig = matplotlib.figure.Figure()
    ax = fig.gca()
    if mode == 'direct':
        ax.scatter(peakmod.xs, peakmod.data, c="k",\
                   edgecolors="k", marker="+",\
                   s=20, linewidth=1  )
    else:
        peakmod.distribution.scatter()
    ax.plot(peakmod.xs, peakmod.prediction(peakmod.xs), color="r")
    ax = hep.visual.adjust_minor_ticks(ax)
    ax.set_xlabel('energy [keV]')
    ax.set_ylabel('events')
    
    # error handling
    errdict = peakmod.errors

    # compose a title
    stripname = get_stripname(channel)
    title = f'det. {detid} {stripname} pt {ptime/1000}$\mu$s'
    ax.set_title(title, loc='right')

    # create a textbox with some output
    infotext  = r"\begin{tabular}{ll}"
    infotext += r"spectra fit:&\\ "
    if mode == 'direct':
        infotext += r"entries: & {:4.2f} \\".format(sum(peakmod.data))
    else:
        infotext += r"entries: & {:4.2f} \\".format(peakmod.distribution.stats.nentries)
    infotext += r"peak (line)  & {:4.2f} $\pm$ {:4.2f} \\ ".format(peakmod.best_fit_params[0], errdict['mu0'])
    infotext += r"FWHM (line)  & {:4.2f} $\pm$ {:4.2f} \\ ".format(peakmod.best_fit_params[1], errdict['fwhm0'])
    if len(peakmod.best_fit_params) > 3:
        infotext += r"peak (shoulder) & {:4.2f} $\pm$ {:4.2f} \\ ".format(peakmod.best_fit_params[3], errdict['mu1'])
        infotext += r"FWHM (shoulder) & {:4.2f} $\pm$ {:4.2f} \\ ".format(peakmod.best_fit_params[4], errdict['fwhm1'])
        infotext += r"$\chi^2/ndf$ & {:4.2f} \\ ".format(peakmod.chi2_ndf)
    infotext += r"\end{tabular}"
    ax.text(0.3, 0.5, infotext,\
            horizontalalignment='center',\
            verticalalignment='center',\
            transform=ax.transAxes,\
            bbox=dict(facecolor='white', alpha=0.7, edgecolor=None),\
           )

    # in case we have one bin per digitizer 
    # channel, cut down on the xrange
    if len(bins) > 1000:
        # ignore the first 10 mV
        # we will not see a peak at zero,
        # but actually we don't care
        # and this helps with the scaling
        ax.set_xlim(left=0, right=110)
        if (ax.get_ylim()[1] > 2*max(energyxs[20:])):
            ax.set_ylim(top= 1.1*max(energyxs[20:]),
                        bottom = 0)

    # return the resolution
    pngfilename = f'det{detid}-stime{ptime}-{stripname}-{metainfo}.png' 
    if savename_prefix:
        pngfilen = savename_prefix + pngfilename
    logger.info(f'Saving {pngfilename} in {plot_dir}... ')
    fig.savefig(os.path.join(plot_dir,pngfilename))
    time.sleep(1) # give the filesystem time to save the figure
                  # not sure why it seems that there are some hickups
    p.close(fig)
    logger.info('..done')
    # copy them to decouple them from the model
    # since some minuit entities are non-picklable
    # which prevents this from running multithreaded
    return_pars = copy(peakmod.best_fit_params)
    return_errs = dict()
    for k in errdict.keys():
        return_errs[k] = errdict[k]
    #@return detid, channel, peakmod.best_fit_params[1], errdict['fwhm0'], pngfilename
    return detid, channel, return_pars, return_errs, pngfilename

