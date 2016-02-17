# -*- coding: utf-8 -*-
"""
Created on Mon Jan 04 13:39:27 2016

@author: Jordan
"""

# imports
import matplotlib.pyplot as plt
import numpy as np    
import scipy as sp


#==========================================================================        
def multiplot(data,Fs=1,col="k",linesep=0):
    
    """
    Plots multiple lines on same graph, spaced apart vertically.
    
    Inputs:
        data = data matrix in column format...numpy matrix/ndarray works best
        Fs = sampling rate...if not provided, default = 1
        col = line color...default = black
        linesep = spacing amount...if not provided, default = max(data)
     """
     
    # create time vector using sampling rate
    timevec = np.array(range(len(data)))/float(Fs)
    
    # get spacing
    if linesep == 0:
        linesep = np.max(data)
    
    # loop and plot
    for i in range(np.size(data)/len(data)):
        plt.plot(timevec,data[:,i] - (linesep * i),color=col)

    ax = plt.gca()
    ax.xaxis.set_tick_params(direction='out')
    ax.yaxis.set_tick_params(direction='out')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
#==========================================================================        


#==========================================================================    
def fillplot(data,Fs=1,errtype='sem',linecol='k',fillcol=[.1,.4,.9],edgecol=0,edgesize=0,opacity=1):
        
    """    
    Plots mean of a data matrix (across columns) and error as shaded bar.
    
    Inputs:
        data = nxm matrix with n = observations, m = variables
        Fs = sampling rate...defualt = 1
        errtype = 'sd','sem', or 'ci'...defualt = 'sem'
        linecol = color of main line...default = 'k'
        fillCol = fill color of error...default = [.1 .4 .9]
        edgecol = color of edge of fill...default = same color as fillCol
        edgesize = width of line of fill...default = 0
        opacity = transparency of fill...default = 1 (opaque)
    Outputs:
        hf = figure handle
    """
    
    # calculate error
    sd = np.std(data,1) # gets std for each data point across columns
    sem = sd/data.shape[1]
    ci = sem*2
    
    # set error to errtype
    if errtype == "sem":
        err = sem
    elif errtype == "sd":
        err = sd
    elif errtype == "ci":
        err = ci
    
    # get mean across columns (trials) and timevec
    meanval = np.mean(data,1)        
    timevec = np.array(range(len(meanval)))/float(Fs)
    
    # set edgecol to fillcol if edgecol == 0
    if edgecol == 0:
        edgecol = fillcol
        
    # plot the data and fill
    hf = plt.plot(timevec,meanval,color=linecol,linewidth=2) 
    plt.fill_between(timevec,meanval+err,meanval-err,
                     facecolor=fillcol,edgecolor=edgecol,
                     linewidth=edgesize,alpha=opacity)
    
    return hf
#==========================================================================        


#==========================================================================    
def raster(data,start,Fs=0,pretime=0.5,posttime=1.0,tickheight=2,tickwidth=1,matshape='c'):
    '''
    Creates a raster plot of an nxm data mtarix.
    
    Inputs:
        data = data matrix
        start = scalar or array of event starting times in seconds
        Fs = sampling rate...default = 0.
            - If specified, assumes spiketimes are in samples and will convert to seconds
        pretime = pretime relative to start for plot...default = 0.5 second
        posttime = posttime relative to start for plot...default = 1 second
        tickheight = height of tickmark...defualt = 2
        tickwidth = width of each tickmark...default = 1
        matshape = 'r' or 'c' (row or colum)...default = 'c'
    '''
    
    # check matrix format
    if matshape == 'r':
        data = data.T # transpose    
    matsize = data.shape
    
    # convert spiketimes to seconds if Fs provided
    if Fs > 0:
        data = data/float(Fs)
    
    for trial in range(matsize[1]):
        spikes = data[:,trial]
        valid = sp.logical_and(spikes >= start-pretime, spikes <= start+posttime) # boolean
        spikes = spikes[valid] - start # extract spikes that are within the pre/post times relative to start tim
        del valid 
        
        if spikes.any():
            plt.plot(spikes,np.ones(spikes.shape)*trial,'k|',markersize=tickheight,mew=tickwidth)
        del spikes
    
    # tweak plot
    plt.ylabel('trials')
    plt.xlabel('time (s)')
    plt.axis([-pretime,posttime,-1,matsize[1]])
    ax = plt.gca()
    
    # remove top/right axes
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.yaxis.set_ticks_position('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.xaxis.set_tick_params(direction='out') 
#==========================================================================        
      

#==========================================================================          
def psth(data,start,Fs=0,bw=0.01,pretime=0.5,posttime=1.0,col='k',matshape='c',kernal='none',plotting=0):   
    '''
    Calculates a peri-stimulus time histogram (psth) of neural activity.
    Use can additionaly specify to use smoothing kernals of the psth histogram.
    
    Inputs:
        data = nxm data matrix of spikes, where n = observations, m = variables
        start = scalar or mx1 array of event start times
        Fs = sampling rate...default = 0.
            -If provided, will convert spiketimes (assumed as samples) to seconds using Fs
        bw = bin width for extracting spikes...default = 0.01 seconds
        pretime = pretime to plot psth relative to start...default = 0.5 seconds
        posttime = post-time to plot psth relative to start...default = 1 second
        col = color of histogram...default = 'k'
        matshape = 'r' or 'c' (row or column)...default = 'c'
        kernal = 'none','gaussian','causal'...default = 'none'
        plotting = 1 or 0...default = 0. 
            -If 1, plots mean psth, variance, and fanofactor in separate figure
    Outputs:
        spikerate = full psth of the data matrix
        meanrate = mean psth across variables/trials
        meanvar = variance of fullrate
        fanofactor = meanvar/meanrate
        smoothrate = smoothed spikerate via optional smoothing filter
    '''  
    
    # check matrix format
    if matshape == 'r':
        data = data.T # transpose    
    
    # convert spiketimes to seconds if Fs provided
    if Fs > 0:
        data = data/float(Fs)    
    
    ntrials = data.shape[1]
    nbins = (posttime+pretime)/bw # number of bins for histogram
    spikerate = np.zeros((nbins,ntrials)) # initiate with zeros
    
    for trial in range(ntrials):
        spikerate[:,trial],edges = np.histogram(data[:,trial],int(nbins),(start-pretime,start+posttime))
        
    # calculate the histogram with the given bin edges
    #count,edges = np.histogram(data,nbins,(start-pretime,start+posttime)) 
    
    
    meanrate = np.mean(spikerate,1) * (1/bw)
    meanvar = np.var(spikerate,1) * (1/bw)
    fanofactor = meanvar/meanrate
    spikerate = spikerate * (1/bw)
    if 
    
    if plotting == 1:
        f, (ax1,ax2,ax3) = plt.subplots(3,1,sharex=True)
        
        ax1.bar(edges[0:nbins]-start,meanrate,bw,color=col,linewidth=0)
        ax1.set_title('Mean spike rate')
        ax1.tick_params(direction='out')
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax1.xaxis.set_ticks_position('bottom')
        ax1.yaxis.set_ticks_position('left')
        ax1.axis([-pretime,posttime,0,max(meanrate)])

        ax2.plot(edges[0:nbins]-start,meanvar,'b')
        ax2.set_title('Variance')
        ax2.tick_params(direction='out')
        ax2.spines['top'].set_visible(False)
        ax2.spines['right'].set_visible(False)
        ax2.xaxis.set_ticks_position('bottom')
        ax2.yaxis.set_ticks_position('left')
        
        ax3.plot(edges[0:nbins]-start,fanofactor,'r')
        ax3.set_title('Fanofactor')
        ax3.tick_params(direction='out')
        ax3.spines['top'].set_visible(False)
        ax3.spines['right'].set_visible(False)
        ax3.xaxis.set_ticks_position('bottom')
        ax3.yaxis.set_ticks_position('left')
        
    return spikerate,meanrate,meanvar,fanofactor
#========================================================================== 


#========================================================================== 
def isi(data,Fs=0,matshape='c',plotting=0):
    '''
    Calculate the inter-spike interval (ISI) of a spike-train array or matrix.
    
    Inputs:
        data = nxm matrix where n = observations (spikes), m = trials/variables
        Fs = sampling rate...default = 0. 
            -If provided, will convert spiketimes to seconds using Fs
        matshape = 'r' or 'c'...default = 'c'
        plotting = 1 or 0...default = 0. 
            -If 1, plots mean isi across trials, variance, and fanofactor
    Outputs:
        fullisi = inter-spike intervals for the entire dataset
        meanisi = mean isi across trials for each time point
        varisi = variance across trials 
        fanofactor = varisi/meanisi
        * note isi is length n-1
    '''
    
    # check matrix format
    if matshape == 'r':
        data = data.T
    
    # check Fs
    if Fs > 0:
        data = data/Fs
    
    # get ISI for entire dataset 
    fullisi = sp.diff(data,axis=0)
    meanisi = fullisi.mean(1)
    varisi = sp.var(fullisi,1)
    fanofactor = varisi/meanisi
    
    if plotting == 1:
        f, (ax1,ax2,ax3) = plt.subplots(3,1,sharex=True)
        
        ax1.plot(meanisi,color='k')
        ax1.set_title('Mean ISI')
        ax1.tick_params(direction='out')
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax1.xaxis.set_ticks_position('bottom')
        ax1.yaxis.set_ticks_position('left')

        ax2.plot(varisi,'b')
        ax2.set_title('Variance')
        ax2.tick_params(direction='out')
        ax2.spines['top'].set_visible(False)
        ax2.spines['right'].set_visible(False)
        ax2.xaxis.set_ticks_position('bottom')
        ax2.yaxis.set_ticks_position('left')
        
        ax3.plot(fanofactor,'r')
        ax3.set_title('Fanofactor')
        ax3.tick_params(direction='out')
        ax3.spines['top'].set_visible(False)
        ax3.spines['right'].set_visible(False)
        ax3.xaxis.set_ticks_position('bottom')
        ax3.yaxis.set_ticks_position('left')
        
    return fullisi,meanisi,varisi,fanofactor
#==========================================================================    


#==========================================================================  
def harmonicfilter(data,Fs,fund=60,amp=0,offset=0,harmonics=0,maxtol=1e-10,maxiter=10,plotting=0):
    '''
    Fits a sinusoid to a 1D signal using "fund" and any amount of harmonics using least squares.
    
    Notes:
        Fits best when the frequency and amplitude of the inline noise is approximately known. 
        Recursively calls itself for the # of harmonics to fit

    Inputs:
        data = 1D signal of noisy data to be filtered
        Fs = sampling rate
        fund = fundamental frequency to fit (default = 60 Hz)
        amp = amplitude of noise (default = std(data))
        offset = offset of the noise (default = mean(data))
        harmoincs = # of harmonics to fit and remove (default = 0)
        maxtol = max error tolerance (default = 1e-10)
        maxiter = maximum # of times to increase maxtol if not converging (default = 10)
        plotting = 1 or 0...if 1, plots the raw and filtered data (default = 0)
    Outputs:
        filtdata = filtered data by removing fundamental and harmonic noises
        fitparams = amplitude, frequency, and offset parameters of the fitted sinusoid (for fundamental only)
    '''

    # import optimize.leastsq from scipy
    from scipy.optimize import leastsq

    # create time vector using Fs and len(data)
    tm = sp.linspace(0,len(data)/Fs,len(data))

    # get default parameters for the fit if not supplied by the user
    if amp == 0:
        amp = sp.std(data)
    if offset == 0:
        offset = sp.mean(data)

    # create sine function and cost function for the fit, and initial parameters 
    sinefunc = lambda p: p[0] * sp.sin(2*sp.pi*tm*p[1]) + p[2]
    costfunc = lambda p: sinefunc(p) - data
    guess = [amp, fund, offset]

    # run least squares and subtract the fit...loop until convergence, mult. error by 10 each time
    success = 0; iterval = 0;
    while success != 1 and iterval <= maxiter:
        results = leastsq(costfunc,guess,ftol=maxtol)
        fitparams = results[0]
        success = results[1]

        print('Not converged...increasing max tolerated error')
        iterval += 1 # add 1 to the iterval to keep track of loop
        maxtol *= 10 # mult by 10 for next iteration if not converged
    
    
    # create a fit from the fit params
    if success == 1:    
        print('Converged!')
        fit = sinefunc(fitparams) # make the fit
        filtdata = data - fit # subtract the fit from the data

        # plot if plotting == 1
        if plotting == 1:
            plt.figure() # raw data + fit
            plt.subplot(2,1,1)
            plt.plot(tm,data,'k')
            plt.plot(tm,fit,'r')
            plt.legend(['raw','fit'])

            plt.subplot(2,2,3) # raw data + filtered data
            plt.plot(tm,data,'k')
            plt.plot(tm,filtdata,'r')
            plt.legend(['raw','filtered'])

            plt.subplot(2,2,4) # difference (filtered-raw)
            plt.plot(tm,filtdata-data)
            plt.legend(['difference (filtered - raw)'])

            plt.suptitle('Freq: ' + str(round(fitparams[1],4)) + ' | Amp: ' + str(round(fitparams[0],4)) + ' | Offset: ' + str(round(fitparams[2],4)))


        # recursively call the function for the number of harmonics supplied
        if harmonics > 0:
            print
            print('Done with fundamental...fitting harmonics')
            for i in range(harmonics):
                harmonics -= 1 # subtract one from the harmonics
                fund *= 2 # multiply fundamental by two
                filtdata,_ = harmonicfilter(filtdata,Fs,fund=fund,plotting=1)
    else:
        print
        print('Could not converge...ending function.')

    # return the filtered data and fit params for the fundamental frequency
    return filtdata, fitparams, fit
#==========================================================================  


#========================================================================== 
def spectrum(data,Fs,frequencies=0,type='pgram',scal='spectrum',plottype='log',col='k',plotting=1,fitting=0):
    '''
    Caclulates the power spectrum of a signal and optionally fits 1/f noise and a best-fit line.

    Inputs:
        data = data array
        Fs = sampling rate
        frequencies = min and max of frequencies of interest (default = 0, or all frequenceis)
        type = 'pgram' or 'welch'...defines the method for spectral decomposition (default = 'pgram')
        scal = 'spectrum' or 'density'...(V**2 vs. V**2/Hz...default = 'spectrum')
        plottype = 'log' or 'linear' (default = log)
        col = line color (defualt = 'k')
        plotting = 0 or 1...defualt = 1
        fitting = 0 or 1...if 1, fits pink noise to the data (default = 0)
    Outputs:
        freqs = frequency array of the sepctral power
        pxx = power spectrum (or density) of the data
        fit = vector of fitted pink noise if fitting = 1
        fitparams = slope + y-intercept of fit line on a log-log plot (1/f^a + e^b on linear plot)
    '''

    # imports
    from scipy import signal
    from scipy.optimize import leastsq
    from statsmodels.nonparametric.smoothers_lowess import lowess

    # compute the spectral power 
    if type == 'pgram':
        freqs,pxx = signal.periodogram(data,Fs,scaling=scal)
    elif type == 'welch':
        freqs,pxx = signal.periodogram(data,Fs,scaling=scal)

    # get only the frequencies of interest
    if frequencies != 0:
        idx = (freqs >= frequencies[0]) * (freqs <= frequencies[1])
        freqs = freqs[idx]
        pxx = pxx[idx]

    # fit pink noise if fitting = 1
    if fitting == 1:

        # perform fit
        linefunc = lambda p: p[0] * sp.log(freqs) + p[1] # mx + b in log space, 1/f^m + e^b in linear space
        costfunc = lambda p: linefunc(p) - sp.log(pxx) # cost function between line and data
        pinkfunc = lambda p: -1*sp.log(freqs) + p # purely pink noise + offset...no change in slope
        pinkcost = lambda p: pinkfunc(p) - sp.log(pxx)

        # make a guess of parameters
        guess = [-1,0]

        # perform the fit using leastsq
        fitparams = leastsq(costfunc,guess)[0]
        fit = linefunc(fitparams) # make 1/f^m fit
        pinkparams = leastsq(pinkcost,0)[0] # guess zero offset
        pinkfit = pinkfunc(pinkparams) # make 1/f fit

        # changing plotting type to 'log' and plotting = 1
        plottype = 'log'
        plotting = 1
    
    else:
        fit = 0
        fitparams = 0


    # plot the result and fit (if fitting = 1)
    if plotting == 1:
        plt.figure()        
        plt.grid(b=True,which='minor',color=[.7,.7,.7],linestyle='--')
        plt.grid(b=True,which='major',color=[.7,.7,.7],linestyle='-')
        
        plt.plot(freqs,pxx,'k')
        if fitting == 1:
            plt.plot(freqs,sp.exp(fit),'b')
            plt.plot(freqs,sp.exp(pinkfit),'m')
            plt.legend(['signal','best fit: ' + str(round(fitparams[0],3)),'1/f'])
        
        ax = plt.gca()
        ax.set_yscale(plottype)
        ax.set_xscale(plottype)


        ax.set_xlabel('Frequency (Hz)')
        if scal == 'spectrum':
            ax.set_ylabel('Power (V^2)')
        elif scal == 'density':
            ax.set_ylabel('Power (V^2 / Hz)')

    return freqs,pxx,sp.exp(fit),fitparams # return fit in linear space 
#========================================================================== 


