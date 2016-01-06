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
        data = data matrix in column format...numpy matrix works best
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
        valid = sp.logical_and(spikes>=start-pretime,spikes<=start+posttime) # boolean
        spikes = spikes[valid] - start; # extract spikes that are within the pre/post times relative to start tim
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
def psth(data,start,Fs=0,bw=0.01,pretime=0.5,posttime=1.0,col='k',matshape='c',plotting=0):   
    '''
    Calculates a peri-stimulus time histogram (psth) of neural activity.
    
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
        plotting = 1 or 0...default = 0. 
            -If 1, plots mean psth, variance, and fanofactor in separate figure
    Outputs:
        spikerate = full psth of the data matrix
        meanrate = mean psth across variables/trials
        meanvar = variance of fullrate
        fanofactor = meanvar/meanrate
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
    
    
    meanrate = np.mean(spikerate,1) * (1/bw) / ntrials
    meanvar = np.var(spikerate,1) * (1/bw) / ntrials
    fanofactor = meanvar/meanrate
    spikerate = spikerate * (1/bw) / ntrials
    
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
    
    if plotting==1:
        f, (ax1,ax2,ax3,ax4) = plt.subplots(3,1,sharex=True)
        
        ax1.plot(meanisi,'k')
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
        
        ax4.hist(meanisi,color='k')
        ax4.set_title('Mean ISI hist')
        ax3.tick_params(direction='out')
        ax4.spines['top'].set_visible(False)
        ax4.spines['right'].set_visible(False)
        ax4.xaxis.set_ticks_position('bottom')
        ax4.yaxis.set_ticks_position('left')
        
    return fullisi,meanisi,varisi,fanofactor
#==========================================================================     
           