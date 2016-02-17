function [BI,Ib,It,Sb,St] = burstindex(varargin)
%
% ----------- [BI,Ib,It,Sb,St] = burstindex(varargin) -----------------
%
%   Calculates the burst index of a rebound-bursting neuron from a CCIV ranging from
%   negative to postivie currents. Assumes that the CCIV protocol sweeps
%   through a large enough range of currents to elicit rebound bursts and
%   tonic firing. 
%
%   Burst index is defined as follows:
%   
%       BI = 1 - (B / (B+T))
%           where 
%               B = abs(Ib)/sqrt(Sb)
%               T = abs(It)/sqrt(St)
%           and 
%               Ib = threshold current for eliciting rebound burst    
%               It = threshold current for eliciting tonic spikes
%               Sb = (# spikes in rebound) / (stim length in seconds)
%               St = (# of tonic spikes) / (stim length in seconds)
%                   * dividing by stim length normalizes the spike counts
%
%   We take the square root of St and Sb to avoid excessive weighting of spike
%   rate in the equation, and instead emphasize the threshold currents. 
%   This is due to the physiology of rebound bursts, which typically plateau 
%   at 4-5 spikes despite increasing hyperpolarizing current, while tonic
%   spikes plateau at much higher counts, and depend on the stim length
%   
%
%               >>> INPUTS >>>
% Optional:
%   stimon = on time of pulse from CCIV protocol in seconds
%   stimoff = off time of pulse from CCIV protocol in seconds
%   Spiketimes = spiketimes from CCIV in column format.
%       units of spiketimes needs to equal units of stimon/stimoff
%   current = current levels used in CCIV protocol
%       * if you provide these 4 optional inputs, this function
%       will calculate the parametrs Ib, It, Sb, and St for you.
%
%   If you already have [Ib,It,Sb,St] parameters calculated, 
%   leave the 1st 4 inputs blank and proceed with the following four inputs. 
%
%   Ib = threshold burst current
%   It = threshold tonic current
%   Sb = # burst spikes / length of pulse (seconds)
%   St = # tonic spikes / length of pulse (seconds)
%
%               <<< OUTPUTS <<<
%   BI = burst index, ranging between [0,1]. Higher values = more "bursty"
%   Ib = threshold burst current (optional)
%   It = threshold tonic current (optional)
%   Sb = normalized # rebound spikes (optional)
%   St = normalized # tonic spikes (optional)
%
% By JMS, 02/12/2016
%
%-------------------------------------------------------------------------------------

% check inputs
if nargin < 4
    error('need at least 4 inputs...refer to document')
end

if nargin>0 && ~isempty(varargin{1}); stimon = varargin{1}; end
if nargin>1 && ~isempty(varargin{2}); stimoff = varargin{2}; end
if nargin>2 && ~isempty(varargin{3}); spiketimes = varargin{3}; end
if nargin>3 && ~isempty(varargin{4}); current = varargin{4}; end
if nargin>4 && ~isempty(varargin{5}); Ib = varargin{5}; end
if nargin>5 && ~isempty(varargin{6}); It = varargin{6}; end
if nargin>6 && ~isempty(varargin{7}); Sb = varargin{7}; end
if nargin>7 && ~isempty(varargin{8}); St = varargin{8}; end

%% calculate equation parameters if presets not provided
if ~exist('Ib','var') % could be any of the four...so long as one is absent, assumes you want to calculate the parameters below
    % get stim duration
	stimdur = stimoff-stimon;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % rebound parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    neg_curr = current < 0; 
    if max(neg_curr) > 0
        
        % get # spikes that occur after the end of the stim for the neg-current sweeps only
        n = sum(spiketimes(:,neg_curr) > stimoff);
        ind = max(find(n > 0));
        
        % get Ib, Sb if ind not empty
        if ~isempty(ind)
            c = current(neg_curr);
            Ib = c(ind); % threshold current
            Sb = n(ind); % # of spikes
        else
            Ib = 0;
            Sb = 1;
        end
    else
        Ib = nan;
        Sb = nan;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % tonic parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    pos_curr = current > 0;
    if max(pos_curr) > 0
        
        % get tonic spikes and find first current needed
        n = sum(spiketimes(:,pos_curr) > stimon & spiketimes(:,pos_curr) < stimoff);
        ind = min(find(n > 0));
        
        if ~isempty(ind)
            c = current(pos_curr);
            It = c(ind); % threshold current
            St = n(ind); % # of spikes
        else
            It = 0;
            St = 1;
        end
    else
        It = nan;
        St = nan;
    end    
end

%% calculate burst index
% correct for any zeros in Sb or St, which will create a NaN value in BI
Sb(Sb==0) = 1;
St(St==0) = 1;

% calculate BI
if ~max(isnan(Ib)) && ~max(isnan(It)) % if either nan, means CCIV current does not cover both negative and positive values
    % take absolute values
    Ib = abs(Ib); It = abs(It);

    % compute B and T 
    B = Ib ./ sqrt(Sb/stimdur); % use ./ notation in case user provides vectors, not scalars
    T = It ./ sqrt(St/stimdur);

    % compute BI
    BI = 1 - (B ./ (B+T));
    
    % correct for extremes
    BI(BI==0) = 2; % only way for BI = 0 is if Ib = 0, It > 0...put 2 as a place holder
    BI(BI==1) = 0; % only way for BI = 1 is if It = 0, Ib > 0
    BI(BI==2) = 1; % now recorrect the first condition
else
    error('CCIV protocol does not cover a large enough current range');
end
