function [StimOn,StimOff,StimThrsh] = findtrigs(trigger_chan,varargin)
% ----------[Lindex,Sindex] = findtrigs(trigger_chan,varargin)-------------
%
% Finds stim indexes from TTL stim channel given user-defined min/max 
% thresholds for stims. 
%
%               >>> INPUTS >>>
% Required:
%   trigger_chan: trigger data
% Optional:
%   Lthresh: min, max for stim (i.e. [.2 .6])
%   Fs: sampling rate...to return stim times in seconds
%   mpd: minimum distance between peaks in samples (default = 0)
%
%               <<< OUTPUTS <<<
%   StimOn: time (in samples) of stim
%   StimOff: time (in samples) of off-times of stim
%   StimIndex: indicies of stims found (useful if stim channel has
%       multiplexed stimuli)
%   StimThrsh: OPTIONAL output thresh for stims (to use in loop if desired)
%
%  by JMS, 09/29/2015.... for seizure detection analysis
% --------------------------------------------------------------------------

if nargin>1; StimThrsh = varargin{1};
else StimThrsh= []; end
if nargin>2; Fs = varargin{2};
else Fs = 0; end
if nargin>3; mpd = varargin{3};
else mpd = 0; end

% ask uesr for stim thresholds if StimTrhsh = []
if isempty(StimThrsh)
    StimThrsh = input('min/max of stim? ');
end

% convert min peak distance to sample points if provided
if Fs>0
    mpd = mpd*Fs;
end

% differentiate the stim to get rising-portion only
stimdiff = diff(trigger_chan); 

% search through stims and extract
if min(StimThrsh)>0
    
    % get stim indicies
    StimIndex = find(stimdiff>StimThrsh(1) & stimdiff<StimThrsh(2));
    
    % get the stim on/off times using stimindex...remove those less than
    % mpd if mpd is provided
    StimOn = StimIndex;
    StimOff = StimIndex;
    if mpd>0
        extra = find(diff(StimIndex)<mpd);
        StimOn(extra+1)=[]; % add 1 to keep first Stim
        StimOff(extra)=[]; % don't add 1 to keep first off stim
    end
else
    StimOn = [];
    StimOff = [];
    StimIndex = [];
end

% make StimIndex = 0 if empty
if isempty(StimIndex)
    disp('No stims found');
    StimIndex = 0;
end

% replace by time-values
if Fs>0
    StimOn = StimOn/Fs;
    StimOff = StimOff/Fs;
end

end



