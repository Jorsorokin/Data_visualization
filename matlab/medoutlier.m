function [idx,vals,filt] = medoutlier(x,sigthresh,pad)
%------------[idx,vals,filt] = medoutlier(x,sigthresh,pad)----------------
%   iterative outlier detection using median-absoulte deviation (MAD)
%   via the equation: 
%       scores = abs(x - median(x)) / mad(x) 
%
%   scores higher than a set threshold will be removed
%
%               >>> INPUTS >>>
% Required:   
%   x = data vector
% Optional:
%   sigthresh = threshold to remove outliers (default = 3)
%   pad = will pad output vector with nans for each removed data point to
%       preserve the length (if output vector "filt" is called)
%
%               <<< OUTPUTS <<<
%   idx = index of outliers
%   vals = values of the outliers
%   filt = new vector with outliers removed (optionally padded with nans)
%
% By JMS, 1/15/2016
%-------------------------------------------------------------------------

% check optionals
if nargin > 1 && ~isempty(varargin{1})
    sigthresh = varargin{1};
else sigthresh = 3; end
if nargin > 2 && ~isempty(varargin{2})
    pad = varargin{1};
else pad = 0; end

if ~ismember(pad,[0,1])
    pad = 0;
    disp('pad must be 0 or 1...reverting back to 0');
end

% create empty variables for outlier loop
idx = [];
valid = 0;
xloop = x; 

% begin the loop
while valid == 0;
    m = nanmedian(xloop);
    ma = mad(xloop);
    scores = abs(xloop-m) / ma;

    % detect outliers...break loop if none found
    idvec = scores > sigthresh;
    bad = find(idvec==1);
    if isempty(bad)
        valid = 1;
        break
    end
    
    % concatenate outlier indexes and reloop after replacing outliers with
    % nans in the 'xloop' variable    
    if ~isrow(bad)
        bad = bad';
    end
    idx = [idx,bad];
    xloop(idx) = nan;
end

% return values using original vector x
if nargout==2
    vals = x(idx);
elseif nargout==3
    vals = x(idx);
    filt = x;
    if pad == 1 %pad with NaN if pad == 1
        filt(idx) = nan;
    else
        filt(idx) = [];
    end
end
        
end
    