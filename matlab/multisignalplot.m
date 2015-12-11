function [fh,fl] = multisignalplot(d,varargin)
% -------- [fh,fl] = multisignalplot(d,varargin) ---------
% 
% Plots all of the signals in a data matrix in one figure, separated on the y-axis by a
% certain value. Note that this function will plot on top of a currently
% active figure or subplot.
%
%               >>> INPUTS >>>
% Required:
%   d = data matrix
% Optional:
%   Fs = sampling rate...if provided, will plot relative to time on x-axis
%   format = 'r','c'...lets the function know if your data is in row or
%       column format...by default assumes 'c'.
%   col = line color (default 'k');
%   maxval = value to separate lines by...default is mean(max(d)); 
%
%               <<< OUTPUTS <<<
%   fh = handle to gca
%   lh = handle to each line separately
%  
% Example:
%   vals = rand(1000,10) + repmat(2*sin(linspace(0,4*pi,1000))',1,10);
%   samplingrate = 500; 
%   [fh,fl] = multisignalplot(vals,samplingrate,[],'r');
%       % plots 10 red sinewaves corrupted with noise
%
% By JMS, 11/04/2015
% ----------------------------------------------------------------------------

% check optionals
if nargin>1; Fs = varargin{1}; 
else Fs = []; end
if nargin>2 && ~isempty(varargin{2}); format = varargin{2};
else format = 'c'; end
if nargin>3 && ~isempty(varargin{3}); col = varargin{3};
else col = 'k'; end
if nargin>4 && ~isempty(varargin{4}); maxval = varargin{4};
else maxval = []; end

% convert to column if row format
if strcmp(format,'r')
    d = d';
end

% get time vec if supplied sampling rate
if ~isempty(Fs)
    time = (1:size(d,1))/Fs;
else
    time = (1:size(d,1));
end

% get the mean maximum value of the data to separate lines by
n = size(d,2); 
if isempty(maxval)
    maxval = mean(max(d)); % mean of maximum value
end

% plot the values
hold on;
for i = 1:n
    fl(i) = plot(time,d(:,i)-maxval*(i-1),col);
end
fh = gca;
set(fh,'box','off','tickdir','out');
axis tight
end