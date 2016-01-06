function [hL,hF] = fillPlot(mat,varargin)
%
%----------------------- [hL,hF] = fillPlot(mat,varargin) ------------------
%
% Takes in a matrix of n x m values and plots the mean across all trials (mean(mat),1) 
% and error boundaries with fill between lower/upper bounds.  
%
%
%                       >>> INPUTS >>>
%   REQUIRED:
%       mat = matrix of data (nxm matrix, with n = num trials/variables, m = samples over time/space)
%           *this should be in row-vector form*
%   OPTIONAL:
%       tm = time vector (default = 1:length(mat))
%       err = standard error ('sem'),standard deviation ('sd'), 95% confidence interval ('ci')
%           (default = 'sem')
%       lineCol = line color of mean value (default = 'k')
%       edgeCol = edge colors of boundaries (default = 'none')
%       faceCol = face color of fill (default = [.3 .7 1])
%       fAlpha = transparency of shade (default = 1 (opaque));
%           * set optinonal inputs = [] to use skip and use default *
%
%                       <<< OUTPUTS <<<
%   OPTIONAL:
%       hL = handle to line 
%       hF = handle to fill 
%
%
%   Example:
%       figure; hold on;
%       data1 = rand(10,1000); 
%       data2 = 0.5*rand(10,1000) + repmat(0.25*sin(linspace(0,4*pi,1000)),10,1);
%       err = 'sd';
%       lineCol = 'b';
%       lineCol2 = 'k';
%       faceCol = [.8 0 .2];
%
%       [hL,hF] = fillPlot(data1,[],err,lineCol);
%       [hL2,hF2] = fillPlot(data2,[],err,lineCol2,[],faceCol,.5);
%           
%----------------------------------------------------------------------------------


% check optional inputs
if nargin<7;fAlpha=1; % default no transparency
else fAlpha=varargin{6};end 
if nargin<6;faceCol=[.3 .7 1]; % default light blue
else faceCol=varargin{5};end
if nargin<5;edgeCol='none'; % default no edge color
else edgeCol=varargin{4};end
if nargin<4;lineCol='k'; % default black line
else lineCol=varargin{3};end
if nargin<3;err='sem'; % default SEM 
else err=varargin{2};end
if nargin<2;tm=1:size(mat,2);% default y-values = samples
else tm=varargin{1};end
if nargin<1;error('Must input matrix');end

% set defaults if empty
if isempty(tm);tm=1:length(mat);end
if isempty(err);err='sem';end
if isempty(lineCol);lineCol='k';end
if isempty(edgeCol);edgeCol='none';end
if isempty(faceCol);faceCol=[.3 .7 1];end
if isempty(fAlpha);fAlpha=1;end 

% descriptive calculations
y = mean(mat);
SD = std(mat);
SEM = SD/sqrt(size(mat,1));
CI = 1.96*SEM; % 95 percent confidence interval

% time vector
if ~isrow(tm)
    tm=tm';
end

% error lines
switch err
    case 'sem'
        low = y-SEM;
        hi = y+SEM;
    case 'sd'
        low = y-SD;
        hi = y+SD;
    case 'ci'
        low = y-CI;
        hi = y+CI;
    otherwise
        low = y-SEM;
        hi = y+SEM;
end

% plot the data
hF = fill([tm fliplr(tm)],[hi fliplr(low)],faceCol);
hold on;
set(hF,'edgecolor',edgeCol,'faceAlpha',fAlpha);
hL = plot(tm,y,'color',lineCol,'linewidth',2);
set(gca,'tickdir','out','layer','top');
set(gcf,'color','w');
end
    