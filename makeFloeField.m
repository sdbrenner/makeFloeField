%% MAKE FLOE FIELD

%% Clean workspace
clear;
clc;
close all;

%% Define boundaries and targets
% EDIT VALUES IN THIS SECTION BASED ON DESIRED TARGETS
% (Try different examples by commenting/uncommenting)

% Define domain window (for FloeDyn)
% in the form: [x1,x2,y1,y2] the x,y extents
x1 = 0; 
y1 = 0; 
x2 = 25e3;
y2 = 25e3;
window = [x1,x2,y1,y2];

% Define boundary within which the floea are placed
% in the form: [X,Y] where X,Y are point pairs that define a polygon
% EXAMPLE 1: boundary equal to domain
bx = [x1,x1,x2,x2,x1];
by = [y1,y2,y2,y1,y1];
boundary = [ bx(:),by(:) ];
%
% EXAMPLE 2: boundary equal to half-domain
% bx = [x1,x1,x2,x2,x1];
% by = [(y2-y1)/2,y2,y2,(y2-y1)/2,(y2-y1)/2];
% boundary = [ bx(:),by(:) ];
%
% EXAMPLE 3: circular boundary
% r = 0.5*min(x2-x1,y2-y1);
% th = linspace(0,2*pi,180);
% [bx,by] = pol2cart(th,r);
% boundary = [ bx(:),by(:) ]+0.5*[x2-x1,y2-y1];


% Calculate area enclosed by boundary 
% (used for calculating target sea ice area)
Abound = polyarea( boundary(:,1), boundary(:,2) );

% Sea ice concentration target (fractional):
targetSIC = 0.5;
targetIceArea = targetSIC*Abound;


% Define FSD size range and properties
minSize = 100;    % minimum floe size
maxSize = 1000;   % maximum floe size
% EXAMPLE 1: 
fsdType = "powerlaw";
fsdParams = {2};
%
% EXAMPLE 2: 
% fsdType = "lognormal"; 
% modeSize = 2*minSize;
% sig = 1+1/3;
% mu = sig^2 + log( modeSize-minSize );
% fsdParams = {mu,sig};
% Use `help getFloeSizes` for more information about FSD inputs

%% Load floe shape inventory

% Load default FloeDyn inventory of floe shapes 
load('floeShapes.mat');

%% Generate FSD and place floes

% GET LIST OF FLOE SIZES
floeSizes = getFloeSizes( fsdType,[minSize,maxSize],fsdParams,targetIceArea  );
% Displace message with FSD/SIC properties
fprintf('FSD generated with %g floes from %2.1f m to %2.1f m\n\n',...
        numel(floeSizes), min(floeSizes), max(floeSizes) );    

% PLACE FLOES
% NOTE: `placeFloes` has a range of optional arguments; 
% Use `help getFloeSizes` for more information.
% EXAMPLE 1: Default options
% [floeOutlines,floeCentroids,iceMask] = placeFloes( boundary, floeSizes, G );
%
% EXAMPLE 2: Periodic boundaries
[floeOutlines,floeCentroids,iceMask] = placeFloes( boundary, floeSizes, G, periodicBCs=true );

%% Visualize results 

% Subset to non-empty values (accounting for "ghost floes" in periodic BCs)
% and convert to 'polyshape' array
ind = find( cellfun( @(C) ~isempty(C), floeOutlines ) );
floeField = cellfun(@(C) polyshape(C,KeepCollinearPoints=true), floeOutlines(ind) );

% Get window aspect ratio
aspectRatio = (y2-y1)/(x2-x1);

% Plot
fH = figure(1); clf;
fH.Position(4) = aspectRatio*fH.Position(3);

ax = gca;
hold on; 
plot( window([1,1,2,2,1]), window([3,4,4,3,3]),'--',...
      Color=0.25*[1,1,1],LineWidth=1.5 ); 
plot( boundary(:,1),boundary(:,2),'--',...
      Color=0.65*[1,1,1],LineWidth=1.5 ); 
plot( floeField ); 
plot( floeCentroids(:,1),floeCentroids(:,2),'k.',MarkerSize=1 );

ax.XLim = [x1,x2] + 0.1*[-1,1]*(x2-x1);
ax.YLim = [y1,y2] + 0.1*[-1,1]*(y2-y1);
ax.Position = [0,0,1,1];
axis off;
daspect([1,1,1]);
drawnow;


%% Make corresponding input file for FloeDyn

fileName = 'floeDynInput.h5';
makeFloeDynInput( fileName,window,floeOutlines,floeCentroids );

