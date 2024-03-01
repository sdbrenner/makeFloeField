function [floeOutlines,floeCentroids,iceMask,successFlag] = placeFloes( boundary, floeSizes, floeInv, opts )
% PLACEFLOES generates initial floe placements for a given floe-size list
% for use as initial conditions for FloeDyn
%   
%   [floeXY,floe_states] = placeFloes( boundary, floeSizes, floeInv )
%
% Floes are randomly selected from the inventory, then randomly rotated and
% placed within the boundary starting from the largest and working down the
% list. The random placement location follows a probability distribution
% that increases with increasing distance away from other floes (which
% requires binarizing the floe field on some specified grid).
% 
% After each placement, the floe position is rejected if it overlaps
% with another floe or the edge of the boundary. To avoid infinite loops, a
% limit is placed on the number of "rejections" allowed for any floe--if
% the limit is reached then placement fails. 
% 
% Convergence is not guaranteed, and is generally difficult for high
% concentrations and/or small size spans.
%
% The variable 'boundary' is a set of [X,Y] values in which the floes are
% placed. Any arbitrary polygon is allowed for boundary; however, for
% non-rectangular boundaries periodic behaviour (allowed through optional
% 'periodicBCs') is not well posed and weird things will happen.
%
%   [...] = placeFloes( boundary, floeSizes, floeInv,'Name','Value' )
%
% Allows for optional name/value pair arguments:
%   'gridRes'           [ default: min(floeSizes)/5 ]
%   'attemptsPerFloe'   [ default: 5e3 ]
%   'floeResolution'    [ default: min(floeSizes)/2 ]
%   'periodicBCs'       [ default: true ]
%   'includeProgress'   [ default: true ]
%
% S.D.Brenner 2022 (updated Feb. 2024)

    %% Parse input data
    arguments
        boundary             (:,2) {mustBeNumeric}
        floeSizes            (:,1) {mustBeNumeric}
        floeInv 
        opts.gridRes         (1,1) {mustBeNumeric} = min(floeSizes)/5;
        opts.attemptsPerFloe (1,1) {mustBeNumeric} = 5e3;
        opts.floeResolution  (1,1) {mustBeNumeric} = min(floeSizes)/2; 
        opts.periodicBCs     (1,1) = false;
        opts.includeProgress (1,1) = true;
    end

    % extract number of floes
    numFloes = numel(floeSizes);
    
    % extract domain information (for creating binary ice mask)
    X = floor( min(boundary(:,1)) );
    W = ceil( range(boundary(:,1)) );
    Y = floor( min(boundary(:,2)) );
    H = ceil( range(boundary(:,2)) );


    %% Set up loop

    % Set up initial floe placement probability matrix   
    nx = ceil(W/(opts.gridRes)) + 1;
    ny = ceil(H/(opts.gridRes)) + 1;
    fx = linspace(X,X+W,nx);
    fy = linspace(Y,Y+H,ny);
    [FX,FY] = meshgrid(fx,fy);
    Pb = binMask( boundary,fx,fy); % boundary mask
    P = Pb;
   
    % Pre-allocate variables
    iceMask = zeros(size(FX));
    floeOutlines = cell(numFloes,9);
    floeCentroids = NaN(numFloes,2);
    
    % Set up loop info
    aptRecord = NaN(1,numFloes);
    updateLabIncrement = 0.1;
    updateLab = updateLabIncrement;
    attemptsPerFloe = opts.attemptsPerFloe;
    successFlag = false;

    %% Loop through floes
    if opts.includeProgress
        disp('placing floes...')
    end
    for k = 1:numFloes

        % Generate possible floe positions
        fsXYCandidates = randP2D( FX,FY,P,attemptsPerFloe );
%         fsXYCandidates = [X,Y] + [W,H].*rand(attemptsPerFloe,2);

        % loop through 'attemptsPerFloe'
        for j = 1:attemptsPerFloe

            % Get random floe from inventory
            floeK = selectAndScale( floeInv,floeSizes(k),opts.floeResolution );

            % Generate random floe position and orientation
            fsXY = fsXYCandidates(j,:); %invTransSample2D( FX,FY,P,1 );
            theta = 2*pi*rand(1);            
            R = [ cos(theta), sin(theta) ;
                 -sin(theta), cos(theta) ];

            % Check that the floe doesn't overlap with any other floes
            floePosTest = (floeK*R)+fsXY;
            allFloes = floeOutlines( cellfun( @(C) ~isempty(C), floeOutlines ) );
            % allFloes = floePos(1:k-1,1); 
            [hitTestResult,boundaryTestResult] = hitTest( floePosTest,allFloes,boundary,opts.periodicBCs );
            % [hitTestResult,boundaryTestResult] = hitTest( floePosTest,floePos(1:k-1),boundary,...
            %                             fsXY,floeSizes(k),floe_states(:,1:k-1),floeSizes(1:k-1),k );

    
            % If floe position is good, accept placement and advance loop
            if ~hitTestResult 
                
                % Put floe into position
                floeOutlines{k,1} = floePosTest;
                floeCentroids(k,:) = fsXY;

                % Put extra floe ghosts into PBC matrix
                if boundaryTestResult && opts.periodicBCs
                    floePosPBC = makeGhostFloes( floePosTest, boundary );
                    numGhosts = numel(floePosPBC);
                    ind = 1+(1:numGhosts);
                    floeOutlines(k,ind) = floePosPBC;
                end
    
                % update sampling probability matrix
                % (probability increases with the distance from existing
                % floes)
                numFloeK = sum(cellfun(@(C) ~isempty(C), floeOutlines(k,:) ));
                for p = 1:numFloeK    
                    floeMask = binMask( floeOutlines{k,p}, fx,fy );
                    iceMask = iceMask+floeMask;
                end
                % Get matrix D of distances to the outside of existing floes
                D = opts.gridRes*bwdist( iceMask );
                P = Pb.*D.^2;
                % decrease the probability within a floeSize/2 distance
                % away from existing floes (minimized overlaps):
                ind =  D<=0.5*floeSizes( min([numFloes,k+1]) );
                P(ind) = 0.25*P(ind);

                % report progress
                aptRecord(k) = j;
                if (k/numFloes)>=updateLab && opts.includeProgress
                    fprintf('%02g%% placed...\n',100*updateLab );
                    updateLab = updateLab+updateLabIncrement;
                end
                % return to main loop
                break;
            elseif hitTestResult &&  j == attemptsPerFloe
                % If at the end of the 'attemptsPerFloe' loop without
                % placement, then give up
                warning('Could not place all floes. Try a wider size span or lower SIC.');
                return;
            end
        end  
    end

    if opts.includeProgress
        fprintf('floe placement succeeded\n\n');
    end
    successFlag = true;
    

end













