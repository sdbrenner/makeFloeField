function [hitTestResult,boundaryTestResult] = hitTest( floePos,allFloes,boundary,periodicBCs )
    
    % start with no hits
    hitTestResult = 0;

    % check if any points are outside boundary
    boundaryTestResult = false;
    if any( ~inpolygon( floePos(:,1),floePos(:,2), ...
                        boundary(:,1),boundary(:,2) ))
        boundaryTestResult = true;
    end

    % Do we care about hitting the boundaries?
    if boundaryTestResult && ~periodicBCs
        hitTestResult = 1;
        return;
    end

    % are there other floes?
    if isempty(allFloes)
        return;
    end
        
    % If boundary is intersected, create "ghost" floes
    floeChecks{1} = floePos;
    if boundaryTestResult && periodicBCs
        floePosPBC = makeGhostFloes( floePos, boundary );
        numGhosts = numel(floePosPBC);
        ind = 1+(1:numGhosts);
        floeChecks(ind) = floePosPBC;
    end

    % check if any points are inside other floes:
    %   - first, check if boundary boxes around floes overlap (to identify
    %     "potential" overlaps, and reduce computation for floes far away)
    %   - then, loop through the floes with overlapping bounding boxes to
    %     check for actual intersections

    % get bounding boxes for all floes
    bboxFun = @(X) double([min(X(:,1)),max(X(:,1)),min(X(:,2)),max(X(:,2))]);
    bboxMat = NaN(numel(allFloes),4);
    for k = 1:numel(allFloes)
        bboxMat(k,:) = bboxFun( allFloes{k} );
    end

    % Loop through test floes (including ghosts)
    for p = 1:numel(floeChecks)
    
        % get bounding box for test floe     
        bboxTest = bboxFun( floeChecks{p} );
    
        % get list of floes whose bounding boxes overlap the test floe
        % (i.e., possible intersecting floes)
        isOverlapping = @(x,xr) max(x(1),xr(:,1)) <= min(x(2),xr(:,2));
        xOverlaps = isOverlapping(bboxTest(:,1:2),bboxMat(:,1:2));
        yOverlaps = isOverlapping(bboxTest(:,3:4),bboxMat(:,3:4));
        possibleHits = find( xOverlaps & yOverlaps );
        % isbetweenVec = @(x,xr) (x>=xr(:,1) & x<=xr(:,2));
        % possibleHits = find( ( isbetweenVec(bboxTest(:,1),bboxMat(:,1:2)) |...
        %                         isbetweenVec(bboxTest(:,2),bboxMat(:,1:2)) ) &...
        %                       ( isbetweenVec(bboxTest(:,3),bboxMat(:,3:4)) |...
        %                         isbetweenVec(bboxTest(:,4),bboxMat(:,3:4)) ) ).';
        % loop through list of possible intersecting floes and check 
        if ~isempty(possibleHits)
            for k = 1:numel(possibleHits) %1:length(allFloes)
                ki = possibleHits(k);
                p1x =  floeChecks{p}(:,1);
                p1y =  floeChecks{p}(:,2);
                p2x = allFloes{ki}(:,1);
                p2y = allFloes{ki}(:,2);
                2+2; %debugging stop
                if any( [inpolygon(p1x,p1y,p2x,p2y) ; inpolygon(p2x,p2y,p1x,p1y)] )
                    hitTestResult = 1;
                    break;
                end
            end
        end

    end


end
