function floeK = selectAndScale( floeShapes,floeSize,floeResolution )
% SELECTANDSCALE select random floe from inventory and rescale

    % Resample outline
    minVerts = 6; % minimum number of vertices to describe floe
    M = length(floeShapes);
    m = randi(M,1);
    T = size(floeShapes{m},1);
    N = max( [ round(5*(floeSize/floeResolution)), minVerts+1] );
    t = 1:T; tp = linspace(1,T,N);
    fx = interp1(t,floeShapes{m}(:,1),tp).';
    fy = interp1(t,floeShapes{m}(:,2),tp).';
    
    % Get centroid location and area of inventory floe
    pgon = polyshape( fx(1:N-1), fy(1:N-1), KeepCollinearPoints=true );
    [fcx,fcy] = centroid(pgon);
    floeArea = area(pgon);

    % Center
    F(:,1) = fx(1:N-1)-fcx;
    F(:,2) = fy(1:N-1)-fcy;
    
    % Scale
    A = (pi*floeSize.^2);
    floeK = F * sqrt( A/floeArea );
end