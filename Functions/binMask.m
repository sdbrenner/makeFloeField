function mask = binMask( XY, fx, fy )
% BINMASK Create binary mask for set of boundary points in XY

    XY = double(XY);

    % Parse domain size/grid data
    Xf = min(fx);
    Yf = min(fy);
    Wf = range(fx);
    Hf = range(fy);
    nx = length(fx);
    ny = length(fy);
    dx = Wf/(nx-1);
    dy = Hf/(ny-1);
    nxy = [nx,ny];

    % Get normalized vertices for floe polygon
    vnX = [ (XY(:,1)-Xf)./dx + 1; NaN ]  ;
    vnY = [ (XY(:,2)-Yf)./dy + 1; NaN ]  ;
    % Create combined [X,Y] vertex matrix and binanarize
    vert = [vnX(:),vnY(:)];

    % Make empty mask grid
    mask = false( nxy(2),nxy(1) );

    % Extract subset of domain and create mask
    % (faster than calling poly2mask on whole domain [?])
    extMin = max( [floor( min(vert));[0,0]] );
    extMax = min( [ceil(max(vert));nxy] );
    nc = extMax-extMin;
    vertN = (vert-extMin);
    % maskN = mpoly2mask(vertN,fliplr(nc));
    maskN = poly2mask( vertN(1:end-1,1),vertN(1:end-1,2),nc(2),nc(1) );

    % Write back to main grid
    xInd = extMin(1)+(1:nc(1));
    yInd = extMin(2)+(1:nc(2));
    mask(yInd,xInd) = maskN;
end