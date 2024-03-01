% function [sampleXY] = invTransSample2D( X,Y,P,N )
function [sampleXY] = randP2D( X,Y,P,N )
    % 2D inverse transform sample, selects a set of N random points that
    % follow the 2D probability distribution P

    % Add small random noise to P matrix:
    noise = 1-0.03*rand(size(P));
    P = double(P).*noise;
    % Generate 1D CDF
    CDF = [0;  cumsum( P(:)/sum(P(:)) ) ]; %first term out of of cumsum is not zero
    x = linspace(0,1,numel(CDF));
    % % normalize CDF
    % CDF = CDF./max(CDF);
    
    % Ensure only unique CDF values (which exist because of P=0 regions)
    udx = ~(P(:)==0);
    CDF = CDF(udx);
    x = x(udx);   

    % Generate random values and invert CDF
    randVals = rand(N,1);  
    outVals = interp1( CDF,x,randVals ); % spans zero to one
    
    % Interpolate to X,Y pair values
    xs = interp1( linspace(0,1,numel(X)), X(:), outVals );
    ys = interp1( linspace(0,1,numel(Y)), Y(:), outVals );
    
    % 'jiggle' the x-values
    dx = mean(diff(X(1,:)));
    xs = xs + (dx)*(rand(N,1)-0.5);
    xs(xs<min(X(:))) = min(X(:));
    xs(xs>max(X(:))) = max(X(:));

    % Create output
    sampleXY = [xs(:),ys(:)];

end