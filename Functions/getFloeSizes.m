function floeSizes = getFloeSizes( fsdType,sizeRange,fsdParams,targetIceArea )
% GETFLOESIZES generates a list of floe sizes following a prescribed
% (truncated) FSD, such that the total surface area of the floes satisfies
% some target value.
%
%   floeSizes = getFloeSizes( fsdType,sizeRange,fsdParams,targetIceArea )
%   
%   where fsdType specifies the FSD function (see below); sizeRange is a
%   2-element vector [rMin,rMax]; fsdParams is a cell array containing
%   relevant parameters for each distrubtion (see below); and targetIceArea
%   is the target total surface area covered by the ice floes
%
%   fsdTypes and corresponding fsdParams are:
%
%       fsdTypes            fsdParams
%       ----------------    ----------------
%       "powerlaw"          {alpha}
%       "lognormal"         {mu,sigma}
%       "uniformdist"       {}
%       "singlefloesize"    {}
%       "custom"            {function handle for FSD PDF}
%
%   The output list `floeSizes` is sorted in descending order.
%
%   S.D.Brenner, 2024

    % Parse input data
    % arguments
    %     fsdType
    %     sizeRange
    %     fsdParams 
    %     totalIceArea
    % end
    
    % Extract min/max floe sizes
    rMin = sizeRange(1);
    rMax = sizeRange(2);
    
    % Check if rMin and rMax are equal:
    if rMin==rMax && ~strcmp(fsdType)
        fsdType = "singlefloesize";
        warning('Minimum and maximum floe sizes are equal. Using "singlefloesize" fsd type');
    end

    % Generate semi-continuous floe size number distribution, N(r), and
    % where:
    %   N(r) is the number of floes *per unit area* of size r, and is
    %   normalized such that: ∫πr^2 N(r)dr = totalIceArea; 
    %   and r is the effective floe radius, r = √(A/π) or r = √A
    A = @(r) pi*r.^2;
    fsdType = lower(fsdType);
    switch fsdType
        case "powerlaw"
            [alf] = deal(fsdParams{:});
            pfun = @(r) r.^-alf;
        case "lognormal"
            [mu,sig] = deal(fsdParams{:});
            pfun = @(r) (1./(r-rMin)).*exp( -(log(r-rMin)-mu).^2/(2*sig^2) );
        case "uniformdist"
            pfun = @(r) 1+0*r;
        case "singlefloesize"
            numFloes = round( targetIceArea/A(rMin) );
            floeSizes = repelem( rMin,numFloes);
            return;
        case "custom"
            pfun = fsdParams{1};
    end
    
    % Generate FSD functions  
    P = @(r) pfun(r) /integral( pfun,rMin,rMax );
    N = @(r) targetIceArea *pfun(r) ./integral( @(r) A(r).*pfun(r),rMin,rMax );

    % Estimate number of floes
    numFloes = round( integral(N,rMin,rMax ) ); 
    
    % Numerically revise estimate (to account for discretization errors):
    % Generate candidate floe sizes
    nF = floor(0.95*numFloes):ceil(1.05*numFloes);
    M = numel(nF);
    iceArea = zeros(1,M);
    for m = 1:M
        % Generate floe sizes by numerically inverting the CCDF
        k = linspace(0,nF(m),nF(m));    
        Fval = 1-k/nF(m);
        rval = invertCCDF( P, Fval, sizeRange );
        % Get total ice area for given distribution of floes
        iceArea(m) = sum( A(rval) );
    end
    % interpolate to find the number of floes that gives an ice area that
    % is closest target value
    numFloes = round( interp1( iceArea, nF, targetIceArea,"pchip",'extrap' ) );


    % Generate floe sizes by numerically inverting the CCDF
    k = linspace(0,numFloes,numFloes); 
    Fval = 1-k/(numFloes);
    r = invertCCDF( P, Fval, sizeRange );
    floeSizes = sort(r,"descend");


    % Re-adjust floe sizes if necessary
    floeSizes = floeSizes*sqrt( targetIceArea/sum(A(floeSizes)) ); 


end


function rval = invertCCDF( P, Fval, sizeRange )
    
    rMin = sizeRange(1);
    rMax = sizeRange(2);
    r = linspace( rMin, rMax, 1e3 );
    Pr = P(r); 
    Pr(isnan(Pr)) = 0; % accounts for lognormal weirdness
    F = 1-cumtrapz(r,Pr);
    % subset to unique points
    [~,udx] = unique(F);
    F = F(udx); r = r(udx);
    % invert
    rval = interp1( F,r,Fval,'linear','extrap' );
    
end