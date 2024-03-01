function floePosPBC = makeGhostFloes( floePos, boundary )
% makeGhostFloes creates copies of floe outlines 

    % Get boundary info
    X = min(boundary(:,1));
    Y = min(boundary(:,2));
    W = range( boundary(:,1) );
    H = range( boundary(:,2) );

    % Identify which boundaries are interesected:
    xIdx = [ any(floePos(:,1)>(X+W)), 0, any(floePos(:,1)<X ) ];
    yIdx = [ any(floePos(:,2)>(Y+H)), 0, any(floePos(:,2)<Y ) ];

    % Define offset matrices defaults
    xOff = W*reshape(xIdx.*repmat(-1:1,3,1) ,1,[]);
    yOff = H*reshape((yIdx.*repmat(-1:1,3,1))',1,[]);
    off = [xOff(:),yOff(:)];
    % Retain only nonzero offsets and no duplicates
    off = off( any(off~=0,2),: );
    off = unique(off,'rows');

    % Generate offset floes
    nr = size( floePos,1);
    floePosPBC = floePos+permute(off,[3,2,1]);
    % Convert to cell array
    numGhosts = size(floePosPBC,3);
    floePosPBC = squeeze(mat2cell(floePosPBC, nr, 2, ones(1,numGhosts) ));

end