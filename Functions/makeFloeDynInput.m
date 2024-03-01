function makeFloeDynInput(fName,window,floeOutlines,floeCentroids)
% MAKEFLOEDYNINPUT generates an input file compatible with FloeDyn from a
% given floe field.
%
%   makeFloeDynInput(fName,window,floeOutlines,floeCentroids)
%   where fName is the full filename of the FloeDyn .hdf5 input file, 
%   window is a 1x4 vector defining the domain of the FloeDyn simulation
%   [x1,x2,y1,y2], and floeOutlines, floeCentroids are outputs from the
%   'placeFloes' function.
%   
%   CAUTION: If the file 'fName' already exists, this function will remove and
%   replace it.
%
%   S.D.Brenner, 2024

    %% Organize inputs

    numFloes = size(floeOutlines,1);
    floeXY = cell(1,numFloes);
    for n = 1:numFloes
        floeXY{n} = floeOutlines{n,1} - floeCentroids(n,:);
    end

    floe_states = zeros(9,numFloes);
    floe_states(1:2,:) = floeCentroids.';
    floe_states(8:9,:) = floeCentroids.';


    %% Save file

    % Remove any existing files with the same name
    if exist(fName,'file'), delete(fName); end

    % Create empty HDF5 file
    h5create(fName,'/floe_states',size(floe_states,1:3),'ChunkSize',size(floe_states,1:3));
    h5create(fName,'/window',length(window) );
    for k = 1:numel(floeXY)
        S = size(floeXY{k});
        floeLabel = sprintf('%u',k-1);
        h5create(fName, ['/floe_shapes/',floeLabel], S );
    end

    % Write data to file
    h5write(fName,'/floe_states',floe_states);
    h5write(fName,'/window',window);
    for k = 1:length(floeXY)
        floeLabel = sprintf('%u',k-1);
        h5write(fName, ['/floe_shapes/',floeLabel], floeXY{k} );
    end
    fprintf('Saved File: %s\n',fName)

end