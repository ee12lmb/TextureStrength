function [ J, strain ] = J_ngrains(input_texture,n,seed)
%J_ngrains calculates the J index for a specified number of grains
%   A number of grains (n) are selected at random from a VPSC input file.
%   From this the J index is then calculated for each strain step present in the
%   input file.

%% Setup & read data

% set up MTEX package
addpath /nfs/see-fs-01_teaching/ee12lmb/project/source/mtex-4.1.3/
startup_mtex;

% add path to read files
addpath /nfs/see-fs-01_teaching/ee12lmb/project/source/dev/readfiles/

% check if input is raw VPSC or texture array
if (ischar(input_texture) == 1)
    
    % if a file path is given then read in
    [textures,ngrains,strain,blocks] = sample_VPSC(input_texture,n,seed);
    
else
    
   % if input is not a file then pass to sample_texture to deal with
   [textures,blocks] = sample_texture(input_texture,n,seed);
   
   % strain information cannot be extracted from inputted texture
   %+but should already be known from previous read_VPSC
   strain = 'Input is texture - strain already extracted'; 
    
end

% Set up symmetry
CS = crystalSymmetry('Pbnm', [4.75, 10.20, 5.98]);
SS = specimenSymmetry('-1');

    
%% Calculate J-index 

% check how many textures are given (i.e. do we need to reference cell
% array or just normal matrix)
if (blocks == 1)
  
    eulers_r = textures*degree; % textures is SAMPLED TEXTURE

    % calculate J index for this block and store
    g = orientation('Euler', eulers_r(1,:), eulers_r(2,:), eulers_r(3,:), ...
        CS, SS, 'Bunge');
    odf = calcODF(g,'HALFWIDTH', 10*degree, 'silent');
    J = textureindex(odf);  
    
    
else % now deal with multiple textures (e.g. cell arrays)
    
    % pull out samples for each time step
    for i = 1:blocks % deal with repeated texture

        eulers_r = textures{i}*degree; % textures is SAMPLED TEXTURE

        % calculate J index for this block and store
        g = orientation('Euler', eulers_r(1,:), eulers_r(2,:), eulers_r(3,:), ...
            CS, SS, 'Bunge');
        odf = calcODF(g,'HALFWIDTH', 10*degree, 'silent');
        J(i) = textureindex(odf);

    end
    
end

%% Build output

% strain vector - * need to update to pull this from file *

% 
% for i = 1:length(J)
%     fprintf('%f\t%f',J(i),strain(i))
% end

end

