function [ J, strain ] = J_ngrains(infile,n,seed)
%J_ngrains calculates the J index for a specified number of grains
%   A number of grains (n) are selected at random from a VPSC input file.
%   From this the J index is then calculated for each strain step present in the
%   input file.

%% Setup & read data file

% set up MTEX package
addpath /nfs/see-fs-01_teaching/ee12lmb/project/source/mtex-4.1.3/
startup_mtex;

% add path to read files
addpath /nfs/see-fs-01_teaching/ee12lmb/project/source/dev/readfiles/

% Sample data file
[textures,ngrains,strain,blocks] = sample_VPSC(infile,n,seed) ;

% Set up symmetry
CS = crystalSymmetry('Pbnm', [4.75, 10.20, 5.98]);
SS = specimenSymmetry('-1');


%% Sample data and calculate J-index 

% pull out samples for each time step
for i = 1:blocks % deal with repeated texture
   
    eulers_r = textures{i}*degree; % textures is SAMPLED TEXTURE
    
    % calculate J index for this block and store
    g = orientation('Euler', eulers_r(1,:), eulers_r(2,:), eulers_r(3,:), ...
        CS, SS, 'Bunge');
    odf = calcODF(g,'HALFWIDTH', 10*degree, 'silent');
    J(i) = textureindex(odf);
    
end

%% Build output

% strain vector - * need to update to pull this from file *

% 
% for i = 1:length(J)
%     fprintf('%f\t%f',J(i),strain(i))
% end

end

