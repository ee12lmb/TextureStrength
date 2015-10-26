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

% Read whole data file
[textures,ngrains,blocks] = read_VPSC(infile) ;

% * ADD INPUT CHECK THAT n

% Set up symmetry
CS = crystalSymmetry('Pbnm', [4.75, 10.20, 5.98]);
SS = specimenSymmetry('-1');


%% Generate random numbers 

rng(seed);                      % set seed to repeat analysis
sample = randperm(ngrains,n);   % random permutation of 1:ngrains, n times
                                % to generate index to pull samples from
                                % array
                                
% * Possible mod: change grains for each time step? *

%% Sample data and calculate J-index 

% pull out samples for each time step
for i = 1:blocks-1 % deal with repeated texture
   
    eulers_r = textures{i}*degree;
    
    for j = 1:n % pull desired number of samples 
        % sample from this block 
        eulers_r_samp(1,j) = eulers_r(1,sample(j));
        eulers_r_samp(2,j) = eulers_r(2,sample(j));
        eulers_r_samp(3,j) = eulers_r(3,sample(j));
    end     
    
    % calculate J index for this block and store
    g = orientation('Euler', eulers_r_samp(1,:), eulers_r_samp(2,:), eulers_r_samp(3,:), ...
        CS, SS, 'Bunge');
    odf = calcODF(g,'HALFWIDTH', 10*degree, 'silent');
    J(i) = textureindex(odf);
    
end

%% Build output

% strain vector - * need to update to pull this from file *
strain = linspace(0.02,0.5,25);

for i = 1:length(J)
    fprintf('%f\t%f',J(i),strain(i))
end

end

