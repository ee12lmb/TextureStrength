function [sample_texture,ngrains,strain,blocks] = sample_VPSC(infile,n,seed)
% SAMPLE_VPSC randomly samples from an input VPSC file
%
%   Input texture (VPSC file path) is read in using read_VPSC. Following
%   this, n random grains (without replacement) are taken from the total of
%   N Euler angle sets. Meta data is also returned.
%
%   Inputs:  infile          - VPSC file path
%            n               - number of grains to sample
%            seed            - set the same to repeat the same random no.s
%
%   Outputs: sample_texture  - nx3 matrix of sampled Euler angles
%            ngrains         - total number of grains in the sample N
%            strain          - vector of strain at each step
%            blocks          - number of strain stesps
%
%   Lewis Bailey - University of Leeds, School of Earth and Environment 
%   2015-16 Undergraduate final year project
%
%   Usage: [sample_texture,ngrains,strain,blocks] = sample_VPSC(infile,n,seed)
%
%   See also: READ_VPSC, SAMPLE_TEXTURE, READ_EBSD, SAMPLE_EBSD


%% Read VPSC & set up sample index

addpath /nfs/see-fs-01_teaching/ee12lmb/project/source/dev/
setup_env;

% Read whole data file
[textures,nxtl,strain,blocks] = read_VPSC(infile);

% check that nxtl is same for all time steps, and assign one scalar
assert((range(nxtl) == 0),'Number of grains not consitent across time steps!')
ngrains = nxtl(1);

% check that n is not greater than number of grains in sample
assert((n <= ngrains), ...
    'Number of samples requested is greater than total number in sample!')

% generate sample index
rng(seed);                      % set seed to repeat analysis
sample = randperm(ngrains,n);   % sample without repeating


%% Extract textures  

% check how many blocks (i.e. do we need cell array?)
if (blocks == 1)  
    
    for j = 1:n   
        sample_texture(1,j) = textures(1,sample(j));
        sample_texture(2,j) = textures(2,sample(j));
        sample_texture(3,j) = textures(3,sample(j));
    end 
    
else % more than one time step requires cell array
    
    for i = 1:blocks
        
        tmp_eulers = textures{i};
        
        for j = 1:n
            tmp_sample(1,j) = tmp_eulers(1,sample(j));
            tmp_sample(2,j) = tmp_eulers(2,sample(j));
            tmp_sample(3,j) = tmp_eulers(3,sample(j));
        end
        
        sample_texture(i) = {tmp_sample}; 
end


end

