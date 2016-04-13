function [ sampled_texture, ngrains ] = sample_EBSD(infile,crystal,n,seed)
%SAMPLE_EBSD Returns a random sample of Euler angles from EBSD file
%
%   Uses read_EBSD to get Euler angles from an EBSD (.ctf) file for a
%   specific crystal phase defined by 'crystal'. This martix of Euler
%   angles is then passed to sample_texture to return a matrix of n
%   randomly sampled sets of Euler angles. Set the seed to allow repitition
%   
%   Inputs:     infile          - input EBSD file (must be .ctf)
%               crystal         - string containing crystal phase name
%               n               - number of grains to sample
%               seed            - sets the random numbers to allow repeatability
%
%   Outputs:    sampled_texture - the sampled texture
%               ngrains         - number of grains in the whole sample
%
%   Lewis Bailey - University of Leeds, School of Earth and Environment 
%   2015-16 Undergraduate final year project
%
%   Usage: [ sampled_texture, ngrains ] = sample_EBSD(infile,crystal,n,seed)
%
%   See also: READ_EBSD, SAMPLE_TEXTURE


%% Setup envrionment
addpath /nfs/see-fs-01_teaching/ee12lmb/project/source/dev/
setup_env;

%% Read in and sample file

% read in entire file
eulers = read_EBSD(infile,crystal);

% pass to sample texture to deal with 
[sampled_texture, ~, ngrains] = sample_texture(eulers,n,seed);


end

