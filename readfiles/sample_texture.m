function [sample_texture,blocks,ngrains] = sample_texture(textures,n,seed)
%SAMPLE_TEXTURE randomly samples from a texutre already loaded into matlab
%
%   SAMPLE_TEXTURE allows resampling without the need to read a texture in
%   again. To repeat a sample, specify the same seed as previously ran.
%   Samples are taken without replacement.
%                                         
%   Inputs:  texture        - array or cell array of euler angles read in by
%                             read_VPSC or read_EBSD
%            n              - number of samples to pull out (cannot be bigger than
%                             the total in the sample).
%            seed           - sets the random number generator to a certain state
%                             to allow for repeatability. 
%
%   Outputs: sample_texture - set of randomly samples Euler angles of
%                             length n
%            blocks         - number of strain steps 
%            ngrains        - number of total grains in sample
%
%   Lewis Bailey - University of Leeds, School of Earth and Environment 
%   2015-16 Undergraduate final year project
%
%   Usage: [sample_texture,blocks,ngrains] = sample_texture(textures,n,seed)
%
%   See also: READ_VPSC, READ_EBSD, SAMPLE_VPSC, SAMPLE_EBSD, GET_INPUTINFO

%% Input checks

% check input texture format
if (iscell(textures) == 1) 
    
    % if texture is cell array then extract number of timesteps (blocks)
    blocks = length(textures);
    
    % extract number of xtls
    for i = 1:blocks
        nxtl(i) = length(textures{i});
    end
    
    % check that the number of grains does not change between time steps
    assert((range(nxtl) == 0), ...
        'Number of grains not consitent across time steps!')
    ngrains = nxtl(1); % assign scalar value for ngrains
    
else
    % if texture is just a matrix then there is only one timestep
    blocks = 1;
    ngrains = length(textures);
end

% check that n is not greater than number of grains in sample
assert((n <= ngrains), ...
    'Number of samples requested is greater than total number in sample!')

%% Sampling 

% generate sample index
rng(seed);                      % set seed to repeat analysis
sample = randperm(ngrains,n);   % sample without repeating


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

