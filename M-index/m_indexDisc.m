function [ m,strain ] = m_indexDisc(input_texture,n,seed)
%M_INDEXDISC calculates the M-index using a discrete method
%   Takes either a VPSC file path or a cell array/matrix of a texture that
%   has already been read in (see read_VPSC).
%
%   Inputs:  input_texture - file path/texture array/texture matrix 
%            n             - number of samples to pull from texture
%            seed          - allows repeatability, i.e. generates the same
%                            'random' samples if set equal to previous run
%
%   Outputs: m             - the M-index as calculated by the discrete
%                            method as outlined by Skemer (REF)
%            strain        - strain vector if input is file path (otherwise
%                            this is known before function call)
%
%   Usage: [ m, strain ] = m_indexDisc(input_texture,n,seed)



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
   [textures,blocks,~] = sample_texture(input_texture,n,seed);
   
   % strain information cannot be extracted from inputted texture
   %+but should already be known from previous read_VPSC
   strain = 'Input is texture - strain already extracted'; 
    
end

% Set up symmetry
CS = crystalSymmetry('Pbnm', [4.75, 10.20, 5.98]);
SS = specimenSymmetry('-1');




%% Calculate misorientation angles

% ** THIS WOULD BE USEFUL TO MODULARISE **
% e.g function discreteMDF

    
end
























end

