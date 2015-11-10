function [ m, strain ] = m_indexCont(input_texture,n,seed)
%M_INDEX returns m index and strain vector for an input texture 
%   Takes either a VPSC file path or a cell array/matrix of a texture that
%   has already been read in (see read_VPSC).
%
%   Inputs:  input_texture - file path/texture array/texture matrix 
%            n             - number of samples to pull from texture
%            seed          - allows repeatability, i.e. generates the same
%                            'random' samples if set equal to previous run
%
%   Outputs: m             - the M-index as calculated by the continuous
%                            function method outlined in Mainprice (REF)
%            strain        - strain vector if input is file path (otherwise
%                            this is known before function call)
%
%   Usage: [ m, strain ] = m_indexCont(input_texture,n,seed)

%% Setup & read data

% set up MTEX package & source all functions
addpath /nfs/see-fs-01_teaching/ee12lmb/project/source
startup_mtex;

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

% get misorientation angle distribution for specified symmetry
[uniform_density,uniform_angles] = calcAngleDistribution(CS);
uniform_density_N = uniform_density/sum(uniform_density); % normalise

%% Calculate M-index 

% check how many textures are given (i.e. do we need to reference cell
% array or just normal matrix)
if (blocks == 1)
  
    % call function to retrive relevant MDF info
    [~, MDFdensity,~] = textureMDF(textures);
    MDF_density = MDFdensity/sum(MDFdensity); % normalise
    
    % calculate M-index
    m = (sum((abs(uniform_density_N - MDF_density))/2));
    
else % now deal with multiple timesteps (e.g. cell arrays)
    
    for i = 1:blocks
        
        % call function to retrive relevant MDF info
        [~, MDFdensity,~] = textureMDF(textures{i});
        MDF_density = MDFdensity/sum(MDFdensity); % normalise
    
        % calculate M-index
        m(i) = (sum((abs(uniform_density_N - MDF_density))/2));

    end
    
end

end

