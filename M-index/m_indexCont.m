function [ m, strain, MDF_density, MDF_angles, uniform_density, uniform_angles ] ...
         = m_index(input_texture,n,seed)
%M_INDEX returns m index for an input texture file (Bunge format).
%   Detailed explanation goes here

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

% get misorientation angle distribution for specified symmetry
[uniform_density,uniform_angles] = calcAngleDistribution(CS);
uniform_density_N = uniform_density/sum(uniform_density); % normalise

%% Calculate M-index 

% check how many textures are given (i.e. do we need to reference cell
% array or just normal matrix)
if (blocks == 1)
  
    eulers_r = textures*degree; % textures is SAMPLED TEXTURE

    % calculate M-index for this block and store
    g = orientation('Euler', eulers_r(1,:), eulers_r(2,:), eulers_r(3,:), ...
        CS, SS, 'Bunge');
    odf = calcODF(g,'HALFWIDTH', 10*degree, 'silent');
    
    % calculate MDF
    MDF = calcMDF(odf);
    
    % find density and normalise
    [MDF_density,MDF_angles] = calcAngleDistribution(MDF,'resolution',1*degree);
    MDF_density_N = MDF_density/sum(MDF_density);
    
    % calculate M-index
    m = (sum((abs(uniform_density_N - MDF_density_N))/2));
    
else % now deal with multiple timesteps (e.g. cell arrays)
    
    % pull out samples for each time step
    for i = 1:blocks % deal with repeated texture

        eulers_r = textures{i}*degree; % textures is SAMPLED TEXTURE

        % calculate M-index for this block and store
        g = orientation('Euler', eulers_r(1,:), eulers_r(2,:), eulers_r(3,:), ...
            CS, SS, 'Bunge');
        odf = calcODF(g,'HALFWIDTH', 10*degree, 'silent');
       
        % calculate MDF
        MDF = calcMDF(odf);
    
        % find density and normalise
        [MDF_density{i},MDF_angles{i}] = calcAngleDistribution(MDF,'resolution',1*degree);
        MDF_density_N = MDF_density{i}/sum(MDF_density{i});
    
        % calculate M-index
        m(i) = (sum((abs(uniform_density_N - MDF_density_N))/2));

    end
    
end

end

