function [ MDF, density, angles] = continuousMDF(input_texture,CS,SS)
%textureMDF returns the MDF, the angles, and angle desnsities for an input
%texture
%   If input texture has more than one timestep (i.e. requires cell array),
%   then the outputs will also be cell arrays.
%
%   Inputs:  input_texture - either VPSC file path, texture cell array or
%                            single texture.
%            CS            - crystal symmetry
%            SS            - sample symmetry
% 
%   Outputs: MDF           - The misorientation distribution function
%            density       - The angle density (for angles in ANGLES)
%            angles        - Array of angles (up to max for cystal symmetry)
%
%   Usage: [ MDF, density, angles ] = textureMDF(input_texture,CS,SS)
%

addpath /nfs/see-fs-01_teaching/ee12lmb/project/source/dev/
setup_env


%% Input checks

% check if input is raw VPSC, texture array or single matrix
if (ischar(input_texture) == 1)
    
    % if a file path is given then read in
    [textures,ngrains,strain,blocks] = read_VPSC(input_texture);
    
elseif (iscell(input_texture) == 1)
    
    % if texture is cell array then extract number of timesteps (blocks)
    blocks = length(input_texture);
        
    % assign to textures for consistancy
    textures = input_texture;
    
else % texture now must be single matrix
    
    % assign to textures for consistancy
    textures = input_texture;
    
    % only one block of data
    blocks = 1;
    
end

%% Calculate outputs

if (blocks == 1)
    
    eulers_r = textures*degree;

    % find ODF
    g = orientation('Euler', eulers_r(1,:), eulers_r(2,:), eulers_r(3,:), ...
        CS, SS, 'Bunge');
    odf = calcODF(g,'HALFWIDTH', 10*degree, 'silent');
    
    % calculate MDF
    MDF = calcMDF(odf);
    
    % find density and normalise
    [density,angles] = calcAngleDistribution(MDF,'resolution',1*degree);
    
else % deal with multiple textures
    
    for i = 1:blocks
        
        eulers_r = textures{i}*degree;

        % find ODF
        g = orientation('Euler', eulers_r(1,:), eulers_r(2,:), eulers_r(3,:), ...
            CS, SS, 'Bunge');
        odf = calcODF(g,'HALFWIDTH', 10*degree, 'silent');

        % calculate MDF
        MDF{i} = calcMDF(odf);

        % find density and normalise
        [density{i},angles{i}] = calcAngleDistribution(MDF{i},'resolution',1*degree);  
        
    end
end

