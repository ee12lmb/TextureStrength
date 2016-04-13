function [ CS ] = lookupSym(crystal)
%LOOKUPSYM finds symmetry for a given crystal
%
%   Add more definitions here to extend the functionality of the toolbox
%
%   Lewis Bailey - University of Leeds, School of Earth and Environment 
%   2015-16 Undergraduate final year project
%
%   See also: GET_INPUTINFO

addpath /nfs/see-fs-01_teaching/ee12lmb/project/source/dev/
setup_env;

%% Find symmetry

switch lower(crystal)
    case 'olivine'
        
        CS = crystalSymmetry('Pbnm', [4.75, 10.20, 5.98]);
        
    case 'post-perovskite'
        
        CS = crystalSymmetry('Cmcm', [2.47 8.09 6.11]);
        
    case 'quartz'
        
%         CS = crystalSymmetry('pointId', 18, [4.913 4.913 5.504],...
%             [90,90,120]*degree,'X||a*', 'Y||b', 'Z||c', 'mineral',...
%             'Quartz','color', 'red');

        % histograms for this agree with Wheeler (2001)
        CS = loadCIF('quartz');
        
        
    otherwise
        
        error('Input crystal not indexed')
        
end


end

