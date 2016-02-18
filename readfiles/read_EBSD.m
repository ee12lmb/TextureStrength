function [ eulers, CS, phase_names ] = read_EBSD(filename,crystal)
%READ_EBSD extracts Euler angles for specific crystals from EBSD
%   Given the path of a .ctf EBSD file, read_EBSD will extract the Euler
%   angles for the phase defined by 'crystal'. It returns them in degrees.
%
%   Input:   filename    - file path to .ctf file
%            crystal     - which phase to extract Euler angles for
%
%   Outputs: eulers      - Nx3 matrix containing Euler angles for each pixel 
%                          in the form (phi_1,Phi,phi_2)
%            CS          - Cell array of crystal symmetries from .ctf
%            phase_names - list of phases present in the file
%
%   Usage: [ eulers, CS, phase_names ] = read_EBSD(filename,crystal)


addpath /nfs/see-fs-01_teaching/ee12lmb/project/source/dev/
setup_env;

%% Check input and extract/check headers

% check input is .ctf
[~,~,ext] = fileparts(filename);
assert((strcmp(ext,'.ctf') == 1),'Must be *.ctf file.');

% pull out header info from *.ctf file (REF AMIECIA)
[CS, phase_names] = get_symmetry(filename);

% check that the requested crystal is contained in the file
check = any(ismember(lower(phase_names),lower(crystal)));
assert((check == 1),'Requested crystal not indexed in this EBSD file.');


%% Load and extract Euler angles for requested crystal

% load the EBSD object
ebsd = loadEBSD(filename,CS,'interface','ctf','convertEuler2SpatialReferenceFrame');

% extract individual orientations
o = ebsd(crystal).orientations;

% extract Euler angles from individual orientation objects
eulers = Euler(o);
eulers = transpose(eulers)/degree;     % transpose for consistency with read_VPSC

end

