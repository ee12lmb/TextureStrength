%% Setup Environment
%
% Sets up relevent paths and sets a global variable that
% allows other functions to know whether the neccessary paths have
% already been added. Initilises MTEX package and adds path to MSAT.
%
% Script will not setup up the environment if the global variable has
% already been set.

%% Add paths and initialise MTEX

global ENVIRONMENT_SET;
if (ENVIRONMENT_SET == 1) % environment is already set up
    return
end

% add path and initialise MTEX package
addpath /nfs/see-fs-01_teaching/ee12lmb/project/source/mtex-4.1.3/
startup_mtex

% add path to MSAT
addpath /nfs/see-fs-01_teaching/ee12lmb/project/source/msat-1.1.1/msat

% add paths to function folders
addpath /nfs/see-fs-01_teaching/ee12lmb/project/source/dev/J-index
addpath /nfs/see-fs-01_teaching/ee12lmb/project/source/dev/M-index
addpath /nfs/see-fs-01_teaching/ee12lmb/project/source/dev/readfiles
addpath /nfs/see-fs-01_teaching/ee12lmb/project/source/dev/analysis

% set global variable for other functions to see
ENVIRONMENT_SET = 1;
