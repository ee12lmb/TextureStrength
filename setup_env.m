%% Setup Environment
%
% Sets up relevent paths and sets a global variable that
% allows other functions to know whether the neccessary paths have
% already been added. Initilises MTEX package.
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

% add paths to function folders
addpath /nfs/see-fs-01_teaching/ee12lmb/project/source/dev/Jstrain
addpath /nfs/see-fs-01_teaching/ee12lmb/project/source/dev/M-index
addpath /nfs/see-fs-01_teaching/ee12lmb/project/source/dev/Ngrains
addpath /nfs/see-fs-01_teaching/ee12lmb/project/source/dev/readfiles

% set global variable for other functions to see
ENVIRONMENT_SET = 1;