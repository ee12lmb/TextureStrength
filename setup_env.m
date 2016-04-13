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

fprintf('\nStarting up texture strength toolbox...\n')
fprintf('-----------------------------------------------------------\n')
fprintf('Lewis Bailey 2015-16 Undergraduate final year project\n\n')
   

fprintf('Initialising MTEX package...\n')
% add path and initialise MTEX package
addpath ../mtex-4.1.3/
startup_mtex

fprintf('Adding path to MSAT package...\n')
% add path to MSAT
addpath ../msat-1.1.1/msat

fprintf('Adding path to toolbox functions...\n\n')
% add paths to function folders
addpath ./J-index
addpath ./M-index
addpath ./readfiles
addpath ./analysis

% list available functions
help toolbox_functions

% set global variable for other functions to see
ENVIRONMENT_SET = 1;
