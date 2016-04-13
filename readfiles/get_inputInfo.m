function [ textures, strain, blocks, input, output ] = get_inputInfo(input,n,seed,crystal)
%GET_INPUTINFO Extracts Euler angles and other info from various texture
%inputs
%
%   GET_INPUTINFO will deal with either a VPSC file path, an EBSD (*.ctf)
%   file path or an previously inputted texture (either cell array or
%   single matrix). In all cases 
%
%   Input:      input    - either VPSC, EBSD or matrix/cell array
%               n        - number of grains to sample
%               seed     - sets the random numbers for repeatability
%               crystal  - crystal type (only relevant for EBSD)
%
%   Output:     textures - a matrix or cell array of Euler angles
%               strain   - either strain vector (VPSC) or string (EBSD/matrix)
%               blocks   - number of strain steps present
%               input    - returns the input string, or 'Matlab matrix'
%               output   - flag for output type in other functions
%
%   Lewis Bailey - University of Leeds, School of Earth and Environment 
%   2015-16 Undergraduate final year project
%
%   Usage: [ textures, strain, blocks, input, output ] = get_inputInfo(input,n,seed,crystal)
%
%   See also: SAMPLE_VPSC, SAMPLE_EBSD, SAMPLE_TEXTURE

%%

% check if input is raw VPSC, EBSD or texture array
if (ischar(input) == 1)
    
    % check file extension to see if EBSD
    [~,~,ext] = fileparts(input);
    
    if (strcmp(ext,'.ctf') == 1)
        
        % if .ctf then call sample_EBSD to read from file
        [ textures, ~ ] = sample_EBSD(input,crystal,n,seed);
        blocks = 1; % only ever one block in .ctf file
        strain = 'Input is EBSD, cannot return strain';
        output = 1;  % cannot output strain

        
    else
        
        % if a file path is given and not EBSD then must be VPSC
        [textures,~,strain,blocks] = sample_VPSC(input,n,seed);
        output = 0; % specify output format for file
    
    end
    
else
    
   % if input is not a file then pass to sample_texture to deal with
   [textures,blocks] = sample_texture(input,n,seed);
   
   % strain information cannot be extracted from inputted texture
   %+but should already be known from previous read_VPSC
   strain = 'Input is texture - strain already extracted'; 
   input = 'Matlab matrix';
   output = 1;
    
end
end

