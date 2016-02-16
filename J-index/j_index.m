function [ J, strain ] = j_index(input_texture,n,seed,varargin)
%J_INDEX calculates the J index for a specified number of grains
%   A number of grains (n) are selected at random from either an input file
%   or a texture array already inputted into matlab (see read_VPSC). The
%   J-index is then calculated for each time step.
%
%   Inputs:  input_texture - file path/texture array/texture matrix
%                        n - number of grains to use
%                     seed - allows repeatability of 'random' numbers
%
%            OPTIONAL: 'filename','[ name out output file ]'
%                          - specify location to dump data to file
%
%   Outputs:             J - either single value/vector of J-index values
%                   strain - if read from file, will be vector of strains
%
%   Usage: [ J, strain ] = j_index(input_texture,n,seed,varargin)

tic;
t = clock;
%% Setup & read data

addpath /nfs/see-fs-01_teaching/ee12lmb/project/source/dev/
setup_env;

% check for optional arguments
iarg = 1;
wantout = 1; % we don't want output unless the 'filename' flag is active

% setup defautl symmetry (olivine)
CS = crystalSymmetry('Pbnm', [4.75, 10.20, 5.98]);
SS = specimenSymmetry('-1');
crystal = 'olivine';

while iarg<(length(varargin))
    switch lower(varargin{iarg})
        case 'outfile'
            
            iarg = iarg + 1; % take next argument as filename 
            outfile = varargin{iarg};
           
            % check that we are not overwriting a file
            check = exist(outfile,'file');
            assert((check == 0),'Output file already exists!')
           
            wantout = 0;  % we do want the output to file
            
        case 'crystal'  % find the appropriate symmetry 
            
            iarg = iarg + 1; % take next argument
            crystal = varargin{iarg};
            CS = lookupSym(crystal);
            
        otherwise
            error('Unknown flag')
    end
    iarg = iarg + 1;
end


% check if input is raw VPSC or texture array
if (ischar(input_texture) == 1)
    
    % if a file path is given then read in
    [textures,~,strain,blocks] = sample_VPSC(input_texture,n,seed);
    output = 0; % specifies the ouput format (as we DO have strain info)
    
else
    
   % if input is not a file then pass to sample_texture to deal with
   [textures,blocks] = sample_texture(input_texture,n,seed);
   
   % strain information cannot be extracted from inputted texture
   %+but should already be known from previous read_VPSC
   strain = 'Input is texture - strain already extracted'; 
   output = 1; % specifies the output format (we DO NOT have strain info)
    
end

% Set up symmetry


    
%% Calculate J-index 

% check how many textures are given (i.e. do we need to reference cell
% array or just normal matrix)
if (blocks == 1)
  
    eulers_r = textures*degree; % textures is SAMPLED TEXTURE

    % calculate J index for this block and store
    g = orientation('Euler', eulers_r(1,:), eulers_r(2,:), eulers_r(3,:), ...
        CS, SS, 'Bunge');
    odf = calcODF(g,'HALFWIDTH', 10*degree, 'silent');
    J = textureindex(odf);  
    
    
else % now deal with multiple timesteps (e.g. cell arrays)
    
    % pull out samples for each time step
    for i = 1:blocks % deal with repeated texture

        eulers_r = textures{i}*degree; % textures is SAMPLED TEXTURE

        % calculate J index for this block and store
        g = orientation('Euler', eulers_r(1,:), eulers_r(2,:), eulers_r(3,:), ...
            CS, SS, 'Bunge');
        odf = calcODF(g,'HALFWIDTH', 10*degree, 'silent');
        J(i) = textureindex(odf);

    end
    
end

%% Build output

time = toc;

if (wantout == 0) % if the filepath has been given as an option
    
    % assume that filepath checked in shell script/matlab can handle this
    fid = fopen(outfile,'a'); % open file for writing (append, so can add headers in shell)

    switch output
        case 0     % our input was a file path so we know strain
            
            % build header            
            fprintf(fid,'J2\t%i\n',length(J)); % code for read_texout 
            fprintf(fid,'+Function:\tj_index\n');
            fprintf(fid,'+Time/date:\t%i:%i %i/%i/%i\n',t(4),t(5),t(3),t(2),t(1));
            fprintf(fid,'+Input file:\t%s\n',input_texture);
            fprintf(fid,'+Crystal:\t%s\n',crystal);
            fprintf(fid,'+Grains:\t%i\n',n);
            fprintf(fid,'+Seed:\t\t%i\n',seed);
            fprintf(fid,'+Time taken(s):\t%f\n',time);
            fprintf(fid,'+Columns:\tStrain,J-index\n');
            fprintf(fid,'Data\n');

              for i = 1:length(J)
                  fprintf(fid,'%10.5f %10.5f\n',strain(i),J(i));
              end
                
        case 1     % our input was inputted texture so we don't know strain
            
            % build header
            fprintf(fid,'J1\t%i\n',length(J)); % code for read_texout 
            fprintf(fid,'+Function:\tj_index\n');
            fprintf(fid,'+Time/date:\t%i:%i %i/%i/%i\n',t(4),t(5),t(3),t(2),t(1));
            fprintf(fid,'+Input file:\tn/a\n');
            fprintf(fid,'+Crystal:\t%s\n',crystal);
            fprintf(fid,'+Grains:\t%i\n',n);
            fprintf(fid,'+Seed\t\t%i\n',seed);
            fprintf(fid,'+Time taken(s):\t%f\n',time);
            fprintf(fid,'+Columns:\tJ-index\n');
            fprintf(fid,'Data\n');

              for i = 1:length(J)
                  fprintf(fid,'%10.5f\n',J(i));
              end
    end
    fclose(fid);
    
    
end

end

