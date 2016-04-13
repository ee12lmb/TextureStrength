function [ J, strain ] = j_index(input_texture,n,seed,varargin)
%J_INDEX calculates the J index for a specified number of grains
%
%   Takes input texture as either a VPSC or EBSD (*.ctf) file path, 
%   or a cell array/matrix of a Euler angles that has already been read in, 
%   randomly sampling n grains (without replacement). 
%
%   Inputs:  input_texture - file path of VPSC or EBSD (*.ctf) data
%                            OR matrix/cell array of euler angles
%            n             - number of grains to pull from texture
%            seed          - allows repeatability, i.e. generates the same
%                            'random' samples if set equal to previous run
%
%   Outputs: J             - the J-index as calculated by the continuous
%                            function method outlined in Mainprice (2014)
%            strain        - strain vector if input is file path (otherwise
%                            this is known before function call)
%
%   Optional arguments...
%
%   j = J_INDEX(input,n,seed,'crystal',CRYSTAL) will calculate the
%   continuous J-index using the theoretical random distribution for a 
%   specified crystal symmetry. CRYSTAL can be; 'olivine','quartz' or
%   'post-perovskite'. The default is 'olivine'.
%
%   j = J_INDEX(input,n,seed,'outfile',PATH) will output to the file
%   specified by PATH (should be a string). The output will include meta
%   data in headers.
%
%   Lewis Bailey - University of Leeds, School of Earth and Environment 
%   2015-16 Undergraduate final year project
%
%   References
%
%   Mainprice, D., Bachmann, F., Hielscher, R. and Schaeben, H. (2014).
%   "Descriptive tools for the analysis of texture projects with large
%   datasets using MTEX: strength, sym- metry and components". Geological
%   Society, London, Special Publications, 409. doi: 10.1144/SP409.8.
%
%   Usage: [ j, strain ] = j_index(input_texture,n,seed,...)
%
%   See also: M_INDEXCONT, M_INDEXDISC, INDEX_REPEAT

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

% determine input type and extract relevant information
[ textures, strain, blocks, input_texture, output ] = get_inputInfo(input_texture,n,seed,crystal);

    
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
            fprintf(fid,'+Input file:\t%s\n',input_texture);
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

