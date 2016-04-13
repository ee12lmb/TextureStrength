function [ m, strain ] = m_indexCont(input_texture,n,seed,varargin)
%M_INDEXCONT returns  continuous M-index and strain vector for an input texture 
%
%   Takes either a VPSC file path or a cell array/matrix of a texture that
%   has already been read in (see read_VPSC). Input can also be an EBSD file
%   (*.ctf). 
%
%   Inputs:  input_texture - file path of VPSC or EBSD (*.ctf) data
%                            OR matrix/cell array of euler angles
%            n             - number of grains to pull from texture
%            seed          - allows repeatability, i.e. generates the same
%                            'random' samples if set equal to previous run
%
%   Outputs: m             - the M-index as calculated by the continuous
%                            function method outlined in Mainprice (REF)
%            strain        - strain vector if input is file path (otherwise
%                            this is known before function call)
%
%   Optional arguments...
%
%   m = M_INDEXCONT(input,n,seed,'crystal',CRYSTAL) will calculate the
%   continuous M-index using the theoretical random distribution for a 
%   specified crystal symmetry. CRYSTAL can be; 'olivine','quartz' or
%   'post-perovskite'. The default is 'olivine'.
%
%   m = M_INDEXCONT(input,n,seed,'outfile',PATH) will output to the file
%   specified by PATH (should be a string). The output will include meta
%   data in headers.
%
%   Lewis Bailey - University of Leeds, School of Earth and Environment 
%   2015-16 Undergraduate final year project
%
%   Usage: [ m, strain ] = m_indexCont(input_texture,n,seed,...)
%
%   See also: J_INDEX, M_INDEXDISC, INDEX_REPEAT

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
    switch varargin{iarg}
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


% get misorientation angle distribution for specified symmetry
[uniform_density,~] = calcAngleDistribution(CS);
uniform_density_N = uniform_density/sum(uniform_density); % normalise

%% Calculate M-index 

% check how many textures are given (i.e. do we need to reference cell
% array or just normal matrix)
if (blocks == 1)
  
    % call function to retrive relevant MDF info
    [~, MDFdensity,~] = continuousMDF(textures,CS,SS);
    MDF_density = MDFdensity/sum(MDFdensity); % normalise
    
    % calculate M-index
    m = (sum((abs(uniform_density_N - MDF_density))/2));
    
else % now deal with multiple timesteps (e.g. cell arrays)
    
    for i = 1:blocks
        % call function to retrive relevant MDF info
        [~, MDFdensity,~] = continuousMDF(textures{i},CS,SS);
        MDF_density = MDFdensity/sum(MDFdensity); % normalise
    
        % calculate M-index
        m(i) = (sum((abs(uniform_density_N - MDF_density))/2));

    end
    
end

%% Build output

time = toc;

if (wantout == 0) % if the filepath has been given as an option
    
    % assume that filepath checked in shell script/matlab can handle this
    fid = fopen(outfile,'a'); % open file for writing (append, so can add headers in shell)

    switch output
        case 0     % our input was a file path so we know strain
            
            fprintf(fid,'MC2\t%i\n',length(m)); % code for read_texout 
            fprintf(fid,'+Function:\tm_indexCont\n');
            fprintf(fid,'+Time/date:\t%i:%i %i/%i/%i\n',t(4),t(5),t(3),t(2),t(1));
            fprintf(fid,'+Input file:\t%s\n',input_texture);
            fprintf(fid,'+Crystal:\t%s\n',crystal);
            fprintf(fid,'+Grains:\t%i\n',n);
            fprintf(fid,'+Seed:\t\t%i\n',seed);
            fprintf(fid,'+Time taken(s):\t%f\n',time);
            fprintf(fid,'+Columns:\tStrain,M-index\n');
            fprintf(fid,'Data\n');

              for i = 1:length(m)
                  fprintf(fid,'%10.5f %10.5f\n',strain(i),m(i));
              end
                
        case 1     % our input was matrix/EBSD so we don't know strain
            
            % build header
            fprintf(fid,'MC1\t%i\n',length(m)); % code for read_texout 
            fprintf(fid,'+Function:\tm_indexCont\n');
            fprintf(fid,'+Time/date:\t%i:%i %i/%i/%i\n',t(4),t(5),t(3),t(2),t(1));
            fprintf(fid,'+Input file:\t%s\n',input_texture);
            fprintf(fid,'+Crystal:\t%s\n',crystal);
            fprintf(fid,'+Grains:\t%i\n',n);
            fprintf(fid,'+Seed:\t\t%i\n',seed);
            fprintf(fid,'+Time taken(s):\t%f\n',time);
            fprintf(fid,'+Columns:\tM-index\n');
            fprintf(fid,'Data\n');

              for i = 1:length(m)
                  fprintf(fid,'%10.5f\n',m(i));
              end
    end
    fclose(fid);
    
    
end
end

