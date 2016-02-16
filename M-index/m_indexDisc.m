function [ m,strain ] = m_indexDisc(input_texture,n,seed,varargin)
%M_INDEXDISC calculates the M-index using a discrete method
%   Takes either a VPSC file path or a cell array/matrix of a texture that
%   has already been read in (see read_VPSC).
%
%   Inputs:  input_texture - file path/texture array/texture matrix 
%            CS            - crystal symmetry
%            n             - number of samples to pull from texture
%            seed          - allows repeatability, i.e. generates the same
%                            'random' samples if set equal to previous run
%
%   Outputs: m             - the M-index as calculated by the discrete
%                            method as outlined by Skemer (REF)
%            strain        - strain vector if input is file path (otherwise
%                            this is known before function call)
%
%   Usage: [ m, strain ] = m_indexDisc(input_texture,n,seed)

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


% check if input is raw VPSC or texture array
if (ischar(input_texture) == 1)
    
    % if a file path is given then read in
    [textures,ngrains,strain,blocks] = sample_VPSC(input_texture,n,seed);
    output = 0; % output format (see bottom of function)
    
else
    
   % if input is not a file then pass to sample_texture to deal with
   [textures,blocks,~] = sample_texture(input_texture,n,seed);
   
   % strain information cannot be extracted from inputted texture
   %+but should already be known from previous read_VPSC
   strain = 'Input is texture - strain already extracted';
   output = 1;
 
end


%% Calculate and bin theoretical random dist

% find theoretic distribution density
[ uniform_density, uniform_angles ] = calcAngleDistribution(CS);

% find max theoretical angle
theta_max = uniform_angles(length(uniform_angles));

% turn to degrees
theta_max = theta_max/degree;
uniform_angles = uniform_angles/degree;

% calculate bin dimensions (here set to one degree - *could take input?*)
bins = linspace(0,theta_max,theta_max+1);

% calcAngleDistribution returns angles and densities as a fraction of 300
% (not sure why, seems arbitrary). To solve, need to sum up density in the
% bins that we have defined

% find number of angles returned from calcAngleDist in our defined bins
Nangles_in_bins = histc(uniform_angles,bins);

% intialise freq sum for all bins
uniform_freq = zeros(1,length(bins));  


j = 1; % initialise loop counter
for i = 1:length(bins) % now loop through and sum relevant densities into each bin

    for null = 1:Nangles_in_bins(i) % sum up correct indicies for this bin
        
        % e.g. if three angles lie in this bin, loop three times and sum
        uniform_freq(i) = uniform_freq(i) + uniform_density(j);

        % increase counter
        j = j +1;

    end 
    
    % have increased height for this bin, now deal with increased bin width
    uniform_freq(i) = uniform_freq(i)/Nangles_in_bins(i);
    
end % end bin loop

% normalise freq
%uniform_freq = uniform_freq/300
uniform_freq = uniform_freq/sum(uniform_freq);
uniform_density = uniform_density/sum(uniform_density);

% figure(1)
% bar(bins,uniform_freq,'histc')
% axis([0 140 0 0.025])



%% Calculate and bin misorientation angles for each input texture

if (blocks == 1)
    
    % find misorientation angle distribution
    [ disorentation, ~ ] = discreteMDF(textures,CS);

    % turn angle to degrees
    disorentation = disorentation/degree;

    % bin angles
    disor_freq = histc(disorentation,bins);

    % normalise 
    disor_freq = disor_freq/sum(disor_freq);

    %% Calculate M-index
    m = (theta_max/length(bins))*sum((abs(uniform_freq - disor_freq))/2);
    
else

    for i = 1:blocks 

        % find misorientation angle distribution
        [ disorentation, ~ ] = discreteMDF(textures{i},CS);

        % turn angle to degrees
        disorentation = disorentation/degree;

        % bin angles
        disor_freq = histc(disorentation,bins);

        % normalise 
        disor_freq = disor_freq/sum(disor_freq);
%         figure(i+1)
%         bar(bins,disor_freq,'histc')
%         axis([0 140 0 0.025]) 

        %% Calculate M-index

        m(i) = (theta_max/length(bins))*sum((abs(uniform_freq - disor_freq))/2);

    end

end


%% Build output to file (if requested *NEEDS WORK*)

time = toc;

if (wantout == 0) % if the filepath has been given as an option
    
    % assume that filepath checked in shell script/matlab can handle this
    fid = fopen(outfile,'w'); % open file for writing - will overwrite

    switch output
        case 0     % our input was a file path so we know strain
            
            fprintf(fid,['MD2\t%i\n',...  % build header in one fprintf call
                        '+Function:\tm_indexDisc\n',...
                        '+Time/date:\t%i:%i %i/%i/%i\n',...
                        '+Input file:\t%s\n',...
                        '+Crystal:\t%s\n'...
                        '+Grains:\t%i\n',...
                        '+Seed:\t\t%i\n',...
                        '+Time taken(s):\t%f\n',...
                        '+Columns:\tStrain,M-index\n',...
                        'Data\n'],...
                        length(m),t(4),t(5),t(3),t(2),t(1),input_texture,crystal,n,seed,time);
            
            fprintf(fid,'%10.5f %10.5f\n',[strain;m]); % dump data to file
       
                
        case 1     % our input was inputted texture so we don't know strain
            
            % build header
            fprintf(fid,['MD1\t%i\n',...
                         '+Function:\tm_indexDisc\n',...
                         '+Time/date:\t%i:%i %i/%i/%i\n',...
                         '+Input file:\tn/a\n',...
                         '+Crystal:\t%s\n'...
                         '+Grains:\t%i\n',...
                         '+Seed:\t\t%i\n','+Time taken(s):\t%f\n',...
                         '+Columns:\tM-index\n','Data\n'],...
                         length(m),t(4),t(5),t(3),t(2),t(1),crystal,n,seed,time);  

            fprintf(fid,'%10.5f\n',m);
 
    end
    fclose(fid);
    
    
end
   
end

