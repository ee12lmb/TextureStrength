function [ strain ] = get_strain(index_input,varargin)
%GET_STRAIN will find expected strain from a given index value
%
%   Index values calculated from EBSD data can be related to strain by
%   interpolating reference curves (indices vs strain) computed for VPSC
%   models (where strain is known). The crystal type and index used are
%   found from the header information in infile. Strain can also be
%   determined for a single index value (e.g. if passed a value within a
%   script) - in this case, the crystal type and index must be specified. 
%
%   Inputs:     index_input - a data file outputted by an index calculation,
%                             that can be read by read_texout 
%                             OR
%                             a variable containing a single index value                             
%
%   Outputs:    strain      - strain value found from interpolating 
%                             reference curves
%
%   s = GET_STRAIN(filepath,'simple-shear') will get the strain value by 
%   interpolating the simple shear VPSC curve.
%
%   Similarly, s = GET_STRAIN(filepath,'axial-compression') will do the same 
%   for the axial compression curve.
%
%   If the call is s = GET_STRAIN(index_value,...), i.e. not a file path
%   but a single index value, optional args are required to specify the
%   index and crystal symmetry by the following flags...
%
%   'index'  : 'j' OR 'mc' OR 'md'
%
%   'crystal': 'olivine' OR 'quartz' OR 'post-perovskite'
%
%   e.g. s = GET_STRAIN(m(1),'crystal','quartz','index','j')
%
%   Lewis Bailey - University of Leeds, School of Earth and Environment 
%   2015-16 Undergraduate final year project
%
%   Usage: [ strain ] = GET_STRAIN(infile,varargin)


%% Setup and read data

addpath /nfs/see-fs-01_teaching/ee12lmb/project/source/dev/
setup_env;

% set default values for optional arguments
wantPlot   = 0;              % by deafult dont plot
strainType = 'simple-shear'; % by default use simple shear curve
crystal    = 'quartz';
indexType  = 'j';

iarg = 1;
while iarg<=(length(varargin))
    switch varargin{iarg}
        
        case 'simple-shear'      % use simple shear reference curve
            strainType = 'simple-shear';
        case 'axial-compression' % use axial compression reference curve  
            strainType = 'axial-compression';  
            
        case 'crystal'
            iarg = iarg + 1;
            crystal = varargin{iarg};
            
        case 'index'
            iarg = iarg + 1;
            indexType = varargin{iarg};
            
        case 'plot'
            wantPlot = 1;
        otherwise
            error('Unknown flag')
    end
    iarg = iarg + 1;
end

if (ischar(index_input) == 1) % input is file path
    % read in output file
    [index, type, header] = read_texout(index_input);

    % check that the index was calculated for EBSD data (*.ctf)
    [~,~,ext] = fileparts(cell2mat(header{1}));
    assert((strcmp(ext,'.ctf') == 1),'Index was not calculated from EBSD (*.ctf) file')

    % determine crystal symmetry & index type
    crystal = cell2mat(header{2});

    switch type

        case 'IR'
            fprintf('Input is index repeat - taking the average index value')
            index = mean(index);        
            indexType = cell2mat(header{4});
        case 'MD1'
            indexType = 'md';  
        case 'MC1'
            indexType = 'mc';
        case 'J1'
            indexType = 'j';

        otherwise
            error('Input file not suitable')

    end
    
else % input must be matrix
    
    % all other info should have been specified in opt args
    assert((ismatrix(index_input) == 1),'Input not recognised')
    assert((length(index_input) == 1),'Input should be single index value')
    index = index_input;
    
end
    

%% Read references curves and interpolate

% add path to reference curves
addpath('~/project/source/dev/analysis/reference_curves/')

% load curve & separate data
fname = sprintf('%s_%s_%s.out',crystal,indexType,strainType);
reference = read_texout(fname);
ref_strain = reference(:,1);
ref_index  = reference(:,2);

% check that the index value is in range of reference curve
if ((index > max(ref_index)) || (index < min(ref_index)))
    warning('Index value %s = %f out of range of %s reference curve',...
                                             indexType,index,strainType)
    strain = NaN;
    return % cannot interp as index is out of range
end

% use linear interpolation to find strain value
%strain = interp1(ref_index,ref_strain,index)

for i = 1:length(reference)-1
   
    if (ref_index(i+1) > index)  % we have found upper bound
        
        % find gradient
        m = (ref_strain(i+1) - ref_strain(i))/(ref_index(i+1) - ref_index(i));
        
        % find x (index) incriment from lower bound 
        dind = index - ref_index(i);
        
        % add on apporpriate amount of strain (dy = m*dx)
        strain = ref_strain(i) + m*dind;
        
        break
        
    end
    
end

if (wantPlot == 1)
    figure('Name','Index-strain interp')
    hold on
    plot(ref_strain,ref_index,'k-');
    plot(strain,index,'rx');
end



end

