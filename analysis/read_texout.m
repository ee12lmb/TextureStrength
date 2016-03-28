function [ data, type, header ] = read_texout(file)
%READ_TEXOUT Reads in data from output of an index calculation
%   Data files outputted by j_index, m_indexCont, m_indexDisc and
%   index_repeat can be read in using read_texout. Header information is
%   also extracted from these files.
%
%   Inputs:     file   - file outputted from other texture function
%
%   Outputs:    data   - array of data read from file, either two columns
%                        where strain is available, or a single column for 
%                        all other cases
%
%               type   - label at the top of output file that represents the
%                        data within the file e.g. MD2 is an output from
%                        m_indexDisc with two columns (strain | index)
%
%               header - meta data extracted from data file
%                    
%   Usage: [ data, type, header ] = read_texout(file)

%% setup environment
addpath /nfs/see-fs-01_teaching/ee12lmb/project/source/dev/
setup_env

%% Read in file 

% open file for reading
fid = fopen(file,'r');

% determine type of file and length of data from first line
info = textscan(fid,'%s %d',1);

% unpack info (bit clunky)
type   = info{1};
type   = type{1};
length = info{2};


switch type
    case 'IR' % reading in index repeat file
        
        for i = 1:14 % capture headers
            header_tmp{i} = fgetl(fid);
        end

        % extract useful header info
        fname     = textscan(header_tmp{4},'%s %s %s');
        index     = textscan(header_tmp{5},'%s %s %s %s %s %s %s');
        xtal      = textscan(header_tmp{8},'%s %s');
        numgrains = textscan(header_tmp{9},'%s %d %s %s %s %s');
        repeat    = textscan(header_tmp{11},'%s %d');
        
        header{1} = fname{3};
        header{2} = xtal{2};
        header{3} = numgrains{2};
        header{4} = index{2};
        header{5} = repeat{2};
        
        % read in data
        data = fscanf(fid, '%g', [length 1]);
        
    case {'MC2','MD2','J2'} % reading in two column index output file
        
        for i = 1:10 % capture headers
            header_tmp{i} = fgetl(fid);
        end
        
        % extract useful header info
        fname     = textscan(header_tmp{4},'%s %s %s');
        xtal      = textscan(header_tmp{5},'%s %s');
        numgrains = textscan(header_tmp{6},'%s %d %s %s %s %s');
        
        header{1} = fname{3};
        header{2} = xtal{2};
        header{3} = numgrains{2};
        
        % read in data
        data = fscanf(fid, '%g %g', [2 length]);
        data = data';
        
    case {'MC1','MD1','J1'} % reading in single column index output file
        
        for i = 1:10
            header_tmp{i} = fgetl(fid);
        end
        
        % extract useful header info
        fname     = textscan(header_tmp{4},'%s %s %s');
        xtal      = textscan(header_tmp{5},'%s %s');
        numgrains = textscan(header_tmp{6},'%s %d %s %s %s %s');
        
        header{1} = fname{3};
        header{2} = xtal{2};
        header{3} = numgrains{2};
        
        % read in data
        data = fscanf(fid, '%g', [length 1]);
        
    otherwise 
        error('Input file type not recognised')
end


end

