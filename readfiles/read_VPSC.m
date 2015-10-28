% MVT_read_VPSC_file - Read euler angles from a VPSC texture file 
%
% Given the name of a VPSC formatted texture file, return an array 
% of Euler angles and the number of crystals in the file. If only one 
% texture exists in the file return Matlab arrays for eulers and nxtl.
% If multiple textures exists return cell arrays (where the indicies 
% map to each texture.
% 
% Usage: 
%     [eulers, nxtl] = MVT_read_VPSC_file(filename)
%
% See also: MVT_write_VPSC_file


% Copyright (c) 2012, Andrew Walker
% All rights reserved.
% 
% Redistribution and use in source and binary forms, 
% with or without modification, are permitted provided 
% that the following conditions are met:
% 
%    * Redistributions of source code must retain the 
%      above copyright notice, this list of conditions 
%      and the following disclaimer.
%    * Redistributions in binary form must reproduce 
%      the above copyright notice, this list of conditions 
%      and the following disclaimer in the documentation 
%      and/or other materials provided with the distribution.
%    * Neither the name of the University of Bristol nor the names 
%      of its contributors may be used to endorse or promote 
%      products derived from this software without specific 
%      prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS 
% AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED 
% WARRANvTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A 
% PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL 
% THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY 
% DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF 
% USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE 
% OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

function [eulers,nxtl,strain,blocks] = readVPSC(filename)
    % Reads data from a VPSC output or input texture file.
    % Must be Bunge convention and in degrees. Returns eulers, a (3,nxtl) 
    % array of Euler angles (Bunge convention, (1,:) = phi1, (2,:) = Phi
    % and (3,:) = phi2) and a scalar nxtl, specifying the number of 
    % crystals / Euler angle triples.

    % Read data from the file
    fid = fopen(filename); % Read - the default
    
    blocks = 0; % Which number texture block are we on?
    
    while (true) % infinite loop (break with conditional)
        
        % get strain header info
        tmp_l = fgetl(fid);
        if (tmp_l == -1)     % if no more lines then break loop 
            break 
        end
        
        tmp_header = sscanf(tmp_l, '%s %s %s %s %g');
        
        fgetl(fid); % Lengths of phase ellipsoid axes - ignore
        fgetl(fid); % Euler angles for phase ellipsoid - ignore
        L = sscanf(fgetl(fid), '%s %d'); % Convention and number of crystals

        % Get hold of header info
        assert((char(L(1))=='B'), ... % Check Euler angle convention
            'Could not read VPSC file - not Bunge format\n');
        tmp_nxtl = L(2); % Number of crystals

        % Read this set of Euler angles...
        E = fscanf(fid, '%g %g %g %g', [4 tmp_nxtl]);
    
        % Build Euler angles array.
        tmp_eulers = zeros(3,tmp_nxtl);
        tmp_eulers(1,:) = E(1,:);
        tmp_eulers(2,:) = E(2,:);
        tmp_eulers(3,:) = E(3,:);
        
        % Bulild output...
        blocks = blocks + 1;
        if blocks == 1
            eulers = tmp_eulers;
            nxtl(blocks) = tmp_nxtl;
            strain(blocks) = tmp_header(17);
        elseif blocks == 2
            % Multiple textures, bundle into a cell array...
            eulers = {eulers};
            eulers(2) = {tmp_eulers};
            nxtl(blocks) = tmp_nxtl;
            strain(blocks) = tmp_header(17);
        else
            eulers(blocks) = {tmp_eulers};
            nxtl(blocks) = tmp_nxtl;
            strain(blocks) = tmp_header(17);
        end
            
        % Put the current file position in a sensible place.
        % This is yucky.
        frewind(fid);
        for i = 1:blocks
            fgetl(fid);
            fgetl(fid);
            fgetl(fid);
            fgetl(fid);
            for j = 1:tmp_nxtl % do we need to handle changing the number of xtals?
                fgetl(fid);
            end
        end
        
    end
    fclose(fid);
   
end