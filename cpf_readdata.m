function [R2, R2sig, nucp, resname, carrier_freq, usef, nodatpts] =...
    cpf_readdata(filename)


%
% R2 is a cell array with a matrix for each carrier_freq. The matrix has
%   one row per residue and one column per nucp, so that the size is 
%   [number of residues x number of nucp].
%
% R2sig is a cell array with the corresponding errors.
%
% nucp is a cell array with a vector for each carrier_freq, 
%   with the CPMG field strengths [nucp_1 nucp_2 nucp_3 ...]
%
% usef is a cell array (same dimenstions as R2) containing boolean flags for
%   R2 values that should be fit. Replace R2 value by -1 to exclude from
%   fit.
%


fid = fopen(filename,'rt');
n = 0;
maxdelay = 0;
currdata = '';
R2sig = [];   % Might not be present
while 1
	tline = fgetl(fid);
    if isempty(tline), tline = fgetl(fid); end  % Step up to next non-empty line
	if ~ischar(tline), break, end
    
    lab = textscan(tline,'%s');
    if isequal(lab{1}{1},'DATA')
        if isequal(lab{1}{2},'CARRIER_FREQ')
            currdata = 'CARRIER_FREQ';
            currcf = 1;
        elseif isequal(lab{1}{2},'NUCPMG')
            currdata = 'NUCPMG';
        elseif isequal(lab{1}{2},'ASSIGNMENTS')
            currdata = 'ASSIGNMENTS';
            assix = 1;
        elseif isequal(lab{1}{2},'R2')
            currdata = 'R2';
            r2ix = 1;
        elseif isequal(lab{1}{2},'R2ERRORS')
            currdata = 'R2ERRORS';
            r2six = 1;
        end
    elseif isequal(lab{1}{1},'@CARRIER_FREQ')
        currcf = str2num(lab{1}{2});
        r2ix = 1;
        r2six = 1;
    else
        if isequal(currdata,'CARRIER_FREQ')
            carrier_freq(currcf) = sscanf(tline,'%f');
            currcf = currcf+1;
        elseif isequal(currdata,'NUCPMG')
            nucp{currcf} = sscanf(tline,'%f')';
        elseif isequal(currdata,'ASSIGNMENTS')
            [tmp, tmps] = strread(tline,'%f %s');
            resname(assix) = tmps;
            assix = assix+1;
        elseif isequal(currdata,'R2')
            tmp = sscanf(tline,'%f');
            R2{currcf}(r2ix,:) = tmp(2:length(tmp));
            % Flag for data that is not -1
            usef{currcf}(r2ix,:) = R2{currcf}(r2ix,:) ~= -1;
            % Number of data points per data set
            nodatpts{currcf}(r2ix,1) = sum(usef{currcf}(r2ix,:));
            r2ix = r2ix+1;
        elseif isequal(currdata,'R2ERRORS')
            tmp = sscanf(tline,'%f');
            R2sig{currcf}(r2six,:) = tmp(2:length(tmp));
            r2six = r2six+1;
        end
    end
end
fclose(fid);
