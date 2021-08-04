function x = binfileload(path,IDname,IDnum,CHnum,varargin)
% x = binfileload(path,IDname,IDnum,CHnum,'NStart',Nstart,'NRead',N,'CHname','CH')
% 'NumSamplesOffset', 'NumSamplesRead', and 'CHname' are optional inputs and
% do not need to be specified.
% This function loads single-precision, little-endian binary files without header information.
% The file name has the format: IDnameIDnum_CHnum.bin, where IDnum and CHnum have %03.0f format. 
% Also accomodates the format: IDnameIDnum_CHnameCHnum.bin with optional
% CHname input. Defaults to '' if not specified.
% Example: % ID001_004.bin
% path: file path, e.g., 'C:\Data'
% IDname: Root test name, e.g., 'ID'
% IDnum: Test number, e.g., 4
% CHnum: Channel number, e.g., 12
% Nstart: number of samples to offset from beginning of file.  Default is beginning of file
% NRead: Number of samples to read.  Default is the entire file
% CHname: Channel identifier. Default is nothing. e.g. 'CH'
%
% Note - the path and IDname variables must be characters (' '), not strings (" ")
%
% Author: Kent Gee, 11/14/13
% Update: Includes provisions for four-digit ID numbers native to AFR
% version 12.16.2.76 and newer. LM 7 JUL 2020
%
% Update: Adds functionality to read files with channel identifiers, not
% backwards compatible if using Nstart and N. Uses input parsing methods to achieve optional
% inputs. LM/MCA

% Parse inputs
p = inputParser;
p.addParameter('NStart',0);
p.addParameter('NRead',inf);
p.addParameter('CHname','');
p.parse(varargin{:});
    
% Assigning optional inputs
Nstart = p.Results.NStart; % number of samples to offset from beginning of file.  Default is beginning of file
N = p.Results.NRead; % Number of samples to read.  Default is the entire file
CHname = p.Results.CHname; % Channel identifier. Default is nothing.

filename=[path,filesep,IDname,sprintf('%03.0f',IDnum),'_',CHname,sprintf('%03.0f',CHnum),'.bin'];

% Error Checking: See if file exists
if ~isfile(filename) % Check to see if the three-digid ID number does not work, then try the four digit ID number.
    
    fourDigitFilename=[path,filesep,IDname,sprintf('%04.0f',IDnum),'_',CHname,sprintf('%03.0f',CHnum),'.bin'];
    
    if ~isfile(fourDigitFilename)
        disp('binfileload: The following file name is not found (with 3 or 4-digit ID number):')
        disp(filename) % display the 3-digit filename (for clarity)
        x = NaN;
        return;
    else
        disp('binfileload: 4-digit ID detected in file, using 4-digit ID number')
        filename = fourDigitFilename;
    end
    
end

% Read in data

fid = fopen(filename,'r');
Nstart = Nstart*4;   % Convert from samples to bytes
fseek(fid,Nstart,'bof');
x = fread(fid,N,'single');
fclose(fid);

end

