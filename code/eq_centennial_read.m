function eq_data=eq_centennial_read(centennial_file,check_plot)
% EQ eq data read
% NAME:
%   eq_centennial_read
% PURPOSE:
%   Read earthquake (EQ) database, a raw text file downloaded from
%   earthquake.usgs.gov/data/centennial/ usually the file
%   http://earthquake.usgs.gov/data/centennial/centennial_Y2K.CAT
%   
%   see the associated readme at
%   http://earthquake.usgs.gov/data/centennial/centennial_README.rtf
%
%   next step: see eq_global_probabilistic
% CALLING SEQUENCE:
%   eq_data=eq_centennial_read(centennial_file,check_plot)
% EXAMPLE:
%   eq_data=eq_centennial_read('',1) % read std file, show epicenters plot
% INPUTS:
% OPTIONAL INPUT PARAMETERS:
%   centennial_file: the filename of the ASCII text file with the raw data
%       if empty, the code tries a default name, if it does not exist, it
%       promptes the user to locate the file
%   check_plot: show a check plot (=1), or not (=0, default)
% OUTPUTS:
%   eq_data, a structure with
%       yyyy(line_i): year
%       mm(line_i)  : month
%       dd(line_i)  : day
%       hr(line_i)  : origin hour
%       min(line_i) : origin min
%       sec(line_i) : origin sec
%       glat(line_i): geographic latitude
%       glon(line_i): geographic longitude
%       dep(line_i) : focal depth
%       mag(line_i) : magnitude
%   also stored as .mat file. Please delete the .mat file manually in order
%   to re-read from the original .txt file 
% MODIFICATION HISTORY:
% David N. Bresch, david.bresch@gmail.com, 20141010, initial
% David N. Bresch, david.bresch@gmail.com, 20141014, some minor improvements
%-

eq_data=[]; % init

global climada_global
if ~climada_init_vars,return;end % init/import global variables

%%if climada_global.verbose_mode,fprintf('*** %s ***\n',mfilename);end % show routine name on stdout

if ~exist('centennial_file','var'),centennial_file='';end
if ~exist('check_plot','var'),check_plot=0;end

% PARAMETERS
%
% where we prompt for file open dialogs
eq_dir=[climada_global.additional_dir filesep 'eq_global' filesep 'data'];
%
% set default value for centennial_file if not given
centennial_file_default=[eq_dir filesep 'epicenters' filesep 'centennial_Y2K.CAT.txt'];
if isempty(centennial_file),centennial_file=centennial_file_default;end
%
% approx. line count of raw data file
n_lines=13541; % only used for waitbar

if ~exist(centennial_file,'file')
    [filename, pathname] = uigetfile(centennial_file_default, 'Open epicenter database:');
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        centennial_file=fullfile(pathname,filename);
    end
end

if ~exist(centennial_file,'file'),fprintf('ERROR: file %s not found\n',centennial_file);return;end
centennial_file_mat=strrep(centennial_file,'.txt','.mat');

if ~exist(centennial_file_mat,'file')
    fid=fopen(centennial_file,'r');
    
    % init
    line_i=0;
    mod_step = 10; % first time estimate after 10 lines, then every 100
    
    fprintf('reading raw data from %s ...\n',centennial_file);
    
    h = waitbar(0.5,'Reading and converting data ...');
    set(h,'Name','Hazard EQ (centennial)');
    
    t0       = clock;
    % read raw data
    while not(feof(fid))
        
        % read one line
        line=fgetl(fid);
        line_i=line_i+1; % incerement
                
        % three lines for index, fllowed by a few sample data lines
        %0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001111111111111111111111111111111111111111111111111111111111111111111111
        %0000000001111111111222222222233333333334444444444555555555566666666667777777777888888888899999999990000000000111111111122222222223333333333444444444455555555556666666666
        %1234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
        % ABE        1900  1  5  19  0  0.00   -3.000 102.000   0.0 274   0 7.0 Ms AN2   0.0          0.0          0.0          0.0          0.0          0.0          0.0
        % EHB  AFEQ  1932 11  2  11  3 24.61  -22.977-112.338  15.0 685  38 6.6 UK G&R   6.8 UK PAS   0.0          0.0          0.0          0.0          0.0          0.0
        
        eq_data.yyyy(line_i)=str2double(line(13:16)); % year
        eq_data.mm(line_i)  =str2double(line(18:19)); % month
        eq_data.dd(line_i)  =str2double(line(21:22)); % day
        eq_data.hr(line_i)  =str2double(line(25:26)); % origin hour
        eq_data.min(line_i) =str2double(line(28:29)); % origin min
        eq_data.sec(line_i) =str2double(line(31:35)); % origin sec
        eq_data.glat(line_i)=str2double(line(37:44)); % geographic latitude
        eq_data.glon(line_i)=str2double(line(45:52)); % geographic longitude
        eq_data.dep(line_i) =str2double(line(53:58)); % focal depth
        eq_data.mag(line_i) =str2double(line(68:70)); % magnitude
                       
        if mod(line_i,mod_step)==0
            mod_step          = 100;
            t_elapsed         = etime(clock,t0)/line_i;
            lines_remaining  = n_lines-line_i;
            t_projected_sec   = t_elapsed*lines_remaining;
            if t_projected_sec<60
                msgstr = sprintf('est. %3.0f sec left (%i/%i lines)',t_projected_sec,   line_i,n_lines);
            else
                msgstr = sprintf('est. %3.1f min left (%i/%i lines)',t_projected_sec/60,line_i,n_lines);
            end
            waitbar(line_i/n_lines,h,msgstr); % update waitbar
        end
        
    end % while not(feof(fid))
    
    fclose(fid);
    if exist('h','var'), close(h), end % close waitbar
    
    % add a serial date number
    eq_data.nodetime_mat=datenum(eq_data.yyyy,eq_data.mm,eq_data.dd,eq_data.hr,eq_data.min,eq_data.sec);
    eq_data.filename=centennial_file;
    eq_data.orig_event_flag=eq_data.yyyy*0+1;
    eq_data.n_epicenters_orig=length(eq_data.yyyy);
    eq_data.ens_size=0; % to indicate only original events
    
    fprintf('%i lines read, stored as %s\n',line_i,centennial_file_mat);
    save(centennial_file_mat,'eq_data');
else
    load(centennial_file_mat);
end

if check_plot
    % raw
    %     climada_plot_world_borders
    %     hold on
    %     plot(eq_data.glon,eq_data.glat,'.r')
    fprintf('preparing epicenter plot...\n');
    climada_circle_plot(exp(eq_data.mag),eq_data.glon,eq_data.glat,'EQ centennial',20);
end % check_plot

return
