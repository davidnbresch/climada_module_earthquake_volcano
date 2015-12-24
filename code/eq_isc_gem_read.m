function [eq_data,isc_gem_file_mat]=eq_isc_gem_read(isc_gem_file,check_plot)
% read ISC-GEM earthquake catalogue
% EQ eq data read
% NAME:
%   isc_gem_read
% PURPOSE:
%   Read earthquake (EQ) database, a .csv file downloaded from
%   www.isc.ac.uk/iscgem/download.php
%   for further information on the earthquake catalogue, see
%   www.isc.ac.uk/iscgem/index.php
%
%   Note: it can sometimes happen that the original epicenters file
%   isc-gem-cat.csv gets currupted. Just donwload it from
%   www.isc.ac.uk/iscgem/download.php (or directly
%   http://colossus.iris.washington.edu/iscgem/download/isc-gem-cat.zip)
%   again. 
%
%   next step: see eq_global_probabilistic
% CALLING SEQUENCE:
%   [eq_data,isc_gem_file_mat]=eq_isc_gem_read(isc_gem_file,check_plot)
% EXAMPLE:
%   eq_data=eq_isc_gem_read('',1) % read std file, show epicenters plot
% INPUTS:
% OPTIONAL INPUT PARAMETERS:
%   isc_gem_file: the filename of the .csv file with the raw data
%       if empty, the code tries a default name, if it does not exist, it
%       prompts the user to locate the file
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
%       mag(line_i) : magnitude (moment magnitude, mw)
%       > also stored as .mat file. Please delete the .mat file manually in order
%       to re-read from the original .csv file 
%   isc_gem_file_mat: the name of the .mat file eq_data is stored to
% MODIFICATION HISTORY:
% Melanie Bieli, melanie.bieli@bluewin.ch, 20141117
% David N. Bresch, david.bresch@gmail.com, 20141210, date conversion revised
% David N. Bresch, david.bresch@gmail.com, 20141223, issue with .csv file noted
%-

eq_data=[]; % initialize output
isc_gem_file_mat=''; % initialize output

% global climada_global
if ~climada_init_vars,return;end % init/import global variables

if ~exist('isc_gem_file','var'),isc_gem_file='';end
if ~exist('check_plot','var'),check_plot=0;end

% PARAMETERS
%
% prompt for file open dialogs
eq_dir=[fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];
%
% set default value for isc_gem_file if not given
isc_gem_file_default=[eq_dir filesep 'epicenters' filesep 'isc-gem-cat.csv'];
if isempty(isc_gem_file),isc_gem_file=isc_gem_file_default;end

if ~exist(isc_gem_file,'file')
    [filename, pathname] = uigetfile(isc_gem_file_default,'Choose ISC-GEM epicenter database:');
    if isequal(filename,0) || isequal(pathname,0)
        fprintf('No ISC-GEM epicenter database selected, consider downloading earthquake_volcano module again\n');
        fprintf('> See https://github.com/davidnbresch/climada_module_earthquake_volcano\n');
        return; % cancel
    else
        isc_gem_file=fullfile(pathname,filename);
    end
end

if ~exist(isc_gem_file,'file'),fprintf('ERROR: file %s not found\n',isc_gem_file);return;end

% Check whether there is already a .mat file containing the epicenter data
[fP,fN]=fileparts(isc_gem_file);
isc_gem_file_mat=[fP filesep fN '.mat'];

if ~climada_check_matfile(isc_gem_file,isc_gem_file_mat)

    % Open and read the ISC-GEM database 
    fid=fopen(isc_gem_file);
   
    % Skip commented lines
    % Not the most elegant way to do this 
    current_line = fgetl(fid);
    number_of_lines = 1;
    previous_line = '';
    while current_line(1) == '#' 
        previous_line = current_line;
        % current_line is a comment line (starting with #)
        current_line = fgetl(fid);
        number_of_lines = number_of_lines+1;    
    end
    
    fclose(fid); 
   
    % previous_line now contains the field names
    % get rid of the "#" at the beginning of previous_line
    previous_line = previous_line(2:end);
    
    % Create a vector containing the field names by splitting
    % previous_lines at the commas
    vector_of_field_names = strsplit(previous_line,',');
    index_date = 0;     % initialize indices
    index_lat = 0;
    index_lon = 0;
    index_depth = 0;
    index_mw = 0;
  
    % parse vector_of_field_names for the desired variables and store their
    % indices; e.g., the index 3 means that the corresponding variable is
    % the third comma separated value in each data line of isc_gem_file
    % TODO: Fehler ausgeben, wenn ein Feld nicht gefunden wird
    for i = 1:length(vector_of_field_names)
        if ~isempty(cell2mat(strfind(vector_of_field_names(i),'date')))
            % vector element i contains string 'date'
            index_date = i;
        elseif ~isempty(cell2mat(strfind(vector_of_field_names(i),'lat')))
            index_lat = i;
        elseif ~isempty(cell2mat(strfind(vector_of_field_names(i),'lon')))
            index_lon = i;
        elseif ~isempty(cell2mat(strfind(vector_of_field_names(i),'depth')))
            index_depth = i;
        elseif ~isempty(cell2mat(strfind(vector_of_field_names(i),'mw')))
            index_mw = i;
        end
    end
    
    % generate format string that will be used as input to textscan
    parse_string = ' %s';
    for i = 1:(length(vector_of_field_names) - 1)
        parse_string = [parse_string, ' %s'];
    end
    
    % read the isc_gem_file (starting after the header lines)
    fid = fopen(isc_gem_file);
    data = textscan(fid, parse_string,'HeaderLines',number_of_lines-1,'Delimiter',',','CollectOutput',0);
    fclose(fid);
    
    % construct eq_data by extracting the variables of interest from data
    formatIn = 'yyyy-mm-dd HH:MM:SS';
    
    % Note: in case an ERROR shows in the next line, download the
    % epicenters file (isc-gem-cat.csv) again.
    datenum_data = datenum(data{index_date},formatIn); 
    eq_data.yyyy=str2num(datestr(datenum_data,'yyyy'))';
    eq_data.mm=str2num(datestr(datenum_data,'mm'))';
    eq_data.dd=str2num(datestr(datenum_data,'dd'))';
    eq_data.hr=str2num(datestr(datenum_data,'HH'))';
    eq_data.min=str2num(datestr(datenum_data,'MM'))';
    eq_data.sec=str2num(datestr(datenum_data,'SS'))';
    
    eq_data.glat = str2double(data{index_lat}');
    eq_data.glon = str2double(data{index_lon}');
    eq_data.dep = str2double(data{index_depth}');
    eq_data.mag = str2double(data{index_mw}');

    % add further information: serial date number, filename etc.
    eq_data.datenum=datenum(eq_data.yyyy,eq_data.mm,eq_data.dd,eq_data.hr,eq_data.min,eq_data.sec);
    eq_data.filename=isc_gem_file;
    eq_data.orig_event_flag=eq_data.yyyy*0+1;
    eq_data.n_epicenters_orig=length(eq_data.yyyy);
    eq_data.ens_size=0; % to indicate only original events
    
    fprintf('saving as %s\n',isc_gem_file_mat);
    save(isc_gem_file_mat,'eq_data');
else
    load(isc_gem_file_mat);
end

if check_plot
    % raw
    %     climada_plot_world_borders
    %     hold on
    %     plot(eq_data.glon,eq_data.glat,'.r')
    fprintf('preparing epicenter plot...\n');
    climada_circle_plot(exp(eq_data.mag),eq_data.glon,eq_data.glat,'Earthquakes ISC-GEM',20);
end % check_plot

return

