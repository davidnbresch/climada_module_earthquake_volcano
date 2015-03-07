function [vq_data,volcano_list_file_mat]=vq_volcano_list_read(volcano_list_file,check_plot)
% read ISC-GEM earthquake catalogue
% EQ eq data read
% NAME:
%   eq_GVP_volcano_list_read
% PURPOSE:
%   Read volcano (VQ) database, either
%
%   GVP_Volcano_List.xls from http://www.volcano.si.edu
%   (more volcanoes, but less information)
%
%   or
%
%   NGDC_SignificantVolcanicEvents.xls from
%   http://www.ngdc.noaa.gov/nndc/struts/results?type_15=Like&query_15=&op_30=eq
%   &v_30=&ge_23=&le_23=&op_29=eq&v_29=&type_16=Like&query_16=&le_17=&ge_18=
%   &le_18=&ge_17=&op_20=eq&v_20=&ge_7=&le_7=&bt_24=&st_24=&type_25=EXACT
%   &query_25=None+Selected&bt_26=&st_26=&type_27=EXACT&query_27=None+Selected
%   &type_12=Exact&query_12=&type_11=Exact&query_11=&t=102557&s=50&d=50
%   (more information, but less volcanoes)
%
%   next step: see vq_global_probabilistic
% CALLING SEQUENCE:
%   [vq_data,vq_data_file_mat]=vq_volcano_list_read(volcano_list_file_mat,check_plot)
% EXAMPLE:
%   [vq_data,vq_data_file_mat]=vq_volcano_list_read % returns NGDC...
% INPUTS:
% OPTIONAL INPUT PARAMETERS:
%   volcano_list_file: the filename of the .xls file with the tab
%       volacano_list with volcano information
%       if empty, the code tries a default name, if it does not exist, it
%       prompts the user to locate the file
%       Two short options are possible, too: 'GVP' or 'NGDC', see above
%       Default is 'NGDC', as this list contains more parameters
%   check_plot: show a check plot (=1), or not (=0, default)
% OUTPUTS:
%   vq_data, a structure with either (for
%       NGDC_SignificantVolcanicEvents.xls)
%       yyyy(line_i): year
%       mm(line_i)  : month
%       dd(line_i)  : day
%       lat(line_i): geographic latitude
%       lon(line_i): geographic longitude
%       elevation_m(line_i) : the elevation of the volcano
%       mag(line_i) : magnitude (moment magnitude, mw)
%       ...and other fields, same names as column headers in Excel
%   or for GVP_Volcano_List.xls
%       lat(i): geographic latitude
%       lon(line_i): geographic longitude
%       elevation_m(line_i) : the elevation of the volcano
%       filename: the filename the data has been read from
%       Volcano_Number(i)
%       Volcano_Name{i}
%       country_name{i}
%       Primary_Volcano_Type{i}
%       Activity_Evidence{i}
%       Last_Known_Eruption{i}
%       Region{i}
%       Subregion{i}
%       Dominant_Rock_Type{i}
%       Tectonic_Setting{i}
%       > also stored as .mat file. Please delete the .mat file manually in order
%       to re-read from the original .xls file
%   volcano_list_file_mat: the name of the .mat file vq_data is stored to
% MODIFICATION HISTORY:
% David N. Bresch, david.bresch@gmail.com, 20150217, initial
% David N. Bresch, david.bresch@gmail.com, 20150223, some trimming
%-

vq_data=[]; % initialize output
volcano_list_file_mat=''; % initialize output

% global climada_global
if ~climada_init_vars,return;end % init/import global variables

if ~exist('volcano_list_file','var'),volcano_list_file='';end
if ~exist('check_plot','var'),check_plot=0;end

% PARAMETERS
%
eq_dir=[fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];
%
% define default value(s) for volcano list file if not given
GVP_volcano_list_file_default =[eq_dir filesep 'volcanoes' filesep 'GVP_Volcano_List.xls'];
NGDC_volcano_list_file_default=[eq_dir filesep 'volcanoes' filesep 'NGDC_SignificantVolcanicEvents.xls'];

if isempty(volcano_list_file),volcano_list_file='NGDC';end

if strcmp(volcano_list_file,'GVP')
    volcano_list_file=GVP_volcano_list_file_default;
elseif strcmp(volcano_list_file,'NGDC')
    volcano_list_file=NGDC_volcano_list_file_default;
end

if ~exist(volcano_list_file,'file')
    volcano_list_file=NGDC_volcano_list_file_default;
    [filename, pathname] = uigetfile(volcano_list_file,'Choose volcano database:');
    if isequal(filename,0) || isequal(pathname,0)
        fprintf('No volcano database selected, consider downloading eq_global module again\n');
        fprintf('> See https://github.com/davidnbresch/climada_module_eq_global\n');
        return; % cancel
    else
        volcano_list_file=fullfile(pathname,filename);
    end
end

if ~exist(volcano_list_file,'file'),fprintf('ERROR: file %s not found\n',volcano_list_file);return;end

% Check whether there is already a .mat file containing the epicenter data
[fP,fN]=fileparts(volcano_list_file);
volcano_list_file_mat=[fP filesep fN '.mat'];

if ~climada_check_matfile(volcano_list_file,volcano_list_file_mat)
    
    vq_data=climada_xlsread('no',volcano_list_file,'volcano_list',1,-999);
    
    % we check for Last_Known_Eruption, otherwise assume the data contains the
    % field Year
    if isfield(vq_data,'Last_Known_Eruption')
        vq_data.Year=vq_data.lon*0+NaN; % init
        for record_i=1:length(vq_data.Last_Known_Eruption)
            record_str=char(vq_data.Last_Known_Eruption{record_i});
            BCE_pos=strfind(record_str,'BCE');
            CE_pos=strfind(record_str,'CE');
            if ~isempty(BCE_pos)
                vq_data.Year(record_i)=-str2num(record_str(1:BCE_pos-1));
            elseif ~isempty(CE_pos)
                vq_data.Year(record_i)=str2num(record_str(1:CE_pos-1));
            end % record_i
        end
    end % Last_Known_Eruption
    
    vq_data.orig_event_flag=vq_data.lon*0+1; % indicate historic eruptions
    vq_data.ens_size=0; % no ensemble, just the original events
    vq_data.datenum=datenum(vq_data.Year,0,0);
    
    if ~isfield(vq_data,'duration_h')
        vq_data.duration_h=vq_data.lon*0+1; % 1 hour
    end

    fprintf('saving as %s\n',volcano_list_file_mat);
    save(volcano_list_file_mat,'vq_data');
else
    load(volcano_list_file_mat);
end

if check_plot
    if isfield(vq_data,'damage')
        fprintf('preparing volcano plot (circle diameter for damage)...\n');
        climada_circle_plot(sqrt(vq_data.damage),vq_data.lon,vq_data.lat,'Volcanoes',20);
        hold on;plot(vq_data.lon,vq_data.lat,'.r')
    else
        climada_plot_world_borders
        hold on
        plot(vq_data.lon,vq_data.lat,'.r')
    end
end % check_plot

end % vq_volcano_list_read

