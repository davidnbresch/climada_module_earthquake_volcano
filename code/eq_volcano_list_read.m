function [vq_data,volcano_list_file_mat]=eq_volcano_list_read(volcano_list_file,check_plot)
% read ISC-GEM earthquake catalogue
% EQ eq data read
% NAME:
%   eq_GVP_volcano_list_read
% PURPOSE:
%   Read volcano (VQ) database, either
%
%   GVP_Volcano_List.xls from http://www.volcano.si.edu
%
%   or
%
%   NGDC_SignificantVolcanicEvents.xls from
%   http://www.ngdc.noaa.gov/nndc/struts/results?type_15=Like&query_15=&op_30=eq
%   &v_30=&ge_23=&le_23=&op_29=eq&v_29=&type_16=Like&query_16=&le_17=&ge_18=
%   &le_18=&ge_17=&op_20=eq&v_20=&ge_7=&le_7=&bt_24=&st_24=&type_25=EXACT
%   &query_25=None+Selected&bt_26=&st_26=&type_27=EXACT&query_27=None+Selected
%   &type_12=Exact&query_12=&type_11=Exact&query_11=&t=102557&s=50&d=50
%
%   next step: see vq_global_probabilistic
% CALLING SEQUENCE:
%   [vq_data,vq_data_file_mat]=eq_volcano_list_read(volcano_list_file_mat,check_plot)
% EXAMPLE:
%   [vq_data,vq_data_file_mat]=eq_volcano_list_read
% INPUTS:
% OPTIONAL INPUT PARAMETERS:
%   volcano_list_file: the filename of the .xls file with the tab
%       volacano_list with volcano information
%       if empty, the code tries a default name, if it does not exist, it
%       prompts the user to locate the file
%   check_plot: show a check plot (=1), or not (=0, default)
% OUTPUTS:
%   vq_data, a structure with
%       yyyy(line_i): year
%       mm(line_i)  : month
%       dd(line_i)  : day
%       lat(line_i): geographic latitude
%       lon(line_i): geographic longitude
%       elevation_m(line_i) : the elevation of the volcano
%       mag(line_i) : magnitude (moment magnitude, mw)
%       > also stored as .mat file. Please delete the .mat file manually in order
%       to re-read from the original .xls file
%   volcano_list_file_mat: the name of the .mat file vq_data is stored to
% MODIFICATION HISTORY:
% David N. Bresch, david.bresch@gmail.com, 20150217, initial
%-

vq_data=[]; % initialize output
vq_data_file_mat=''; % initialize output

% global climada_global
if ~climada_init_vars,return;end % init/import global variables

if ~exist('volcano_list_file','var'),volcano_list_file='';end
if ~exist('check_plot','var'),check_plot=0;end

% PARAMETERS
%
% prompt for file open dialogs
eq_dir=[fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];
%
% set default value for isc_gem_file if not given
volcano_list_file_default=[eq_dir filesep 'volcanoes' filesep 'NGDC_SignificantVolcanicEvents.xls'];
if isempty(volcano_list_file),volcano_list_file=volcano_list_file_default;end

if ~exist(volcano_list_file,'file')
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

end % eq_volcano_list_read

