function eq_data=eq_signigeq_read(signigeq_file,check_plot)
% EQ eq data read
% NAME:
%   eq_signigeq_read
% PURPOSE:
%   Read earthquake (EQ) database, an .xlsx file downloaded from
%   www.ngdc.noaa.gov/nndc/struts/form?t=101650&s=1&d=1
%
%   for MATLAB to be able to read the .xlsx file, open it once in Excel,
%   convert columns SECOND, LATITUDE and LONGITUDE to Number format and
%   save as signigeq_CLEAN.xlsx (otherwise you will likely get errors of
%   the kind "Error using xlsread, Invalid zip file ... "
%
%   next step: see eq_global_probabilistic
% CALLING SEQUENCE:
%   eq_data=eq_signigeq_read(signigeq_file,check_plot)
% EXAMPLE:
%   eq_data=eq_signigeq_read('',1) % read std file, show epicenters plot
% INPUTS:
% OPTIONAL INPUT PARAMETERS:
%   signigeq_file: the filename of the .xlsx file with the raw data
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
%   to re-read from the original .xlsx file 
% MODIFICATION HISTORY:
% David N. Bresch, david.bresch@gmail.com, 20141015, initial
%-

eq_data=[]; % init

%global climada_global
if ~climada_init_vars,return;end % init/import global variables

%%if climada_global.verbose_mode,fprintf('*** %s ***\n',mfilename);end % show routine name on stdout

if ~exist('signigeq_file','var'),signigeq_file='';end
if ~exist('check_plot','var'),check_plot=0;end

% PARAMETERS
%
% where we prompt for file open dialogs
eq_dir=[fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];
%
% set default value for signigeq_file if not given
signigeq_file_default=[eq_dir filesep 'epicenters' filesep 'signigeq_CLEAN.xlsx'];
if isempty(signigeq_file),signigeq_file=signigeq_file_default;end
%
% approx. row count of data file
n_lines=5829; % only used for waitbar, an approx number is ok

if ~exist(signigeq_file,'file')
    [filename, pathname] = uigetfile(signigeq_file_default, 'Open epicenter database:');
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        signigeq_file=fullfile(pathname,filename);
    end
end

if ~exist(signigeq_file,'file'),fprintf('ERROR: file %s not found\n',signigeq_file);return;end
signigeq_file_mat=strrep(signigeq_file,'.xlsx','.mat');

if ~exist(signigeq_file_mat,'file')
    
    fprintf('reading raw data from %s ...\n',signigeq_file);
    eq_data_raw=climada_xlsread(0,signigeq_file,'results1',0);
    
    % eq_data_raw has the fields:
    %     YEAR
    %     MONTH
    %     DAY
    %     HOUR
    %     MINUTE
    %     SECOND
    %     FOCAL_DEPTH
    %     EQ_MAG_MW
    %     EQ_MAG_MS
    %     EQ_MAG_MB
    %     EQ_MAG_ML
    %     EQ_MAG_MFA
    %     EQ_MAG_UNK
    %     INTENSITY
    %     COUNTRY
    %     STATE
    %     LOCATION_NAME
    %     LATITUDE
    %     LONGITUDE
    %     REGION_CODE
    %     DEATHS
    %     DEATHS_DESCRIPTION
    %     MISSING
    %     MISSING_DESCRIPTION
    %     INJURIES
    %     INJURIES_DESCRIPTION
    %     DAMAGE_MILLIONS_DOLLARS
    %     DAMAGE_DESCRIPTION
    %     HOUSES_DESTROYED
    %     HOUSES_DESTROYED_DESCRIPTION
    %     HOUSES_DAMAGED
    %     HOUSES_DAMAGED_DESCRIPTION
    %     TOTAL_DEATHS
    %     TOTAL_DEATHS_DESCRIPTION
    %     TOTAL_MISSING
    %     TOTAL_MISSING_DESCRIPTION
    %     TOTAL_INJURIES
    %     TOTAL_INJURIES_DESCRIPTION
    %     TOTAL_DAMAGE_MILLIONS_DOLLARS
    %     TOTAL_DAMAGE_DESCRIPTION
    %     TOTAL_HOUSES_DESTROYED
    %     TOTAL_HOUSES_DESTROYED_DESCRIPTION
    %     TOTAL_HOUSES_DAMAGED
    %     TOTAL_HOUSES_DAMAGED_DESCRIPTION
    
    eq_data.yyyy=eq_data_raw.YEAR';
    eq_data.mm=eq_data_raw.MONTH';
    eq_data.dd=eq_data_raw.DAY';
    eq_data.hr=eq_data_raw.HOUR';
    eq_data.min=eq_data_raw.MINUTE';
    eq_data.sec=eq_data_raw.SECOND';
    eq_data.glat=eq_data_raw.LATITUDE'; % geographic latitude
    eq_data.glon=eq_data_raw.LONGITUDE'; % geographic longitude
    eq_data.dep=eq_data_raw.FOCAL_DEPTH'; % focal depth
    eq_data.mag=eq_data_raw.INTENSITY';
    
    % the brute way, not worth the effort to loop over all elements of
    % eq_data.* etc.    
    if ~isnumeric(eq_data.sec)
        fprintf('WARNING: fixing cell issue for SECOND (sec)\n');
        temp_data=zeros(1,length(eq_data.sec));
        for record_i=1:length(eq_data.sec)
            if isnumeric(eq_data.sec{record_i})
            temp_data(record_i)=eq_data.sec{record_i};
            end
        end % record_i
        eq_data=rmfield(eq_data,'sec');
        eq_data.sec=temp_data;
    end
    
    if ~isnumeric(eq_data.glat)
        fprintf('WARNING: fixing cell issue for LATITUDE (glat)\n');
        temp_data=zeros(1,length(eq_data.glat));
        for record_i=1:length(eq_data.glat)
            if isnumeric(eq_data.glat{record_i})
                temp_data(record_i)=eq_data.glat{record_i};
            end
        end % record_i
        eq_data=rmfield(eq_data,'glat');
        eq_data.glat=temp_data;
    end
    
    if ~isnumeric(eq_data.glon)
        fprintf('WARNING: fixing cell issue for LONGITUDE (glon)\n');
        temp_data=zeros(1,length(eq_data.glon));
        for record_i=1:length(eq_data.glon)
            if isnumeric(eq_data.glon{record_i})
                temp_data(record_i)=eq_data.glon{record_i};
            end
        end % record_i
        eq_data=rmfield(eq_data,'glon');
        eq_data.glon=temp_data;
    end
    
    % add a serial date number
    eq_data.nodetime_mat=datenum(eq_data.yyyy,eq_data.mm,eq_data.dd,eq_data.hr,eq_data.min,eq_data.sec);
    eq_data.filename=signigeq_file;
    eq_data.orig_event_flag=eq_data.yyyy*0+1;
    eq_data.n_epicenters_orig=length(eq_data.yyyy);
    eq_data.ens_size=0; % to indicate only original events
    
    fprintf('%i records read, stored as %s\n',length(eq_data.yyyy),signigeq_file_mat);
    save(signigeq_file_mat,'eq_data');
else
    load(signigeq_file_mat);
end


if check_plot
    % %raw version (fast)
    % climada_plot_world_borders
    % hold on
    % plot(eq_data.glon,eq_data.glat,'.r')
    fprintf('preparing epicenter plot...\n');
    climada_circle_plot(exp(eq_data.mag),eq_data.glon,eq_data.glat,'EQ signigeq',20);
end % check_plot

return
