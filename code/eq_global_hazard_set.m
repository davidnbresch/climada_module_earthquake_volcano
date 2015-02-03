function hazard=eq_global_hazard_set(eq_data,hazard_set_file,centroids,TEST_epicenter_preselection,correction,a1,a2,a3,a4,b)
% Create earthquake hazard set
% MODULE
% eq_global
% NAME:
%   eq_global_hazard_set
% PURPOSE:
%   generate an earthqake (EQ) hazard event set, starting from epicenters
%   calculating attenuation (distance and depth) to convert Richter
%   magnitude to a modified Mercalli intensity (MMI)
%
%   previous step: see eq_global_probabilistic or eq_isc_gem_read (or
%   eq_centennial_read, or eq_signieq_read)
% CALLING SEQUENCE:
%   hazard=eq_global_hazard_set(eq_data,hazard_set_file,centroids)
% EXAMPLE:
%   hazard=eq_global_hazard_set(eq_global_probabilistic(eq_isc_gem_read))
% INPUTS:
% OPTIONAL INPUT PARAMETERS:
%   eq_data: a structure with EQ epicenters, see eq_isc_gem_read
%       > promted for if not given
%   hazard_set_file: the name and path of the hazard set file
%       > promted for if not given
%   centroids: the variable grid centroids (see climada_centroids_read)
%       a structure with
%           Longitude(1,:): the longitudes
%           Latitude(1,:): the latitudes
%           centroid_ID(1,:): a unique ID for each centroid, simplest: 1:length(Longitude)
%       or a .mat-file which contains a centroids struct (saved by
%       climada_centroids_read) or the filename of an Excel file (the original
%       input to climada_centroids_read) which holds the centroids, in
%       which case climada_centroids_read is called.
%       > promted for .mat or .xls filename if not given
%       NOTE: if you then select Cancel, a regular default grid is used
%       (TEST mode), see hard-wired definition in code (a rectangular area in California)
%   TEST_epicenter_preselection: whether we show the epicenters selected
%       for processing (=1) or not (=0, default)
%       if =2, STOP after plot of preselection (to check preselection only)
%   correction,a1,a2,a3,a4,b: parameters defining the attenuation function. See
%       eq_global-master/data/system/attenuation_parameters.xlsx to use
%       parameters for specific regions; otherwise the function
%       eq_global_attenuation will use default values that represent a "global
%       average attenuation function"
% OUTPUTS:
%   hazard: a struct, the hazard event set, more for tests, since the
%       hazard set is stored as hazard_set_file, see code
%       lon(centroid_i): the longitude of each centroid
%       lat(centroid_i): the latitude of each centroid
%       centroid_ID(centroid_i): a unique ID for each centroid
%       peril_ID: just an ID identifying the peril, e.g. 'EQ' for
%       earthquake
%       comment: a free comment, normally containing the time the hazard
%           event set has been generated
%       orig_years: the original years the set is based upon
%       orig_event_count: the original events
%       event_count: the total number of events in the set, including all
%           probabilistic ones, hence event_count>=orig_event_count
%       orig_event_flag(event_i): a flag for each event, whether it's an original
%           (1) or probabilistic (0) one
%       event_ID: a unique ID for each event
%       date: the creation date of the set
%       arr(event_i,centroid_i),sparse: the hazard intensity of event_i at
%           centroid_i
%       frequency(event_i): the frequency of each event
%       matrix_density: the density of the sparse array hazard.intensity
%       filename: the filename of the hazard event set (if passed as a
%           struct, this is often useful)
%
%   simple check for hazard content: hist(full(hazard.intensity(find(hazard.intensity))))
% MODIFICATION HISTORY:
% David N. Bresch, david.bresch@gmail.com, 20141012
% Melanie Bieli, melanie.bieli@bluewin.ch, 20141223, added correction,a1,a2,a3,a4
% Melanie Bieli, melanie.bieli@bluewin.ch, 20150120, hazard_arr_density=0.01
%-

hazard=[]; % init

global climada_global
if ~climada_init_vars,return;end % init/import global variables

%if climada_global.verbose_mode,fprintf('*** %s ***\n',mfilename);end % show routine name on stdout

% poor man's version to check arguments
if ~exist('eq_data','var'),eq_data=[];end
if ~exist('hazard_set_file','var'),hazard_set_file=[];end
if ~exist('centroids','var'),centroids=[];end
if ~exist('TEST_epicenter_preselection','var'),TEST_epicenter_preselection=0;end

if ~exist('correction','var') || isempty(correction), correction = 0;  end
if ~exist('a1','var') || isempty(a1), a1 = 1.7;         end
if ~exist('a2','var') || isempty(a2), a2 = 1.5;         end
if ~exist('a3','var') || isempty(a3), a3 = 1.1726;      end
if ~exist('a4','var') || isempty(a4), a4 = 0.00106;     end
if ~exist('b','var')  || isempty(b),  b = 0;            end

% PARAMETERS
%
% where we prompt for file open dialogs
eq_dir=[fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];
%
% since we store the hazard as sparse array, we need an a-priory estimation
% of its density
%hazard_arr_density=0.03; % 3% sparse hazard array density (estimated)
hazard_arr_density=0.01; % 1% sparse hazard array density (estimated)
%
% define the reference year for this hazard set
hazard_reference_year = climada_global.present_reference_year; % does not really matter for EQ
%
% the earthquake epicenter preselection margin, i.e. to widen the box
% around the centroids by EPM degrees.
EPM=2; % in degrees, default=2 (approx. 200km)

% prompt for eq_data if not given
if isempty(eq_data) % local GUI
    eq_data_dir = [eq_dir filesep '*.mat'];
    [filename, pathname] = uigetfile(eq_data_dir,'Select epicenters data:');
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        eq_data=fullfile(pathname,filename);
    end
end
if ~isstruct(eq_data) % load, if filename given
    eq_data_file=eq_data;eq_data=[];
    load(eq_data_file);
end

% prompt for hazard_set_file if not given
if isempty(hazard_set_file) && TEST_epicenter_preselection<2 % local GUI
    hazard_set_file      = [climada_global.data_dir filesep 'hazards' filesep 'EQXX_hazard.mat'];
    [filename, pathname] = uiputfile(hazard_set_file, 'Save EQ hazard set as:');
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        hazard_set_file = fullfile(pathname,filename);
    end
end

% prompt for centroids if not given
if isempty(centroids) % local GUI
    centroids_default    = [climada_global.system_dir filesep '*.mat'];
    [filename, pathname] = uigetfile({'*.mat;*.xls'},'Select centroids (.mat or .xls):',centroids_default);
    if isequal(filename,0) || isequal(pathname,0)
        % TEST centroids
        fprintf('WARNING: Special mode, TEST centroids grid (California) created in %s\n',mfilename);
        ii=0;
        for lon_i=-130:1:-110
            for lat_i=25:1:45
                ii=ii+1;
                centroids.lon(ii)=lon_i;
                centroids.lat(ii)=lat_i;
            end
        end
        centroids.centroid_ID=1:length(centroids.lon);
    else
        centroids_file=fullfile(pathname,filename);
        [~,~,fE]=fileparts(centroids_file);
        if strcmp(fE,'.xls')
            fprintf('reading centroids from %s\n',centroids_file);
            centroids=climada_centroids_read(centroids_file);
        else
            centroids=centroids_file;
        end
        
    end
end

if ~isstruct(centroids) % load, if filename given
    centroids_file=centroids;centroids=[];
    fprintf('centroids read from %s\n',centroids_file);
    load(centroids_file); % contains centrois as a variable
end

% figure which epicenters can affect the region
centroids_rect = [min(centroids.lon)-EPM max(centroids.lon)+EPM min(centroids.lat)-EPM max(centroids.lat)+EPM];
centroids_edges_x = [centroids_rect(1), centroids_rect(1), centroids_rect(2), centroids_rect(2)];
centroids_edges_y = [centroids_rect(3), centroids_rect(4), centroids_rect(4), centroids_rect(3)];
in_centroids_poly = inpolygon(eq_data.glon,eq_data.glat,centroids_edges_x,centroids_edges_y);
if TEST_epicenter_preselection>0
    climada_plot_world_borders; hold on
    plot(eq_data.glon,eq_data.glat,'.b')
    plot(eq_data.glon(in_centroids_poly),eq_data.glat(in_centroids_poly),'.g')
    if TEST_epicenter_preselection>1
        fprintf('STOP after epicenter preselection\n');
        return
    end
end

min_year   = min(eq_data.yyyy);
max_year   = max(eq_data.yyyy);
orig_years = max_year - min_year+1;

% fill the hazard structure
hazard.reference_year   = hazard_reference_year;
hazard.lon              = centroids.lon;
hazard.lat              = centroids.lat;
hazard.centroid_ID      = centroids.centroid_ID;
hazard.event_count      = length(eq_data.mag);
hazard.event_ID         = 1:hazard.event_count;
hazard.orig_years       = orig_years;
hazard.orig_event_count = eq_data.n_epicenters_orig;
hazard.orig_event_flag  = zeros(1,hazard.event_count);
hazard.orig_event_flag(1:hazard.orig_event_count)=1;
hazard.yyyy             = eq_data.yyyy;
hazard.mm               = eq_data.mm;
hazard.dd               = eq_data.dd;
hazard.nodetime_mat     = eq_data.datenum;

% allocate the hazard array (sparse, to manage memory)
fprintf('%s: spalloc(%i,%i,%i)\n',mfilename,hazard.event_count,length(hazard.lon),...
    ceil(hazard.event_count*length(hazard.lon)*hazard_arr_density));
hazard.intensity = spalloc(hazard.event_count,length(hazard.lon),...
    ceil(hazard.event_count*length(hazard.lon)*hazard_arr_density));

t0       = clock;
n_events = hazard.event_count;
n_events_eff=sum(in_centroids_poly);

msgstr   = sprintf('processing %i (of globally %i) epicenters',n_events_eff,n_events);
if climada_global.waitbar
    fprintf('%s (updating waitbar with estimation of time remaining every 100th epicenter)\n',msgstr);
    h        = waitbar(0,msgstr);
    set(h,'Name','Hazard EQ: shaking intensity (MMI)');
    mod_step = 10; % first time estimate after 10 events, then every 100
else
    fprintf('%s (waitbar suppressed)\n',msgstr);
    mod_step=n_events+10;
end

event_i_eff=0; % since we only process a subset
for event_i=1:n_events
    
    if in_centroids_poly(event_i)
        
        event_i_eff=event_i_eff+1;
        
        hazard.intensity(event_i,:)=eq_global_attenuation(eq_data.glat(event_i),...
            eq_data.glon(event_i),eq_data.mag(event_i),...
            centroids, 0, eq_data.dep(event_i), correction, a1,a2,a3,a4,b);
        
        if mod(event_i_eff,mod_step)==0
            mod_step          = 100;
            t_elapsed   = etime(clock,t0)/event_i_eff;
            events_remaining  = n_events_eff-event_i_eff;
            t_projected_sec   = t_elapsed*events_remaining;
            if t_projected_sec<60
                msgstr = sprintf('est. %3.0f sec left (%i/%i epicenters)',t_projected_sec,   event_i_eff,n_events_eff);
            else
                msgstr = sprintf('est. %3.1f min left (%i/%i epicenters)',t_projected_sec/60,event_i_eff,n_events_eff);
            end
            waitbar(event_i_eff/n_events_eff,h,msgstr); % update waitbar
        end
        
    end % in_centroids_poly
    
end %event_i
if exist('h','var'),close(h);end % dispose waitbar

t_elapsed = etime(clock,t0);
msgstr    = sprintf('processing %i epicenters took %3.2f min (%3.4f sec/event)',n_events_eff,t_elapsed/60,t_elapsed/n_events_eff);
fprintf('%s\n',msgstr);

% number of derived tracks per original one
event_frequency = 1/(orig_years*(eq_data.ens_size+1));
fprintf('%i original, %i probabilistic years, event frequency = %f\n',orig_years,orig_years*(eq_data.ens_size+1),event_frequency);

% not transposed, just regular
hazard.frequency         = ones(1,hazard.event_count)*event_frequency;
hazard.matrix_density    = nnz(hazard.intensity)/numel(hazard.intensity);
hazard.eq_comment        = msgstr;
hazard.peril_ID          = 'EQ';
hazard.filename          = hazard_set_file;
hazard.comment           = sprintf('EQ hazard event set, generated %s',datestr(now));
hazard.date              = datestr(now);

fprintf('saving EQ hazard set as %s\n',hazard_set_file);
save(hazard_set_file,'hazard','-v7.3'); % see note on next line:
% Warning: Variable 'hazard' cannot be saved to a MAT-file whose version is older than 7.3. To save this variable, use the -v7.3 switch.
% to avoid this warning, the switch is used. david's comment: only shows for large hazard sets, seems to be due to huge size of hazard
return
