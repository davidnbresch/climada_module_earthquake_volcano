function hazard=vq_global_hazard_set(vq_data,hazard_set_file,centroids,TEST_volcano_preselection)
% Create volcanoe hazard set
% MODULE
% eq_global
% NAME:
%   vq_global_hazard_set
% PURPOSE:
%   generate a volcano (VQ) hazard event set, starting from volcano data
%   calculating ash thickness (distance and depth).
%
%   previous step: see vq_volcano_list_read
% CALLING SEQUENCE:
%   hazard=vq_global_hazard_set(vq_data,hazard_set_file,centroids)
% EXAMPLE:
%   hazard=vq_global_hazard_set(vq_volcano_list_read)
%   hazard=vq_global_hazard_set(vq_volcano_list_read,'',[],1) % TEST
% INPUTS:
% OPTIONAL INPUT PARAMETERS:
%   vq_data: a structure with volcano information, see vq_volcano_list_read
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
%   TEST_volcano_preselection: if =1, check and show centroids and volcano
%       positions, if =2 STOP after check, Default=0, no check
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
% David N. Bresch, david.bresch@gmail.com, 20150302, initial
%-

hazard=[]; % init

global climada_global
if ~climada_init_vars,return;end % init/import global variables

%if climada_global.verbose_mode,fprintf('*** %s ***\n',mfilename);end % show routine name on stdout

% poor man's version to check arguments
if ~exist('vq_data','var'),vq_data=[];end
if ~exist('hazard_set_file','var'),hazard_set_file=[];end
if ~exist('centroids','var'),centroids=[];end
if ~exist('TEST_volcano_preselection','var'),TEST_volcano_preselection=0;end

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
% define radius around centroids we search for volcanoes
EPM=2; % in degrees, i.e. =2 means about 220 km

% prompt for vq_data if not given
if isempty(vq_data) % local GUI
    vq_data_dir = [eq_dir filesep '*.mat'];
    [filename, pathname] = uigetfile(vq_data_dir,'Select volcanoes data:');
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        vq_data=fullfile(pathname,filename);
    end
end
if ~isstruct(vq_data) % load, if filename given
    vq_data_file=vq_data;vq_data=[];
    load(vq_data_file);
end

% prompt for hazard_set_file if not given
if isempty(hazard_set_file) && TEST_volcano_preselection<2 % local GUI
    hazard_set_file      = [climada_global.data_dir filesep 'hazards' filesep 'VQXX_hazard.mat'];
    [filename, pathname] = uiputfile(hazard_set_file, 'Save VQ hazard set as:');
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
        fprintf('WARNING: Special mode, TEST centroids grid (Naples) created in %s\n',mfilename);
        lon=14.426;lat=40.821; % Vesuvius
        [X,Y]=meshgrid(lon-1:.01:lon+1,lat-2:.01:lat+1); % 2D grid
        centroids.lon=reshape(X,numel(X),1);centroids.lat=reshape(Y,numel(X),1); % 2D grid
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
in_centroids_poly = inpolygon(vq_data.lon,vq_data.lat,centroids_edges_x,centroids_edges_y);

% skip small events (for TESTs)
in_centroids_poly(vq_data.Cloud_height_km<10)=0;

if TEST_volcano_preselection>0
    climada_plot_world_borders; hold on
    plot(vq_data.lon,vq_data.lat,'.b','MarkerSize',2)
    plot(vq_data.lon(in_centroids_poly),vq_data.lat(in_centroids_poly),'og','MarkerSize',3)
    %plot(vq_data.lon(in_centroids_poly),vq_data.lat(in_centroids_poly),'xg')
    plot(centroids.lon,centroids.lat,'.r','MarkerSize',1)
    axis equal
    fprintf('%i (of %i) volcanoes in range\n',sum(in_centroids_poly),length(vq_data.lon));
    in_centroids_pos=find(in_centroids_poly);
    for i=1:length(in_centroids_pos)
        fprintf('%s: %f (damage %f, cloud height %f km)\n',...
            char(vq_data.Name(in_centroids_pos(i))),...
            vq_data.VEI(in_centroids_pos(i)),...
            vq_data.damage(in_centroids_pos(i)),...
            vq_data.Cloud_height_km(in_centroids_pos(i)));
    end %i
    if TEST_volcano_preselection>1
        fprintf('STOP after epicenter preselection\n');
        return
    end
end

if isfield(vq_data,'Year')
    min_year   = min(vq_data.Year);
    max_year   = max(vq_data.Year);
else
    fprintf('Warning: no Year information, assuming dummy 0..2000\n');
    min_year=0;
    max_year=2000;
end
orig_years = max_year - min_year+1;

% fill the hazard structure
hazard.reference_year   = hazard_reference_year;
hazard.lon              = centroids.lon;
hazard.lat              = centroids.lat;
hazard.centroid_ID      = centroids.centroid_ID;
hazard.event_count      = length(vq_data.lon);
hazard.event_ID         = 1:hazard.event_count;
hazard.orig_years       = orig_years;
hazard.orig_event_count = sum(vq_data.orig_event_flag);
hazard.orig_event_flag  = vq_data.orig_event_flag;
hazard.yyyy             = vq_data.Year;
hazard.datenum          = vq_data.datenum;
hazard.vq_data_filename = vq_data.filename;


% allocate the hazard array (sparse, to manage memory)
% fprintf('%s: spalloc(%i,%i,%i)\n',mfilename,hazard.event_count,length(hazard.lon),...
%     ceil(hazard.event_count*length(hazard.lon)*hazard_arr_density));
hazard.intensity = spalloc(hazard.event_count,length(hazard.lon),...
    ceil(hazard.event_count*length(hazard.lon)*hazard_arr_density));

% init eruption parameters
if isfield(vq_data,'Cloud_height_km')
    H=vq_data.Cloud_height_km; % H, eruptive column height in km
else
    H=vq_data.lon*0+14; % default, 14 km
end
if isfield(vq_data,'U_vel_kmh')
    % wind velocity U (in km/h, indicative 50-100 km/h most often)
    U_vel=vq_data.U_vel_kmh; % km/h
else
    U_vel=vq_data.lon*0+50; % default, 50 km/h
end
if isfield(vq_data,'U_phi_rad')
    % direction of wind in radian (-pi..pi), counterclockwise from (dx=1,dy=0), i.e. from direction East, West=pi, North=pi/2
    U_phi=vq_data.U_phi_rad; % in rad
else
    U_phi=vq_data.lon*0; % default, wind eastward (westerlies dominate)
end
if isfield(vq_data,'duration_h')
    % duration of the high-intensity phase of the eruption in hours
    tau=vq_data.duration_h; % in hours
else
    tau=vq_data.lon*0+1; % 1h
end

t0       = clock;
n_events = hazard.event_count;
n_events_eff=sum(in_centroids_poly);

msgstr   = sprintf('processing %i (of globally %i) volcanoes',n_events_eff,n_events);
mod_step = 10; % first time estimate after 10 volcanoes, then every 100
if climada_global.waitbar
    fprintf('%s (updating waitbar with estimation of time remaining every 100th volcano)\n',msgstr);
    h        = waitbar(0,msgstr);
    set(h,'Name','Hazard VQ: thephra/ash thickness (cm)');
else
    fprintf('%s (waitbar suppressed)\n',msgstr);
    format_str='%s';
end

event_i_eff=0; % since we only process a subset
for event_i=1:n_events
    
    if in_centroids_poly(event_i)
        
        event_i_eff=event_i_eff+1;
        
        hazard.intensity(event_i,:)=sparse(vq_tephra_field_cm(centroids,...
            vq_data.lon(event_i),vq_data.lat(event_i),...
            H(event_i),U_vel(event_i),U_phi(event_i),tau(event_i)));
        
        % progress update
        if mod(event_i_eff,mod_step)==0
            mod_step          = 100;
            t_elapsed   = etime(clock,t0)/event_i_eff;
            events_remaining  = n_events_eff-event_i_eff;
            t_projected_sec   = t_elapsed*events_remaining;
            if t_projected_sec<60
                msgstr = sprintf('est. %3.0f sec left (%i/%i volcanoes)',t_projected_sec,   event_i_eff,n_events_eff);
            else
                msgstr = sprintf('est. %3.1f min left (%i/%i volcanoes)',t_projected_sec/60,event_i_eff,n_events_eff);
            end
            if climada_global.waitbar
                waitbar(event_i_eff/n_events_eff,h,msgstr); % update waitbar
            else
                fprintf(format_str,msgstr);
                format_str=[repmat('\b',1,length(msgstr)) '%s'];
            end
        end
        
    end % in_centroids_poly
    
end %event_i
if climada_global.waitbar
    close(h) % dispose waitbar
else
    fprintf(format_str,''); % move carriage to begin of line
end

t_elapsed = etime(clock,t0);
msgstr    = sprintf('processing %i volcanoes took %3.2f min (%3.4f sec/event)',n_events_eff,t_elapsed/60,t_elapsed/n_events_eff);
fprintf('%s\n',msgstr);

hazard

% number of probabilistic eruptions per original one
event_frequency = 1/(hazard.orig_years*(hazard.orig_event_count/hazard.event_count));
fprintf('%i original, %i probabilistic years, event frequency = %f\n',orig_years,orig_years*(vq_data.ens_size+1),event_frequency);

% not transposed, just regular
hazard.frequency         = ones(1,hazard.event_count)*event_frequency;
hazard.matrix_density    = nnz(hazard.intensity)/numel(hazard.intensity);
hazard.eq_comment        = msgstr;
hazard.peril_ID          = 'VQ';
hazard.filename          = hazard_set_file;
hazard.comment           = sprintf('VQ hazard event set, generated %s',datestr(now));
hazard.date              = datestr(now);
hazard.units             = 'cm';

fprintf('saving VQ hazard set as %s\n',hazard_set_file);
save(hazard_set_file,'hazard');
%save(hazard_set_file,'hazard','-v7.3'); % see note on next line:
% Warning: Variable 'hazard' cannot be saved to a MAT-file whose version is older than 7.3. To save this variable, use the -v7.3 switch.
% to avoid this warning, the switch is used. david's comment: only shows for large hazard sets, seems to be due to huge size of hazard

end % vq_global_hazard_set