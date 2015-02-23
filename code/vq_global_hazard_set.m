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
%       positions, Default=0, no check
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
% Melanie Bieli, melanie.bieli@bluewin.ch, 20150223, initial
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
        ii=0;
        for lon_i=13:.025:16
            for lat_i=38:.025:42
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
in_centroids_poly = inpolygon(vq_data.lon,vq_data.lat,centroids_edges_x,centroids_edges_y);
if TEST_volcano_preselection>0
    climada_plot_world_borders; hold on
    plot(vq_data.lon,vq_data.lat,'.b')
    plot(vq_data.lon(in_centroids_poly),vq_data.lat(in_centroids_poly),'og')
    plot(vq_data.lon(in_centroids_poly),vq_data.lat(in_centroids_poly),'xg')
    plot(centroids.lon,centroids.lat,'.r')
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
hazard.orig_event_count = vq_data.n_volcanoes_orig;
hazard.orig_event_flag  = zeros(1,hazard.event_count);
hazard.orig_event_flag(1:hazard.orig_event_count)=1;
hazard.yyyy             = vq_data.Year;
% hazard.mm               = vq_data.mm;
% hazard.dd               = vq_data.dd;
hazard.nodetime_mat     = vq_data.datenum;

% allocate the hazard array (sparse, to manage memory)
% fprintf('%s: spalloc(%i,%i,%i)\n',mfilename,hazard.event_count,length(hazard.lon),...
%     ceil(hazard.event_count*length(hazard.lon)*hazard_arr_density));
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

cos_vq_data_lat = cos(vq_data.lat/180*pi);

event_i_eff=0; % since we only process a subset
for event_i=1:n_events
    
    if in_centroids_poly(event_i)
        
        event_i_eff=event_i_eff+1;
                
        % according to A. O. Gonzalez-Mellado and S. De la Cruz-Reyna, 2010: A
        % simple semi-empirical approach to model thickness of ash-deposits for
        % different eruption scenarios. Nat. Hazards Earth Syst. Sci., 10,
        % 2241?2257, 2010, direct:
        % www.nat-hazards-earth-syst-sci.net/10/2241/2010/nhess-10-2241-2010.pdf
        
        if isfield(vq_data,'Cloud_height_km')
            H=vq_data.Cloud_height_km(event_i); % H, eruptive column height in km
        else
            H=14; % default, 14 km
        end
        rho=1100; % density of the falling material (kgm-3, indicative value 1100),
        alpha_param=2.535-0.051*H; % the rate at which the deposit thickness decays with distance
        Htropopause=15.5; % tropopause height in km
        tau=24; % duration of the high-intensity phase of the eruption in hours
        
        U=100; % wind velocity U (in km/h, indicative 50-100 km/h most often), duration of the high-intensity phase of the eruption ? (in hours, e.g. Pinatubo 1-5h) and D (in km2/h) as follows (to distinguish between event that do and do not penetrate the stratosphere):
        phi=pi/4; % angle (in radian, North is 0)
        
        if H<Htropopause
            D=-4.189*H+114.407;
        else % above)
            D=52.822*H-770.17;
        end
        
        for centroid_i=1:length(centroids.lon)
            
            % distance to eruption center im km
            r=sqrt(((centroids.lon(centroid_i)-vq_data.lon(event_i))*cos_vq_data_lat(event_i))^2+...
                (centroids.lat(centroid_i)-vq_data.lat(event_i))^2)*111.12; % km           
            
            T=14.28*H^4*tau^(alpha_param/2) / ...
                ( (2*D)^(1-alpha_param/2)*rho ) * ...
                exp(-U/(2*D)*r*(1-cos(phi)))*r^-alpha_param;
            
            hazard.intensity(event_i,centroid_i)=T;
            
        end % centroid_i
        
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
            waitbar(event_i_eff/n_events_eff,h,msgstr); % update waitbar
        end
        
    end % in_centroids_poly
    
end %event_i
if exist('h','var'),close(h);end % dispose waitbar

t_elapsed = etime(clock,t0);
msgstr    = sprintf('processing %i volcanoes took %3.2f min (%3.4f sec/event)',n_events_eff,t_elapsed/60,t_elapsed/n_events_eff);
fprintf('%s\n',msgstr);

% number of derived tracks per original one
event_frequency = 1/(orig_years*(vq_data.ens_size+1));
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
return
