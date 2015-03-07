function vq_data=vq_global_probabilistic(vq_data,ens_size,check_plot)
% generate probabilistic hazard set
% MODULE:
%   eq_global
% NAME:
%   vq_global_probabilistic
% PURPOSE:
%   Given a volcano (VQ) database, create a probabilistic version by
%   wiggling duration, cloud height, wind direction and speed...
%
%   First, add wind speed and direction (from climatolog) to historic
%       eruptions. For those historic eruptions with no known month, add mean
%       wind.
%   Second, generate 11 times the events and add the 'missing' months winds
%       For historic eruptions without wind, use January wind for historic
%   Third, generate probabilistic eruptions
%
%   previous step: see vq_volcano_list_read
%   next step: see vq_global_hazard_set
% CALLING SEQUENCE:
%   vq_data=vq_global_probabilistic(vq_volcano_list_read,ens_size,check_plot)
% EXAMPLE:
%   vq_data=vq_global_probabilistic(vq_volcano_list_read,9,1)
% INPUTS:
%   vq_data: an volcano (VQ) database, see vq_volcano_list_read
% OPTIONAL INPUT PARAMETERS:
%   ens_size: ensemble size, the number of 'copies' of the original
%       database we create, default=9 (means 10x as original data)
%       if =0, the data is only enriched by wind climatology and no
%       probabilistic eruptions are generated (for historic events with no
%       month indicated, the annual mean wind is used).
%       if =1, the data is 'replicated' 11 additional times, to reflect the
%       wind climatology each month.
%       if >1, it's really replicating other parameters, mainly Cloud_height_km
%       parameters are sampled ens_size-1 times, the final event count is
%       thus n_orig_events*12*(ens_size-1)
%   check_plot: show a check plot (=1), or not (=0, default)
% OUTPUTS:
%   vq_data, same as output of vq_volcano_list_read (see there), just for
%       (many) more eruptions, plus U_vel_kmh an U_vel_rad added
% MODIFICATION HISTORY:
% David N. Bresch, david.bresch@gmail.com, 20150307, initial
%-

%global climada_global
if ~climada_init_vars,return;end % init/import global variables

%%if climada_global.verbose_mode,fprintf('*** %s ***\n',mfilename);end % show routine name on stdout

if ~exist('vq_data','var'),fprintf('vq_data required, STOPPED\n');end
if ~exist('ens_size','var'),ens_size=[];end
if ~exist('check_plot','var'),check_plot=0;end

module_data_dir=[fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];

% PARAMETERS
%
% probabilistic scattering parameters (see also code below for details)
Cloud_height_amp=1.5; % i.e. max cloud can be 0..Cloud_height_amp times historic cloud height
duration_amp=2; % i.e. max duration can be 0..duration_amp times historic duration
%
% the pressure level we take the wind from
pressure_level=500; % in mb, like 1000 (surface) or 500 (about 5km)
%
% whether we use the same seed each time we call this (=1) or a true new seed (=0)
% If =1, use always the same seed, in order to allow for reproduceability in exercises
force_seed_0=1; % default=1 (for reproduceability)
%
% the number of created epicenters per original epicenter
if isempty(ens_size),ens_size=9;end
%
% the netCDF file with wind climatology (2.5x2.5 deg, 17 levels, monthly)
wind_clim_file=[module_data_dir filesep 'volcanoes' filesep 'NCEP_wind_data.nc'];


% get wind climatology (NCAR)
% ---------------------------
% as the tephra (ash thickness) calculation later needs wind speed and
% direction, we add it here


[fP,fN]=fileparts(wind_clim_file);
wind_clim_file_mat=[fP filesep fN '.mat'];
if ~climada_check_matfile(wind_clim_file,wind_clim_file_mat) % Check whether there is already a .mat file
    
    %FINFO = ncinfo(wind_clim_file);
    wind_clim.u=ncread(wind_clim_file,'u'); % u(lon,lat,P_level,month)
    wind_clim.v=ncread(wind_clim_file,'v');
    wind_clim.speed=ncread(wind_clim_file,'speed'); % m/s
    wind_clim.X=ncread(wind_clim_file,'X'); % lon
    wind_clim.Y=ncread(wind_clim_file,'Y'); % lat
    wind_clim.P=ncread(wind_clim_file,'P');
    %T=ncread(wind_clim_file,'T')
    
    fprintf('saving wind climatology in %s\n',wind_clim_file_mat);
    save(wind_clim_file_mat,'wind_clim');
else
    load(wind_clim_file_mat);
end

n_eruptions=length(vq_data.Year);

vq_data.U_vel_kmh=vq_data.lon*0; % init
vq_data.U_phi_rad=vq_data.lon*0; % init
wind_clim_pos.xi=vq_data.lon*0; % init
wind_clim_pos.yi=vq_data.lon*0; % init
wind_clim_pos.p_i=vq_data.lon*0; % init

% step 1: just add wind
% ---------------------

for eruption_i=1:n_eruptions
    
    % find position in wind climatology
    [~,xi]=min(abs(wind_clim.X-vq_data.lon(eruption_i)));
    [~,yi]=min(abs(wind_clim.Y-vq_data.lat(eruption_i)));
    [~,p_i]=min(abs(wind_clim.P-pressure_level));
    mi=vq_data.Month(eruption_i);
    
    if mi==0,mi=1:12;end % if month unknown, get all
    
    % obtain wind speed and direction climatology for month of eruption
    U_vel_kmh=squeeze(wind_clim.speed(xi,yi,p_i,mi))*3.6; % km/h
    % four quadrant arctangent in radion (-pi<=atan2(Y,X)<= pi, in degree: ./pi*180:
    U_phi_rad=squeeze(atan2(wind_clim.v(xi,yi,p_i,mi),wind_clim.u(xi,yi,p_i,mi)));
    
    vq_data.U_vel_kmh(eruption_i)=mean(U_vel_kmh); % mean of one scalar is faster then asking each time
    vq_data.U_phi_rad(eruption_i)=mean(U_phi_rad);
    
    % store indices (for speedup below)
    wind_clim_pos.xi(eruption_i)=xi;
    wind_clim_pos.yi(eruption_i)=yi;
    wind_clim_pos.p_i(eruption_i)=p_i;
    
end % eruption_i

if ens_size==0,return;end % only wind added

% step 2: generate one eruption a month (add 11 times, all except original month)
% -------------------------------------

n_orig_eruptions=n_eruptions; % store

% replicate all data 12 times
field_names=fieldnames(vq_data);
for field_i=1:length(field_names)
    if ~strcmp(char(field_names{field_i}),'filename')
        vq_data.(field_names{field_i})=repmat(vq_data.(field_names{field_i}),12,1);
    end
end % field_i
vq_data.ens_size=12-1;

vq_data.orig_event_flag(n_orig_eruptions+1:end)=0;
    
for eruption_i=1:n_orig_eruptions % still loop over all original eruptions
    
    xi=wind_clim_pos.xi(eruption_i);
    yi=wind_clim_pos.yi(eruption_i);
    p_i=wind_clim_pos.p_i(eruption_i);
    % obtain wind speed and direction climatology for all months
    U_vel_kmh=squeeze(wind_clim.speed(xi,yi,p_i,:))*3.6; % km/h
    % four quadrant arctangent in radion (-pi<=atan2(Y,X)<= pi, in degree: ./pi*180:
    U_phi_rad=squeeze(atan2(wind_clim.v(xi,yi,p_i,:),wind_clim.u(xi,yi,p_i,:)));
    
    if vq_data.Month(eruption_i)==0
        % we replace the mean with January
        % hence the historic eruption has January wind
        vq_data.U_vel_kmh(eruption_i)=U_vel_kmh(1); % mean of one scalar is faster then asking each time
        vq_data.U_phi_rad(eruption_i)=U_phi_rad(1);
        mi=2:12;
        U_vel_kmh=U_vel_kmh(mi);
        U_phi_rad=U_phi_rad(mi);
    else
        % historic eruption known by month, all other months
        % this way, the historic eruption has the correct monthyl wind
        mi=1:12;
        U_vel_kmh=U_vel_kmh(mi~=vq_data.Month(eruption_i));
        U_phi_rad=U_phi_rad(mi~=vq_data.Month(eruption_i));
        mi=mi(mi~=vq_data.Month(eruption_i)); % also reduce mi
    end
    target_i=eruption_i+(1:length(mi))*n_orig_eruptions;
    
    %fprintf('%i %i %i %i\n',eruption_i,n_eruptions,max(target_i),length(mi));
        
    vq_data.U_vel_kmh(target_i)=U_vel_kmh;
    vq_data.U_phi_rad(target_i)=U_phi_rad;
    
end % eruption_i

if ens_size==1;return;end % we only added months

% step 3: generate random cloud height etc.
% -----------------------------------------

n_month_eruptions=length(vq_data.lon); % store

% replicate all data ens_size-1 times
field_names=fieldnames(vq_data);
for field_i=1:length(field_names)
    if ~strcmp(char(field_names{field_i}),'filename')
        vq_data.(field_names{field_i})=repmat(vq_data.(field_names{field_i}),ens_size,1);
    end
end % field_i
vq_data.ens_size=length(vq_data.lon)/n_orig_eruptions-1;

vq_data.orig_event_flag(n_orig_eruptions+1:end)=0;

% generate random starting points for ensemble members
if force_seed_0,rand('seed',0),end % always the same seed, in order to allow for reproduceability in exercises

% Note: rand generates uniformly distributed random numbers ]0..1[
% -1..1: 2*(rand(1,(ens_size-1)*n_month_eruptions)-0.5);
% 0..1:  rand(1,(ens_size-1)*n_month_eruptions);
    
amp=Cloud_height_amp*rand(1,(ens_size-1)*n_month_eruptions)';
vq_data.Cloud_height_km(n_month_eruptions+1:end)=...
    amp.*vq_data.Cloud_height_km(n_month_eruptions+1:end);

amp=duration_amp*rand(1,(ens_size-1)*n_month_eruptions)';
vq_data.duration_h(n_month_eruptions+1:end)=...
    amp.*vq_data.duration_h(n_month_eruptions+1:end);

fprintf('adding %i derived eruptions to each of the %i original eruptions -> %i events\n',...
    12*ens_size,n_orig_eruptions,length(vq_data.lon));

if check_plot
    
    % compare distributions
    
    figure('Name','Eruption properties');
    subplot(2,2,1)
    X=0:2:ceil(max(vq_data.Cloud_height_km)); % the bin centers
    [nelements]      = hist(vq_data.Cloud_height_km(1:n_orig_eruptions),X);
    [nelements_prob] = hist(vq_data.Cloud_height_km,X);
    bar(X,[nelements*12*ens_size;nelements_prob]');legend('original','probabilistic')
    title(strrep('Cloud_height_km','_',' '))
    hold off
    
    subplot(2,2,2)
    X=0:0.25:ceil(max(vq_data.duration_h)); % the bin centers
    [nelements]      = hist(vq_data.duration_h(1:n_orig_eruptions),X);
    [nelements_prob] = hist(vq_data.duration_h,X);
    bar(X,[nelements*12*ens_size;nelements_prob]');legend('original','probabilistic')
    title(strrep('duration_h','_',' '))
    hold off
    
    subplot(2,2,3)
    X=0:1:ceil(max(vq_data.U_vel_kmh)); % the bin centers
    [nelements]      = hist(vq_data.U_vel_kmh(1:n_orig_eruptions),X);
    [nelements_prob] = hist(vq_data.U_vel_kmh,X);
    bar(X,[nelements*12*ens_size;nelements_prob]');legend('original','probabilistic')
    title(strrep('U_vel_kmh','_',' '))
    hold off
    
    subplot(2,2,4)
    X=-180:30:180; % the bin centers
    [nelements]      = hist(vq_data.U_phi_rad(1:n_orig_eruptions)/pi*180,X);
    [nelements_prob] = hist(vq_data.U_phi_rad/pi*180,X);
    bar(X,[nelements*12*ens_size;nelements_prob]');legend('original','probabilistic')
    title(strrep('U_phi_rad','_',' '))
    hold off
    
    set(gcf,'Color',[1 1 1]); % white background
    
end % check_plot

return
