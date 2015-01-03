function eq_data=eq_global_probabilistic(eq_data_in,ens_size,check_plot)
% EQ eq data probabilistic
% NAME:
%   eq_global_probabilistic
% PURPOSE:
%   Given an earthquake(EQ) database, create a probabilistic version by
%   wiggling location, depth and magnitudes of epicenters
%
%   previous step: see eq_isc_gem_read, eq_centennial_read, or eq_signigeq_read
%   next step: see eq_global_hazard_set
% CALLING SEQUENCE:
%   eq_data=eq_global_probabilistic(eq_data_in,check_plot)
% EXAMPLE:
%   eq_data=eq_global_probabilistic(eq_centennial_read,9,1)
% INPUTS:
%   eq_data_in: an EQ database, see eq_centennial_read
% OPTIONAL INPUT PARAMETERS:
%   ens_size: ensemble size, the number of 'copies' of the original
%       database we create, default=9 (means 10x as original data)
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
% MODIFICATION HISTORY:
% David N. Bresch, david.bresch@gmail.com, 20141010, initial
%-

eq_data=[]; % init

%global climada_global
if ~climada_init_vars,return;end % init/import global variables

%%if climada_global.verbose_mode,fprintf('*** %s ***\n',mfilename);end % show routine name on stdout

if ~exist('eq_data_in','var'),fprintf('eq_data_in required, STOPPED\n');end
if ~exist('ens_size','var'),ens_size=[];end
if ~exist('check_plot','var'),check_plot=0;end

% PARAMETERS
%
% whether we use the same seed each time we call this (=1) or a true new seed (=0)
% If =1, use always the same seed, in order to allow for reproduceability in exercises
force_seed_0=1; % default=1 (for reproduceability)
%
% the number of created epicenters per original epicenter
if isempty(ens_size),ens_size=9;end
%
% ensemble generation parameters
geo_amp = 1.0; % amplitude of max random location shift [degree]
dep_mul = 0.2; % depth multiplier (0.1 means depth of derived events will be in the range +/-10% of original)
mag_mul = 0.5; % magnitude multiplier (works analog to depth)

n_epicenters=length(eq_data_in.mag);

% generate random starting points for ensemble members
if force_seed_0,rand('seed',0),end % always the same seed, in order to allow for reproduceability in exercises
dlon =   geo_amp*2*(rand(1,ens_size*n_epicenters)-0.5); % rand: uniformly distributed random numbers
dlat =   geo_amp*2*(rand(1,ens_size*n_epicenters)-0.5);
ddep = 1+dep_mul*2*(rand(1,ens_size*n_epicenters)-0.5);
dmag = 1+mag_mul*2*(rand(1,ens_size*n_epicenters)-0.5);

eq_data.yyyy=eq_data_in.yyyy;
eq_data.mm=eq_data_in.mm;
eq_data.dd=eq_data_in.dd;
eq_data.hr=eq_data_in.hr;
eq_data.min=eq_data_in.min;
eq_data.sec=eq_data_in.sec;
eq_data.datenum=eq_data_in.datenum;
eq_data.n_epicenters_orig=n_epicenters;
eq_data.ens_size=ens_size;
% allocate memory
n_elements=fix((ens_size+1)*n_epicenters);
eq_data.glat=zeros(1,n_elements);
eq_data.glon=zeros(1,n_elements);
eq_data.dep =zeros(1,n_elements);
eq_data.mag =zeros(1,n_elements);
eq_data.orig_event_flag =zeros(1,n_elements); % will remain 0 for prob events
% copy historic
eq_data.glat(1:n_epicenters)=eq_data_in.glat;
eq_data.glon(1:n_epicenters)=eq_data_in.glon;
eq_data.dep(1:n_epicenters)=eq_data_in.dep;
eq_data.mag(1:n_epicenters)=eq_data_in.mag;
eq_data.orig_event_flag(1:n_epicenters)=eq_data_in.orig_event_flag;

fprintf('adding %i derived epicenters to each of the %i original epicenters\n',ens_size,n_epicenters);

% allow for distribution abovde minmal magnitude
min_mag=min(eq_data_in.mag);
eq_data_in.mag=eq_data_in.mag-min_mag;
eq_data.mag=eq_data.mag-min_mag;

for ens_i=1:ens_size
    eq_data.glat(ens_i*n_epicenters+1:(ens_i+1)*n_epicenters)=eq_data_in.glat+dlat((ens_i-1)*n_epicenters+1:ens_i*n_epicenters);
    eq_data.glon(ens_i*n_epicenters+1:(ens_i+1)*n_epicenters)=eq_data_in.glon+dlon((ens_i-1)*n_epicenters+1:ens_i*n_epicenters);
    eq_data.dep(ens_i*n_epicenters+1:(ens_i+1)*n_epicenters)=eq_data_in.dep.*ddep((ens_i-1)*n_epicenters+1:ens_i*n_epicenters);
    eq_data.mag(ens_i*n_epicenters+1:(ens_i+1)*n_epicenters)=eq_data_in.mag.*dmag((ens_i-1)*n_epicenters+1:ens_i*n_epicenters);
    % and replicate the time/date info
    eq_data.yyyy(ens_i*n_epicenters+1:(ens_i+1)*n_epicenters)=eq_data_in.yyyy;
    eq_data.mm(ens_i*n_epicenters+1:(ens_i+1)*n_epicenters)=eq_data_in.mm;
    eq_data.dd(ens_i*n_epicenters+1:(ens_i+1)*n_epicenters)=eq_data_in.dd;
    eq_data.hr(ens_i*n_epicenters+1:(ens_i+1)*n_epicenters)=eq_data_in.hr;
    eq_data.min(ens_i*n_epicenters+1:(ens_i+1)*n_epicenters)=eq_data_in.min;
    eq_data.sec(ens_i*n_epicenters+1:(ens_i+1)*n_epicenters)=eq_data_in.sec;
    eq_data.datenum(ens_i*n_epicenters+1:(ens_i+1)*n_epicenters)=eq_data_in.datenum;
end % ens_i

% set magnitude back
eq_data_in.mag=eq_data_in.mag+min_mag; % just for safety
eq_data.mag=eq_data.mag+min_mag;

try
    epicenters_file_prob_mat=strrep(eq_data_in.filename,'.txt','_prob.mat');
    fprintf('probabilistic set stored as %s\n',epicenters_file_prob_mat);
    save(epicenters_file_prob_mat,'eq_data');
catch
    fprintf('WARNING: failed to save probabilistic epicenters\n');
end

if check_plot
    
    % compare distributions
    
    figure('Name','Epicenter properties');
    subplot(2,2,1)
    X=-90:10:90; % the bin centers
    [nelements] = hist(eq_data_in.glat,X);
    [nelements_prob] = hist(eq_data.glat,X);
    bar(X,[nelements*(ens_size+1);nelements_prob]');legend('original','probabilistic'),title('latitude')
    hold off
    
    subplot(2,2,2)
    X=-180:20:180; % the bin centers
    [nelements] = hist(eq_data_in.glon,X);
    [nelements_prob] = hist(eq_data.glon,X);
    bar(X,[nelements*(ens_size+1);nelements_prob]');legend('original','probabilistic'),title('longitude')
    hold off
    
    subplot(2,2,3)
    X=0:40:800; % the bin centers
    [nelements] = hist(eq_data_in.dep,X);
    [nelements_prob] = hist(eq_data.dep,X);
    bar(X,[nelements*(ens_size+1);nelements_prob]');legend('original','probabilistic'),title('depth')
    hold off
    
    subplot(2,2,4)
    X=5:0.25:10; % the bin centers
    [nelements] = hist(eq_data_in.mag,X);
    [nelements_prob] = hist(eq_data.mag,X);
    bar(X,[nelements*(ens_size+1);nelements_prob]');legend('original','probabilistic'),title('magnitude')
    xlabel(sprintf('%i samples',ens_size));
    hold off
    
    set(gcf,'Color',[1 1 1]); % white background
    
    figure('Name','Epicenter locations (zoom in)');
    climada_plot_world_borders
    hold on
    plot(eq_data.glon,eq_data.glat,'or')
    plot(eq_data_in.glon,eq_data_in.glat,'.b')
    
end % check_plot

return
