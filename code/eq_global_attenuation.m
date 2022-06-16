function intensity_at_centroids = eq_global_attenuation(glat,glon,mag,centroids,check_plot,dep,correction,a1,a2,a3,a4,b)
% Calculate the attenuation of an earthquake's intensity with
% increasing distance from the epicenter
% MODULE:
% eq_global
% NAME:
%   eq_global_attenuation
% PURPOSE:
%   given an earthquake (EQ) epicenter, defined by glat,glon,dep,mag,
%   calculate the resulting intensity at a set of centroids, using an attenuation
%   function. Distance is on sphere iro correction for latitude, but no
%   great circle distance (for speedup reasons, other errors are bigger).
%
%   In order to speed up the calculation, no consistency checks are done,
%   code assumes proper glat,glon,dep,mag and reasonable centroids input
%
%   For basic theory, see e.g. https://en.wikipedia.org/wiki/Peak_ground_acceleration
%
%   To test attenuation functions, see eq_global_attenuation_TEST
%   then edit between *********** lines beow
%
%   mainly called from: eq_global_hazard_set
% CALLING SEQUENCE:
%   intensity_at_centroids = eq_global_attenuation(glat,glon,dep,mag,centroids,check_plot,a1,a2,a3,a4)
% EXAMPLE:
%   for an earthquake in Vancouver, Canada:
%   glat=49.25;glon=-123.10;
%   intensity_at_centroids = eq_global_attenuation(glat,glon,5,6,centroids,0,1,2.4,1.5,1.13,0.0063)
%   The parameters a1=2.4, a2=1.5, a3=1.13, a4=0.0063 describe the
%   attenuation function for Western Canada; see
%   eq_global-master/data/system/attenuation_parameters.xlsx
% INPUTS:
%   glat: latitude of epicenter [in degrees]
%   glon: longitude of epicenter [in degrees]
%   dep: depth [km] of epicenter
%   mag: magnitude (Richter) of epicenter
%   centroids: a structure with the centroids information
%       centroids.lat: the latitude of the centroids
%       centroids.lon: the longitude of the centroids
% OPTIONAL INPUT PARAMETERS:
%   correction,a1,a2,a3,a4,b: parameters defining the attenuation function. 
%   See eq_global-master/data/system/attenuation_parameters.xlsx to use
%   parameters for specific regions; otherwise default values representing
%   a "global average attenuation function" will be used
%   check_plot: =1, if a check plot shall be drawn (default=0)
% OUTPUTS:
%   intensity_at_centroids: the MMI at centroids
% RESTRICTIONS:
%   code does not quality or even consistency checks, optimized for speed.
% MODIFICATION HISTORY:
% David N. Bresch, david.bresch@gmail.com, 20141013
% Melanie Bieli, melanie.bieli@bluewin.ch, 20141106
%-

% initialize output
intensity_at_centroids = []; 

% global climada_global
% if ~climada_init_vars,return;end % init/import global variables

% default values for attenuation parameters: dep, correction, a1, a2, a3, a4
if ~exist('dep','var') || isempty(dep), dep = 0; end
if ~exist('correction','var') || isempty(correction), correction = 0; end
if ~exist('a1','var') || isempty(a1), a1 = 1.7; end
if ~exist('a2','var') || isempty(a2), a2 = 1.5; end
if ~exist('a3','var') || isempty(a3), a3 = 1.1726; end
if ~exist('a4','var') || isempty(a4), a4 = 0.00106; end

if ~exist('glat','var'),return; end
if ~exist('glon','var'),return; end
if ~exist('dep','var'),return; end
if ~exist('mag','var'),return; end
if ~exist('centroids','var'),return; end
if ~exist('check_plot','var'), check_plot= 0; end

intensity_at_centroids = centroids.lon*0; % init output

% PARAMETERS
%
% in order to speed up, we do not calculate intensity too far away from epicenter
% note that 1 deg (at equator) is approx 111.12km
max_centroids_dist = 150; % max distance to epicenter we calculate intensity (in km)

% Great circle distance of centroids from epicenter using the spherical law 
% of cosines (see http://en.wikipedia.org/wiki/Great-circle_distance)
% Calculates the radial distance between the epicenter and all centroids
R_mean_Earth = 6371;    % the Earth's mean radius (in km)
glat_rad = deg2rad(glat);
glon_rad = deg2rad(glon);
raddist = acos(sin(glat_rad)*sin(deg2rad(centroids.lat)) + cos(glat_rad)*cos(deg2rad(centroids.lat)).*cos(deg2rad(centroids.lon)-glon_rad));
if raddist < 0 
   raddist = raddist + pi;
end
centroids_dist = R_mean_Earth * raddist;    

% Only look at centroids within the defined maximum distance from epicenter
eff_centroids=find(centroids_dist<(max_centroids_dist)); 
    
for centroid_ii=1:length(eff_centroids) %
    %now loop over all centroids close enough
    
    centroid_i=eff_centroids(centroid_ii); % index in centroids
    
    R = centroids_dist(centroid_i); % epicentral distance [km]

    % ********************************************************************
    % calling the function MMI_attenuation_calc to calculate the intensity at the
    % centroids
    intensity_at_centroids(centroid_i) = MMI_attenuation_calc(mag, R, dep, correction, a1, a2, a3, a4, b);    
    % ********************************************************************
    
end % centroid_ii

if check_plot
    
    fprintf('max MMI: %f\n',max(intensity_at_centroids));
    fprintf('preparing footprint plot\n')
    % create gridded values
    [X, Y, gridded_VALUE] = climada_gridded_VALUE(intensity_at_centroids,centroids);
    %fprintf('gridded values %f\n',gridded_VALUE);
    gridded_max = max(max(gridded_VALUE));
    gridded_max_round = gridded_max; % 90
    contourf(X, Y, full(gridded_VALUE),...
        2:gridded_max_round,'edgecolor','none')
    hold on
    climada_plot_world_borders(0.7)
    
    plot(centroids.lon, centroids.lat, '+r','MarkerSize',0.8,'linewidth',0.1)
    
    axis equal
    axis([min(centroids.lon) max(centroids.lon) ...
        min(centroids.lat)  max(centroids.lat)]);
    caxis([2 gridded_max_round])
    colorbar
end

return

