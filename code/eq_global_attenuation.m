function intensity_at_centroids = eq_global_attenuation(glat,glon,dep,mag,centroids,check_plot,a1,a2,a3,a4)
% eq attenuation calculation footprint
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
%   intensity_at_centroids = eq_global_attenuation(glat,glon,5,6,centroids,1, 2.4,1.5,1.13,0.0063)
%   The parameters a1=2.4, a2=1.5, a3=1.13, a4=0.0063 describe the
%   attenuation function for Western Canada; see
%   eq_global-master/data/system/attenuation_parameters.xlsx
% INPUTS:
%   glat: latitude of epicenter [in degrees]
%   glon: longitude of epicenter [in degrees]
%   dep: depth [km] of epicenter
%   mag: magnitude (Richter) of epicenter
%   centroids: a structure with the centroids information
%       centroids.Latitude: the latitude of the centroids
%       centroids.Longitude: the longitude of the centroids
% OPTIONAL INPUT PARAMETERS:
%   a1,a2,a3,a4: parameters defining the attenuation function. See
%   eq_global-master/data/system/attenuation_parameters.xlsx to use
%   parameters for specific regions; otherwise default values representing
%   a "global average attenuation function" will be used
%   check_plot: =1, if a check plot shall be drawn (default=0)
% OUTPUTS:
%   intensity_at_centroids: the MMI at centroids
% RESTRICTIONS:
%   code does not quality or even consistency checks, optimized for speed.
% MODIFICATION HISTORY:
% Melanie Bieli, melanie.bieli@bluewin.ch, 20141106
% David N. Bresch, david.bresch@gmail.com, 20141013
%-

intensity_at_centroids = []; % init output

%% default values for attenuation parameters a1, a2, a3, a4
if ~exist('a1','var') || isempty(a1), a1 = 1.67;       end
if ~exist('a2','var') || isempty(a2), a2 = 1.67;       end
if ~exist('a3','var') || isempty(a3), a3 = 1.3;        end
if ~exist('a4','var') || isempty(a4), a4 = 0.0026;     end

%global climada_global % currently not used
if ~climada_init_vars, return; end

if ~exist('glat','var'),return; end
if ~exist('glon','var'),return; end
if ~exist('dep','var'),return; end
if ~exist('mag','var'),return; end
if ~exist('centroids','var'),return; end
if ~exist('check_plot','var'), check_plot= 0; end

intensity_at_centroids = centroids.Longitude*0; % init output

% PARAMETERS
%
% in order to speed up, we do not calculate intensity too far away from epicenter
% note that 1 deg (at equator) is approx 111.12km
max_centroids_dist = 400; % max distance to epicenter we calculate intensity (in km)

% Great circle distance of centroids from epicenter using the spherical law 
% of cosines (see http://en.wikipedia.org/wiki/Great-circle_distance)
% Calculates the radial distance between the epicenter and all centroids
R_mean_Earth = 6371;    % the Earth's mean radius (in km)
glat_rad = degtorad(glat);
glon_rad = degtorad(glon);
raddist = acos(sin(glat_rad)*sin(degtorad(centroids.Latitude)) + cos(glat_rad)*cos(degtorad(centroids.Latitude)).*cos(degtorad(centroids.Longitude)-glon_rad));
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
    % calling the function simple_eq_MMI to calculate the intensity at the
    % centroids
    intensity_at_centroids(centroid_i) = simple_eq_MMI(mag, R, a1, a2, a3, a4);    
    % ********************************************************************
    
end % centroid_ii

if check_plot
    
    fprintf('max MMI: %f\n',max(intensity_at_centroids));
    fprintf('preparing footprint plot\n')
    % create gridded values
    [X, Y, gridded_VALUE] = climada_gridded_VALUE(intensity_at_centroids,centroids);
    fprintf('gridded values %f\n',gridded_VALUE);
    gridded_max = max(max(gridded_VALUE));
    gridded_max_round = gridded_max; % 90
    contourf(X, Y, full(gridded_VALUE),...
        2:gridded_max_round,'edgecolor','none')
    hold on
    climada_plot_world_borders(0.7)
    
    plot(centroids.Longitude, centroids.Latitude, '+r','MarkerSize',0.8,'linewidth',0.1)
    
    axis equal
    axis([min(centroids.Longitude) max(centroids.Longitude) ...
        min(centroids.Latitude)  max(centroids.Latitude)]);
    caxis([2 gridded_max_round])
    colorbar
end

return

