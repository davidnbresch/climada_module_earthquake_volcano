function intensity_at_centroids = eq_global_attenuation(glat,glon,dep,mag,centroids,check_plot)
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
%   mainly called from: eq_hazard_set
% CALLING SEQUENCE:
%   intensity_at_centroids = eq_global_attenuation(glat,glon,dep,mag,centroids,check_plot)
% EXAMPLE:
%   centroids.Longitude=10;centroids.Latitude=45;
%   intensity_at_centroids = eq_global_attenuation(24,10,5,6,centroids,1)
% INPUTS:
%   glat: latitude of epicenter [in degrees]
%   glon: longitude of epicenter [in degrees]
%   dep: depth [km] of epicenter
%   mag: magnitude (Richter) of epicenter
%   centroids: a structure with the centroids information
%       centroids.Latitude: the latitude of the centroids
%       centroids.Longitude: the longitude of the centroids
% OPTIONAL INPUT PARAMETERS:
%   check_plot: =1, if a check plot shall be drawn (default=0)
% OUTPUTS:
%   intensity_at_centroids: the MMI at centroids
% RESTRICTIONS:
%   code does not quality or even consistency checks, optimized for speed.
% MODIFICATION HISTORY:
% David N. Bresch, david.bresch@gmail.com, 20141013
%-

intensity_at_centroids = []; % init output

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
max_centroids_dist=2; % max disatance in degrees to epicenter we calculate intensity
%
% calculate square of distance of all centroids from epicenter
centroids_dist2=((glon-centroids.Longitude).*cos(glat/180*pi)).^2+(glat-centroids.Latitude).^2; % in deg^2

eff_centroids=find(centroids_dist2<(max_centroids_dist^2)); % still in deg^2
centroids_dist2(centroids_dist2<=1)=1; % all within 1 km shows center intensity

for centroid_ii=1:length(eff_centroids) % now loop over all centroids close enough
    
    centroid_i=eff_centroids(centroid_ii); % index in centroids
    
    R = centroids_dist2(centroid_i); % hypocentral distance (km), here still square of
    
    if R>1 % only take square root if >1 for speedup
        R = sqrt(R)*111.12; % now in km
    end
    
    % ********************************************************************
    % Po-Shen Lin and Chyi-Tyi Lee, 2008: Ground-Motion Attenuation
    % Relationships for Subduction-Zone Earthquakes in Northeastern Taiwan,
    % Bulletin of the Seismological Society of America, Vol. 98, No. 1, pp.
    % 220?240, February 2008, doi: 10.1785/0120060002
    % formulas 13 and 14 (p 226) and table 3/4 (p 227/228)
    % Zt= 0 for interface earthquakes, and Zt=1 for intraslab earthquakes
    M=mag;H=dep; % common names in paper, Zt ignored
    lnPGA_soil = -0.900 +1.000*M -1.900*log(R +0.99178*exp(0.52632*M)) +0.0040*H; % +0.310*Zt; % soil sites
    
    intensity_at_centroids(centroid_i) = exp(lnPGA_soil);
    % ********************************************************************
    
end % centroid_ii

if check_plot
    
    fprintf('max MMI: %f\n',max(intensity_at_centroids));
    fprintf('preparing footprint plot\n')
    
    % create gridded values
    [X, Y, gridded_VALUE] = climada_gridded_VALUE(intensity_at_centroids,centroids);
    gridded_max = max(max(gridded_VALUE));
    gridded_max_round = gridded_max; % 90
    
    contourf(X, Y, full(gridded_VALUE),...
        0:10:gridded_max_round,'edgecolor','none')
    hold on
    climada_plot_world_borders(0.7)
    
    plot(centroids.Longitude, centroids.Latitude, '+r','MarkerSize',0.8,'linewidth',0.1)
    
    axis equal
    axis([min(centroids.Longitude) max(centroids.Longitude) ...
        min(centroids.Latitude)  max(centroids.Latitude)]);
    caxis([0 gridded_max_round])
    colorbar
end

return

