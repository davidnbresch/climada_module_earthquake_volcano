function T=vq_tephra_field_cm(centroids,lon,lat,H,U_vel,U_phi,tau,check_plot)
% Create volcanoe hazard set
% MODULE
% eq_global
% NAME:
%   vq_tephra_field_cm
% PURPOSE:
%   Calculate the tephra (ash) field thickness for one eruption at
%   centroids. The eroption is characterized by (obviously) the coordinate
%   of the eruption center, the eruptive column height, the prevailing wind
%   and the duration of the eruption.
%
%   according to A. O. Gonzalez-Mellado and S. De la Cruz-Reyna, 2010: A
%   simple semi-empirical approach to model thickness of ash-deposits for
%   different eruption scenarios. Nat. Hazards Earth Syst. Sci., 10,
%   2241?2257, 2010, direct:
%   www.nat-hazards-earth-syst-sci.net/10/2241/2010/nhess-10-2241-2010.pdf
%
%   Variable parameters are distance r between centroid(i) and eruption
%   center and phi, angle between the wind direction and the vector from
%   the ash emission center to centroid(i).
%
%   Some PARAMETERS set in code, see below, also defaults
%
%   previous step: see vq_global_hazard_set (usually colled from there)
% CALLING SEQUENCE:
%   T=vq_tephra_field_cm(centroids,lon,lat,H,U_vel,U_phi,tau)
% EXAMPLE:
%   lon=14.426;lat=40.821; % Vesuvius
%   centroids.lon=lon-1:.01:lon+1;centroids.lat=centroids.lon*0+lat; % 1D cross-section
%   T=vq_tephra_field_cm(centroids,lon,lat,15,50,0,1);
%   plot((centroids.lon-lon)*100,T);hold on;plot([0 0],[0 max(T)],'-r');
%   xlabel('approx km');ylabel('ash thickness [cm]')
%
%   [X,Y]=meshgrid(lon-1:.01:lon+1,lat-1:.01:lat+1); % 2D grid
%   centroids.lon=reshape(X,numel(X),1);centroids.lat=reshape(Y,numel(X),1); % 2D grid
%   T=vq_tephra_field_cm(centroids,lon,lat,15,50,0,1,1); % last '1' for plot
%
%   % As in Fig. 11 of cited paper: 30 June 1997 Popocatepetl eruption
%   T=vq_tephra_field_cm(centroids,lon,lat,8,80,312/360*pi,0.58,1); % last '1' for plot
%
% INPUTS:
%   centroids: a structure with lon(i) and lat(i) of all the points we'd
%       like to obtain the tephra field thickness for
%   lon: the longitude (scalar) of the eruption center
%   lat: the latitude of the eruption center
% OPTIONAL INPUT PARAMETERS:
%   H: eruptive column height in km
%   U_vel: wind velocity in km/h, indicative 50-100 km/h most often.
%   U_phi: direction of wind in radian (-pi..pi), counterclockwise from
%       (dx=1,dy=0), i.e. from direction East, West=pi, North=pi/2
%   tau: duration of the high-intensity phase of the eruption in hours
%       (e.g. Pinatubo 1-5h) default=8h
%   check_plot: =1, show check plot, =0 not (default)
% OUTPUTS:
%   T: the height of the tephra field (in cm) at each centroid
% MODIFICATION HISTORY:
% David N. Bresch, david.bresch@gmail.com, 20150302, initial
%-

T=[]; % init (in case of return)

if ~exist('centroids','var'),return;end
if ~exist('lon','var'),return;end
if ~exist('lat','var'),return;end
if ~exist('H','var'),    H=14;end
if ~exist('U_vel','var'),U_vel=50;end
if ~exist('U_phi','var'),U_phi=0;end
if ~exist('tau','var'),tau=8;end
if ~exist('check_plot','var'),check_plot=0;end

T=centroids.lon*0; % init (with zero)

% PARAMETERS
%
rho=1100; % density of the falling material (density of the uncompressed fresh deposit)
%           (kgm-3, indicative value 1100)
Htropopause=15.5; % tropopause height in km
alpha_param=2.535-0.051*H; % the rate at which the deposit thickness decays with distance

if H<Htropopause
    D=-4.189*H+114.407;
else % above)
    D=52.822*H-770.17;
end

% distance to eruption center im km
dx=(centroids.lon-lon)*cos(lat/180*pi);
dy=centroids.lat-lat;

r2=dx.^2+dy.^2; % in degree^2

% find the centroids close enough, speeds up calculation
valid_centroid_pos=find(r2<2); % means within order of sqrt(2) degrees, ~140km

% see further below for explicit form

r2=r2(valid_centroid_pos);
dx=dx(valid_centroid_pos);
dy=dy(valid_centroid_pos);

r=sqrt(r2)*111.12; % km
r(r<1)=1; % area 1km around crater the same
phi0=atan2(dy,dx); % four quadrant arctangent in radion (-pi<=atan2(Y,X)<= pi, in degree: ./pi*180)
% measured counterclockwise from dx=1,dy0, i.e. atan2(0,1)=0, atan2(1,1)=pi/4, atan2(1,0)=pi/2

phi=phi0-U_phi;

T(valid_centroid_pos)=14.28*H^4*tau^(alpha_param/2) ./ ...
    ( (2*D)^(1-alpha_param/2)*rho ) * ...
    exp(-U_vel/(2*D).*r.*(1-cos(phi))).*r.^(-alpha_param);

% below in explicit form
%
% centroid_count=length(valid_centroid_pos);
%
% for centroid_ii=1:centroid_count % now loop over all valid centroids
%
%     centroid_i=valid_centroid_pos(centroid_ii);
%
%     r=sqrt(r2(centroid_i))*111.12; % km
%     if r<1;r=1;end % area 1km around crater the same
%
%     phi0=atan2(dy(centroid_i),dx(centroid_i)); % four quadrant arctangent in radion (-pi<=atan2(Y,X)<= pi, in degree: ./pi*180)
%     % measured counterclockwise from dx=1,dy0, i.e. atan2(0,1)=0, atan2(1,1)=pi/4, atan2(1,0)=pi/2
%
%     phi=phi0-U_phi;
%
%     T(centroid_i)=14.28*H^4*tau^(alpha_param/2) / ...
%         ( (2*D)^(1-alpha_param/2)*rho ) * ...
%         exp(-U_vel/(2*D)*r*(1-cos(phi)))*r^(-alpha_param);
%
% end % centroid_i

if check_plot
    
    % the 2D - plot
    subplot(2,2,1)
    %climada_circle_plot(values,lon,lat,title_str,circle_diam,circle_format,marker_size,marker_format,overlay_plot,axis_range,max_value)
    climada_circle_plot(T,centroids.lon,centroids.lat,'',50,'',1,'',0);
    hold on
    plot(lon,lat,'og');plot(lon,lat,'Xg');
    axis([min(centroids.lon) max(centroids.lon) min(centroids.lat)-1 max(centroids.lat)]); % set axis for good zoom
    
    % and the cross-sections
    subplot(2,2,3)
    mlat=min(abs(centroids.lat-lat));
    pos=find(centroids.lat-lat==mlat);
    plot((centroids.lon(pos)-lon)*111.12*cos(lat/180*pi),T(pos));hold on;plot([0 0],[0 max(T(pos))],'-r');
    xlabel('longitude, approx km');ylabel('ash thickness [cm]')
    subplot(2,2,2)
    mlon=min(abs(centroids.lon-lon));
    pos=find(centroids.lon-lon==mlon);
    plot(T(pos),(centroids.lat(pos)-lat)*111.12);hold on;plot([0 max(T(pos))],[0 0],'-r');
    ylabel('latitude, approx km');xlabel('ash thickness [cm]')
    
    % histogram
    subplot(2,2,4)
    
    hist(T(T>eps));xlabel('ash thickness [cm]');ylabel('# centroids')
    
end

end % vq_tephra_field