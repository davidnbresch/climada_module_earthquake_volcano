function plot_attenuation_parameters(magnitude, attenuation_file)
% MODULE:
% eq_global
% NAME:
%   plot_attenuation_parameters
% PURPOSE:
%   plot the attenuation functions defined by the sets of parameters in
%   climada_module_eq_global/data/system/attenuation_parameters.xlsx, for a
%   given magnitude
% CALLING SEQUENCE:
%   plot_attenuation_parameters(magnitude, filename)
% EXAMPLE:
%   plot_attenuation_parameters(8)
%   plot_attenuation_parameters
% INPUT:
%   magnitude: Richter magnitude at epicenter. If empty, 8 is used as
%   default value
%   attenuation_file: Name of the file that contains the attenuation 
%   parameters. 
%   If empty, the code tries a default name, if it does not exist, it
%   prompts the user to locate the file
% OUTPUTS:
%   A plot showing the attenuation of intensity with distance for each set
%   of parameters and for the magnitude specified in the input
% COMMENT:
%   The parameters a1, a2, a3, and a4 describe an attenuantion of the type
%   MMI(dist) = a1 + a2 * mag - a3 * log(dist) - a4 * dist, where
%       mag is the Richter magnitude at the epicenter
%       dist is the distance from the epicenter (in km)
%       MMI is the Mercalli Modified Intensity at distance dist from the
%       epicenter
%
% MODIFICATION HISTORY:
% Melanie Bieli, melanie.bieli@bluewin.ch, 20141209, initial
%-

%global climada_global
if ~climada_init_vars,return;end % init/import global variables

if ~exist('magnitude','var'),magnitude=8;end
if ~exist('attenuation_file','var'),attenuation_file='';end

% set default value for attenuation_file if not given
eq_dir=[fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];
attenuation_file_default=[eq_dir filesep 'system' filesep 'attenuation_parameters.xlsx'];
if isempty(attenuation_file),attenuation_file=attenuation_file_default;end

if ~exist(attenuation_file_default,'file')
    [filename, pathname] = uigetfile(fullfile([eq_dir filesep 'system']), 'Open file containing the attenuation parameters:');
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        attenuation_file=fullfile(pathname,filename);
    end
end

% read attenuation parameters from .xlsx file
attenuation_data = climada_xlsread('no',attenuation_file);
a1_vector = attenuation_data.A1;
a2_vector = attenuation_data.A2;
a3_vector = attenuation_data.A3;
a4_vector = attenuation_data.A4;

% set colormap 
number_of_parameter_sets = min([length(a1_vector),length(a2_vector),length(a3_vector),length(a4_vector)]);
mycolormap = colormap(hsv(number_of_parameter_sets+3));

% calculate MMI up to a distance of 200 km from the epicenter, for each set
% of parameters a1,a2,a3,a4
% MMI does not exceed a maximum value defined by the equation 
% I_0 = 1.5*(mag - 1), where I_0 is the epicentral intensity (in MMI)
% and mag the Richter magnitude of the earthquake
% Source: Y-X. Hu, S-C. Liu, W. Dong: Earthquake Engineering,

distance=1:200;
MMI = zeros(1,200);
maximum_MMI = 1.5*(magnitude-1);
for parameter_set_i = 1:number_of_parameter_sets
    for dist = 1:200
        MMI(dist) = a1_vector(parameter_set_i) + a2_vector(parameter_set_i)* magnitude - a3_vector(parameter_set_i) * log(dist) - a4_vector(parameter_set_i) * dist;
        if MMI(dist) > maximum_MMI, MMI(dist) = maximum_MMI; end
    end
    % Settings to plot the curve used for the 'World average' as a bold
    % black line
    scope_string = char(attenuation_data.geographical_scope(parameter_set_i));
    if strcmp(scope_string,'World1')
      plot(distance,MMI,'k','Linewidth',4);
    else plot(distance,MMI,'Color',mycolormap(parameter_set_i,:),'Linewidth',1.5);
    end
    
    if parameter_set_i < number_of_parameter_sets
        hold on;
    else
        hold off;
    end
end

attenuation_legend = attenuation_data.geographical_scope(1:number_of_parameter_sets);       
xlabel('Distance from epicenter (km)');
ylabel('Intensity (Modified Mercalli)');
axis([0,200,0,13])
title('Attenuation functions');
legend(attenuation_legend);
        