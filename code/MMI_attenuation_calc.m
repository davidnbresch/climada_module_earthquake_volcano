function MMI = MMI_attenuation_calc(mag, dist, dep, correction, a1, a2, a3, a4)
% MODULE:
% eq_global
% NAME:
%   simple_eq_MMI
% PURPOSE:
%   calculate MMI at a given distance from the epicenter 
% CALLING SEQUENCE:
%   MMI=MMI_attenuation_calc(mag, dist, correction,a1, a2, a3, a4)
% EXAMPLE:
%   MMI=MMI_attenuation_calc(7.44,10.0,1.67,1.67,1.3,0.0026);
%   MMI=MMI_attenuation_calc(8,15.5)   % using global average values for 
%   correction and the attenuation parameters a1, a2, a3, a4
% INPUTS:
%   mag:    magnitude
%   dist:   distance to epicentre [in km] 
% OPTIONAL INPUTS:
%   correction, a1,a2,a3,a4: parameters defining an attenuation function
%   of the type MMI = a1 + a2 * mag - a3 * log(dist+correction) - a4 * dist
%   Have a look at eq_global-master/data/system/attenuation_parameters.xlsx 
%   for a compilation of parameters for selected geographic regions
% OUTPUT:
%   MMI:    the MMI at distance dist from the epicenter (Modified Mercalli
%   Intensity)
% MODIFICATION HISTORY:
% Melanie Bieli, melanie.bieli@bluewin.ch, 20141116
% Melanie Bieli, melanie.bieli@bluewin.ch, 20141209, added correction parameter
% Melanie Bieli, melanie.bieli@bluewin.ch, 20151018, minor clean-up, renaming

%% default values for attenuation parameters: dep, correction, a1, a2, a3, a4
if ~exist('dep','var') || isempty(dep), dep = 0; end
if ~exist('correction','var') || isempty(correction), correction = 0; end
if ~exist('a1','var') || isempty(a1), a1 = 1.7; end
if ~exist('a2','var') || isempty(a2), a2 = 1.5; end
if ~exist('a3','var') || isempty(a3), a3 = 1.1726; end
if ~exist('a4','var') || isempty(a4), a4 = 0.00106; end

%% compute the attenuation equation 
% MMI does not exceed a maximum value defined by the equation 
% I_0 = 1.5*(mag - 1), where I_0 is the epicentral intensity (in MMI)
% and mag the Richter magnitude of the earthquake
% Source: Y-X. Hu, S-C. Liu, W. Dong: Earthquake Engineering

maximum_MMI = 1.5*(mag-1);
MMI = a1 + a2 * mag - a3 * log(dist+correction) - a4 * dist;
if MMI > maximum_MMI, MMI = maximum_MMI; end

end




        
 
 

