function MMI = simple_eq_MMI(mag, dist, a1, a2, a3, a4)
% MMI = simple_eq_MMI(mag, dist, a1, a2, a3, a4)
% NAME:
%   simple_eq_MMI
% PURPOSE:
%   calculate MMI at a given distance from the epicenter
%
% CALLING SEQUENCE:
%   MMI=simple_eq_MMI(mag, dist, a1, a2, a3, a4)
% EXAMPLE:
%   MMI=simple_eq_MMI(7.44,10.0,1.67,1.67,1.3,0.0026);
%   MMI=simple_eq_MMI(8,15.5)   % using global average values for the
%   attenuation parameters a1, a2, a3, a4
% INPUTS:
%   mag:    magnitude
%   dist:   distance to epicentre [in km] 
%   a1, a2, a3, a4: attenuation parameters defining an attenuation function
%   of the type MMI = a1 + a2 * mag - a3 * log(dist) - a4 * dist
%   Have a look at eq_global-master/data/system/attenuation_parameters.xlsx 
%   for a compilation of parameters for selected geographic regions
% OUTPUT:
%   MMI:    the MMI at distance dist from the epicenter (Modified Mercalli
%   Intensity)
%
% MODIFICATION HISTORY:
% Melanie Bieli, melanie.bieli@bluewin.ch, 20141116

%% default values for attenuation parameters a1, a2, a3, a4
if ~exist('a1','var') || isempty(a1), a1 = 1.67;       end
if ~exist('a2','var') || isempty(a2), a2 = 1.67;       end
if ~exist('a3','var') || isempty(a3), a3 = 1.3;        end
if ~exist('a4','var') || isempty(a4), a4 = 0.0026;     end

%% compute the attenuation equation 
MMI = a1 + a2 * mag - a3 * log(dist) - a4 * dist;
% make sure MMI does not exceed 13.5
if MMI > 13.5 
    MMI = 13.5;
end




        
 
 

