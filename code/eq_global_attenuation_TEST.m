function eq_global_attenuation_TEST
% eq attenuation test
% NAME:
%   eq_global_attenuation_TEST
% PURPOSE:
%   Test different attenuation functions
%
%   all set in code to allow for maximum flexibility
%
% CALLING SEQUENCE:
%   eq_global_attenuation_TEST
% EXAMPLE:
%   eq_global_attenuation_TEST
% INPUTS:
% OPTIONAL INPUT PARAMETERS:
% OUTPUTS:
%   plots attenuation function
% MODIFICATION HISTORY:
% David N. Bresch, david.bresch@gmail.com, 20141016
%-

R=1:1:500; % distance from epicenter, in km

legend_str={};legend_entry_i=1;cr=0.1;cg=0.5;cb=0.1; % init

for M=3:2:9 % magnitude (Richter, I presume)
    for H=0:0 % focal depth, in km
        
        % Po-Shen Lin and Chyi-Tyi Lee, 2008: Ground-Motion Attenuation
        % Relationships for Subduction-Zone Earthquakes in Northeastern Taiwan,
        % Bulletin of the Seismological Society of America, Vol. 98, No. 1, pp.
        % 220?240, February 2008, doi: 10.1785/0120060002
        % formulas 13 and 14 (p 226) and table 3/4 (p 227/228)
        Zt=0.5; % Zt= 0 for interface earthquakes, and Zt=1 for intraslab earthquakes
        lnPGA_rock = -2.500 +1.205*M -1.905*log(R +0.51552*exp(0.63255*M)) +0.0075*H +0.275*Zt; % rock sites
        lnPGA_soil = -0.900 +1.000*M -1.900*log(R +0.99178*exp(0.52632*M)) +0.0040*H +0.310*Zt; % soil sites
        % where PGA is the geometric mean value (acceleration of gravity) of the
        % horizontal PGA, M is the moment magnitude, R is the hypocentral distance,
        % H is the focal depth (kilometers), and Zt indicates the subduction zone
        % earthquake type.
        
        loglog(R,exp(lnPGA_rock),'--','Color',[cr cg cb],'MarkerSize',1); hold on
        cr=cr+0.05;
        loglog(R,exp(lnPGA_soil),'Color',[cr cg cb],'MarkerSize',1);
        cb=cb+0.05;
        legend_str{legend_entry_i}=sprintf('rock M%i H%i',M,H);legend_entry_i=legend_entry_i+1;
        legend_str{legend_entry_i}=sprintf('soil M%i H%i',M,H);legend_entry_i=legend_entry_i+1;
        
    end % H
end % M

legend(legend_str);
set(gcf,'Color',[1 1 1]);
title('Po-Shen Lin and Chyi-Tyi Lee, 2008');
grid on

return


% other TESTs below

% convert Richter to MMI
% if dep>1
%     MMI = 0.44 + 1.70 * mag - 0.0048*dep - 2.73*log(dep);
% else
%     MMI = 0.44 + 1.70 * mag;
% end


mag=6;
dep=0; % to keep it simple
D=1:1:500;

M=mag;H=dep;R=D; % often these names are used

%   for the attenuation function, see www.seismo.ethz.ch/static/gshap/earift/report.html
%   (another option might be
%   www.geo.arizona.edu/gsat/1887eq/other/Bakun_2006.pdf)
%
%   C. B. Crouse (1991) Ground?Motion Attenuation Equations for Earthquakes
%   on the Cascadia Subduction Zone. Earthquake Spectra: May 1991, Vol. 7,
%   No. 2, pp. 201-236. doi: http://dx.doi.org/10.1193/1.1585626
%

lna_Jonathatn  = 3.024 + 1.030*mag - 1.351*log(D) - 0.0008*D; %(Jonathan, 1996)
lna_Twesigomwe = 2.832 + 0.866*mag -       log(D) - 0.0025*D; %(Twesigomwe, 1997)
log10a = 0.41*mag - log10(D + 0.032*10^(0.41*mag)) - 0.0034*D + 1.3; % log10
lnPGA_Crouse =        6.36 + 1.76*mag - 2.73*log(D + 1.58*exp(0.608*mag)) + 0.00916*dep; % C. B. Crouse, 1991

% Po-Shen Lin and Chyi-Tyi Lee, 2008: Ground-Motion Attenuation
% Relationships for Subduction-Zone Earthquakes in Northeastern Taiwan,
% Bulletin of the Seismological Society of America, Vol. 98, No. 1, pp.
% 220?240, February 2008, doi: 10.1785/0120060002
% formulas 13 and 14 (p 226) and table 3/4 (p 227/228)
Zt=0.5; % Zt= 0 for interface earthquakes, and Zt=1 for intraslab earthquakes
lnPGA_rock = -2.500 +1.205*M -1.905*log(R +0.51552*exp(0.63255*M)) +0.0075*H +0.275*Zt; % rock sites
lnPGA_soil = -0.900 +1.000*M -1.900*log(R +0.99178*exp(0.52632*M)) +0.0040*H +0.310*Zt; % soil sites
% where PGA is the geometric mean value (acceleration of gravity) of the
% horizontal PGA, M is the moment magnitude, R is the hypocentral distance,
% H is the focal depth (kilometers), and Zt indicates the subduction zone
% earthquake type.

loglog(D,exp(lna_Jonathatn),'-b');hold on
loglog(D,exp(lna_Twesigomwe),'-r');
loglog(D,10.^log10a,'-g');
loglog(D,exp(lnPGA_Crouse),'-m');
loglog(D,exp(lnPGA_rock),':r');
loglog(D,exp(lnPGA_soil),':g');
legend('Jonathan','Twesigomwe','log10','Crouse','Lin&Lee 2008 rock','Lin&Lee 2008 soil')
% compare with www.seismo.ethz.ch/static/gshap/earift/figure4.jpeg
