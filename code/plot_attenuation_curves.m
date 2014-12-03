function plot_attenuation_curves(I_epi) 

% DESCRIPTION: plot_attenuation plots the regional attenuation functions 
% (decrease of an earthquake's intensity (Modified Mercalli Intensity MMI) 
% with distance R (in km) from epicenter) that are currently implemented
% y axis: I(R)(intensity at distance R to the epicenter, in Modified 
% Mercalli Intensity (MMI)
% x axis: R, i.e. the distance to the epicenter (in km)
%
% SYNTAX:
%    plot_attenuation(I_epi)
%
% EXAMPLE:
%    plot_attenuation(6)
% INPUT:
%   I_epi: Intensity at epicenter (on Modified Mercalli Scale, i.e., I_epi 
%   should be in the range [1,12]
%  
%  
% MODIFICATION HISTORY:
% Melanie Bieli, melanie.bieli@bluewin.ch, 20141107

% Defining the different attenuation functions for specific geographic 
% regions
    
    I_0 = I_epi;    % intensity at epicenter
    
    WesternUS_Chen = [];
    Western_China = [];
    WesternUS_Murphy_OBrien = [];
    Japan = [];
    SouthernEurope = [];
    Himalaya = [];
    India = [];
    Turkey = [];
    
    i = 1;
    
    for R=1:250
        % Seismic Hazard and Risk Analysis: A Simplified Approach
        % CHEN Yong, CHEN Qifu, LIU Jie, CHEN Ling, LI Juan
        % Science Press, Beijing, 2002
        % p. 80 f, using the conversion between earthquake magnitude m and
        % epicenter intensity I_0 proposed by Gutenberg and Richter, 1956:
        % m = 2/3*I_0 + 1
        WesternUS_Chen(i) = I_0 + 3.2 - 0.00106*R - 2.7*log10(R);
        Western_China(i) = 7.181 + 1.0253*I_0 - 2.109*log(R+25);
        
        % Murphy and O'Brien, 1977:
        % log(PGA) = 0.14*MMI + 0.24*mag - 0.68*log10(R) + beta_k
        % using log(PGA) = 0.25*I + 0.25
        % beta_k refers to a specific region k
        % beta_WesternUS = 0.60
        % beta_Japan = 0.69
        % beta_SouthernEurope = 0.88
        WesternUS_Murphy_OBrien(i) = 1.4545*I_0 - 0.0909 - 6.1818*log10(R) + 0.60;
        Japan(i) = 1.4545*I_0 - 0.0909 - 6.1818*log10(R) + 0.69;
        SouthernEurope(i) = 1.4545*I_0 - 0.0909 - 6.1818*log10(R) + 0.88;
        
        % Ghosh, G.K. and Mahajan, A.K. (2011). Interpretation of intensity
        % attenuation relation of 1905 Kangra earthquake with epicentral
        % distance and magnitude in the Northwest Himalayan region
        % Volume 77, Issue 6, pp 511-520
        % DOI: 10.1007/s12040-012-0261-z
        Himalaya(i) = 4.1660 + 0.8733*I_0 - 0.0017*R - 0.9598*log(R);
        
        % Szeliga, W. et al. (2010). Intensity, Magnitude, Location, and
        % Attenuation in India for Felt Earthquakes since 1762
        % Bulletin of the Seismological Society of America;
        % doi: 10.1785/0120080329
        India(i) = 6.6300 + 0.7067*I_0 - 0.0010*R - 3.37*log10(R);
        
        % Erdik and Oner, (1982).
        Turkey(i) = 6.00 + 1.3867*I_0 - 0.98*log(R);
        
        i = i+1;
    end
    
    R = 1:250;
    fig = figure;
    h = plot(R, WesternUS_Chen, 'c', R, WesternUS_Murphy_OBrien, 'g', R, Western_China, 'r', R, Japan, 'b', R, SouthernEurope, 'k', R, Himalaya, 'm', R, India, 'b', R, Turkey, 'y');
    set(h,'LineWidth',1.5)
    xlabel('Distance from epicenter (km)');     
    ylabel('Intensity (Modified Mercalli)');
    title('Regional attenuation relations');
    legend('Western U.S. (Chen)','Western U.S. (Murphy & OBrien)', 'Western China', 'Japan', 'Southern Europe', 'Himalaya', 'India', 'Turkey');

end

