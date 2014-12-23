function single_eq_event = climada_get_single_event(eq_data, event_i)

% NAME:
%   climada_get_single_event
% PURPOSE:
%   Extract a single event from a set of hazards
% PREVIOUS STEP:
%   Generation of a hazard set, either by reading a set of historic hazard
%   events (see eq_isc_gem_read) or by creating a probabilistic hazard
%   event set based on a historic event set (see eq_global_probabilistic)
% INPUTS:
%   eq_data, a structure with the following fields:
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
%       country(line_i): country
%   event_i: number of the event to be extracted from eq_data. If not
%   given, event_i is set to 1
% OUTPUT:
%   single_event: a structure with the same fields as eq_data, but 
%   containing only the data for the (event_i)th hazard in eq_data
%
% MODIFICATION HISTORY:
% Melanie Bieli, melanie.bieli@bluewin.ch, 20141213, initial
%-

% check input arguments
if ~exist('eq_data','var') || ~isstruct(eq_data)
    fprintf('error: climada_get_single_event requires input argument eq_data of type struct')
end
if ~exist('event_i', 'var'), event_i = 1; end
if mod(event_i,1) ~= 0 
    fprintf('error: input event_i for climada_get_single_event has to be an integer.');
    return;
end

% initialize struct single_event 
fields = fieldnames(eq_data);
single_eq_event = struct(fields{1},eq_data.(fields{1})(event_i));
% extract the (event_i)th event from eq_data
for i = 2 : numel(fields)
    if size(eq_data.(fields{i}),1) > 1
        single_eq_event.(fields{i}) = eq_data.(fields{i})(event_i);
    else
        single_eq_event.(fields{i}) = eq_data.(fields{i})(1);
    end
end

% set some of the fields manually 
% TO DO: This could be solved more elegantly (i.e. in a way that needs no 
% manual adjustment, but for the time being I leave it at that
single_eq_event.n_epicenters_orig = 1;
single_eq_event.country = char(eq_data.country{event_i});
single_eq_event.filename = eq_data.filename;
end

 