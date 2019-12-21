function y = isoctave()
% Test whether the function is running in OCTAVE or not.
% Retuns 1 (true) if yes, 0 (false) otherwise.
%
% _____________________________________
% Anderson Winkler
% FMRIB / University of Oxford
% Jan/2012

persistent isoct;
if isempty(isoct),
    isoct = exist('OCTAVE_VERSION','builtin') ~= 0;
end;
y = isoct;
