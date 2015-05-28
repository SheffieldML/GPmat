function cmpndNoiseDisplay(noise, varargin)


% CMPNDNOISEDISPLAY Display parameters of the CMPND noise.
% FORMAT
% DESC displays the parameters of the compound
% noise and the noise type to the console.
% ARG noise : the noise to display.
%
% FORMAT does the same as above, but indents the display according
% to the amount specified.
% ARG noise : the noise to display.
% ARG spacing : how many spaces to indent the display of the noise by.
%
% SEEALSO : cmpndNoiseParamInit, modelDisplay, noiseDisplay
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005

% NOISE


for i = 1:length(noise.comp)
  noiseDisplay(noise.comp{i}, varargin{:});
end