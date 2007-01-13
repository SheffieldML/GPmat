function noiseDisplay(noise, varargin)

% NOISEDISPLAY Display the parameters of the noise model.
% FORMAT
% DESC display the type of noise model and any associated
% parameters.
% ARG noise : noise model to display.
% ARG spacing : any spacing to place before the noise model.
% 
% SEEALSO modelDisplay
%
% COPYRIGHT : Neil D. Lawrence, 2005

% NOISE

fhandle = str2func([noise.type 'NoiseDisplay']);
fhandle(noise, varargin{:});