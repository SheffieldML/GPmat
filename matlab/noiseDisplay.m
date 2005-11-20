function noiseDisplay(noise)

% NOISEDISPLAY Display the parameters of the noise model.

% NOISE

fhandle = str2func([noise.type 'NoiseDisplay']);
fhandle(noise);