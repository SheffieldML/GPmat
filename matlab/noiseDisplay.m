function noiseDisplay(noise)

% NOISEDISPLAY Display the parameters of the noise model.

% IVM

feval([noise.type 'NoiseDisplay'], noise)