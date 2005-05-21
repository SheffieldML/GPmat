function noiseDisplay(noise)

% NOISEDISPLAY Display the parameters of the noise model.

% NOISE

feval([noise.type 'NoiseDisplay'], noise)