function noiseDisplay(noise)

% KERNDISPLAY Display the parameters of the noise model.

feval([noise.type 'NoiseDisplay'], noise)