#!/usr/bin/perl

# FOr creating files for new noise models. Run as 
# perl noiseGenerator.pl {noise model short name} {noise model long name}
use warnings;
use strict;


my $prefix = $ARGV[0];
my $longName = $ARGV[1];
my $ucPrefix = uc($prefix);


my %files = 
("${prefix}NoiseParamInit.m",
"function noise = ${prefix}NoiseParamInit(noise, y)\n\n% ${ucPrefix}NOISEPARAMINIT ${longName}'s parameter initialisation.\n\n% IVM\n\n",
"${prefix}NoiseGradVals.m", 
"function [dlnZ_dmu, dlnZ_dvs] = ${prefix}NoiseGradVals(noise, mu, varsigma, y)\n\n% ${ucPrefix}NOISEGRADVALS Gradient wrt mu and varsigma of log-likelihood for ${longName}.\n\n% IVM\n\n",
"${prefix}NoiseGradientParam.m", 
"function g = ${prefix}NoiseGradientParam(noise, mu, varsigma, y)\n\n% ${ucPrefix}NOISEGRADIENTPARAM Gradient of the ${longName}'s parameters.\n\n% IVM\n\n",
"${prefix}Likelihood.m", 
"function L = ${prefix}Likelihood(noise, mu, varsigma, y)\n\n% ${ucPrefix}LIKELIHOOD Likelihood of data under ${longName}.\n\n% IVM",
"${prefix}LogLikelihood.m",
"function L = ${prefix}LogLikelihood(noise, mu, varsigma, y)\n\n% ${ucPrefix}LOGLIKELIHOOD Log-likelihood of data under ${longName}.\n\n% IVM\n\nL = sum(sum(log(${prefix}Likelihood(noise, mu, varsigma, y))));",
"${prefix}NoiseExpandParam.m",
"function noise = ${prefix}NoiseExpandParam(noise, params)\n\n% ${ucPrefix}NOISEEXPANDPARAM Expand ${longName}'s structure from param vector.\n\n% IVM\n\nnoise.bias = params(1:end);",
"${prefix}NoiseOut.m", 
"function y = ${prefix}NoiseOut(noise, mu, varsigma)\n\n% ${ucPrefix}NOISEOUT Ouput from ${longName}.\n\n% IVM\n\n",
"${prefix}NoiseExtractParam.m",
"function [params, names] = ${prefix}NoiseExtractParam(noise)\n\n% ${ucPrefix}NOISEEXTRACTPARAM Extract parameters from ${longName}.\n\n% IVM\n\nparams = [noise.bias];",
"${prefix}NoiseDisplay.m",
"function ${prefix}NoiseDisplay(noise)\n\n% ${ucPrefix}NOISEDISPLAY Display  parameters from ${longName}.\n\n% IVM\n\n");

my @fileNames = keys %files; 

foreach my $fileName (@fileNames)
{
  if(-e $fileName)
  {
    print "$fileName exists\n";
  }
  else
  {
    open(FILEOUT, ">$fileName") or die "Can't open file $fileName";
    print FILEOUT $files{$fileName};
    close(FILEOUT);
  }
}
