#!/usr/bin/perl

# FOr creating files for new kernels use kernelGenerator.pl {kernelName}
use warnings;
use strict;


my $prefix = $ARGV[0];
my $longName = $ARGV[1];
my $ucPrefix = uc($prefix);


my %files = 
("${prefix}NoiseParamInit.m",
"function model = ${prefix}NoiseParamInit(model)\n\n% ${ucPrefix}NOISEPARAMINIT ${longName}'s parameter initialisation.\n\n% IVM\n\nmodel.noise.",
"${prefix}GradX.m", 
"function g = ${prefix}GradX(X, Y, model, prior)\n\n% ${ucPrefix}GRADX Gradient wrt x of log-likelihood for ${longName}.\n\n% IVM\n\nif size(X, 1) > 1\n  error('This function only takes one data-point');\nend",
"${prefix}GradientParam.m", 
"function g = ${prefix}GradientParam(model, params)\n\n% ${ucPrefix}GRADIENTPARAM Gradient of the ${longName}'s parameters.\n\n% IVM\n\n",
"${prefix}Likelihood.m", 
"function L = ${prefix}Likelihood(X, Y, model)\n\n% ${ucPrefix}LIKELIHOOD Likelihood of data under ${longName}.\n\n% IVM",
"${prefix}LogLikelihood.m",
"function L = ${prefix}LogLikelihood(X, Y, model)\n\n% ${ucPrefix}LOGLIKELIHOOD Log-likelihood of data under ${longName}.\n\n% IVM\n\nL = sum(sum(log(${prefix}Likelihood(X, Y, model))));",
"${prefix}NoiseExpandParam.m",
"function noise = ${prefix}NoiseExpandParam(params, noise)\n\n% ${ucPrefix}NOISEEXPANDPARAM Expand probit noise structure from param vector.\n\n% IVM\n\nnoise.bias = params(1:end);",
"${prefix}NoiseExtractParam.m",
"function params = ${prefix}NoiseExtractParam(noise)\n\n% ${ucPrefix}NOISEEXTRACTPARAM Extract parameters from ${longName}.\n\n% IVM\n\nparams = [noise.bias];","${prefix}UpdateParams.m",
"function model = ${prefix}UpdateParams(model, index)\n\n% ${ucPrefix}UPDATEPARAMS Update parameters for ${longName}.\n\n% IVM\n\n",
"${prefix}UpdateSites.m", 
"function model = ${prefix}UpdateSites(model, index)\n\n% ${ucPrefix}UPDATESITES Update site parameters for ${longName}.\n\n% IVM\n\nif nargin < 2\n  error('No site index specified');\nend\n\nmodel = ${prefix}UpdateParams(model, index);\nmodel = updateSites(model, index);");

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
