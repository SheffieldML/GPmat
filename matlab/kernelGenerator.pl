#!/usr/bin/perl

# FOr creating files for new kernels use kernelGenerator.pl {kernelName}
use warnings;
use strict;


my $prefix = $ARGV[0];
my $longName = $ARGV[1];
my $ucPrefix = uc($prefix);


my %files = 
("${prefix}KernCompute.m", 
"function k = ${prefix}KernCompute(x, kern, x2)\n\n% ${ucPrefix}KERNCOMPUTE Compute the $longName kernel given the parameters and X.\n\n% IVM\n\nk =",
"${prefix}KernDiagCompute.m", 
"function k = ${prefix}KernDiagCompute(x, kern)\n\n% ${ucPrefix}KERNDIAGCOMPUTE Compute diagonal of $longName kernel.\n\n% IVM\n\nk = ",
"${prefix}KernExpandParam.m", 
"function kern = ${prefix}KernExpandParam(params, kern)\n\n% ${ucPrefix}KERNEXPANDPARAM Create kernel structure from ${longName}'s parameters.\n\n% IVM\n\nkern.",
"${prefix}KernExtractParam.m",
"function params = ${prefix}KernExtractParam(kern)\n\n% ${ucPrefix}KERNEXTRACTPARAM Extract parameters from $longName kernel structure.\n\n% IVM\n\nparams = ",
"${prefix}KernGradient.m",
"function g = ${prefix}KernGradient(kern, x, covGrad)\n\n% ${ucPrefix}KERNGRADIENT Gradient of $longName kernel's parameters.\n\n% IVM\n",
"${prefix}KernParamInit.m", 
"function kern = ${prefix}KernParamInit(kern)\n\n% ${ucPrefix}KERNPARAMINIT $longName kernel parameter initialisation.\n\n% IVM\n\nkern.",
"${prefix}KernGradX.m", 
"function gX = ${prefix}KernGradX(x, kern, X2)\n\n% ${ucPrefix}KERNGRADX Gradient of $longName kernel with respect to a point x.\n\n% IVM\n\ngX =",
"${prefix}KernDiagGradX.m", 
"function gX = ${prefix}KernDiagGradX(x, kern)\n\n% ${ucPrefix}KERNDIAGGRADX Gradient of $longName kernel's diagonal with respect to a point x.\n\n% IVM\n\ngX =");

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
