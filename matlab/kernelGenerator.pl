#!/usr/bin/perl

# FOr creating files for new kernels. Run as 
# perl kernelGenerator.pl {kernel short name} {kernel long name}
use warnings;
use strict;


my $prefix = $ARGV[0];
my $longName = $ARGV[1];
my $ucPrefix = uc($prefix);


my %files = 
("${prefix}KernCompute.m", 
"function k = ${prefix}KernCompute(kern, x, x2)\n\n% ${ucPrefix}KERNCOMPUTE Compute the $longName kernel given the parameters and X.\n\n% IVM\n\nk =",
"${prefix}KernDiagCompute.m", 
"function k = ${prefix}KernDiagCompute(kern, x)\n\n% ${ucPrefix}KERNDIAGCOMPUTE Compute diagonal of $longName kernel.\n\n% IVM\n\nk = ",
"${prefix}KernExpandParam.m", 
"function kern = ${prefix}KernExpandParam(kern, params)\n\n% ${ucPrefix}KERNEXPANDPARAM Create kernel structure from ${longName}'s parameters.\n\n% IVM\n\nkern.",
"${prefix}KernExtractParam.m",
"function params = ${prefix}KernExtractParam(kern)\n\n% ${ucPrefix}KERNEXTRACTPARAM Extract parameters from $longName kernel structure.\n\n% IVM\n\nparams = ",
"${prefix}KernGradient.m",
"function g = ${prefix}KernGradient(kern, x, covGrad)\n\n% ${ucPrefix}KERNGRADIENT Gradient of $longName kernel's parameters.\n\n% IVM\n",
"${prefix}KernParamInit.m", 
"function kern = ${prefix}KernParamInit(kern)\n\n% ${ucPrefix}KERNPARAMINIT $longName kernel parameter initialisation.\n\n% IVM\n\nkern.",
"${prefix}KernGradX.m", 
"function gX = ${prefix}KernGradX(kern, x, X2)\n\n% ${ucPrefix}KERNGRADX Gradient of $longName kernel with respect to a point x.\n\n% IVM\n\ngX =",
"${prefix}KernDiagGradX.m", 
"function gX = ${prefix}KernDiagGradX(kern, x)\n\n% ${ucPrefix}KERNDIAGGRADX Gradient of $longName kernel's diagonal with respect to a point x.\n\n% IVM\n\ngX =",
"${prefix}KernDisplay.m", 
"function ${prefix}KernDisplay(kern)\n\n% ${ucPrefix}KERNDISPLAY Display parameters of $longName kernel.\n\n% IVM\n\ngX =");

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
