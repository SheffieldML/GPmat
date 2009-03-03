function meanFunction = lfmwhiteMeanCreate(q, d, varargin)

% LFMWHITEMEANFUNCTIONCREATE creates the mean function for a multi output
% GP model based in the LFM-WHITE kernel (second order ODE excited by a
% white noise input process). The outputs of the model are generated
% according to
%
%       mean_q = B_q/D_q
%
% where mean_q is an output constant corresponding to the mean of the 
% output function q, B_q is basal transcription and D_q is the spring 
% constant.
%
% FORMAT
% DESC returns a structure for the mean function for the multiple output 
% Gaussian process model that uses the LFM-WHITE kernel.
% RETURN model : the structure for the multigp model
% ARG q : input dimension size.
% ARG d : output dimension size.
% ARG options : contains the options for the MEAN of the MULTIGP model.
%
% SEE ALSO: lfmwhiteKernParamInit, lfmwhiteKernCompute
%
% COPYRIGHT : Mauricio Alvarez and Neil D. Lawrence, 2008
%
% MODIFICATIONS : David Luengo, 2009

% MULTIGP

if q > 1
  error('LFM-WHITE MEAN FUNCTION only valid for one-D input.')
end

meanFunction.type = 'lfmwhite';
meanFunction.basal  = ones(d,1);
meanFunction.spring = ones(d,1);
meanFunction.transforms.index = d+1:2*d;
meanFunction.transforms.type = optimiDefaultConstraint('positive');
% Only the parameters of basal rates are counted. The springs are already
% counted in the kernel
meanFunction.nParams = 2*d; 
