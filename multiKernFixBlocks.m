function kern = multiKernFixBlocks(kern, blocks)

% MULTIKERNFIXBLOCKS
% FORMAT
% DESC marks a number of blocks fixed to allow caching their values
% ARG kern : the structure containing the kernel.
% ARG blocks : vector of blocks to fix
% RETURN kern : the structure containing the kernel.
%
% SEEALSO : multiKernCreate, multiKernCacheBlock
%
% COPYRIGHT : Antti Honkela, 2008

% KERN

try,
  uuid = char(java.util.UUID.randomUUID);
  uuid = ['a' uuid(~(uuid=='-'))];
catch,
  uuid = sprintf('%.30f', now);
  uuid = ['a' uuid(~(uuid=='.'))];
end

kern.fixedBlocks = zeros(1, kern.numBlocks);
kern.fixedBlocks(blocks) = 1;
kern.uuid = uuid;

multiKernCacheBlock(kern);
