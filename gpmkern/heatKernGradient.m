function g = heatKernGradient(kern, x, varargin)

% HEATKERNGRADIENT Gradient of HEAT kernel's parameters.
% FORMAT
% DESC computes the gradient of functions with respect to the heat
% kernel's parameters. As well as the kernel structure and the
% input positions, the user provides a matrix PARTIAL which gives
% the partial derivatives of the function with respect to the
% relevant elements of the kernel matrix.
% ARG kern : the kernel structure for which the gradients are being
% computed.
% ARG x : the input locations for which the gradients are being
% computed.
% ARG partial : matrix of partial derivatives of the function of
% interest with respect to the kernel matrix. The argument takes
% the form of a square matrix of dimension  numData, where numData is
% the number of rows in X.
% RETURN g : gradients of the function of interest with respect to
% the kernel parameters. The ordering of the vector should match
% that provided by the function kernExtractParam.
%
% FORMAT
% DESC computes the derivatives as above, but input locations are
% now provided in two matrices associated with rows and columns of
% the kernel matrix.
% ARG kern : the kernel structure for which the gradients are being
% computed.
% ARG x1 : the input locations associated with the rows of the
% kernel matrix.
% ARG x2 : the input locations associated with the columns of the
% kernel matrix.
% ARG partial : matrix of partial derivatives of the function of
% interest with respect to the kernel matrix. The matrix should
% have the same number of rows as X1 and the same number of columns
% as X2 has rows.
% RETURN g : gradients of the function of interest with respect to
% the kernel parameters.
%
% SEEALSO heatKernParamInit, kernGradient, heatKernDiagGradient, kernGradX
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% KERN

% if length(varargin)<2
%     covGrad = varargin{1};
%     if size(x, 2) ~= 2
%         error('Input can only have two columns');
%     end
%     % Split the domain into time domain and spatial domain and account for
%     % missing values. If there are no missing values the computation of the
%     % kernel is a pointwise prodruct, otherwise it is a kronecker product.
%     t1 = x(x(:,1)~=Inf,1);
%     s1 = x(x(:,2)~=Inf,2);
%     if (length(t1) == length(s1))
%         ut1 = unique(t1);
%         us1 = unique(s1);
%         if (length(ut1)*length(us1) == length(t1))
%             t1 = ut1; s1 = us1;
%             isPointwise = false;
%             sK = zeros(length(t1)*length(s1));
%         else
%             isPointwise = true;
%             sK = zeros(length(t1));
%         end
%     else
%         isPointwise = false;
%         sK = zeros(length(t1)*length(s1));
%     end
%
%     % Although this is done in heatKernExpandParam.m, we do it here again as a
%     % precaution.
%     kern.sim.inverseWidth = kern.inverseWidthTime;
%     kernTempo = kern;
%     sigmax = sqrt(2/kern.inverseWidthSpace);
%     lengthX = kern.lengthX;
%     nterms = kern.nTerms;
%     decay = kern.decay;
%     diff = kern.diffusion;
%
%     w = ((1:nterms)*(pi/lengthX))';
%     gamma = sqrt(-1)*w;
%     beta = decay + diff*(w.^2);
%     z1 = sigmax*gamma/2;
%     z2 = lengthX/sigmax + z1;
%     wz1 = wofzPoppe(sqrt(-1)*z1);
%     wz2 = wofzPoppe(sqrt(-1)*z2);
%     cK = 4/(lengthX^2);
%
%     g1 = zeros(1,5);
%     g2 = zeros(1,5);
%
%     if isPointwise
%         for i=1:nterms
%             kern.sim.decay = beta(i);
%             Kt = simXsimKernCompute(kern.sim, kern.sim, t1, t1);
%             Ks = sheatKernCompute(sigmax, lengthX, s1, s1, w, gamma, wz1, wz2, i, i);
%             covGradt = covGrad.*Ks;
%             covGrads = covGrad.*Kt;
%             [g1A, g2A] = simXsimKernGradient(kern.sim, kernTempo.sim, t1, t1, covGradt);
%             g1(1) = g1(1) + g1A(1);   % Decay for first kernel
%             g1(2) = g1(2) + (w(i)^2)*g1A(1); % Diffusion rate for the first kernel
%             g1(3) = g1(3) + g1A(2);          % Inverse width for time
%             g2(1) = g2(1) + g2A(1);          % Decay for second kernel
%             g2(2) = g2(2) + (w(i)^2)*g2A(1); % Diffusion rate for the second kernel
%             g1B = sheatKernGradient(sigmax, lengthX, s1, s1, w, gamma, wz1, wz2, i, i, covGrads);
%             g1B = -(1/sqrt(2*kern.inverseWidthSpace^3))*g1B; % Transforms to the derivative of the inverse width
%             g1(4) = g1(4) + g1B;
%             sK = sK + Kt.*Ks;
%             for j=1:i-1
%                 if (mod(i+j,2)==0)
%                     kern.sim.decay = beta(i);
%                     kernTempo.sim.decay = beta(j);
%                     Kt = simXsimKernCompute(kern.sim, kernTempo.sim, t1, t1);
%                     Ks = sheatKernCompute(sigmax, lengthX, s1, s1, w, gamma, wz1, wz2, i, j);
%                     covGradt = covGrad.*Ks;
%                     covGrads = covGrad.*Kt;
%                     [g1A, g2A] = simXsimKernGradient(kern.sim, kernTempo.sim, t1, t1, covGradt);
%                     g1(1) = g1(1) + 2*g1A(1);   % Decay for first kernel
%                     g1(2) = g1(2) + 2*(w(i)^2)*g1A(1); % Diffusion rate for the first kernel
%                     g1(3) = g1(3) + 2*g1A(2);          % Inverse width for time
%                     g2(1) = g2(1) + 2*g2A(1);          % Decay for second kernel
%                     g2(2) = g2(2) + 2*(w(j)^2)*g2A(1); % Diffusion rate for the second kernel
%                     g1B = sheatKernGradient(sigmax, lengthX, s1, s1, w, gamma, wz1, wz2, i, j, covGrads);
%                     g1B = -(1/sqrt(2*kern.inverseWidthSpace^3))*g1B; % Transforms to the derivative of the inverse width
%                     g1(4) = g1(4) + 2*g1B;
%                     sK = sK + Kt.*Ks +  (Kt').*(Ks');
%                 end
%             end
%         end
%     else
%         covGradt = zeros(length(t1));
%         covGrads = zeros(length(s1));
%         for i=1:nterms
%             kern.sim.decay = beta(i);
%             Kt = simXsimKernCompute(kern.sim, kern.sim, t1, t1);
%             Ks = sheatKernCompute(sigmax, lengthX, s1, s1, w, gamma, wz1, wz2, i, i);
%             % These loops might slow everything
%             startOne = 1;
%             endOne = 0;
%             for k=1:length(t1)
%                 endOne = endOne + length(s1);
%                 startTwo = 1;
%                 endTwo = 0;
%                 for l=1:length(t1)
%                     endTwo = endTwo + length(s1);
%                     covGradt(k,l) = sum(sum(covGrad(startOne:endOne, startTwo:endTwo).*Ks));
%                     startTwo = endTwo + 1;
%                 end
%                 startOne = endOne + 1;
%             end
%             for k=1:length(s1)
%                 indRows = (k:length(s1):(length(t1)*length(s1)))';
%                 for l=1:length(s1)
%                     indCols = l:length(s1):(length(t1)*length(s1));
%                     covGrads(k,l) = sum(sum(covGrad(indRows, indCols).*Kt));
%                 end
%             end
%             [g1A, g2A] = simXsimKernGradient(kern.sim, kernTempo.sim, t1, t1, covGradt);
%             g1(1) = g1(1) + g1A(1);   % Decay for first kernel
%             g1(2) = g1(2) + (w(i)^2)*g1A(1); % Diffusion rate for the first kernel
%             g1(3) = g1(3) + g1A(2);          % Inverse width for time
%             g2(1) = g2(1) + g2A(1);          % Decay for second kernel
%             g2(2) = g2(2) + (w(i)^2)*g2A(1); % Diffusion rate for the second kernel
%             g1B = sheatKernGradient(sigmax, lengthX, s1, s1, w, gamma, wz1, wz2, i, i, covGrads);
%             g1B = -(1/sqrt(2*kern.inverseWidthSpace^3))*g1B; % Transforms to the derivative of the inverse width
%             g1(4) = g1(4) + g1B;
%             sK = sK + kron(Kt, Ks);
%             for j=1:i-1
%                 if (mod(i+j,2)==0)
%                     kern.sim.decay = beta(i);
%                     kernTempo.sim.decay = beta(j);
%                     Kt = simXsimKernCompute(kern.sim, kernTempo.sim, t1, t1);
%                     Ks = sheatKernCompute(sigmax, lengthX, s1, s1, w, gamma, wz1, wz2, i, i);
%                     % These loops might slow everything
%                     startOne = 1;
%                     endOne = 0;
%                     for k=1:length(t1)
%                         endOne = endOne + length(s1);
%                         startTwo = 1;
%                         endTwo = 0;
%                         for l=1:length(t1)
%                             endTwo = endTwo + length(s1);
%                             covGradt(k,l) = sum(sum(covGrad(startOne:endOne, startTwo:endTwo).*Ks));
%                             startTwo = endTwo + 1;
%                         end
%                         startOne = endOne + 1;
%                     end
%                     for k=1:length(s1)
%                         indRows = (k:length(s1):(length(t1)*length(s1)))';
%                         for l=1:length(s1)
%                             indCols = l:length(s1):(length(t1)*length(s1));
%                             covGrads(k,l) = sum(sum(covGrad(indRows, indCols).*Kt));
%                         end
%                     end
%                     [g1A, g2A] = simXsimKernGradient(kern.sim, kernTempo.sim, t1, t1, covGradt);
%                     g1(1) = g1(1) + g1A(1) + g2A(1);   % Decay for first kernel
%                     g1(2) = g1(2) + (w(i)^2)*g1A(1) + (w(j)^2)*g2A(1); % Diffusion rate for the first kernel
%                     g1(3) = g1(3) + 2*g1A(2);          % Inverse width for time
%                     g2(1) = g2(1) + g2A(1) + g1A(1);          % Decay for second kernel
%                     g2(2) = g2(2) + (w(j)^2)*g2A(1) + (w(i)^2)*g1A(1); % Diffusion rate for the second kernel
%                     g1B = sheatKernGradient(sigmax, lengthX, s1, s1, w, gamma, wz1, wz2, i, j, covGrads);
%                     g1B = -(1/sqrt(2*kern.inverseWidthSpace^3))*g1B; % Transforms to the derivative of the inverse width
%                     g1(4) = g1(4) + 2*g1B;
%                     sK = sK + kron(Kt, Ks) + kron(Kt', Ks');
%                 end
%             end
%         end
%     end
%     g1(1:4) = kern.sensitivity^2*cK*g1(1:4);
%     g2(1:4) = kern.sensitivity^2*cK*g2(1:4);
%     g1(5) = kern.sensitivity*cK*(sum(sum(covGrad.*sK)));
%     g2(5) = kern.sensitivity*cK*(sum(sum(covGrad.*sK)));
%
%     g = real(g1 + g2);
%
% else
%
if length(varargin)<2
    x2 = x;
end
[g1, g2] = heatXheatKernGradient(kern, kern, x, x2, varargin{end});
g = g1 + g2;
%end


