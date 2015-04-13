seqLen = 600;
nStates = 3;
if false
  % Generate random transitions.
  lt = log(rand(nStates));
  lt = log(exp(lt)./repmat(sum(exp(lt), 2), 1, 3));
  priorProb = ones(nStates, 1)/nStates;
else
  % This is a non-ergodic chain (like used in speech).
  lt = [log(199/200) -log(200) -inf; -inf log(199/200) -log(200); -inf -inf 0];
  priorProb = [1 0 0]';
end
r = rand(seqLen, 1);
for j = 1:nStates
  s(1) = 1;
  if r(1)>cumsum(priorProb(1:j))
    s(1) = j;
  end
end
for i = 1:seqLen-1
  s(i+1) = 1;
  for j=1:nStates-1
    if r(i+1)>cumsum(exp(lt(s(i), 1:j)))+eps
      s(i+1) = j+1;
    end
  end
end

% Simple observation model, gaussians with means given by the state.
sd = 1; % standard deviation.
for i = 1:seqLen
  x(i) = randn(1)*sd + s(i);
end
for i = 1:seqLen
  for j = 1:nStates
    ll(j, i) = -.5*log(2*pi) - log(sd) -.5/(sd*sd)*(x(i) - j)^2;
  end
end

shat = viterbiAlign(lt, ll);
disp('Transition')
disp(exp(lt))

disp(['Got ' num2str(sum(shat==s)) ' out of ' num2str(seqLen) ' right.']);

