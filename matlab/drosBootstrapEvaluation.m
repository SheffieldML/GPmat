function [res, accs] = drosBootstrapEvaluation2(drosexp, drosinsitu, rankings, validation, repeats, threshold, do_isfilter)

res = zeros([length(rankings), length(rankings), length(threshold)]);
accs = zeros([repeats, length(rankings), length(threshold)]);

baseids = find(~isnan(validation));

if do_isfilter,
  baseids = is_positive_ids(baseids, drosinsitu, drosexp);
end

[B, I, J] = unique(drosexp.fbgns(baseids), 'first');

N = length(B);
val = validation(baseids(I));

for k=1:length(rankings),
  [myr, I2, J2] = intersect(baseids, rankings{k});
  [dummy, I3] = sort(J2);
  r{k} = J(I2(I3));
end

for j=1:length(threshold),
  for n=1:repeats,
    I = randint(N, 1, [1, N]);
    inds = accumarray(I, 1, [N, 1]);
    
    myacc = zeros(1, length(rankings));
    for k=1:length(rankings),
      J = cumsum(inds(r{k}));
      t = find(J >= threshold(j), 1);
      myinds = zeros(size(inds));
      myinds(r{k}(1:t-1)) = inds(r{k}(1:t-1));
      myinds(r{k}(t)) = inds(r{k}(t)) - J(t) + threshold(j);
      myacc(k) = sum(myinds .* val) ./ sum(myinds);
    end
    accs(n, :, j) = myacc;
  end

  for k=1:length(rankings),
    for l=k+1:length(rankings),
      res(k, l, j) = mean(accs(:, k, j) >= accs(:, l, j));
      res(l, k, j) = mean(accs(:, k, j) < accs(:, l, j));
    end
  end
end


function r = randint(rows, cols, range),

%r = ((range(1)):(range(2)))';
r = floor(range(2)*rand(rows, cols)) + range(1);
