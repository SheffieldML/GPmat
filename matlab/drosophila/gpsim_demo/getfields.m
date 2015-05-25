function r = getfields(s, fields),

sz = size(s.(fields{1}));

%assert(sz(1) == 1);

r = zeros(length(fields), sz(2));
for k = 1:length(fields),
  r(k, :) = s.(fields{k});
end
