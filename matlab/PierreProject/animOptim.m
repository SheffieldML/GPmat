function [] = animOptim(model, xvec, sizePath)

clf;
hold on;
X = reshape(xvec(1:2*sizePath), sizePath, 2);
Yopt = xvec(2*sizePath+1:2*sizePath+62);
Xopt = fgplvmoptimisepoint(model, X(sizePath,:), Yopt, false, 200);
plot(Xopt(1), Xopt(2), 'ks');
plot(model.X(:,1), model.X(:,2), 'g.')
plot(X(:,1),X(:,2),'ro');

drawnow;



end

