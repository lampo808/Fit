function [pars, errs] = test_GaussianParametersErrors(M, N, noise, fixed, Aeq, beq)
% Compare the estimated errors on fit parameters with real ones
% Input:
%  - <M>: number of points per curve
%  - <N>: number of simulated curves
%  - <noise>: Gaussian noise amplitude
%  - <fixed>: vector of fixed parameters
%  - <Aeq>: matrix for equality constraints
%  - <beq>: vector for equality constraints

model = @(x, p) p(1).*exp(-(x-p(2)).^2/(2*p(3).^2)) + p(4);

xData = linspace(-5, 5, M);

f = Fit();
f.setModel(model, 4);
f.setStart([1, 0, 1, 0]);
if nargin > 3
    f.fixParameters(fixed)
end
if nargin > 4
    f.setEqualityConstraints(Aeq, beq);
end
for i=1:N
    yData = model(xData, [1, 0, 1, 0]) + noise.*randn(size(xData));
    f.setData(xData, yData, ones(size(xData))./noise^2);
    f.fit();
    
    pars(i,:) = f.getFittedParameters();
    pars_errors(i,:) = f.getParametersErrors();
end

errs = mean(pars_errors, 1);
errs = [errs; std(pars, 1)];
pars = mean(pars, 1);

end
