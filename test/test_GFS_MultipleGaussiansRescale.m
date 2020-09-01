function [pars, fit_pars, fit_errs] = test_GFS_MultipleGaussiansRescale(...
    offset, scaling, seed)

if nargin < 3
    rng(1)  % Set a seed for the random number generation (for reproducibility)
else
    rng(seed)
end

Gaussian = @(x, A, x0, s) A*exp(-(x-x0).^2/(2*(s*1e10).^2));
model = @(x, p) Gaussian(x, p(1), p(2), p(3)) + Gaussian(x, p(4), p(5), p(6)) + ...
    Gaussian(x, p(7), p(8), p(9));  % Sum of three Gaussians

% Parameters for the simulated dataset
% Five datasets
x0_1 = 2e6;  % Common positions and widths
x0_2 = 2e6 + 1;
x0_3 = 2e6 + 6;
s_1 = 1e-10;
s_2 = 2e-10;
s_3 = 0.5e-10;
for i=1:5
    pars(i,:) = [1 + tanh((i-2)*2), x0_1, s_1, ...
        0.3*(1 + tanh(3-i)), x0_2, s_2, ...
        0.5*(1 + tanh(i-3)), x0_3, s_3];
end

start = pars;
start(:,[1,4,7]) = start(:,[1,4,7])/2;

N = 50;  % Points per curve
noise = 0.05;  % Gaussian noise on the data

% Generate data
for i=1:size(pars, 1)
    xData{i} = linspace(2e6 - 5, 2e6 + 8, N);
    yData{i} = model(xData{i}, pars(i,:)) + noise*randn(size(xData{i}));
end

gf = GlobalFitSimple();
gf.setData(xData, yData);
gf.setModel(model, 9, [0 1 1 0 1 1 0 1 1])
% Set start point (for all data).
gf.setStart(start, offset, scaling)

gf.fit(1)

fit_pars = gf.getFittedParameters();
fit_errs = gf.getParamersErrors();

if any(imag(fit_errs) ~= 0, 'all')
    disp('Complex errors')
end

% % Evaluate and plot the data and fit - for debug
% figure()
% hold on
% for i=1:size(pars, 1)
%     plot(xData{i}, yData{i}, '.')
%     yData{i} = model(xData{i}, fit_pars(i,:));
%     plot(xData{i}, yData{i}, '-')
% end
% hold off

end