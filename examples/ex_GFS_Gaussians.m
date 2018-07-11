% Example for the global fit class (GlobalFitSimple)

% Fit a set of Gaussians whose amplitude is varying but whose
% position and width are fixed

clear all
close all

rng(1)  % Set a seed for the random number generation (for reproducibility)

toColumn = @(x) x(:);  % Helper function: the data must be in column vectors

Gaussian = @(x, A, x0, s) A*exp(-(x-x0).^2/(2*s.^2));
model = @(x, p) Gaussian(x, p(1), p(2), p(3)) + Gaussian(x, p(4), p(5), p(6)) + ...
    Gaussian(x, p(7), p(8), p(9));  % Sum of three Gaussians

% Parameters for the simulated dataset
% Five datasets (one row per dataset): exponentials with the same decay
% time, but different amplitudes

x0_1 = 2;  % Common positions and widths
x0_2 = 3;
x0_3 = 6;
s_1 = 1;
s_2 = 2;
s_3 = 0.5;
for i=1:5
    pars(i,:) = [1 + tanh((i-2)*2), x0_1, s_1, ...
        0.3*(1 + tanh(3-i)), x0_2, s_2, ...
        0.5*(1 + tanh(i-3)), x0_3, s_3];
end

N = 50;  % Points per curve
noise = 0.05;  % Poisson noise on the data

% Generate and plot data
figure()
hold on
for i=1:size(pars, 1)
    xData{i} = toColumn(linspace(-1, 8, N));
    yData{i} = toColumn(model(xData{i}, pars(i,:))) + ...
        noise/10*poissrnd(10, size(xData{i}));
    plot(xData{i}, yData{i}, '.')
end
hold off

gf = GlobalFitSimple();  % Instantiate the class
gf.setData(xData, yData);  % Set the data to fit

% 9 -> number of parameters. The last array tells whether a parameter is local or global
gf.setModel(model, 9, [0 1 1 0 1 1 0 1 1])
% Set start point (for all data). randn moves them from the true value
gf.setStart(pars.*(1+0.1*randn(size(pars))))
gf.fit()  % Run the fit!
fit_pars = gf.getFittedParameters();  % Retrieve the fitted parameters...
fit_errs = gf.getParamersErrors();  % ...and their errors

% Evaluate and plot the fit
hold on
for i=1:size(pars, 1)
    yData{i} = toColumn(model(xData{i}, fit_pars(i,:)));
    plot(xData{i}, yData{i}, '-')
end
hold off

% Compare simulated and fitted populations
figure()
hold on
plot(pars(:,1))
plot(fit_pars(:,1), 'o')
plot(pars(:,4))
plot(fit_pars(:,4), 'o')
plot(pars(:,7))
plot(fit_pars(:,7), 'o')
hold off
legend('First population (true)', 'First population (fitted)', ...
    'Second population (true)', 'Second population (fitted)',...
    'Third population (true)', 'Third population (fitted)');