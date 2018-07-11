% Example for the global fit class (GlobalFitSimple)

% Fit a set of multiexponential decays where the decay times are common
% between the datasets

clear all
close all

rng(1)  % Set a seed for the random number generation (for reproducibility)

toColumn = @(x) x(:);  % Helper function: the data must be in column vectors

model = @(x, p) p(1)*exp(-x/p(2)) + p(3)*exp(-x/p(4));  % Bi-exponential decay

% Parameters for the simulated dataset
% Seven datasets (one row per dataset): exponentials with the same decay
% time, but different amplitudes

t1 = 150;  % Common decay times
t2 = 720;
for i=1:7
    pars(i,:) = [230*(1 + tanh(i-4)), t1, 120*(1 + tanh(3-i)), t2];
end

N = 120;  % Points per curve

% Generate and plot data
figure()
hold on
for i=1:size(pars, 1)
    xData{i} = toColumn(linspace(0, 5*t2, N));
    % Poisson sampling of the point ordinate (like photon counting)
    yData{i} = toColumn(poissrnd(model(xData{i}, pars(i,:))));
    plot(xData{i}, yData{i}, '.')
end
hold off

gf = GlobalFitSimple();  % Instantiate the class
gf.setData(xData, yData);  % Set the data to fit

% 4 -> number of parameters. The last array tells whether a parameter is local or global
gf.setModel(model, 4, [0 1 0 1])
% Set start point (for all data). randn moves them from the true value
gf.setStart(pars.*(1+0.2*randn(size(pars))))
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
plot(pars(:,3))
plot(fit_pars(:,3), 'o')
hold off
legend('First population (true)', 'First population (fitted)', ...
    'Second population (true)', 'Second population (fitted)');