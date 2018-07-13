% Example for the global fit class (GlobalFitSimple)

% Consider multiple datasets with a linear dependence on the independent
% variable, all with the same slope but different vertical positions.
% Perform a global fit of the parameters of the linear function.

clear all
close all

rng(1)  % Set a seed for the random number generation (for reproducibility)

model = @(x, p) p(1)*x + p(2);  % Linear function

% Parameters for the simulated datasets
% Four datasets (one row per dataset): shifted lines
pars = [0.3, -0.2;
        0.3, 2;
        0.3, 3.4;
        0.3, 1.7];

N = 20;  % Points per line
noise = 0.05;  % Absolute noise level

% Generate and plot data
figure()
hold on
for i=1:size(pars, 1)
    xData{i} = linspace(-2, 5, N);
    yData{i} = model(xData{i}, pars(i,:)) + noise*randn(size(xData{i}));
    plot(xData{i}, yData{i}, '.')
end
hold off


gf = GlobalFitSimple();  % Instantiate the class
gf.setData(xData, yData);  % Set the data to fit
% 2 -> number of parameters. The last array tells whether a parameter is local or global
gf.setModel(model, 2, [1 0])
gf.setStart(pars)  % Set start point (for all data)
gf.fit()  % Run the fit!
fit_pars = gf.getFittedParameters()  % Retrieve the fitted parameters...
fit_errs = gf.getParamersErrors()  % ...and their errors

% Evaluate and plot the fit
hold on
for i=1:size(pars, 1)
    yData{i} = model(xData{i}, fit_pars(i,:));
    plot(xData{i}, yData{i}, '-')
end
hold off