% Example for the global fit class (GlobalFitSimple)

% Fit a set of Gaussians with different position but the same amplitude and
% width.

clear all
close all

rng(1)  % Set a seed for the random number generation (for reproducibility)

toColumn = @(x) x(:);  % Helper function: the data must be in column vectors

model = @(x, p) p(1)*exp(-(x-p(2)).^2/(2*p(3).^2));  % Gaussian function

% Parameters for the simulated datasets
% Four datasets (one row per dataset): shifted Gaussians
pars = [1, 0, 0.5; 
        1, 0.5, 0.5;
        1, 1.88, 0.5;
        1, 5.5, 0.5];

N = 100;  % Points per curve
noise = 0.05;  % Absolute noise level

% Generate and plot data
figure()
hold on
for i=1:size(pars, 1)
    x0 = pars(i,2);
    s = pars(i,3);
    xData{i} = toColumn(linspace(x0-6*s, x0+6*s, N));
    yData{i} = toColumn(model(xData{i}, pars(i,:))) + noise*randn(size(xData{i}));
    plot(xData{i}, yData{i}, '.')
end
hold off

% Define the marices for upper and lower bounds, and for the fixed
% parameters. NaN means no bound/fix
ub = [2, 10, 2; 
      NaN, 10, NaN;
      NaN, 10, NaN;
      NaN, 5, NaN];
lb = [-1, -1, 0; 
      NaN, 0, NaN;
      NaN, 0, NaN;
      NaN, 0, NaN];
fixed = [NaN, 0, 0.7; 
         NaN, NaN, NaN;
         NaN, NaN, NaN;
         NaN, NaN, NaN];

gf = GlobalFitSimple();  % Instantiate the class
gf.setData(xData, yData);  % Set the data to fit

% 3 -> number of parameters. The last array tells whether a parameter is local or global
gf.setModel(model, 3, [1 0 1]) 
gf.setStart(pars)  % Set start point (for all data)

% (optional) Set upper/lower bounds, and fixed parameters. Uncomment to see
% the effect of the constraints on the fit
% gf.setUb(ub);
% gf.setLb(lb);
% gf.fixParameters(fixed);
gf.fit()  % Run the fit! 
fit_pars = gf.getFittedParameters()  % Retrieve the fitted parameters...
fit_errs = gf.getParamersErrors()  % ...and their errors

% Evaluate and plot the fit
hold on
for i=1:size(pars, 1)
    yData{i} = toColumn(model(xData{i}, fit_pars(i,:)));
    plot(xData{i}, yData{i}, '-')
end
hold off