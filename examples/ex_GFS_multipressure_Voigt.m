% Example for the global fit class (GlobalFitSimple)

% Fit a set of Voigt profiles that represent the same absorption line
% measured at different pressures.

clear all
close all

addpath('./fadf')

rng(1)  % Set a seed for the random number generation (for reproducibility)

% The model represents a pressure-boradened line
% p(1) -> Line amplitude
% p(2) -> pressure (not a fitting parameter, must be held fixed)
% p(3) -> Doppler width
% p(4) -> line center frequency (at zero pressure)
% p(5) -> pressure shift coefficient
% p(6) -> pressure broadening coefficient
model = @(x, p) pressureDependentVoigtProfile(x, p(1), p(2), p(3), p(4), p(5), p(6));

% Pressures [atm]
p = [0.1, 1, 10, 20, 50, 200, 500]/760;

% Parameters for the simulated datasets
for i=1:length(p)
    pars(i,:) = [1+0.2*randn(1,1), p(i), 0.7, 0, 0.225, 2.2];
end

N = 300;  % Points per curve
noise = 0.005;  % Absolute noise level

% Generate and plot data
figure()
hold on
for i=1:size(pars, 1)
    xData{i} = linspace(-5, 5, N);
    yData{i} = model(xData{i}, pars(i,:)) + noise*randn(size(xData{i}));
    plot(xData{i}, yData{i}, '.')
end
hold off

% Define the marices for upper and lower bounds, and for the fixed
% parameters. NaN means no bound/fix
for i=1:length(p)
    ub(i,:) = [1.5, Inf, 0.8, 0.1, 0.5, 5];
    lb(i,:) = [0.1, -Inf, 0.6, -0.1, 0.1, 1];
    fixed(i,:) = [NaN, p(i), NaN, NaN, NaN, NaN];
end

gf = GlobalFitSimple();  % Instantiate the class
gf.setData(xData, yData);  % Set the data to fit

% Set the model
gf.setModel(model, 6, [0 0 1 1 1 1])
gf.setStart(pars.*(1+0.1*randn(size(pars))))  % Set start point

% Set upper/lower bounds, and fixed parameters. For this kind of complicate
% fitting, it's necessary to set bounds to guide the fitting routine
gf.setUb(ub);
gf.setLb(lb);
gf.fixParameters(fixed);
gf.fit()  % Run the fit!
fit_pars = gf.getFittedParameters();  % Retrieve the fitted parameters...
fit_errs = gf.getParamersErrors();  % ...and their errors

% Print the difference between the original and the retrieved parameters
disp(pars - fit_pars);

% Evaluate and plot the fit
hold on
for i=1:size(pars, 1)
    yData{i} = model(xData{i}, fit_pars(i,:));
    plot(xData{i}, yData{i}, '-')
end
hold off

function y = pressureDependentVoigtProfile(x, A, p, s, x0, d0, g0)
x_p = x0 + p.*d0;  % Pressure-dependent shift
g_p = p.*g0;  % Pressure-dependent Lorentzian width

y = A.*voigt(x - x_p, s, g_p);
end