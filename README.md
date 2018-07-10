# Fit: a class for constrained fitting in MATLAB

MATLAB offers a plethora of methods to perform data fitting, based on different minimizers whith different syntaxes and capabilities. The most general one, `fmincon`, give the user much flexibility in terms of linear constraints and minimization algorithm, at the expense of ease of use.

`Fit` is a MATLAB class that exposes an user-friendly interface to perfor data fitting operations (or function minimization) using `fmincon` in a simple, straightforward way.
Only few lines of code are needed to perform a simple data fit:
```MATLAB
% Generate data
x = linspace(0, 10, 100);
y = sin(3*x) + 4;

model = @(x, par) sin(par(1)*x) + par(2);

% Initialize the fitting
F = Fit();
F.setData(x, y);
F.setModel(model, 2);
F.setStart([2.8, 4.5]);

% Fit and retrieve parameters!
F.fit()
pars = F.getFittedParameters();
```
while more advanced options can be added easily via other user-friendly methods:
```MATLAB
...
fit = Fit();
fit.setData(xData, yData);
fit.setWeights(1./yData);  % Weights on data points
fit.setModel(model, 8);
fit.setIRF(IRF);  % Fit the data taking automagically into account the system IRF
fit.setStart(fitCond(1,:));
fit.setLb(fitCond(2,:));  % Lower/upper bound support
fit.setUb(fitCond(3,:));

[yFit, xFit] = fit.fitEval(xData);  % Evaluate the best fitting model
fitParams = fit.getFittedParameters();
chi2 = fit.getChiSquare()/numel(time);  % Get the chi square
parErrors = fit.getParametersErrors();  % Retrieve the parameters uncertainties
...
```
<!-- TODO: put examples here -->


<!-- ## Content

 * Fit.m
 * fshift.m -->