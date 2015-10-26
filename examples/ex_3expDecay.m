%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example of the usage of the class Fit for a 3-exponential decay fit.
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Generate triexponential decay
time = linspace(0, 2.44*3000, 3001);
decay_par = [30000, 350, 8000, 1000, 1000, 3000, 0];

decay = decay_par(1)*exp(-time/decay_par(2)) +...
    decay_par(3)*exp(-time/decay_par(4)) + decay_par(5)*exp(-time/decay_par(6));

decay = poissrnd(circshift(decay, [0, 2]));  % Poissonian noise on the data points

h = figure();
semilogy(time, decay);  % Plot of the data

%% Fit

xData = time(5:end);  % Start the fitting after the rising edge
yData = decay(5:end);

model = @(x, par) par(1)*exp(-x/par(2)) + par(3)*exp(-x/par(4)) + ...
    par(5)*exp(-x/par(6)) + par(7);  % Fitting model

fit = Fit();               % Istantiate the class
fit.setData(xData, yData);     % Set the experimental data
fit.setWeights(1./yData)       % Set the weights
fit.setModel(model, 7);        % Set the fitting model
fit.setStart([50000, 400, 5000, 1000, 1000, 2000, 0]);  % Start point for the fit
fit.setLb([0, 200, 3000, 500, 500, 1500, 0]);           % Lower bounds
fit.setUb([1e6, 500, 10000, 1500, 5000, 5000, 1e2]);    % Upper bounds

fit.fit();  % Perform the fit

fitParams = fit.getFittedParameters();  % Get the fitted parameters
chi2 = fit.getChiSquare()/numel(time);  %  and calculate the reduced chi square
% Notice that the fitted decay amplitudes are 

disp('-----------------')
disp(['Parameters: ', num2str(fitParams)]);
disp(['Reduced chi square: ', num2str(chi2)]);

yFit = fit.fitEval(xData);  % Fit evaluation
figure(h)
hold all
plot(xData, yFit)           % Superimpose the fitted curve on the exp. data
hold off