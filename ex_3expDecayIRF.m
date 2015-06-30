%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example of the usage of the class Fit for a
%  3-exponential decay fit with IRF convolution.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define a periodic convoultion with pseudo-normalization
defconv = @(data, IRF) ifft(fft(data) .* fft(IRF))*sum(data)/ ...
    sum(ifft(fft(data) .* fft(IRF)));

time = linspace(0, 2.44*3000, 3001);
decay_par = [30000, 350, 8000, 1000, 1000, 3000, 0];

%% Generate triexponential decay
decay = decay_par(1)*exp(-time/decay_par(2)) + ...
    decay_par(3)*exp(-time/decay_par(4)) + ...
    decay_par(5)*exp(-time/decay_par(6)) + ...
    decay_par(7);

% Gaussian IRF
FWHM = 70;
sigma = FWHM/2.35;

IRF = exp(-(time-2.44*round(3001/2)).^2/(2*sigma.^2));
IRF = circshift(IRF, [1, round(-3001*3/8)]);  % IRF at 1/8th of the total time
decay = defconv(decay, IRF);

decay = poissrnd(circshift(decay, [0, 2]));

%% Fit
xData = time;  % Due to periodic boundary conditions, the fit is
yData = decay; %  performed on the whole data set

h = figure();
semilogy(xData, yData);

model = @(x, par) par(1)*exp(-x/par(2)) + par(3)*exp(-x/par(4)) + ...
    par(5)*exp(-x/par(6)) + par(7);  % Fitting model

fit = Fit();                % Istantiate the class
fit.setData(xData, yData);  % Set the experimental data
fit.setWeights(1./yData)    % Set the weights on the data points
fit.setModel(model, 8);     % First parameter is automatically the shift, if performing convolution
fit.setIRF(IRF);            % Set the IRF: automatically enable convolution
fit.setStart([-200, 50000, 400, 5000, 1000, 1000, 2000, 0]);
fit.setLb([-3000, 0, 200, 3000, 500, 500, 1500, 0]);
fit.setUb([3000, 1e6, 500, 10000, 1500, 5000, 5000, 1e2]);

% Plot the model evaluated on the start point
[yFit, xFit] = fit.fitEval(xData);
figure(h)
hold all
plot(xFit, yFit)
hold off

fit.fit();  % Perform the fit

fitParams = fit.getFittedParameters();  % Get the fitted parameters
chi2 = fit.getChiSquare()/numel(time);  %  and calculate the reduceh chi square

disp('-----------------')
disp(['Shift: ', num2str(fitParams(1))])
disp(['Parameters: ',num2str(fitParams(2:end))]);
disp(['Reduced chi square: ', num2str(chi2)]);

[yFit, xFit] = fit.fitEval(xData);  % Fit evaluation
figure(h)
hold all
plot(xFit, yFit)                    % Fit plotting
hold off
