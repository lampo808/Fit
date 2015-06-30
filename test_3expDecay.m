function [fitParams, chi2] = test_3expDecay(display, time, decay_par, fitCond, fitPoints)
if nargin < 1 || isempty(display)
    display = 1;
end
if nargin < 2 || isempty(time)
    time = linspace(0, 2.44*3000, 3001);
end
if nargin < 3 || isempty(decay_par)
    decay_par = [30000, 350, 8000, 1000, 1000, 3000, 0];
end
if nargin < 4 || isempty(fitCond)
    fitCond = [50000, 400, 5000, 1000, 1000, 2000, 0;
        0, 200, 3000, 500, 500, 1500, 0;
        1e6, 500, 10000, 1500, 5000, 5000, 1e2];
end
if nargin < 5 || isempty(fitPoints)
    fitPoints = 5:length(time);  % The first 3 points are the leading edge
end

% Generate triexponential decay
decay = decay_par(1)*exp(-time/decay_par(2)) +...
    decay_par(3)*exp(-time/decay_par(4)) + decay_par(5)*exp(-time/decay_par(6));

decay = poissrnd(circshift(decay, [0, 2]));

if display
    h = figure();
    semilogy(time, decay);
end

% Fit

xData = time(fitPoints);
yData = decay(fitPoints);

model = @(x, par) par(1)*exp(-x/par(2)) + par(3)*exp(-x/par(4)) + ...
    par(5)*exp(-x/par(6)) + par(7);

fit = Fit();
fit.setData(xData, yData);
fit.setWeights(1./yData)
fit.setModel(model, 7);
fit.setStart(fitCond(1,:));
fit.setLb(fitCond(2,:));
fit.setUb(fitCond(3,:));

fit.fit();

fitParams = fit.getFittedParameters();
chi2 = fit.getChiSquare()/numel(time);

if display
    disp('-----------------')
    disp(['Parameters: ', num2str(fitParams)]);
    disp(['Reduced chi square: ', num2str(chi2)]);
end

if display
    yFit = fit.fitEval(xData);
    figure(h)
    hold all
    plot(xData, yFit)
    hold off
end