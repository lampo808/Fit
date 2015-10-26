function [fitParams, chi2] = test_3expDecayIRF(display, time, decay_par, fitCond, IRF_par, fitPoints)
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
    fitCond = [-200, 50000, 400, 5000, 1000, 1000, 2000, 0;  % First parameter is the shift
        -3000, 0, 200, 3000, 500, 500, 1500, 0;
        3000, 1e6, 500, 10000, 1500, 5000, 5000, 1e2];
end
if nargin < 5 || isempty(IRF_par)
    IRF_par = 70;
end
if nargin < 6 || isempty(fitPoints)
    fitPoints = 1:length(time);
end

% Periodic convoultion with pseudo-normalization
defconv = @(data, IRF) ifft(fft(data) .* fft(IRF))*sum(data)/ ...
    sum(ifft(fft(data) .* fft(IRF)));

% Generate triexponential decay
decay = decay_par(1)*exp(-time/decay_par(2)) + ...
    decay_par(3)*exp(-time/decay_par(4)) + ...
    decay_par(5)*exp(-time/decay_par(6)) + ...
    decay_par(7);

% Gaussian IRF
FWHM = IRF_par;
sigma = FWHM/2.35;

IRF = exp(-(time-2.44*round(3001/2)).^2/(2*sigma.^2));
IRF = circshift(IRF, [1, round(-3001*3/8)]);  % IRF at 1/8th of the total time
decay = defconv(decay, IRF);

decay = poissrnd(circshift(decay, [0, 2]));

% Fit
xData = time(fitPoints);
yData = decay(fitPoints);

if display
    h = figure();
    semilogy(xData, yData);
end
    
model = @(x, par) par(1)*exp(-x/par(2)) + par(3)*exp(-x/par(4)) + ...
    par(5)*exp(-x/par(6)) + par(7);

fit = Fit();
fit.setData(xData, yData);
fit.setWeights(1./yData)
fit.setModel(model, 8);  % First parameter is the shift
fit.setIRF(IRF);
fit.setStart(fitCond(1,:));
fit.setLb(fitCond(2,:));
fit.setUb(fitCond(3,:));

if display
    [yFit, xFit] = fit.fitEval(xData);
    figure(h)
    hold all
    plot(xFit, yFit)
    hold off
end

fit.fit();

fitParams = fit.getFittedParameters();
chi2 = fit.getChiSquare()/numel(time);

if display
    disp('-----------------')
    disp(['Shift: ', num2str(fitParams(1))])
    disp(['Parameters: ',num2str(fitParams(2:end))]);
    disp(['Reduced chi square: ', num2str(chi2)]);
end

if display
    [yFit, xFit] = fit.fitEval(xData);
    figure(h)
    hold all
    plot(xFit, yFit)
    hold off
end