classdef TestFit < matlab.unittest.TestCase

    methods (Test)
        % Wrong input to functions
        function wrongInput(tc)
            F = Fit();

            tc.assertError(@()F.setData([0, 1]), ?MException, 'FIT:setDataArgumentsNumber');
        end

        % Missing fit conditions
        % TBW

        % Test fit
        function fitSine(tc)
            x = linspace(0, 10, 100);
            y = sin(3*x) + 4;

            model = @(x, par) sin(par(1)*x) + par(2);

            F = Fit();
            F.setData(x, y);
            F.setModel(model, 2);
            F.setStart([2.8, 4.5]);

            F.fit()
            pars = F.getFittedParameters();
            tc.assertEqual(pars, [3,4], 'AbsTol', 1e-6)
        end

        function fitSineNoise(tc)
            x = linspace(0, 10, 100);
            y = sin(3*x) + 4 + randn(size(x))*0.1;

            model = @(x, par) sin(par(1)*x) + par(2);

            F = Fit();
            F.setData(x, y);
            F.setModel(model, 2);
            F.setStart([2.8, 4.5]);

            F.fit()
            pars = F.getFittedParameters();
            tc.assertEqual(pars, [3,4], 'AbsTol', 5e-2)
        end

        function decay3exp_1(tc)
            decay_par = [30000, 350, 8000, 1000, 1000, 3000, 100];  % Decay params.
            start = [50000, 400, 5000, 1000, 1000, 2000 0];       % Start point
            lb = [0, 200, 3000, 500, 500, 1500, 0];               % Lower bound
            ub = [1e6, 500, 10000, 1500, 5000, 5000, 1e2];         % Upper bound

            [fitParameters, chi2] = test_3expDecay(0, [], decay_par, ...
                [start; lb; ub]);

            tc.assertEqual(fitParameters([2, 4]), decay_par([2, 4]), 'RelTol', 1e-1);
            tc.assertLessThanOrEqual(chi2, 1.5);
        end

        function decay3expIRF_1(tc)
            decay_par = [30000, 350, 8000, 1000, 1000, 3000, 100];  % Decay params.
            start = [-200, 50000, 400, 5000, 1000, 1000, 2000 0];   % Start point
            lb = [-3000, 0, 200, 3000, 500, 500, 1500, 0];          % Lower bound
            ub = [3000, 1e6, 500, 10000, 1500, 5000, 5000, 1e2];    % Upper bound

            [fitParameters, chi2] = test_3expDecayIRF(0, [], decay_par, ...
                [start; lb; ub]);

            tc.assertEqual(fitParameters([3, 5]), decay_par([2, 4]), 'RelTol', 1e-1);
            tc.assertLessThanOrEqual(chi2, 1.5);
        end

        function fitSineFixedParameter_1(tc)
            x = linspace(0, 10, 100);
            y = sin(3*x) + 4 + randn(size(x))*0.1;

            model = @(x, par) sin(par(1)*x) + par(2);

            F = Fit();
            F.setData(x, y);
            F.setModel(model, 2);
            F.setStart([2.8, 4.5]);
            F.fixParameter(2, 5);

            F.fit()
            pars = F.getFittedParameters();

            tc.assertEqual(pars(2), 5, 'AbsTol', 1e-6)
        end

        function GaussianRescale_1(tc)
            pars = [1e-3, 1e6, 0.1, 1e-10];

            Gaussian = @(x, A, x0, s) A*exp(-(x-x0).^2/(2*s.^2));
            model = @(x, p) Gaussian(x, p(1), p(2), p(3)) + 1e6*p(4);  % Gaussian with offset

            x = linspace(1e6-1, 1e6+1, 100);
            y = model(x, pars) + 1e-4*randn(size(x));

%             plot(x, y)

            % Fit with offset and scaling
            F = Fit();
            F.setData(x, y);
            F.setModel(model, 4);
            F.setStart(pars + [1e-4, -0.02, 0.2, 1e-10], ... % Start point
                [0, 1e6, 0, 0], ...                          % Offset
                [1e-3, 0.1, 0.1, 1e-10]);                    % Scaling

            F.fit()
            fit_pars = F.getFittedParameters();
            fit_errs = F.getParametersErrors();
            yFit = F.fitEval(x);
%             hold on
%             plot(x, yFit)
%             hold off

            % Fit with manual offset and scaling
            scaling = [1e-3, 0.1, 0.1, 1e-10];
            offset = [0, 1e6, 0, 0];

            offsetScaledModel = @(x, p) model(x, p.*scaling + offset);
            F.setModel(offsetScaledModel, 4);
            F.setStart((pars + [1e-4, -0.02, 0.2, 1e-10] - offset)./scaling);

            F.fit()
            fit_pars_2 = F.getFittedParameters();
            fit_pars_2 = fit_pars_2 .* scaling + offset;
            fit_errs_2 = F.getParametersErrors();
            fit_errs_2 = fit_errs_2 .* scaling;
            yFit2 = F.fitEval(x);
%             hold on
%             plot(x, yFit2)
%             hold off

            for i=1:numel(fit_pars)
                tc.assertEqual(fit_pars(i), fit_pars_2(i), 'AbsTol', 1e-6);
                tc.assertEqual(fit_errs(i), fit_errs_2(i), 'AbsTol', 1e-6);
            end
        end

        % TODO:
        % - [x] check errors on parameters
        % - [x] check errors on parameters with fixed parameters
        % - [x] check errors on parameters with equality constraints
        % - [ ] same with global fit
        function GaussianParametersErrors1(tc)
            % Check errors on parameters
            rng(1)

            [pars, errs] = test_GaussianParametersErrors(100, 50, 1e-2);
            tc.assertLessThan(abs(pars - [1 0 1 0]), 2e-3);
            tc.assertLessThan(abs(diff(errs, 1)), 2e-3);
        end

        function GaussianParametersErrors2(tc)
            % Check errors when there are fixed parameters
            rng(1)

            [pars, errs] = test_GaussianParametersErrors(10, 500, 1e-2, ...
                [1, NaN, 1, 0]);
            tc.assertLessThan(abs(pars - [1 0 1 0]), 1e-3);
            tc.assertLessThan(abs(diff(errs, 1)), 1e-3);

            % Check errors when there are equality constraints
            rng(1)

            % In this case, use constraints to fix parameters
            Aeq = [1 0 0 0; 0 0 1 0; 0 0 0 1];
            beq = [1, 1, 0];

            % Note that parameters errors are evaluated on the uncostrained
            % model, so are generally wrong when there are constraints (but
            % not when fixing parameters)
            [pars, errs] = test_GaussianParametersErrors(10, 500, 1e-2, ...
                [NaN, NaN, NaN, NaN], Aeq, beq);

            [~, errs_free] = test_GaussianParametersErrors(10, 500, 1e-2);
            tc.assertLessThan(abs(pars - [1 0 1 0]), 1e-3);
            tc.assertLessThan(abs(errs(1,:) - errs_free(1,:)), 3e-3);
        end

    end  % Methods
end