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

    end
end