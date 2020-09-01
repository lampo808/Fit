classdef TestGFS < matlab.unittest.TestCase

    methods (Test)

        function GFSLinearFit(tc)
            rng(1)  % Set a seed for reproducibility

            model = @(x, p) p(1)*x + p(2);  % Linear function

            % Parameters for the simulated datasets
            % Four datasets (one row per dataset): shifted lines
            pars = [0.3, -0.2;
                    0.3, 2;
                    0.3, 3.4;
                    0.3, 1.7];

            N = 20;  % Points per line
            noise = 0.01;  % Absolute noise level

            % Generate data

            for i=1:size(pars, 1)
                xData{i} = linspace(-2, 5, N);
                yData{i} = model(xData{i}, pars(i,:)) + noise*randn(size(xData{i}));
            end

            gf = GlobalFitSimple();
            gf.setData(xData, yData);
            gf.setModel(model, 2, [1 0])
            gf.setStart(pars)
            gf.fit()
            fit_pars = gf.getFittedParameters();
            % Errors are not used, this is to test that no error is raised
            fit_errs = gf.getParamersErrors();

            for i=1:length(pars(:))
                tc.assertEqual(pars(i), fit_pars(i), 'RelTol', 0.1)
            end
        end

        function GFSMultipleGaussiansFree(tc)
            [pars, fit_pars] = test_GFS_MultipleGaussians();

            for i=1:length(pars(:))
                tc.assertEqual(pars(i), fit_pars(i), 'AbsTol', 0.1)
            end
        end

%         pars =

%   Columns 1 through 8

%     0.0360    2.0000    1.0000    0.5892    3.0000    2.0000    0.0180    6.0000
%     1.0000    2.0000    1.0000    0.5285    3.0000    2.0000    0.1192    6.0000
%     1.9640    2.0000    1.0000    0.3000    3.0000    2.0000    0.5000    6.0000
%     1.9993    2.0000    1.0000    0.0715    3.0000    2.0000    0.8808    6.0000
%     2.0000    2.0000    1.0000    0.0108    3.0000    2.0000    0.9820    6.0000

%   Column 9

%     0.5000
%     0.5000
%     0.5000
%     0.5000
%     0.5000
        function GFSMultipleGaussiansBound(tc)
            % Set lower and upper bounds far from the real parameters
            %  and check that the fit is fine
            for i=1:5
                lb(i,:) = [0, 1, 0, 0, 2, 1, 0, 5, 0];
                ub(i,:) = [3, 3, 2, 1, 4, 3, 1, 7, 1];
            end

            [pars, fit_pars] = test_GFS_MultipleGaussians(lb, ub);

            for i=1:length(pars(:))
                tc.assertEqual(pars(i), fit_pars(i), 'AbsTol', 0.1)
                tc.assertGreaterThanOrEqual(pars(i), lb(i))
                tc.assertLessThanOrEqual(pars(i), ub(i))
            end
        end

        function GFSMultipleGaussiansBoundErrors(tc)
            % Check errors on parameters
            for i=1:5
                lb(i,:) = [0, 1, 0, 0, 2, 1, 0, 5, 0];
                ub(i,:) = [3, 3, 2, 1, 4, 3, 1, 7, 1];
            end

            for i=1:50
                [~, fit_pars(:,:,i), fit_errs(:,:,i)] = ...
                    test_GFS_MultipleGaussians(lb, ub, [], i);
            end

            % Check that the errors from the fit agree with repeated
            % simulations with a reasonably large range - there is a
            % discrepancy on some parameters, hence the large range 0.3-3
            tc.assertLessThan(abs(mean(fit_errs, 3) - std(fit_pars, [], 3))./...
                mean(fit_errs, 3), 2);
            tc.assertLessThan(std(fit_pars, [], 3)./mean(fit_errs, 3), 3)
            tc.assertGreaterThan(std(fit_pars, [], 3)./mean(fit_errs, 3), 0.3)
        end

        function GFSMultipleGaussiansFixed(tc)
            % Fix some parameters to the start value and check that
            %  1. the fit converges
            %  2. the fixed parameters are actually fixed

            % Get the initial parameters
            pars = test_GFS_MultipleGaussians();

            for i=1:5
                lb(i,:) = [i/2-1, 1, 0, 0, 2, 1, 0, 5, 0];
                ub(i,:) = [3, 3, 2, 1, 4, 3, 1, 7, 1];
                fixed(i,:) = [NaN, NaN, pars(i,3), NaN, NaN, NaN, NaN, NaN, NaN];
            end

            [pars, fit_pars] = test_GFS_MultipleGaussians(lb, ub, fixed);

            tc.assertTrue(all(all(fit_pars > lb)))
            tc.assertTrue(all(all(fit_pars < ub)))

            for i=1:length(pars(:))
                tc.assertEqual(pars(i), fit_pars(i), 'AbsTol', 0.1)
                if ~isnan(fixed(i))
                    tc.assertEqual(fit_pars(i), fixed(i), 'AbsTol', 1e-6)
                end
            end
        end
        
        function GFSMultipleGaussiansScaling(tc)
            % Generate badly scaled data, and check that with offset and
            % scaling added the fit converges (it wouldn't converge
            % normally)

            for i=1:5
                offset(i,:) = [0, 2e6, 0, ...
                    0, 2e6, 0, ...
                    0, 2e6, 0];
                scaling(i,:) = [1, 1, 1e-10, ...
                    1, 1, 1e-10, ...
                    1, 1, 1e-10];
            end

            [pars, fit_pars] = test_GFS_MultipleGaussiansRescale(...
                offset, scaling);

            for i=1:numel(pars)
                tc.assertEqual(pars(i), fit_pars(i), 'AbsTol', 0.1)
            end
        end

    end  % Methods
end


