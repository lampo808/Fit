classdef Fit < handle
    % Fit: generic class for constrained fitting using fmincon()
    % v. 0.3.0
    %
    % Changelog
    %   30/06/20 - 0.3.0: added offset and scaling of parameters through
    %       additional parameters passed to setStart()
    %   26/06/20 - 0.2.6: initial step to enable Levenberg-Marquardt minimization:
    %       minFunXXX return also the gradient. Bug fixes and Hessian saved
    %   22/01/20 - 0.2.5: fixed errors on parameters when there is external
    %       chi2 function
    %   14/01/20: changed the way of fixing parameters, now without using
    %       equality constraints
    %   13/01/20: fixed errors on fitted parameters, plus slight cleanup
    %   10/10/19: fixed parHistory not cleared on new fit; added the
    %       possibility to ignore fit warnings
    %   15/11/18: the weights are automatically accounted in the
    %       calculation of the residuals; getChiSquare() returns either the
    %       full or reduced chi square
    %   01/03/16: added fit error estimation via Hessian estimation
    %   15/02/16: introduced dataMask; pretty much everything tested in
    %       real life
    %   04/06/15: added parameter and chi square history
    %   29/05/15: added setChiSquare()
    %   21/05/15 - 0.1: first draft, not complete, but working
    %

    %TODO:
    %  - [ ] expand hessian and covariance matrix
    %  - [ ] offset and scaling does not work with equality/inequality
    %  constraints
    %  - [ ] thorough test suite
    %      - [ ] test that the residuals are calculated the right way
    %      - [ ] test that the chi square is right even before the fit
    %  - [ ] doc: functions
    %  - [ ] doc: class description and examples
    %  - [ ] test + doc: write example functions used for testing purposes
    %  - [ ] possibility not to have periodic convolution
    %  - [x] if IRF is longer than xData, lengthen xData but fit only on
    %  "real" points
    %  - [x] make the convoluted fit work even with 2+ xData intervals
    %  - [ ] optimize the convolution embedding the shift in the class
    %  - [x] add history of parameters
    %  - [x] if fit is stopped through Ctrl+C, fitPar_ get the last value
    %  of the history
    %  - [x] user can provide chi-square function instead of model
    %  - [x] user-friendly global fit
    %  - [ ] set good default options for all the minimizers and algorithms
    %  - [ ] enable the use of Trust-region-reflective by automatically
    %           supplying gradient
    %  - [ ] At the moment, fitEval is limited to the IRF points. Extend.
    %  - [ ] Add parameters offset for big parameters with small variations
    %  - [ ] Add parameters multiplier to equalize the size of the parameters

    properties (GetAccess = public, SetAccess = private)
        xData_       = [];  % Data
        yData_       = [];
        weights_     = [];  % Weights

        dataMask_    = [];  % Mask to subselect fit data

        IRF_         = [];  % Instrument Response Function
        xIRF_        = [];

        fitLength_   = [];  % Number of fit points if IRF is passed

        model_       = [];  % Fit model
        nparam_      = [];  % Number of parameters
        opt_         = [];  % Optimization options

        chi2Func_    = [];  % User-defined chi-square function

        start_       = [];  % Start point
        offset_      = [];  % Offset of the parameters
        scaling_     = [];  % Scaling of the parameters
        ub_          = [];  % Upper bound
        lb_          = [];  % Lower bound
        Aeq_         = [];  % Equality matrix
        beq_         = [];  % Equality vector
        A_           = [];  % Inequality matrix
        b_           = [];  % Inequality vector
        fixed_       = [];  % Vector flag of fixed parameters

        fitPar_      = [];  % Output fit parameters
        hessian_     = [];  % Hessian of the chi2 after the fit
        covmat_      = [];  % Covariance matrix
        fitParError_ = [];  % Errors on the fit parameters
        chi2_        = [];  % Chi square of the fit

        parHistory_  = [];  % History of the parameters during minimization
        chi2History_ = [];  % History of the chi square values

        fitStatus_   = 0;   %   0  --> fit not yet performed
                            % 100  --> fit running
                            % else --> exit flag of the completed fit

        ignoreFitWarnings_ = 0;  % If left 0, returns NaN fit parameters upon warning

    end  % Properties

    % List of methods:  (" " are optional parameters)
    %  - Fit("xData, yData", "weights")        --> constructor
    %  - setData(xData, yData, "weights")      --> set data for the fitting
    %  - setDataMask(mask)                     --> set a mask on experimental data
    %  - setDataMaskByIndices(ind)             --> set a mask given the indices
    %  - setWeights(weights)                   --> set y weights
    %  - setIRF(irf)                           --> set IRF
    %  - setModel(@model, npars, "reset")      --> set fit model
    %  - setChiSquare(@chi2, npars, "reset")   --> set chi-square function
    %  - setStart(sp)                          --> set start point
    %  - getStart()                            --> get start point
    %  - setUb(ub)                             --> set upper bounds
    %  - setLb(lb)                             --> set lower bounds
    %  - setEqualityConstraints(Aeq, beq)      --> set equality matrix and vector
    %  - setInequalityConstraints(A, b)        --> set inequality matrix and vector
    %  - addEqualityConstraints(Arows, bels)   --> add equality constraint
    %  - addInequalityConstraints(Arows, bels) --> add inequality constraint
    %  - fixParameter(i, bel)                  --> equality constraint for parameter i
    %  - setParUb(i, ub)                       --> upper bound for parameter i
    %  - setParLb(i, lb)                       --> lower bound for parameter i
    %  - removeEqualityConstraints("i")        --> remove all/"i-th" equality constraint
    %  - removeInequalityConstraint("i")       --> remove all/"i-th" inequality constraint
    %  - setFminconOptions(opt)                --> set fmincon options
    %  - getFminconOptions()                   --> get fmincon options
    %  - fit("calculateErrors")                                 --> perform the fitting "and calculate the errors on fit parameters"
    %  - getChiSquare()                        --> get the chi square
    %  - getFittedParameters("i")              --> get all/"the i-th" fit parameter
    %  - getParametersErrors("i")              --> get all/"the i-th" fit parameter error
    %  - getResiduals()                        --> get fit residuals
    %  - fitEval("x")                          --> evaluate fit function on xData_/"x"
    %  - getParHistory()                       --> get the parameter history
    %  - getChi2History()                      --> get the chi square history
    %  - resetHistory()                        --> reset all the histories
    methods

        function F = Fit(xData, yData, weights)  % Constructor
            if nargin == 1
                msgID = 'FIT:Constructor_argumentsNumber';
                msg = 'Wrong number of arguments passed to constructor. Zero, two or three allowed.';
                exception = MException(msgID, msg);

                throw(exception);
            end
            if nargin == 2
                F.setData(xData, yData);
            end
            if nargin == 3
                F.setData(xData, yData, weights);
            end

            F.opt_ = optimoptions(@fmincon);
            F.opt_.MaxFunEvals = 5000;
            F.opt_.TolX = 1e-12;
            F.opt_.Display = 'none';
            F.opt_.OutputFcn = @F.updateHistory;
        end


        function setData(F, xData, yData, weights)  % Set data for the fitting
            if nargin < 3
                msgID = 'FIT:setData_argumentsNumber';
                msg = 'Wrong number of arguments passed to setData; expects at least two.';
                exception = MException(msgID, msg);

                throw(exception);
            end

            if length(xData) ~= length(yData)
                msgID = 'FIT:setData_argumentsLength';
                msg = 'xData and yData must have the same length.';
                exception = MException(msgID, msg);

                throw(exception);
            end

            if isempty(xData)
                msgID = 'FIT:setData_argumentsLength';
                msg = 'xData and yData cannot be empty.';
                exception = MException(msgID, msg);

                throw(exception);
            end

            F.IRF_ = [];

            F.xData_ = xData(:);
            F.yData_ = yData(:);

            F.setDataMask(ones(size(F.xData_)));

            if nargin == 4
                F.setWeights(weights);
            else
                F.setWeights(ones(size(F.xData_)));
            end
        end


        function setDataMask(F, mask)
            if isempty(F.xData_) || isempty(F.yData_)
                msgID = 'FIT:setDataMask_dataMissing';
                msg = 'A mask can be set only if xData and yData are already set.';
                exception = MException(msgID, msg);

                throw(exception);
            end
            if size(F.xData_) ~= size(mask)
                msgID = 'FIT:setDataMask_maskLength';
                msg = 'mask must be as long as xData and yData.';
                exception = MException(msgID, msg);

                throw(exception);
            end

            F.dataMask_ = (mask ~= 0);
            F.fitLength_ = sum(F.dataMask_);
        end


        function setDataMaskByIndices(F, ind)
            if any(ind <= 0)
                msgID = 'FIT:setDataMaskByIndices_nonPositiveIndex';
                msg = 'Indices must be strictly positive.';
                exception = MException(msgID, msg);

                throw(exception);
            end
            if any(ind > length(F.xData_))
                msgID = 'FIT:setDataMaskByIndices_indexTooBig';
                msg = 'Indices cannot be greater than length(xData).';
                exception = MException(msgID, msg);

                throw(exception);
            end

            mask = zeros(size(F.xData_));
            mask(ind) = 1;

            F.setDataMask(mask);
        end


        function setWeights(F, weights)  % Set y weights
            if isempty(F.xData_) || isempty(F.yData_)
                msgID = 'FIT:setWeights_dataMissing';
                msg = 'weights can be set only if xData and yData are already set.';
                exception = MException(msgID, msg);

                throw(exception);
            end
            if length(F.xData_) ~= length(weights)
                msgID = 'FIT:setWeights_weightsLengths';
                msg = 'weights must be as long as xData and yData.';
                exception = MException(msgID, msg);

                throw(exception);
            end
            if any(weights <= 0)
                msgID = 'FIT:setWeights_weightsValues';
                msg = 'weights must be strictly positive.';
                exception = MException(msgID, msg);

                throw(exception);
            end
            if any(isinf(weights))
                msgID = 'FIT:setWeights_weightsValues';
                msg = 'weights must be finite.';
                exception = MException(msgID, msg);

                throw(exception);
            end

            F.weights_ = weights(:);
        end


        function setIRF(F, irf)  % set IRF
            % setIRF(irf)
            % Set the Instrument Response Function the data will be
            %  convoluted with. setIRF() must be called after every change
            %  of data and/or model(i.e. every call of setData or setModel
            %  resets the IRF).
            %
            % Input:
            %  <irf>: y-data of the IRF. The x coordinates are made to
            %    correspond to xData. If length(irf) ~= length(xData), IRF
            %    is either cut or zero-padded to match the length of xData.
            %
            if isempty(F.xData_) || isempty(F.yData_)
                msgID = 'FIT:setIRF_dataNotSet';
                msg = 'xData and yData must be set prior to setting the IRF.';
                exception = MException(msgID, msg);

                throw(exception);
            end

            l  = length(F.xData_);
            lI = length(irf);
            deltaL = lI - l;
            if deltaL > 0  % IRF longer than xData
                F.IRF_  = irf(:);
                dx = mean(diff(F.xData_));  % x step
                missing = F.xData_(l) + cumsum(ones(deltaL, 1) * dx);
                F.xIRF_ = [F.xData_; missing];
            else
                F.IRF_  = [irf(:); zeros(l-length(irf))];
                F.xIRF_ = F.xData_;
            end

            F.fitLength_ = length(F.xData_);
        end


        function setModel(F, model, npars, reset)  % Set fit model
            %TODO: documentation: if using convolution with IRF, the first
            %parameter is _always_ the shift
            if nargin < 3
                msgID = 'FIT:setModel_nargin';
                msg = 'Pass at least two parameters: model and number of parameters.';
                exception = MException(msgID, msg);

                throw(exception);
            end
            if nargin < 4
                reset = 1;
            end

            if reset  % Reset old parameters
                F.start_       = [];
                F.ub_          = [];
                F.lb_          = [];
                F.Aeq_         = [];
                F.beq_         = [];
                F.A_           = [];
                F.b_           = [];
                F.fixed_       = false(1, npars);
            end

            F.fitPar_      = [];  % Fitted parameters are always reset
            F.fitParError_ = [];
            F.chi2_        = [];
            F.chi2Func_    = [];
            F.parHistory_  = [];

            F.model_ = model;
            F.nparam_ = npars;
        end


        function setChiSquare(F, chi2, npars, reset)
            %TODO: documentation: if using convolution with IRF, the first
            %parameter is _always_ the shift
            if nargin < 3
                msgID = 'FIT:setChiSquare_nargin';
                msg = 'Pass at least two parameters: model and number of parameters.';
                exception = MException(msgID, msg);

                throw(exception);
            end
            if nargin < 4
                reset = 1;
            end

            if reset  % Reset old parameters
                F.start_       = [];
                F.ub_          = [];
                F.lb_          = [];
                F.Aeq_         = [];
                F.beq_         = [];
                F.A_           = [];
                F.b_           = [];
                F.fixed_       = false(1, npars);
            end

            F.fitPar_      = [];  % Fitted parameters are always reset
            F.fitParError_ = [];
            F.chi2_        = [];
            F.model_       = [];

            F.IRF_         = [];  % If chi-square function is passed, the model
                                  %  evaluation is not performed by Fit

            F.chi2Func_ = chi2;
            F.nparam_ = npars;
        end


        function setStart(F, sp, offset, scaling)  % Set start point
            if (nargin < 3) || isempty(offset)
                offset = zeros(size(sp));
            end
            if (nargin < 4) || isempty(scaling)
                scaling = ones(size(sp));
            end

            if length(offset) ~= F.nparam_
                msgID = 'FIT:setStart_numberParams';
                msg = 'The number of elements of offset vector is different from the number of parameters.';
                exception = MException(msgID, msg);

                throw(exception);
            else
                F.offset_ = offset(:)';
            end

            if length(scaling) ~= F.nparam_
                msgID = 'FIT:setStart_numberParams';
                msg = 'The number of elements of scaling vector is different from the number of parameters.';
                exception = MException(msgID, msg);

                throw(exception);
            else
                F.scaling_ = scaling(:)';
            end

            if length(sp) ~= F.nparam_
                msgID = 'FIT:setStart_numberParams';
                msg = 'The number of elements of startpoint vector is different from the number of parameters.';
                exception = MException(msgID, msg);

                throw(exception);
            else
                F.start_ = sp(:)';
            end
        end


        function sp = getStart(F)
            sp = F.start_;
        end

        function offset = getOffset(F)
            offset = F.offset_;
        end

        function scaling = getScaling(F)
            scaling = F.scaling_;
        end

        function setUb(F, ub)  % Set upper bounds
            if length(ub) ~= F.nparam_
                msgID = 'FIT:setUb_numberParams';
                msg = 'The number of elements of upper bounds vector is different from the number of parameters.';
                exception = MException(msgID, msg);

                throw(exception);
            else
                F.ub_ = ub(:)';
            end
        end


        function setLb(F, lb)  % Set lower bounds
             if length(lb) ~= F.nparam_
                msgID = 'FIT:setLb_numberParams';
                msg = 'The number of elements of lower bounds vector is different from the number of parameters.';
                exception = MException(msgID, msg);

                throw(exception);
            else
                F.lb_ = lb(:)';
            end
        end


        function setEqualityConstraints(F, Aeq, beq)  % Set equality matrix and vector
            s = size(Aeq);
            if s(2) ~= F.nparam_
                msgID = 'FIT:setEqualityConstraints_AeqSize';
                msg = 'The number of columns of Aeq is different from the number of parameters.';
                exception = MException(msgID, msg);

                throw(exception);
            end
            if s(1) ~= length(beq)
                msgID = 'FIT:setEqualityConstraints_AeqSize';
                msg = 'The number of rows of Aeq is different from the length of beq.';
                exception = MException(msgID, msg);

                throw(exception);
            end

            F.Aeq_ = Aeq;
            F.beq_ = beq(:);
        end


        function setInequalityConstraints(F, A, b)  % Set inequality matrix and vector
            s = size(A);
            if s(2) ~= F.nparam_
                msgID = 'FIT:setInequalityConstraints_ASize';
                msg = 'The number of columns of A is different from the number of parameters.';
                exception = MException(msgID, msg);

                throw(exception);
            end
            if s(1) ~= length(b)
                msgID = 'FIT:setInequalityConstraints_ASize';
                msg = 'The number of rows of A is different from the length of b.';
                exception = MException(msgID, msg);

                throw(exception);
            end

            F.A_ = A;
            F.b_ = b(:);
        end


        function addEqualityConstraints(F, Arows, bels)  % Add equality constraint
            s = size(Arows);
            if s(2) ~= F.nparam_
                msgID = 'FIT:addEqualityConstraints_ArowsSize';
                msg = 'The number of columns of Arows is different from the number of parameters.';
                exception = MException(msgID, msg);

                throw(exception);
            end
            if s(1) ~= length(bels)
                msgID = 'FIT:addEqualityConstraints_ArowsSize';
                msg = 'The number of rows of Arows is different from the length of bels.';
                exception = MException(msgID, msg);

                throw(exception);
            end

            F.Aeq_ = [F.Aeq_; Arows];
            F.beq_ = [F.beq_; bels(:)];
        end


        function addInequalityConstraints(F, Arows, bels)  % Add inequality constraint
            s = size(Arows);
            if s(2) ~= F.nparam_
                msgID = 'FIT:addInequalityConstraints_ArowsSize';
                msg = 'The number of columns of Arows is different from the number of parameters.';
                exception = MException(msgID, msg);

                throw(exception);
            end
            if s(1) ~= length(bels)
                msgID = 'FIT:addInequalityConstraints_ArowsSize';
                msg = 'The number of rows of Arows is different from the length of bels.';
                exception = MException(msgID, msg);

                throw(exception);
            end

            F.A_ = [[F.A_]; Arows];
            F.b_ = [[F.b_]; bels(:)];
        end

        function fixParameters(F, fixedValues)
            % fixParameters(F, fixedValues)
            % Set fixed values for multiple parameters at once
            %
            % Input:
            %   <fixedValues>: vector of fixed values.
            %              Long as the number of parameters, NaNs identify parameters whose values are not fixed,
            %              otherwise the value identify the parameter's fixed value
                for i=1:length(fixedValues)
                    if ~isnan(fixedValues(i))
                        F.fixParameter(i, fixedValues(i));
                    end
                end
            end

        function fixParameter(F, i, value)
            % fixParameter(i, value)
            % Set a fixed value <value> for parameter number <i>.
            %
            % Input:
            %  <i>: number of the parameter whose parameter will be fixed.
            %  <value>: value of the parameter.

            if (i < 1) || (i > F.nparam_)
                msgID = 'FIT:fixParameter_iValue';
                msg = 'i is less than 1 or greater than number of parameters.';
                exception = MException(msgID, msg);

                throw(exception);
            end

            %%% Old code
%             % Set lower bound, upper bound and constraint equal to <value>
%             F.setParUb(i, value);
%             F.setParLb(i, value);
%
%             vect = zeros(1, F.nparam_);
%             vect(i) = 1;
%             F.addEqualityConstraints(vect(:)', value);
            %%% New code
            F.fixed_(i) = true;

            % Set initial point
            F.start_(i) = value;
        end


        function setParUb(F, i, ub)
            % setParUb(i, ub)
            % Set an upper bound <ub> for parameter number <i>.
            %
            % Input:
            %  <i>: number of the parameter whose upper bound will be set.
            %  <ub>: value of the upper bound.

            if (i < 1) || (i > F.nparam_)
                msgID = 'FIT:setParUb_iValue';
                msg = 'i is less than 1 or greater than number of parameters.';
                exception = MException(msgID, msg);

                throw(exception);
            end

            F.ub_(i) = ub;
        end


        function setParLb(F, i, lb)
            % setParLb(i, lb)
            % Set a lower bound <lb> for parameter number <i>.
            %
            % Input:
            %  <i>: number of the parameter whose lower bound will be set.
            %  <lb>: value of the lower bound.

            if (i < 1) || (i > F.nparam_)
                msgID = 'FIT:setParLb_iValue';
                msg = 'i is less than 1 or greater than number of parameters.';
                exception = MException(msgID, msg);

                throw(exception);
            end

            F.lb_(i) = lb;
        end


        function removeEqualityConstraints(F, i)
            %%% TODO: write method

            if nargin == 2
            else
            end
        end


        function removeInequalityConstraint(F, i)
            %%% TODO: write method
            if nargin == 2
            else
            end
        end


        function setFminconOptions(F, opt)
            F.opt_ = opt;
        end


        function opt = getFminconOptions(F)
            opt = F.opt_;
        end


        function fit(F, calculateErrors, ignoreFitWarnings)
            if nargin < 2
                calculateErrors = 1;
            end
            if nargin < 3
                F.ignoreFitWarnings_ = 0;
            else
                F.ignoreFitWarnings_ = ignoreFitWarnings;
            end
            %%% TODO: check conditions
            if isempty(F.model_) && isempty(F.chi2Func_)
                msgID = 'FIT:fit_noModelSet';
                msg = 'A model is needed to perform the fit.';
                exception = MException(msgID, msg);

                throw(exception);
            end
            if isempty(F.xData_) || isempty(F.yData_)
                msgID = 'FIT:fit_noData';
                msg = 'No data to fit.';
                exception = MException(msgID, msg);

                throw(exception);
            end
            if isempty(F.start_)
                msgID = 'FIT:fit_noStartPoint';
                msg = 'A start point is needed to initialize the fit.';
                exception = MException(msgID, msg);

                throw(exception);
            end

            F.fitStatus_ = 100;  % Fit started
            cleanupObj = onCleanup(@F.fitInterrupted);

            if isempty(F.chi2Func_)
                if isempty(F.IRF_)
                    minFun = @(par) F.minFun(F.expandFixedPars(par));
                else
                    minFun = @(par) F.minFunConv(F.expandFixedPars(par));
                end
            else
                minFun = @(par) F.extChi2Fun(F.expandFixedPars(par));
            end

            % If no constraint is set, fall back to fminunc
            if (isempty(F.A_) && isempty(F.b_) && isempty(F.Aeq_) && ....
                isempty(F.beq_) && isempty(F.lb_) && isempty(F.ub_))
                [fitted, chi2, exitflag] = fminunc(minFun, ...
                        F.reduceFixedPars(F.start_), ...
                        optimoptions('fminunc', F.opt_));
            else
                if any(F.fixed_)
                    [fitted, chi2, exitflag] = fmincon(minFun, ...
                                F.reduceFixedPars(F.start_), [], [], ...
                                [], [], F.reduceFixedPars(F.lb_), ...
                                F.reduceFixedPars(F.ub_), [], F.opt_);
                else
                    % At the moment, fixing parameters is incompatible with
                    % linear constraints/inequalities
                    [fitted, chi2, exitflag] = fmincon(minFun, ...
                                F.rescalePars(F.start_), ...
                                F.A_, F.b_, ...   %%% TODO: rescale A and b
                                F.Aeq_,F.beq_, ...
                                F.rescalePars(F.lb_), F.rescalePars(F.ub_), [], F.opt_);
                end
            end

            F.fitPar_ = F.expandFixedPars(fitted(:)');
            F.chi2_ = chi2;

            if calculateErrors
                %%%%%%%%% IMPORTANT!!!! %%%%%%%%%
                % NOTE: the errors are adjusted for chi2 = 1; this is the
                % best estimation if the errors on the fit parameters are
                % unknow, otherwise use the correct weights (= 1/sigma^2)
                % and divide the estimated errors by sqrt(F.getChiSquare())
                % (or use bootstrap)
                hess = hessian(minFun, fitted);
%                 for i=1:size(hess,1)  %%% TODO: expand hessian first
%                     for j=1:size(hess,2)
%                         hess(i,j) = hess(i,j)/(F.scaling_(i) * F.scaling_(j));
%                     end
%                 end

                if isempty(F.model_)  % Number of data points not available
                    % In this case the error must be manually corrected for
                    % the value of the reduced chi square
                    err = sqrt(2*diag(inv(hess)));
                else
                    % The factor 2 has been checked experimentally
                    err = sqrt(2*diag(inv(hess))*F.getChiSquare(1));
                end

                F.covmat_ = inv(hess);
                F.hessian_ = hess;  % TODO: expand the hessian and covmat
                F.fitParError_ = F.expandFixedErrors(err(:)');
            end

            if ~F.ignoreFitWarnings_
                switch exitflag
                    case 1
                        % Fit converged
                    case 0
                        warning('Fit stopped: exceeded MaxIter or MaxFunEval.');
                        F.fitPar_ = NaN(size(F.fitPar_));
                    case -1
                        warning('Fit stopped: stopped by output function or plot function');
                        F.fitPar_ = NaN(size(F.fitPar_));
                    case -2
                        warning('Fit stopped: no feasible point found.');
                         F.fitPar_ = NaN(size(F.fitPar_));
                    case 2
                        warning('Fit stopped: change in x better than TolX.');
                    case -3
                        warning('Fit stopped: objective function less than ObjectiveLimit.');
                    otherwise
                        warning(['Fit stopped: flag ', num2str(exitflag), ...
                            ' (generally not a problem if the flag is positive)']);
                end
            end
            F.fitStatus_ = exitflag;
        end


        function chi2 = getChiSquare(F, reduced)
            if nargin < 2
                reduced = 0;
            end

            if (reduced && isempty(F.chi2Func_))
                chi2 = F.chi2_/(F.fitLength_ - F.nparam_ + sum(F.fixed_) + ...
                size(F.beq_, 2));
            else
                chi2 = F.chi2_;
            end
        end


        function pars = getFittedParameters(F, i)
            if isempty(F.fitPar_)
                msgID = 'FIT:getFittedParameters_fitNotPerformed';
                msg = 'Fit not yet performed or not converged.';
                exception = MException(msgID, msg);

                throw(exception);
            end
            if nargin == 2
                if i <= length(F.fitPar_)
                    pars = F.fitPar_(i);
                else
                    msgID = 'FIT:getFittedParameters_wrongIndex';
                    msg = 'i greater than number of parameters.';
                    exception = MException(msgID, msg);

                    throw(exception);
                end
            else
                pars = F.fitPar_;
            end
        end


        function ers = getParametersErrors(F, i)
        % Return the errors estimated on a non-constrained Hessian
        % NOTE: the errors are normalized to the chi square

            if isempty(F.fitPar_)
                msgID = 'FIT:getParametersErrors_fitNotPerformed';
                msg = 'Fit not yet performed or not converged.';
                exception = MException(msgID, msg);

                throw(exception);
            elseif isempty(F.fitParError_)
                msgID = 'FIT:getParametersErrors_ErrorsNotCalculated';
                msg = 'Fit has been performed without error calculation.';
                exception = MException(msgID, msg);

                throw(exception);
            end

            if nargin == 2
                if i <= length(F.fitParError_)
                    ers = F.fitParError_(i);
                else
                    msgID = 'FIT:getParametersErrors_wrongIndex';
                    msg = 'i greater than number of parameters.';
                    exception = MException(msgID, msg);

                    throw(exception);
                end
            else
                ers = F.fitParError_;
            end
        end


        function res = getResiduals(F)
            if isempty(F.fitPar_)
                msgID = 'FIT:getResiduals_fitNotPerformed';
                msg = 'Before evaluating residuals the fit must be performed.';
                exception = MException(msgID, msg);

                throw(exception);
            end

            if isempty(F.chi2Func_)
                if isempty(F.IRF_)
                    yFit = F.fitEval();
                    res = F.residuals(yFit);
                else
                    yFit = F.fitEval();
                    res = F.residuals(yFit(1:F.fitLength_));
                end
            else
                msgID = 'FIT:getResiduals_externalChi2Eval';
                msg = 'Residuals cannot be computed when an external chi-square function is provided.';
                exception = MException(msgID, msg);

                throw(exception);
            end
        end


        function [y, xIRF] = fitEval(F, x)
            if isempty(F.fitPar_)
                par = F.start_;
%                 warning('Fit not performed; evaluation on start point');
            else
                par = F.fitPar_;
            end

            if isempty(F.IRF_)
                if nargin == 2
                    y = F.model_(x(:), par);
                else
                    y = F.fitEval(F.xData_);
                end

                xIRF = [];
            else  % if the IRF is passed, fit can only be evaluated on xIRF_
                yFit = F.model_(F.xIRF_, par(2:end));
                % Convolve fit points with IRF and normalize
                yFit_c = F.altConv(yFit, F.IRF_, par(1));
                y = yFit_c * sum(yFit)/sum(yFit_c);

                xIRF = F.xIRF_;
            end
        end


        function history = getParHistory(F)
            history = F.parHistory_;
        end


        function history = getChi2History(F)
            history = F.chi2History_;
        end


        function resetHistory(F)
            F.parHistory_ = [];
            F.chi2History_ = [];
        end

    end  % Public methods

    methods (Access = private)

        function c = altConv(F, x, y, s)  % Redefine convolution as product of FT
            c = real(ifft(fft(fshift(x, s)) .* fft(y)));
        end

        function res = residuals(F, yFit)
            res = (F.yData_ - yFit) .* F.dataMask_ .* sqrt(F.weights_);
        end


        function chi2 = chisquare(F, res)
            chi2 = sum(abs(res).^2);
        end


% TODO: check the use of xDataMask on all minFun variants
        function [chi2] = minFun(F, par)
            yFit = F.model_(F.xData_, par);
            res = F.residuals(yFit);
            chi2 = F.chisquare(res);
        end

        function [chi2, grad] = minFunGrad(F, par)
            chi2 = F.minFun(par);

            chi2f = @(p) F.chisquare(F.residuals(F.model_...
                (F.xData_, F.expandFixedPars(p))));
            grad = gradest(chi2f, F.reduceFixedPars(par));
        end

        function chi2 = minFunConv(F, par)
            yFit = F.model_(F.xIRF_, par(2:end));
            % Convolve fit points with IRF and normalize
            yFit_c = F.altConv(yFit, F.IRF_, par(1));
            yFit_c = yFit_c * sum(yFit)/sum(yFit_c);
            res = F.residuals(yFit_c(1:F.fitLength_));
            chi2 = F.chisquare(res);
        end

        function [chi2, grad] = minFunConvGrad(F, par)
            chi2 = F.minFunConv(par);

            grad = gradest(chi2f, F.reduceFixedPars(par));
            function chi2 = chi2f(p)
                par = F.expandFixedPars(p);
                yFit = F.model_(F.xIRF_, par(2:end));
                % Convolve fit points with IRF and normalize
                yFit_c = F.altConv(yFit, F.IRF_, par(1));
                yFit_c = yFit_c * sum(yFit)/sum(yFit_c);
                res = F.residuals(yFit_c(1:F.fitLength_));
                chi2 = F.chisquare(res);
            end
        end

        function [chi2, grad] = extChi2Fun(F, par)
            xDataMasked = F.xData_(F.dataMask_);
            chi2 = F.chi2Func_(xDataMasked, par);
        end

        function [chi2, grad] = extChi2FunGrad(F, par)
            xDataMasked = F.xData_(F.dataMask_);
            chi2 = F.extChi2Fun(par);

            chi2f = @(p) F.chi2Func_(xDataMasked, F.expandFixedPars(p));
            grad = gradest(chi2f, F.reduceFixedPars(par));
        end

        function stop = updateHistory(F, x, optimValues, state)
            stop = false;

            F.parHistory_ = [F.parHistory_; x(:)'];
            F.chi2History_ = [F.chi2History_; optimValues.fval];
        end

        function fitInterrupted(F)
            if F.fitStatus_ == 100
                F.fitPar_ = F.parHistory_(end,:);
                disp('')
                disp('-------------------------------------------')
                disp('Last parameters value:')
                disp(F.fitPar_)
                disp('-------------------------------------------')
                disp('')
            end
        end

        function reducedPars = reduceFixedPars(F, fullPars)
            fullPars = F.rescalePars(fullPars);
            reducedPars = fullPars(~F.fixed_);
        end

        function fullPars = expandFixedPars(F, reducedPars)
            fullPars = F.rescalePars(F.start_);
            fullPars(~F.fixed_) = reducedPars;
            fullPars = F.unscalePars(fullPars);
        end

        function rescaledPars = rescalePars(F, unscaledPars)
            rescaledPars = (unscaledPars - F.offset_) ./ F.scaling_;
        end

        function unscaledPars = unscalePars(F, rescaledPars)
            unscaledPars = rescaledPars .* F.scaling_ + F.offset_;
        end

        function fullErrors = expandFixedErrors(F, reducedErrors)
            fullErrors = zeros(size(F.start_));
            fullErrors(~F.fixed_) = reducedErrors;
            fullErrors = fullErrors .* F.scaling_;
        end

    end  % Private methods

end  % Class