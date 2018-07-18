classdef GlobalFitSimple < handle
    % GlobalFitSimple: class for global fitting of data with the same equation
    %  and single or global parameters
    % v. 0.1

    % Changelog:
    %   05/07/2018 - 0.1: first version

    properties (GetAccess = public, SetAccess = private)
    xData_       = {};  % Data
    yData_       = {};
    weights_     = {};  % Weights

    model_       = [];  % Fit model
    nParam_      = [];  % Number of parameters
    opt_         = [];  % Optimization options

    start_       = [];  % Start point
    ub_          = [];  % Upper bound
    lb_          = [];  % Lower bound
    Aeq_         = [];  % Equality matrix
    beq_         = [];  % Equality vector
    A_           = [];  % Inequality matrix
    b_           = [];  % Inequality vector

    fitPar_      = [];  % Output fit parameters
    fitParError_ = [];  % Errors on the fit parameters
    chi2_        = [];  % Chi square of the fit

    parHistory_  = [];  % History of the parameters during minimization
    chi2History_ = [];  % History of the chi square values

    fitStatus_   = 0;   %   0  --> fit not yet performed
                        % 100  --> fit running
                        % else --> exit flag of the completed fit

    F_ = [];  % Instance of the Fit class used to do the fitting
    chi2Fun_ = [];  % Chi square function
    nDataSets_ = 0;  % Number of data sets to fit
    parsAreGlobal_ = [];  % Which parameters are global
    nFlattenedParams_ = 0;  % Number of parameters in the array passed to Fit
    nGlobalPars_ = 0;  % Number of global parameters
    iGlobalPars_ = [];  % Indices of the global parameters
    nLocalPars_ = 0;  % Number of local parameters
    iLocalPars_ = [];  % Indices of the local parameters
    flattenedPars_ = [];  % Last flattened parameters
    unflattenedPars_ = [];  % Last flattened parameters
    end  % Properties

    % List of methods:  (" " are optional parameters)
    %  - Fit("xData, yData", "weights")        --> constructor
    %  - setData(xData, yData, "weights")      --> set data for the fitting

    methods
        function GFS = GlobalFitSimple()
            % TODO
        end

        function setData(GFS, xData, yData)
            % TODO: check that they are the same size etc.
            % TODO: add support for weights
            GFS.xData_ = xData;
            GFS.yData_ = yData;
            GFS.nDataSets_ = length(xData);

            % Ensure the data sets are column vectors
            for i=1:GFS.nDataSets_
                GFS.xData_{i} = GFS.xData_{i}(:);
                GFS.yData_{i} = GFS.yData_{i}(:);
            end

            GFS.F_ = Fit([1], [1]);

        end

        function setModel(GFS, model, nPars, parsAreGlobal)
            GFS.model_ = model;
            GFS.chi2Fun_ = @GFS.calculateChiSquare;

            GFS.nParam_ = nPars;
            GFS.parsAreGlobal_ = parsAreGlobal;  % TODO: check boolean vector
            GFS.nGlobalPars_ = sum(GFS.parsAreGlobal_);
            GFS.iGlobalPars_ = find(GFS.parsAreGlobal_);  % Indices of the global parameters
            GFS.nLocalPars_ = sum(~GFS.parsAreGlobal_);
            GFS.iLocalPars_ = find(~GFS.parsAreGlobal_);  % Indices of the local parameters

            GFS.nFlattenedParams_ = GFS.nGlobalPars_ + GFS.nLocalPars_*GFS.nDataSets_;

            GFS.F_.setChiSquare(GFS.chi2Fun_, GFS.nFlattenedParams_);
        end

        function setStart(GFS, start)
            s = size(start);
            if (s(1) ~= GFS.nDataSets_)
                msgID = 'GFS:setStart_wrongArraySize';
                msg = 'The size of the startpoint array does not match the number of datasets.';
                exception = MException(msgID, msg);

                throw(exception);
            end
            if (s(2) ~= GFS.nParam_)
                msgID = 'GFS:setStart_wrongArraySize';
                msg = 'The size of the startpoint array does not match the number of declared parameters.';
                exception = MException(msgID, msg);

                throw(exception);
            end

            GFS.start_ = start;
            GFS.unflattenedPars_ = start;
            GFS.flattenedPars_ = GFS.flattenParameters(start);

            GFS.F_.setStart(GFS.flattenedPars_);
        end

        function setUb(GFS, ub)  % Set upper bounds
            if size(ub) ~= size(GFS.start_)
                msgID = 'GFS:setUb_arraySize';
                msg = 'The size of the upper bounds matrix is different from that of the start point.';
                exception = MException(msgID, msg);

                throw(exception);
            end

            GFS.ub_ = ub;
            GFS.F_.setUb(GFS.flattenParameters(ub))
        end

        function setLb(GFS, lb)  % Set lower bounds
            if size(lb) ~= size(GFS.start_)
                msgID = 'GFS:setLb_arraySize';
                msg = 'The size of the lower bounds matrix is different from that of the start point.';
                exception = MException(msgID, msg);

                throw(exception);
            end

            GFS.lb_ = lb;
            GFS.F_.setLb(GFS.flattenParameters(lb))
        end

        function fixParameters(GFS, fixedValues)
        % NaN identify parameters whose values are not fixed
            fixed = GFS.flattenParameters(fixedValues);
            for i=1:length(fixed)
                if ~isnan(fixed(i))
                    GFS.F_.fixParameter(i, fixed(i));
                end
            end
        end

        function fit(GFS, calculateErrors)
            if nargin < 2
                calculateErrors = 1;
            end
            GFS.F_.fit(calculateErrors);

            % Get fitted parameters and unflatten them
            pars = GFS.F_.getFittedParameters();
            GFS.fitPar_ = GFS.unflattenParameters(pars);

            if calculateErrors
                ers = GFS.F_.getParametersErrors();
                GFS.fitParError_ = GFS.unflattenParameters(ers);
            end
        end

        function pars = getFittedParameters(GFS)
            if isempty(GFS.fitPar_)
                msgID = 'GFS:getFittedParameters_fitNotPerformed';
                msg = 'Fit not yet performed or not converged.';
                exception = MException(msgID, msg);

                throw(exception);
            end
            pars = GFS.fitPar_;
        end

        function ers = getParamersErrors(GFS)
            if isempty(GFS.fitPar_)
                msgID = 'GFS:getFittedParameters_fitNotPerformed';
                msg = 'Fit not yet performed or not converged.';
                exception = MException(msgID, msg);

                throw(exception);
            end
            if isempty(GFS.fitParError_)
                msgID = 'GFS:getParamersErrors_errorsNotCalculated';
                msg = 'Errors have not been calculated after fit.';
                exception = MException(msgID, msg);

                throw(exception);
            end
            ers = GFS.fitParError_;
        end

        function y = fitEval(GFS, x)
            if nargin < 2
                x = GFS.xData_;
            end

            if isempty(GFS.fitPar_)
                par = GFS.start_;
%                 warning('Fit not performed; evaluation on start point');
            else
                par = GFS.fitPar_;
            end

            for i=1:GFS.nDataSets_
                y{i} = GFS.tocolumn(GFS.model_(x{i}, par(i,:)));
            end
        end

    end  % Methods

    methods (Access = private)

        function chi2 = calculateChiSquare(GFS, unused_x, pars)
            chi2 = 0;
            GFS.unflattenedPars_ = GFS.unflattenParameters(pars);

            l = GFS.nGlobalPars_+1;
            for i=1:GFS.nDataSets_
                res = GFS.yData_{i} - GFS.model_(GFS.xData_{i}, GFS.unflattenedPars_(i,:));
                chi2 = chi2 + sum(res.^2);  % TODO: add weights
            end
        end

        function flattened = flattenParameters(GFS, unflattened)
            flattened(1:GFS.nGlobalPars_) = unflattened(1, GFS.iGlobalPars_);

            l = GFS.nGlobalPars_ + 1;
            for i=1:GFS.nDataSets_
                flattened(l:l+GFS.nLocalPars_-1) = unflattened(i, GFS.iLocalPars_);
                l = l + GFS.nLocalPars_;
            end
        end

        function unflattened = unflattenParameters(GFS, flattened)
            l = GFS.nGlobalPars_ + 1;
            for i=1:GFS.nDataSets_
                unflattened(i, GFS.iGlobalPars_) = flattened(1:GFS.nGlobalPars_);
                unflattened(i, GFS.iLocalPars_) = flattened(l:l+GFS.nLocalPars_-1);
                l = l + GFS.nLocalPars_;
            end
        end

        function col = tocolumn(GFS, v)
            col = v(:);
        end
    end  % Private methods

end  % Class