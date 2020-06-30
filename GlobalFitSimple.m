classdef GlobalFitSimple < handle
    % GlobalFitSimple: class for global fitting of data with the same equation
    %  and single or global parameters
    % v. 0.2.0

        % TODO:
    %  -[x] Evaluate chi square
    %  -[x] Fix getChiSquare in class Fit if the model is not provided
    %  -[x] Normalize the fit errors for the reduced chi square _here_
    %  -[x] Take into account the constraints in the nDof calculation
    %  -[ ] Update list of methods
    %  -[x] Test errors in GFS and Fit

    % Changelog:
    %   26/06/20 - 0.2.0: working well, limited functionality. Some functions
    %       need to be called from the underlying Fit instance
    %   22/01/20 - 0.1.2: fixed errors on parameters
    %   10/10/19 - 0.1.1: fixed constructor
    %   05/07/18 - 0.1: first version

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

    F_                = [];  % Instance of the Fit class used to do the fitting
    chi2Fun_          = [];  % Chi square function
    nDataSets_        = 0;   % Number of data sets to fit
    parsAreGlobal_    = [];  % Which parameters are global
    nFlattenedParams_ = 0;   % Number of parameters in the array passed to Fit
    nGlobalPars_      = 0;   % Number of global parameters
    iGlobalPars_      = [];  % Indices of the global parameters
    nLocalPars_       = 0;   % Number of local parameters
    iLocalPars_       = [];  % Indices of the local parameters
    flattenedPars_    = [];  % Last flattened parameters
    unflattenedPars_  = [];  % Last unflattened parameters

    end  % Properties

    % List of methods:  (" " are optional parameters)
    %  - Fit("xData, yData", "weights")        --> constructor
    %  - setData(xData, yData, "weights")      --> set data for the fitting

    methods
        function GFS = GlobalFitSimple()
            GFS.F_ = Fit([1], [1]);
        end

        function setData(GFS, xData, yData, weights)
            % TODO: check that they are the same size etc.
            % TODO: check support for weights
            if nargin < 4
                weights = cell(size(xData));
                for i=1:length(xData)
                    weights{i} = ones(size(xData{i}));
                end
            end
            GFS.xData_ = xData;
            GFS.yData_ = yData;
            GFS.weights_ = weights;
            GFS.nDataSets_ = length(xData);

            % Ensure the data sets are column vectors
            for i=1:GFS.nDataSets_
                GFS.xData_{i} = GFS.xData_{i}(:);
                GFS.yData_{i} = GFS.yData_{i}(:);
                GFS.weights_{i} = GFS.weights_{i}(:);
            end

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

        function fit(GFS, calculateErrors, ignoreFitWarnings)
            if nargin < 2
                calculateErrors = 1;
            end
            if nargin < 3
                ignoreFitWarnings = 0;
            end
            GFS.F_.fit(calculateErrors, ignoreFitWarnings);

            % Get fitted parameters and unflatten them
            pars = GFS.F_.getFittedParameters();
            GFS.fitPar_ = GFS.unflattenParameters(pars);

            GFS.chi2_ = GFS.F_.getChiSquare();

            if calculateErrors
                ers = GFS.F_.getParametersErrors();
                GFS.fitParError_ = GFS.unflattenParameters(ers)*...
                    sqrt(GFS.getChiSquare(1));
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

        function res = getResiduals(GFS)  % TODO: symmetrize this function wrt Fit
            if isempty(GFS.fitPar_)       %  (and maybe add a private function residuals())
                msgID = 'GFS:getResiduals_fitNotPerformed';
                msg = 'Before evaluating residuals the fit must be performed.';
                exception = MException(msgID, msg);

                throw(exception);
            end

            yFit = GFS.fitEval();
            for i=1:GFS.nDataSets_
                res{i} = (GFS.yData_{i} - yFit{i}) .* sqrt(GFS.weights_{i});
            end
        end

        function chi2 = getChiSquare(GFS, reduced)
             if nargin < 2
                reduced = 0;
             end

             if reduced
                nPoints = 0;
                for i=1:numel(GFS.xData_)
                    nPoints = nPoints + numel(GFS.xData_{i});
                end
                chi2 = GFS.chi2_/(nPoints - GFS.nFlattenedParams_ ...
                    + sum(~GFS.F_.fixed_));
            else
                chi2 = GFS.chi2_;
            end
        end

    end  % Methods

    methods (Access = public)

        function chi2 = calculateChiSquare(GFS, ~, pars)
            chi2 = 0;
            GFS.unflattenedPars_ = GFS.unflattenParameters(pars);

            l = GFS.nGlobalPars_ + 1;
            for i=1:GFS.nDataSets_
                res = GFS.yData_{i} - GFS.model_(GFS.xData_{i}, GFS.unflattenedPars_(i,:));
                chi2 = chi2 + sum(res.^2 .* GFS.weights_{i});  % TODO: check weights
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