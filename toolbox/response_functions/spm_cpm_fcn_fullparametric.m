function varargout = spm_cpm_fcn_fullparametric(P, M, U, varargin)
    % Full parametric model (has parameters for all dimensions of signals 1 D)
    %
    % D - Degrees of visual angle used as units
    % C - Constrained PRF width
    %
    % -------------------------------------------------------------------------
    % FORMAT [y,Z] = spm_prf_fcn_gaussian(P,M,U)
    % Return the BOLD and neuronal predicted timeseries
    %
    % P         parameters
    % M,U       model, inputs
    %
    % y         fMRI time series
    % Z         neuronal response
    % -------------------------------------------------------------------------
    % FORMAT P = spm_prf_fcn_gaussian(P,M,U,'get_parameters')
    % Return the given parameters corrected for display
    %
    % P         parameters
    % M,U       model, inputs
    % -------------------------------------------------------------------------
    % FORMAT S = spm_prf_fcn_template(P,M,U,'get_summary')
    % Summarises the pRF with simple (Gaussian) parameters x,y,width,beta
    %
    % S         structure with fields x,y,width,beta
    % M,U       model, inputs
    % -------------------------------------------------------------------------
    % FORMAT P = spm_prf_fcn_gaussian(P,M,U,'get_response',xy)
    % Return the instantaneous response of the PRF at coordinates xy
    %
    % P         parameters
    % M,U       model, inputs
    % xy        [2xn] vector of coordinates to evaluate the PRF
    %
    % ---------------------------------------------------------------------
    % Copyright (C) 2016 Peter Zeidman
    % This program is free software: you can redistribute it and/or modify
    % it under the terms of the GNU General Public License as published by
    % the Free Software Foundation, either version 3 of the License, or
    % (at your option) any later version.
    %
    % This program is distributed in the hope that it will be useful,
    % but WITHOUT ANY WARRANTY; without even the implied warranty of
    % MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    % GNU General Public License for more details.
    %
    % You should have received a copy of the GNU General Public License
    % along with this program.  If not, see <http://www.gnu.org/licenses/>.
    % ---------------------------------------------------------------------

    if nargin < 4
        % TODO
        % NOTE Not correcting parameters for integration!

        [z, varargout{1}, varargout{2}] = integrate_model(P, M, U);
    else
        action = varargin{1};

        switch action
            case 'get_parameters'
                % Gets the parameters for display

                % Corrected haemodynamic parameters
                P.transit = exp(P.transit);
                P.decay   = exp(P.decay);
                P.epsilon = exp(P.epsilon);

                varargout{1} = P;

            case 'latent_parameters'

                P.transit = log(P.transit);
                P.decay   = log(P.decay);
                P.epsilon = log(P.epsilon);

                varargout{1} = P;

            case 'get_summary'
                % Gets x,y,sigma,beta

                P = correct_parameters(P, M);
                varargout{1} = P;

            case 'is_above_threshold'
                Cp = varargin{2};
                v  = varargin{3};
                if length(varargin) >= 4
                    alpha = varargin{4};
                else
                    alpha = 0.95;
                end
                [varargout{1}, varargout{2}] = is_above_threshold(P, M, Cp, v, alpha);
            case 'get_response'
                error('"get_response", not implemented for this type of function!');
                % P = correct_parameters(P,M);
                % varargout{1} = get_response(P,M,U, U(1).gridpoints);
            case 'get_priors'
                [pE, pC] = get_priors(M);
                varargout{1} = pE;
                varargout{2} = pC;
            case 'glm_initialize'
                y = varargin{2};
                pE = glm_initialize(P, M, U, y);
                varargout{1} = pE;
            otherwise
                error('Unknown action');
        end
    end

    % -------------------------------------------------------------------------
function [z, y, Z] = integrate_model(P, M, U)
    % Integrates the neuronal and BOLD models

    n  = length(U);            % Number of volumes = inputs
    % T0 = M.T0;                 % Microtime offset

    nbins = round(U(1).nbins); % Ensure nbins is ineger...

    w = get_response(P, M, U);

    try
        mc = U(1).microtime;
    catch
        mc = true;
    end

    if mc
        z = zeros(1, nbins);
        for t = 1:n
            % Microtime index for this volume
            ind = U(t).ind;
            resp = U(t).signals1D(:)' * w(:);
            z(ind) = z(ind) + sum(resp);
        end
    else
        z = zeros(n, 1);
        for t = 1:n
            % Microtime index for this volume
            resp = U(t).signals1D(:)' * w(:);
            z(t) = sum(resp);
        end
    end

    % z = z + P.beta0;

    if nargout > 1
        % Integrate BOLD model (calls spm_prf_fx and spm_prf_gx)
        Z.u = z';
        Z.dt = M.dt;
        y = spm_int_sparse(P, M, Z);
    end

    % -------------------------------------------------------------------------
function w = get_response(P, M, U)
    % Get PRF response to a stimulus covering the given pixels
    %
    % P  - parameters
    % xy - [2xn] vector of pixel coordinates
    % is_polar - true if the xy coordinates are polar

    w = zeros(1,  M.options.params);

    for ii = 1:M.options.params
        w(1, ii) = P.(['beta' num2str(ii)]);
    end

    % -------------------------------------------------------------------------
function P = correct_parameters(P, M)
    % Prepare neuronal parameters for use e.g. exponentiate log parameters
    %
    % P - parameters
    % M - model
    error('Not implemented');
    % Constrain PRF centre

    % -------------------------------------------------------------------------
function [tf, Pp] = is_above_threshold(Pp, M, Cp, v, alpha)
    % Evaluate whether the model with parameters P and covariance Cp passes an
    % arbitrary threshold for display
    %
    % P     - parameters
    % M     - model
    % Cp    - covariance matrix
    % alpha - threshold

    pC = M.pC{v};
    pE = M.pE{v};
    Cp = Cp{v};
    Ep = Pp{v};

    rE = spm_vec(pE);
    rC = spm_vec(pC);

    % Indices of parameters to keep in nested model
    np = length(rE);
    q  = zeros(1, np);
    q(end - 2:end) = 1;

    % Remove others
    rC(q ~= 1) = 0;

    % BMR
    F = spm_log_evidence(Ep, Cp, pE, diag(spm_vec(pC)), rE, diag(rC));

    % Convert to post. prob - from spm_dcm_compare
    F    = [F 0];
    i    = F < (max(F) - 32);
    Pp    = F;
    Pp(i) = max(F) - 32;
    Pp    = Pp - min(Pp);
    Pp    = exp(Pp);
    Pp    = Pp / sum(Pp);

    tf = Pp(2) > alpha;
    Pp = Pp(2);

    % -------------------------------------------------------------------------
function [pE, pC] = get_priors(M)

    % Get the neuronal priors

    params = M.options.params;

    pE = {};
    pC = {};

    for fn = 1:params
        pE.(['beta' num2str(fn)]) = 0;
        pC.(['beta' num2str(fn)]) = 1 / params;
    end
    %
    % pE.beta0    = 0;        pC.beta0    = 1;

    % -------------------------------------------------------------------------
function P = glm_initialize(P, M, U, y)

    error('Not yet implemented!');
    % Coordinates to test
    params = M.options.params;

    x = {};

    for ii = 1:length(params)
        x{ii} = [-2:0.2:2]; % This is simply in unconstrained, latent space!
    end

    mu_dim = {1:length(params)};

    % Candidates for mu
    [mu_dim{1:length(params)}] = ndgrid(x{:});

    mu = [];
    for ii = 1:length(params)
        mu = [mu, mu_dim{ii}(:)];
    end

    % Model space

    nm = size(mu, 1);
    sse = zeros(1, nm);

    for i = 1:nm

        P2         = P;
        for pn = 1:length(params)
            P2.(['mu_' params{pn}]) = mu(i, pn);
        end

        P2.beta    = 0;

        P2 = correct_parameters(P2, M);

        % fprintf('M%d  dist: %2.2f angle %2.2f\n',i,P2.dist,P2.angle);

        % Get neuronal response
        z = integrate_model(P2, M, U);

        % Convolve with canonical HRF
        XU.u = z';
        XU.dur = 0;
        XU.dt = U(1).dt;
        XU.name{1} = 'gPRF';
        xBF.name   = 'hrf';
        xBF.order  = 1;
        xBF.length = 32;
        xBF.dt     = U(1).dt;
        xBF        = spm_get_bf(xBF);
        z = spm_Volterra(XU, xBF.bf, 1);

        % Downsample regressor
        fMRI_T = spm_get_defaults('stats.fmri.t');
        ns     = M.ns;
        ind    = (0:(ns - 1)) * fMRI_T + M.T0;
        z      = z(ind, :);

        X = [z ones(ns, 1)];

        % Invert
        beta   = pinv(X) * y;
        e      = y - X * beta;
        sse(i) = e' * e;
    end

    [i, idx] = min(sse);

    for pn = 1:length(params)
        P.(['mu_' params{pn}]) = mu(idx, pn);
    end
