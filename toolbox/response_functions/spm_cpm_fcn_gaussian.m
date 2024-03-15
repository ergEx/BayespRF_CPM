function varargout = spm_cpm_fcn_gaussian(P, M, U, varargin)
    % Gaussian cognitive field model, with eucldiean coordinates and hard
    % constraints. Supports n-dimensional parameter spaces, where the field for
    % each parameter has a location and width parameter.
    %
    % mu_xx - Degrees of visual angle used as units
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
    % xy        [cpxn] vector of coordinates to evaluate the PRF
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
    % Adjusted for CPM: Simon Steinkamp 2023
    %

    if nargin < 4
        % TODO
        % NOTE Not correcting parameters for integration!
        P = correct_parameters(P, M);
        [z, varargout{1}, varargout{2}] = integrate_model(P, M, U);
    else
        action = varargin{1};

        switch action
            case 'get_parameters'
                % Gets the parameters for display

                % Correct neuronal parameters
                P = correct_parameters(P, M);

                % Corrected haemodynamic parameters
                P.transit = exp(P.transit);
                P.decay   = exp(P.decay);
                P.epsilon = exp(P.epsilon);

                varargout{1} = P;

            case 'latent_parameters'

                params = M.options.params;

                for pn = 1:length(params)
                    P.(['mu_' params{pn}]) = latent_parameter(P.(['mu_' params{pn}]), ...
                                                              M.pmin.(['mu_' params{pn}]), ...
                                                              M.pmax.(['mu_' params{pn}]));
                    P.(['width_' params{pn}]) = latent_parameter(P.(['width_' params{pn}]), ...
                                                                 M.pmin.(['width_' params{pn}]), ...
                                                                 M.pmax.(['width_' params{pn}]));
                end

                % Corrected haemodynamic parameters
                P.beta    = log(P.beta);

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
                [varargout{1}, varargout{2}, varargout{3}] = is_above_threshold(P, M, Cp, v, alpha);
            case 'get_response'
                P = correct_parameters(P, M);
                varargout{1} = get_response(P, M, U, U(1).gridpoints);
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
    coords = U(1).gridpoints;

    W = get_response(P, M, U, coords) .* P.beta;

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

            if max(ind) < nbins
                resp = U(t).signals1D' .* W;
                z(ind) = z(ind) + sum(resp);
            end
        end
    else
        z = zeros(n, 1);
        for t = 1:n
            % Microtime index for this volume
            resp = U(t).signals1D' .* W;
            z(t) = sum(resp);
        end
    end

    if nargout > 1
        % Integrate BOLD model (calls spm_prf_fx and spm_prf_gx)
        Z.u = z';
        Z.dt = M.dt;
        y = spm_int_sparse(P, M, Z);
    end

    % -------------------------------------------------------------------------
function x = get_response(P, M, U, coords)
    % Get PRF response to a stimulus covering the given pixels
    %
    % P  - parameters
    % xy - [2xn] vector of pixel coordinates
    % is_polar - true if the xy coordinates are polar

    params = M.options.params;

    mus = zeros(1, length(params));
    widths = zeros(1, length(params));

    for fn = 1:length(params)
        mus(fn) = P.(['mu_' params{fn}]);
        widths(fn) = P.(['width_' params{fn}]);
    end

    % Std dev -> covariance matrix
    sigma = diag(widths);
    sigma = sigma .* sigma;

    % Normalized Gaussian
    x = spm_mvNpdf(coords', mus, sigma);

    x = x ./ sum(x + eps); % Make it a PMF

    if ~isfinite(sum(x))
        warnstring = [];
        for sig = sigma
            warnstring = [warnstring sprintf(' %4.5f', sig)];
        end
        warning(['spm_mvNpdf sums to 0 with sigmas =', warnstring]);
    end

    % -------------------------------------------------------------------------
function x_c = constrain_parameter(x, d_min, d_max)
    % Convert latent variable x to parameter x_c where d_min <= x <= d_max
    % Gives uniform priors with mean 0 and variance 1
    %
    % x     - real-valued parameters
    % d_min - target minimum value
    % d_max - target maximum value
    %
    % x_c - parameters scaled between d_min and d_max

    d_range = d_max - d_min;
    x_c      = (d_range .* spm_Ncdf(x, 0, 1)) + d_min;

    % -------------------------------------------------------------------------
function x_c = latent_parameter(x, d_min, d_max)
    % Convert parameter x to latent parameter x_c where d_min <= x <= d_max
    % Gives uniform priors with mean 0 and variance 1
    %
    % x     - real-valued parameters
    % d_min - target minimum value
    % d_max - target maximum value
    %
    % x_c - scaled parameter back as normal distributed

    d_range = d_max - d_min;
    x_c = norminv(((x - d_min) ./ (d_range)) + eps);

    % -------------------------------------------------------------------------
function P = correct_parameters(P, M)
    % Prepare neuronal parameters for use e.g. exponentiate log parameters
    %
    % P - parameters
    % M - model

    % Constrain PRF centre
    params = M.options.params;

    for pn = 1:length(params)
        P.(['mu_' params{pn}]) = constrain_parameter(P.(['mu_' params{pn}]), ...
                                                     M.pmin.(['mu_' params{pn}]), ...
                                                     M.pmax.(['mu_' params{pn}]));
        P.(['width_' params{pn}]) = constrain_parameter(P.(['width_' params{pn}]), ...
                                                        M.pmin.(['width_' params{pn}]), ...
                                                        M.pmax.(['width_' params{pn}]));
    end

    % Convert log beta -> beta
    P.beta    = exp(P.beta);

    % -------------------------------------------------------------------------
function [tf, Pp, F0] = is_above_threshold(Pp, M, Cp, v, alpha)
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
    q(end - 3:end) = 1;

    % Remove others
    rC(q ~= 1) = 0;

    % BMR
    F0 = spm_log_evidence(Ep, Cp, pE, diag(spm_vec(pC)), rE, diag(rC));

    % Convert to post. prob - from spm_dcm_compare
    F    = [F0 0];
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

    for fn = 1:length(params)
        pE.(['mu_' params{fn}]) = 0;
        pE.(['width_' params{fn}]) = 0;
        pC.(['mu_' params{fn}]) = 1;
        pC.(['width_' params{fn}]) = 1;
    end

    pE.beta    = -2;
    pC.beta    = 4;

    % -------------------------------------------------------------------------
function P = glm_initialize(P, M, U, y)

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

        P2.beta  = 1;
        P2 = correct_parameters(P2, M);

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
        TR = (length(z) * U(1).dt) / M.ns; % Recover repetition time
        % Recover sampling factor:
        sampT =  round(TR / U(1).dt);

        ind    = (0:(ns - 1)) * sampT + M.T0;
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
