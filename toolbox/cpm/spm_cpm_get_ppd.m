function [s_pe, s_ep, prior_pd, post_pd] = spm_cpm_get_ppd(PRF, coords, idx, nsamp)
    % Sample from the prior (or posterior) of a pRF model and run samples
    % through the likelihood, to give the prior (or posterior) predictive
    % density (PPD).
    %
    % PRF   - estimated PRF
    % xy    - coordinates in stimulus space to sample
    % idx   - index of the PRF within the PRF structure
    % nsamp - the number of samples to draw
    %
    % Returns:
    % s_pe      - samples from the prior [parameters x samples]
    % s_ep      - samples from the posterior [parameters x samples]
    % prior_pd - average sampled PRF response (prior predictive density)
    % prior_pd - average sampled PRF response (posterior predictive density)
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
    % Changes
    % 2023.07.05 - By Simon R. Steinkamp removed visual displays for s_peed up

    if nargin < 4
        nsamp = 1;
    end

    % Sample from priors and posteriors
    % -------------------------------------------------------------------------

    pE_struct = PRF.M.pE{idx};

    % Get priors and posteriors
    pE = spm_vec(PRF.M.pE{idx});
    pC = spm_vec(PRF.M.pC{idx});
    Ep = spm_vec(PRF.Ep{idx});
    Cp = full(PRF.Cp{idx});

    if isvector(pC)
        pC = diag(pC);
    end

    % Sample
    s_pe = spm_normrnd(pE, pC, nsamp);
    s_ep = spm_normrnd(Ep, Cp, nsamp);

    if nargout < 3
        return
    end

    % Integrate model with sampled parameters
    % -------------------------------------------------------------------------

    % Create progress window

    prior_pd = 0;
    post_pd  = 0;

    parfor i = 1:nsamp

        % Integrate model using prior sample
        s_pe_struct = spm_unvec(s_pe(:, i), pE_struct);

        g = feval(PRF.M.IS, s_pe_struct, PRF.M, PRF.U, 'get_response', coords);
        prior_pd = prior_pd + g;

        % Integrate model using posterior sample
        s_ep_struct = spm_unvec(s_ep(:, i), pE_struct);

        g = feval(PRF.M.IS, s_ep_struct, PRF.M, PRF.U, 'get_response', coords);
        post_pd = post_pd + g;
    end

    % Average
    prior_pd = prior_pd .* (1 / nsamp);
    post_pd  = post_pd .* (1 / nsamp);
