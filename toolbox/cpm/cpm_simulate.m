function [y, x]  = cpm_simulate(P, PRF, noise_sd, is_latent)
    % % cpm_simulate
    %
    % P - parameters of PRF model
    % PRF    - subject's estimated pRF model
    % noise_sd     - standard deviation of noise to add to outputs
    % is_latent - If the submitted parameters are in latent space already (i.e.
    %             location of the underlying Gaussian distribution), if not, tries
    %             to transform parameters into latent space.
    %
    % Returns:
    % y    - simulated signal, with added noise
    % x     - the model's latent states
    %
    % For details on why one might want to disable conditional dependencies,
    % see https://en.wikibooks.org/wiki/SPM/Bayesian_Parameter_Averaging_(BPA)
    %
    % ---------------------------------------------------------------------
    % Copyright (C) 2023 Simon R. Steinkamp
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
    arguments
        P
        PRF
        noise_sd = 0
        is_latent = false
    end

    P0 = PRF.M.pE{1};

    params = fieldnames(P);

    if is_latent

        Ptemp = P0;

        for fn = 1:length(params)
            Ptemp.(params{fn}) =  P.(params{fn});
        end

    else

        Ptemp =  feval(PRF.M.IS, P0, PRF.M, PRF.U, 'get_parameters');

        for fn = 1:length(params)
            Ptemp.(params{fn}) =  P.(params{fn});
        end

        Ptemp = feval(PRF.M.IS, Ptemp, PRF.M, PRF.U, 'latent_parameters');

    end

    [y, x] = feval(PRF.M.IS, Ptemp, PRF.M, PRF.U);

    y = y + normrnd(0, noise_sd, size(y));

end
