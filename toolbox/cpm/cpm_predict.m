function [y, x] = cpm_predict(prf, idx, from_posterior)

    % % cpm_predict, helper function, predicts signal from prf for given voxel.
    %
    % prf - a prf structure
    % idx - [optional, default=1] voxel index from which parameters should be
    %        predicted
    % from_posterior [optional, default=true] whether to predict from posterior
    %          or prior
    %
    % Returns:
    % true_P    - transformed parameters in real space.
    %
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
        prf
        idx = 1
        from_posterior = true
    end

    if from_posterior
        [y, x] = feval(prf.M.IS, prf.Ep{idx}, prf.M, prf.U);
    else
        [y, x] = feval(prf.M.IS, prf.M.pE{idx}, prf.M, prf.U);
    end
