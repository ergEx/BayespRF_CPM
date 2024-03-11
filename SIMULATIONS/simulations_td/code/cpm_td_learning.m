function [trajectory, trials] = cpm_td_learning(freeparams, fixedparams, data)
    % Simulation configs, file for handling of paths to toolboxes.
    % ---------------------------------------------------------------------
    % Copyright (C) 2023 Simon R. Steinkamp, Iyadh Chaker
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

    if isfield(data, 'trials')
        trials = data.trials;
    elseif isfield(data, 'stimuli')
        trials = data;
    end

    gamma = fixedparams.gamma;
    lambda = fixedparams.lambda;

    % Check for the presence of parameters and assign values to eta
    if isfield(freeparams, 'eta')
        % If eta parameter exists in freeparams, assign its value to eta
        eta = freeparams.eta;
    elseif isfield(fixedparams, 'eta')
        % If eta parameter exists in fixedparams, assign its value to eta
        eta = fixedparams.eta;
    else
        % Assuming linear utility.
        eta = 0;
    end

    % Check for the presence of parameters and assign values to alpha
    if isfield(freeparams, 'alpha')
        % If alpha parameter exists, assign its value to alpha
        alpha = freeparams.alpha;
    elseif isfield(freeparams, 'tau')
        % If tau parameter exists, compute alpha using logit_inv function
        alpha = logit_inv(freeparams.tau);
    elseif isfield(freeparams, 'tauneg') && isfield(freeparams, 'taupos')
        % If both tauneg and taupos parameters exist, compute alpha using logit_inv function
        alpha = [logit_inv(freeparams.tauneg), logit_inv(freeparams.taupos)];
    elseif isfield(freeparams, 'alphaneg') && isfield(freeparams, 'alphapos')
        % If both alphaneg and alphapos parameters exist, assign their values to alpha
        alpha = [freeparams.alphaneg, freeparams.alphapos];
    else
        alpha = 0;
    end

    if ~isfield(fixedparams, 'return')
        fixedparams.return = 'rpe';
    end

    if length(alpha) == 2
        alpha_neg = alpha(1);
        alpha_pos = alpha(2);
    else
        alpha_neg = alpha;
        alpha_pos = alpha;
    end

    % Recovering run parameters from stimuli.
    trial_n = length(trials);
    no_weights = size(trials(1).stimuli, 1);

    W = zeros(no_weights, 1);

    abs_max_diff = -inf;
    % transform trials wealth to reward
    for tn = 1:trial_n
        tmp_wealth = trials(tn).wealth;
        tmp_wealth = isoelastic_utility(tmp_wealth, eta);
        tmp_wealth = [0, diff(tmp_wealth)];
        trials(tn).reward = tmp_wealth;
        trial_steps = length(tmp_wealth);
        trials(tn).W = nan(no_weights, trial_steps);
        trials(tn).RPE = nan(trial_steps, 1);
        trials(tn).V = nan(trial_steps, 1);
        trials(tn).Z = nan(no_weights, trial_steps);

        if max(abs(tmp_wealth)) > abs_max_diff
            abs_max_diff = max(abs(tmp_wealth));
        end
    end

    for tn = 1:trial_n

        timestep_n = length(trials(tn).reward);

        Z = zeros(no_weights, timestep_n);

        for t = 2:timestep_n

            V_t = trials(tn).stimuli(:, t - 1)' * W;
            V_tp1 = trials(tn).stimuli(:, t)' * W;

            delta_t = (trials(tn).reward(t)) + gamma .* V_tp1 - V_t;

            Z(:, t) = gamma .* lambda .* Z(:, t - 1) + trials(tn).stimuli(:, t);

            if delta_t >= 0
                W = W + (alpha_pos .* delta_t .* Z(:, t - 1));
            else
                W = W + (alpha_neg .* delta_t .* Z(:, t - 1));
            end

            trials(tn).RPE(t) = delta_t;
            trials(tn).V(t) = V_t;
            trials(tn).W(:, t - 1) = W;
            trials(tn).Z(:, t) = Z(:, t - 1);

        end
    end

    switch fixedparams.return
        case 'rpe'
            trajectory = cat(1, trials.RPE);
        case 'value'
            trajectory = cat(1, trials.V);
    end

    trajectory(isnan(trajectory)) = 0;
    trajectory = trajectory ./ max(abs(trajectory));

end

function out = isoelastic_utility(c, eta)
    % ISOELASTIC_UTILITY, isoelastic utility function, for wealth (c) > 0 and real parameter
    % eta.
    if eta == 1
        out = log(c);
    else
        out = ((c.^(1 - eta) - 1) ./ (1 - eta));
    end
end

function out = logit_inv(a)
    % LOGIT_INV Summary calculates the inverse logit of a
    out = exp(a) ./ (1 + exp(a));
end
