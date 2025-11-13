function [trajectory, trials] = cpm_decay_learning(freeparams, fixedparams, data)
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


    % Check for the presence of parameters and assign values to alpha
    if isfield(freeparams, 'lambda')
        % If alpha parameter exists, assign its value to alpha
        lambda = freeparams.lambda;
    elseif isfield(freeparams, 'xi')
        % If tau parameter exists, compute alpha using logit_inv function
        lambda = logit_inv(freeparams.xi);
    else
        lambda = 0;
    end

    if ~isfield(fixedparams, 'return')
        fixedparams.return = 'rpe';
    end

    if ~isfield(fixedparams, 'stimpos')
        fixedparams.stimpos = 3;
    end
 
    if ~isfield(fixedparams, 'rewardpos')
        fixedparams.rewardpos = 4;
    end

            
    % Recovering run parameters from stimuli.
    trial_n = length(trials);
    no_weights = size(trials(1).stimuli, 1);

    W = zeros(no_weights, 1);

    abs_max_diff = -inf;
    % transform trials wealth to reward
    for tn = 1:trial_n
        tmp_wealth = trials(tn).wealth;
        tmp_wealth = isoelastic_utility(tmp_wealth, 0);
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

            V_t = trials(tn).stimuli(:, fixedparams.stimpos)' * W;
            Z = trials(tn).stimuli(:, fixedparams.stimpos);

            reward = trials(tn).reward(fixedparams.rewardpos);

            reward_update = Z * reward;

            delta_t = reward - V_t;
            
            
            W =  lambda * W + reward_update;

            trials(tn).RPE(fixedparams.rewardpos) = delta_t;
            trials(tn).RPE(fixedparams.stimpos) = V_t;
            
            trials(tn).V(fixedparams.stimpos) = V_t;
            trials(tn).W(:, fixedparams.stimpos) = W;
            trials(tn).Z(:, fixedparams.stimpos) = Z;

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
