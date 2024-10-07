function [csc_trials, stimulus_ids] = cpm_trials_to_csc(trials, n_timesteps, ...
                                                        csc_idx, ...
                                                        relevant_stimuli, ...
                                                        reward_timestep, ...
                                                        normalize)
    % estimate_sd_for_snr Estimates the SD for white noise, that added to the
    % pure signal in VOI should create the desired SNR in DB
    %   INPUTS:
    %       -trials: Structure with an entry for each trial.
    %        Each trial should have the fields: onsets, wealth stimulus
    %       -n_timesteps: how many time-steps the CSC representation should have
    %       -csc_index: where to place the onsets of the trial structure in the csc
    %        representation
    %       -relevant_stimuli: which stimuli in the trial structure get a label
    %       -csc_idx: place of events from trials struct, ordered and relative to
    %                  timesteps
    %       -relevant stimuli: relative to stimuli in trials struct.
    %       -reward_timestep: relative to n_timesteps, at what timestep the reward occurs
    %   OUTPUTS:
    %       -csc_trials: Trial structure, augmented to have a CSC stimulus representation
    %       -stimulus_ids: Name for the stimulus ids.
    % ---------------------------------------------------------------------
    % (c) Copyright 2023 Simon R. Steinkamp

    arguments
        trials
        n_timesteps
        csc_idx
        relevant_stimuli
        reward_timestep
        normalize = false
    end

    stim_idx = 1:length(csc_idx);

    % There cannot be more indices than onsets.
    assert(length(csc_idx) == length(trials(1).onsets), ...
           "Csc idx out of bounds");
    % There cannot be more indices than stimuli
    assert(length(relevant_stimuli) <= length(trials(1).onsets), ...
           "Stimuli list out of bounds");
    % Indices have to be in the range of n_timesteps
    assert(max(csc_idx) <= n_timesteps, "csc to long");
    assert(max(relevant_stimuli) <= length(trials(1).onsets), ...
           "relevant stims too long");
    assert(reward_timestep <= n_timesteps, "reward out of bounds");

    n_trials = length(trials);

    % First step identify number of stimuli:
    stimulus_ids = cat(1, trials(:).stimuli);

    % Define labels
    id_labels = unique(stimulus_ids(:, relevant_stimuli));

    trial_stimulus_template = zeros(length(id_labels) * n_timesteps, ...
                                    n_timesteps);

    timestep_index = repmat(1:n_timesteps, 1, length(id_labels));
    timestep_index = timestep_index(:);
    stimulus_index = repmat(1:length(id_labels), n_timesteps, 1);
    stimulus_index = stimulus_index(:);

    stimulus_map = containers.Map(id_labels, num2cell(1:length(id_labels)));

    %%
    stimulus_ids = {};

    for nn = 1:length(stimulus_index)
        stimulus_ids{nn} = id_labels{stimulus_index(nn)};
    end

    %% Create empty trials:
    for tn = n_trials:-1:1
        csc_trials(tn).stimuli = trial_stimulus_template;
        csc_trials(tn).onsets = nan(1, n_timesteps);
        csc_trials(tn).wealth = nan(1, n_timesteps);
    end

    for tn = 1:n_trials

        for ts = stim_idx
            csc_trials(tn).onsets(csc_idx(ts)) = trials(tn).onsets(ts);
            csc_trials(tn).wealth(csc_idx(ts)) = trials(tn).wealth(ts);

            if ismember(ts, relevant_stimuli)
                time_idx = timestep_index >= csc_idx(ts);
                map_idx = stimulus_index == stimulus_map(trials(tn).stimuli{ts});
                idxs = find(time_idx .* map_idx);
                csc_trials(tn).stimuli(idxs, csc_idx(ts):end) = eye(length(idxs));
            end

        end

        if normalize
            csc_trials(tn).stimuli = (csc_trials(tn).stimuli ./ ...
                                      max(sum(csc_trials(tn).stimuli, 1)));
        end

        % Clean up stimuli - i.e. not include anything after reward:
        csc_trials(tn).stimuli(:, reward_timestep + 1:end) = 0;
        % Hacky solution
        csc_trials(tn).wealth(1:reward_timestep - 1) = ...
          nanmean(csc_trials(tn).wealth(1:reward_timestep - 1));
        csc_trials(tn).wealth(reward_timestep:end) = ...
          nanmean(csc_trials(tn).wealth(reward_timestep:end));
        csc_trials(tn).trial_no = zeros(size(csc_trials(tn).wealth)) + tn;

    end
