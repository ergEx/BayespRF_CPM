function [trials, wealth, onsets, stimulus_id] = cpm_events_to_trials(input_file, ...
                                                                      maxtime, ...
                                                                      fileOptions)
    % estimate_sd_for_snr Estimates the SD for white noise, that added to the
    % pure signal in VOI should create the desired SNR in DB
    %   INPUTS:
    %       -input_file: A table of an ideally BIDS compatible events.tsv file
    %       -maxtime: If the events file should be truncated based on the `onsets`
    %       column.
    %       -fileOptions: Specifies how the different fields or entries
    %       of the events.tsv should be interpreted.
    %   OUTPUTS:
    %       -trials: Trial structure
    %       -wealth: Wealth at each trial
    %       -onsets: a vector of onsets
    %       -stimulus_ids: Names of the included stimuli
    % ---------------------------------------------------------------------
    % (c) Copyright 2023 Simon R. Steinkamp
    arguments
        input_file table
        maxtime double = inf
        fileOptions.response_cue string = 'Response_window_onset'
        fileOptions.response_stimulus string = 'Stimulus_onset'
        fileOptions.response_wealth string = 'New_wealth_onset'
        fileOptions.old_wealth string = 'old_wealth'
        fileOptions.new_wealth string = 'new_wealth'
        fileOptions.stimulus_id string = 'gFactIncr'
    end

    % Read in BIDS table in .tsv, specifying NaNs.

    % Extract trials of interest
    response_window = input_file(input_file.event_type == fileOptions.response_cue, :);
    stimulus_onset = input_file(input_file.event_type == fileOptions.response_stimulus, :);
    wealth_update = input_file(input_file.event_type == fileOptions.response_wealth, :);

    % Assertion, that extract files are of same length!
    assert (height(wealth_update) == height(response_window) && ...
            height(wealth_update) == height(stimulus_onset));

    % Create a matrix as trial by information
    % Stack onsets for each trial (in order!)

    onsets = [response_window.onset, stimulus_onset.onset, wealth_update.onset];

    old_wealth = response_window.(fileOptions.old_wealth);
    new_wealth = wealth_update.(fileOptions.new_wealth);

    wealth = [old_wealth, old_wealth, new_wealth];

    stimulus_id = [response_window.event_type, ...
                   cellstr(num2str(stimulus_onset.(fileOptions.stimulus_id))), ...
                   wealth_update.event_type];

    n_trials = length(wealth);

    for tn = n_trials:-1:1
        if max(onsets(tn)) < maxtime
            trials(tn).wealth = wealth(tn, :);
            trials(tn).onsets = onsets(tn, :);
            trials(tn).stimuli = stimulus_id(tn, :);
        end
    end

end
