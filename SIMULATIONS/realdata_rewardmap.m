function [] = realdata_rewardmap(redo)

% Script to perform analysis of Niv et al. 2012.
% ---------------------------------------------------------------------
% Copyright (C) 2025 Simon R. Steinkamp
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

%% Adding necessary paths:
configs = simulation_configs();

addpath(genpath('../toolbox'));
addpath(configs.spm_path);
addpath(genpath(configs.brain_slicer)); % TODO Add to configs
addpath(genpath(configs.vba_path));
addpath('simulations_td/code');

vis_template = '../../toolboxes/tpl-MNI152NLin2009cAsym_res-02_T1w.nii.gz';
datapath =configs.realdata;

spm fmri;
close all;

out_folder = 'realdata_rewardMap/';

mkdir(out_folder);
mkdir('realdata_rewardMap/models_init');

model_names = {"td", "rs"};
% ==========================================================================
% PREPARATION
% ==========================================================================
logit_inv = @(x) 1 ./ (1 + exp(-x));
% Loading some example data:

% TE of middle image.
TE = 31.75 ./ 1000; % From paper

rewardmap_participants = configs.participants
% This is where the PRF will be saved, we just set it to be here.
fixedparams = struct('gamma', 0.99, 'lambda', 1.0);

% Setting options for inversion
invert_options = struct('use_parfor', true, ...
    'init', 'GLM_P', ...
    'nograph', true);
%
resolution = 41;

regs = {"lnacc", "rnacc"};

PRFS = {};
% Load participant data:
for sub_idx = 1 : length(rewardmap_participants)

    participant = rewardmap_participants{sub_idx};
    %% Add all to configs
    behavior =  readtable(sprintf(fullfile(datapath, "sub-%s/ses-02/func/sub-%s_ses-02_task-risksensitive_acq-te14ipat2mb2me3pf78_dir-AP_run-04_events.tsv"), participant, participant), ...
        "Delimiter","\t", "FileType","text");
    reward_onset = table2array(behavior(strcmp(behavior.event_type, 'reward'), 'onset'));
    reward = str2double(table2array(behavior(strcmp(behavior.event_type, 'reward'), 'reward')));
    stim = table2array(behavior(strcmp(behavior.event_type, 'reward'), 'current_location'));
    stim_onset = table2array(behavior(strcmp(behavior.event_type, 'selection'), 'onset'));

    spm_mat = sprintf("realdata_rewardMap/preprocessed_data/sub-%s/SPM.mat", participant);
    SPM = load(spm_mat);
    SPM = SPM.SPM;
    SPM.swd = 'realdata_rewardMap/models_init';

    n_trials = length(stim_onset);
    % For simplicity, we can re-use the CSC representation scripts I implemented for
    % the simulation study, but for this we need a "trial" structure first.
    trials = {};
    cc = 1;
    for nt = 1:n_trials
        if ~isnan(reward(nt)) % We remove trials without a response
            trials(cc).onsets = [stim_onset(nt), reward_onset(nt)];
            trials(cc).stimuli = {num2str(stim(nt)), 'reward'};
            trials(cc).wealth = [0, reward(nt)];
            cc = cc + 1;
        end
    end
    %%
    % Create a CSC representation from the trials struct. We use a very simple one
    % containing only the onset of the stimulus and the reward.
    csc_trials = cpm_trials_to_csc(trials, 4, [2, 3], [1], 3);
    % ==========================================================================
    %% The reinforcement learning algorithm
    % ==========================================================================
    % The BayesPRF requires an SPM struct, but only a few fields from there, which
    % we generate here:
    % Create a dummy VOI
    for data_idx = 1 : 2
        single_participant = sprintf("realdata_rewardMap/preprocessed_data/sub-%s/VOI_%s_1.mat", participant, regs{data_idx});

        VOI = load(single_participant);
        dt = SPM.xBF.dt;
        % Creating a data struct
        data = {};
        data.ons = cat(2, csc_trials.onsets);
        data.dur = zeros(size(data.ons)) + 2 .* dt;
        data.dt = zeros(size(data.ons)) + dt;
        data.trials = csc_trials;

        %% PRF
        for mi = 1:length(model_names)
            if strcmp(model_names{mi}, 'rs')
                name = ['sub-' num2str(sub_idx), '_data-', ...
                    num2str(data_idx),  '_rs'];
            elseif strcmp(model_names{mi}, 'td')
                name = ['sub-' num2str(sub_idx), '_data-', ...
                    num2str(data_idx) '_td'];
            end

            if ~ isfile( fullfile(SPM.swd, ['PRF_' name '.mat'])) | redo
                if strcmp(model_names{mi}, 'rs')
                    grid = struct('taupos', [-8.1259, 8.1259, resolution], ...
                        'tauneg', [-8.1259, 8.1259, resolution]);
                    params = {{'taupos', 'tauneg'}};
                    U_prf = cpm_precompute(@cpm_td_learning, ...
                        grid, fixedparams, data, 'U_prf', true);
                    U_prf = cpm_set_constraints(U_prf, 'taupos', [-4.0, 4.0]);
                    U_prf = cpm_set_constraints(U_prf, 'tauneg', [-4.0, 4.0]);
                elseif strcmp(model_names{mi}, 'td')
                    grid = struct('tau', [-8.1259, 8.1259, resolution]);
                    params = {{'tau'}};
                    U_prf = cpm_precompute(@cpm_td_learning, ...
                        grid, fixedparams, data, 'U_prf', true);
                    U_prf = cpm_set_constraints(U_prf, ...
                        'tau', ...
                        [-4.0, 4.0]);
                end

                options = struct('model', 'spm_cpm_fcn_gaussian_normal_beta', ...
                    'name', name, ...
                    'params', params, ...
                    'TE', TE, ...
                    'B0', 3, ...
                    'voxel_wise', true, ...
                    'avg_sess', false);

                PRF = spm_prf_analyse('specify', SPM, VOI, U_prf, options);
                PRF.M.noprint = 1;
                PRFn = spm_prf_analyse('estimate', PRF,  invert_options);
                save_wrapper(PRFn, fullfile(SPM.swd, ['PRF_' name '.mat']));

            else
                tmp = load( fullfile(SPM.swd, ['PRF_' name '.mat']));
                PRFn = tmp.PRF;
            end
            PRFS{sub_idx, data_idx, mi} = PRFn;
        end
    end
end

%%
pp_left_acc = zeros(size(PRFS,1), size(PRFS,3), length(PRFS{1, 1, 1}.F));
pp_right_acc = zeros(size(PRFS,1), size(PRFS,3), length(PRFS{1, 2, 1}.F));

for part = 1 : size(PRFS, 1)
    for mi = 1 : size(PRFS, 3)
        for di = 1 : size(PRFS, 2)
            for fi = 1 : length(PRFS{part, di, mi}.F)
                if di == 1
                    pp_left_acc(part, mi, fi) = PRFS{part, di, mi}.Pp{fi}.beta;
                else
                    pp_right_acc(part, mi, fi) = PRFS{part, di, mi}.Pp{fi}.beta;
                end
            end
        end
    end
end

%% Best subject
[~, b_idx ] = max(squeeze(mean(mean(pp_left_acc, 2), 3)) + squeeze(mean(mean(pp_right_acc, 2), 3)));

%% Get alpha
params_left_acc_td = zeros(2, size(PRFS,1), length(PRFS{1, 1, 1}.F));
params_right_acc_td = zeros(2, size(PRFS,1), length(PRFS{1, 2, 1}.F));
params_left_acc_rs = zeros(3, size(PRFS,1),  length(PRFS{1, 1, 1}.F));
params_right_acc_rs = zeros(3, size(PRFS,1),  length(PRFS{1, 2, 1}.F));

for part = 1 : size(PRFS, 1)
    for mi = 1 : size(PRFS, 3)
        for di = 1 : size(PRFS, 2)
            for fi = 1 : length(PRFS{part, di, mi}.F)

                real_params = cpm_get_true_parameters(PRFS{part, di, mi}, fi);

                if di == 2
                    if mi == 1
                        params_right_acc_td(2, part, fi) = logit_inv(real_params.mu_tau);
                        params_right_acc_td(1, part, fi) =real_params.beta;
                    else
                        params_right_acc_rs(2, part,  fi) = logit_inv(real_params.mu_taupos);
                        params_right_acc_rs(3, part,  fi) = logit_inv(real_params.mu_tauneg);
                        params_right_acc_rs(1, part,  fi) = logit_inv(real_params.beta);
                    end
                else
                    if mi == 1
                        params_left_acc_td(2, part,  fi) = logit_inv(real_params.mu_tau);
                        params_left_acc_td(1, part,  fi) =real_params.beta;
                    else
                        params_left_acc_rs(2, part,  fi) = logit_inv(real_params.mu_taupos);
                        params_left_acc_rs(3, part,  fi) = logit_inv(real_params.mu_tauneg);
                        params_left_acc_rs(1, part,  fi) = logit_inv(real_params.beta);
                    end
                end
            end
        end
    end
end

%%
learning_rate_left = squeeze(params_left_acc_td(2, b_idx, :));
optimism_left = squeeze(params_left_acc_rs(2, b_idx, :) ./ (params_left_acc_rs(2, b_idx, :)  + params_left_acc_rs(3, b_idx, :) ));
lravg_left = squeeze((params_left_acc_rs(2, b_idx, :)  + params_left_acc_rs(3, b_idx, :)) / 2);

learning_rate_right = squeeze(params_right_acc_rs(2, b_idx, :));
optimism_right =squeeze( params_right_acc_rs(2, b_idx, :) ./ (params_right_acc_rs(2, b_idx, :)  + params_right_acc_rs(3, b_idx, :) ));
lravg_right = squeeze((params_right_acc_rs(2, b_idx, :)  + params_right_acc_rs(3, b_idx, :)) / 2);

%%
threshold = 0.75;
optimism_right(pp_right_acc(b_idx, 2, :) < threshold) = nan;
learning_rate_right(pp_right_acc(b_idx, 1, :) < threshold)  = nan;
learning_rate_left(pp_left_acc(b_idx, 1, :) < threshold)  = nan;
optimism_left(pp_left_acc(b_idx, 2, :) < threshold) = nan;
lgBF_td_vs_rs_left = PRFS{b_idx, 1, 1}.F - PRFS{b_idx, 1, 2}.F;
lgBF_td_vs_rs_right = PRFS{b_idx, 2, 1}.F - PRFS{b_idx, 2, 2}.F;
log_bf_to_map = [lgBF_td_vs_rs_left'; lgBF_td_vs_rs_right'];

cpm_write_to_nifti(SPM.xY.P(1, :), [PRFS{1, 1, 1}.xY.XYZmm, PRFS{1, 2, 1}.xY.XYZmm], [optimism_left; optimism_right], [out_folder, 'optimism_best_sub.nii'])
cpm_write_to_nifti(SPM.xY.P(1, :), [PRFS{1, 1, 1}.xY.XYZmm, PRFS{1, 2, 1}.xY.XYZmm], [learning_rate_left; learning_rate_right], [out_folder, 'lr_best_sub.nii'])
cpm_write_to_nifti(SPM.xY.P(1, :), [PRFS{1, 1, 1}.xY.XYZmm, PRFS{1, 2, 1}.xY.XYZmm], log_bf_to_map, [out_folder, 'bf_best_sub.nii'])

%%

%%
slicer_args = {
        'size', 'h50', ...
    'slices', [34 : 37], ...% print one row with 8 slices equally spaced
    'colorMode','w', ...
    'innerMargins', [0.005, 0], ...
    'margins', [0, 0, 0.1, 0], ...
        'cbLocation','east',... % colorbar location can be south or east
            'zoom', {[48 - 15, 48 + 15], [48 - 15, 48 + 15]}, ...
    'resolution', 300, ...
    'minClusterSize',{0,0},
};

bf_lims = max([abs(min(log_bf_to_map)), max(log_bf_to_map)]);

slicer({vis_template,[out_folder ,'bf_best_sub.nii']},...
    'limits',{[],[-bf_lims, bf_lims]},...
    'colormaps', {1, 48}, ...
    'labels',{[],['\tau^*']},... % when a layer's label is empty no colorbar will be printed.
    'title','Bayes Factor map: TD vs RS',...
    'output', [out_folder 'bf.pdf'], ...
    slicer_args{:}) % two background modalities: black ('k') or white ('w')

slicer({vis_template,[out_folder, 'optimism_best_sub.nii']},...
    'limits',{[],[0, 1]},...
    'colormaps', {1, 86}, ...
    'labels',{[],['\tau^*']},... % when a layer's label is empty no colorbar will be printed.
    'title','Learning rate asymmetry',...
    'output', [out_folder, 'tau.pdf'], ...
    slicer_args{:}) % two background modalities: black ('k') or white ('w')

slicer({vis_template,[out_folder, 'lr_best_sub.nii']},...
    'limits',{[],[0 ,1]},...
    'colormaps', {1, 86}, ...
    'labels',{[],['\alpha']},... % when a layer's label is empty no colorbar will be printed.
    'title','Learning rate',...
    'output', [out_folder, 'lr.pdf'], ...
    slicer_args{:}) % two background modalities: black ('k') or white ('w')

%%
end

function save_wrapper(PRF, out)
save(out, 'PRF');
end
