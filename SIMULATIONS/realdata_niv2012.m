function [] = realdata_niv()
    % Script to perform analysis of Niv et al. 2012.
    % ---------------------------------------------------------------------
    % Copyright (C) 2024 Simon R. Steinkamp
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
    addpath(genpath(configs.vba_path));
    addpath('simulations_td/code');

    spm fmri;
    close all;

    mkdir('realdata_niv2012');
    mkdir('realdata_niv2012/models');
    % data url for Wilson & Niv 2015
    if ~isfolder('realdata_niv2012/ROI_data')
        data_url = 'https://doi.org/10.1371/journal.pcbi.1004237.s003';
        save_dir = 'realdata_niv2012/data.zip';
        websave(save_dir, data_url);
        unzip(save_dir, 'realdata_niv2012');
    end
    % ==========================================================================
    % PREPARATION
    %% ==========================================================================
    logit_inv = @(x) 1 ./ (1 + exp(-x));
    % Loading some example data:
    naccLeft = load('realdata_niv2012/ROI_data/lNAC_AllEventsData.mat');
    naccRight = load('realdata_niv2012/ROI_data/rNAC_AllEventsData.mat');

    dt = naccLeft.dt;
    TE = 30 ./ 1000; % From paper
    TR = 2.0; %

    SPM = {};
    SPM.xY.RT = TR;
    % This is where the PRF will be saved, we just set it to be here.
    SPM.swd = 'realdata_niv2012/models';
    fixedparams = struct('gamma', 0.99, 'lambda', 1.0);

    % Setting options for inversion
    invert_options = struct('use_parfor', false, ...
                            'init', 'None', ...
                            'nograph', true);
    %%
    resolution = 41;
    subjects = [];

    % Select the participants that have data
    for ii = 1:length(naccLeft.Data)
        if ~isempty(naccLeft.Data{ii})
            subjects = [subjects, ii];
        end
    end

    model_names = {'td', 'rs'};

    run_estimation = {[1, 2, 3]};

    F = zeros(length(run_estimation), 2, length(subjects),   length(model_names));
    fit_r = zeros(length(run_estimation), 2, length(subjects), length(model_names));
    fit_p = zeros(length(run_estimation), 2, length(subjects), length(model_names));

    fit_orig_r = zeros(size(fit_p));
    fit_orig_td = zeros(size(fit_p));
    fit_orig_rstd = zeros(size(fit_p));

    mse_rstd = zeros(size(fit_p));
    mse_td = zeros(size(fit_p));
    mse_fit = zeros(size(fit_p));

    learning_rates = zeros(length(run_estimation), 2, length(subjects), 3);

    PRFs = {};

    for data_idx = 1:2
        if data_idx == 1
            dataset = naccLeft;
        elseif data_idx == 2
            dataset = naccRight;
        end
        for sub_idx = 1:length(subjects)

            for re = 1:length(run_estimation)

                sub_data = dataset.Data{subjects(sub_idx)}; % Extract the data struct
                % There appear to be three runs of data, let's just use one first:
                % we have basically 1 voxel.

                [cleanData, stimOns, rewardOns, rewards, ...
                 chosen, rstdsig, tdsig, n_trials] = prep_data(sub_data, ...
                                                               run_estimation{re}, ...
                                                               dt, TR);
                % For simplicity, we can re-use the CSC representation scripts I implemented for
                % the simulation study, but for this we need a "trial" structure first.
                trials = {};
                cc = 1;
                for nt = 1:n_trials
                    if ~isnan(chosen(nt)) % We remove trials without a response
                        trials(cc).onsets = [stimOns(nt), rewardOns(nt)];
                        trials(cc).stimuli = {num2str(chosen(nt)), 'reward'};
                        trials(cc).wealth = [0, rewards(nt)];
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
                VOI = {};
                VOI.Y = cleanData;
                VOI.xY.y = cleanData;
                VOI.xY.XYZmm = zeros(3, 1);

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
                                num2str(data_idx),  '_run-', num2str(re), '_rs'];
                        grid = struct('taupos', [-8.1259, 8.1259, resolution], ...
                                      'tauneg', [-8.1259, 8.1259, resolution]);
                        params = {{'taupos', 'tauneg'}};
                        U_prf = cpm_precompute(@cpm_td_learning, ...
                                               grid, fixedparams, data, 'U_prf', true);
                        U_prf = cpm_set_constraints(U_prf, 'taupos', [-4.0, 4.0]);
                        U_prf = cpm_set_constraints(U_prf, 'tauneg', [-4.0, 4.0]);
                    elseif strcmp(model_names{mi}, 'td')
                        name = ['sub-' num2str(sub_idx), '_data-', ...
                                num2str(data_idx) '_run-' num2str(re) '_td'];
                        grid = struct('tau', [-8.1259, 8.1259, resolution]);
                        params = {{'tau'}};
                        U_prf = cpm_precompute(@cpm_td_learning, ...
                                               grid, fixedparams, data, 'U_prf', true);
                        U_prf = cpm_set_constraints(U_prf, ...
                                                    'tau', ...
                                                    [-4.0, 4.0]);
                    end

                    options = struct('model', 'spm_cpm_fcn_gaussian', ...
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

                    F(re, data_idx, sub_idx, mi) = PRFn.F;
                    PRFs{re, data_idx, sub_idx, mi} = PRFn;
                    ypred = cpm_predict(PRFn, 1);
                    [rtmp, ptmp] = corr(PRFn.Y.y, ypred);
                    fit_r(re, data_idx, sub_idx, mi) = rtmp;
                    fit_p(re, data_idx, sub_idx, mi) = ptmp;

                    fit_orig_td(re, data_idx, sub_idx, mi) = corr(PRFn.Y.y, tdsig');
                    fit_orig_rstd(re, data_idx, sub_idx, mi) = corr(PRFn.Y.y, rstdsig');
                    fit_orig_r(re, data_idx, sub_idx, mi) = corr(rstdsig', ypred);

                    mse_rstd(re, data_idx, sub_idx, mi) = mean((PRFn.Y.y - rstdsig').^2);
                    mse_td(re, data_idx, sub_idx, mi) = mean((PRFn.Y.y - tdsig').^2);
                    mse_fit(re, data_idx, sub_idx, mi) = mean((PRFn.Y.y - ypred).^2);

                    if strcmp(model_names{mi}, 'td')
                        true_params = cpm_get_true_parameters(PRFn, 1);
                        learning_rates(re, data_idx, sub_idx, 1) = logit_inv(true_params.mu_tau);
                    elseif strcmp(model_names{mi}, 'rs')
                        true_params = cpm_get_true_parameters(PRFn, 1);
                        learning_rates(re, data_idx, sub_idx, 2) = logit_inv(true_params.mu_tauneg);
                        learning_rates(re, data_idx, sub_idx, 3) = logit_inv(true_params.mu_taupos);
                    end
                end
            end
        end
    end

    %% Model comparisons
    ex_probs = zeros(4, 2, 2);
    for ridx = 1:length(run_estimation)
        for didx = 1:2
            [~, o] = VBA_groupBMC(squeeze(F(ridx, didx, :, :))', struct('DisplayWin', 0));
            ex_probs(ridx, didx, :) = o.ep;
        end
    end

    %%
    taus = zeros(length(run_estimation), 2, length(subjects));
    for ridx = 1:length(run_estimation)
        for didx = 1:2
            for sidx = 1:length(subjects)
                taus(ridx, didx, sidx) = (squeeze(learning_rates(ridx, didx, sidx, 2)) ./ ...
                                          sum(squeeze(learning_rates(ridx, didx, sidx, 2:3))));
            end
        end
    end

    %%

    fig1 = figure('Position', [0, 0, 2000, 800]);
    subplot(2, 3, 1);
    bar([squeeze(ex_probs(1, :, :))]); % , squeeze(ex_probs_sig(1,: ,:))]);
    legend({'left NAcc', 'right NAcc'}); % , 'left NAcc (only sig.)', 'right NAcc (only sig.)'})
    xticklabels({'classic TD', 'risk-sensitive TD'});
    ylabel('exceedance probability');
    title('Model comparison');
    subplot(2, 3, 2);

    scatter(squeeze(fit_r(1, 1, :, 1)), squeeze(fit_orig_td(1, 1, :, 1)), [], 'blue');
    xlabel('r CPM');
    ylabel('r regressors');
    hold on;
    scatter(squeeze(fit_r(1, 2, :, 1)), squeeze(fit_orig_td(1, 2, :, 1)), [], 'red');
    title('Classic TD model');
    h1 = lsline();
    h1(1).Color = 'blue';
    h1(2).Color = 'red';
    legend({'left NAcc', 'right NAcc'});

    subplot(2, 3, 3);

    scatter(squeeze(fit_r(1, 2, :, 1)), squeeze(fit_orig_rstd(1, 2, :, 1)), [], 'blue');
    xlabel('r CPM');
    ylabel('r regressors');
    hold on;
    scatter(squeeze(fit_r(1, 2, :, 2)), squeeze(fit_orig_rstd(1, 2, :, 2)), [], 'red');
    title('risk-sensitive TD model');
    h1 = lsline();
    h1(1).Color = 'blue';
    h1(2).Color = 'red';
    legend({'left NAcc', 'right NAcc'});

    subplot(2, 3, 4);
    scatter(ones(length(subjects), 1)', squeeze(taus(1, 1, :)));
    hold on;
    scatter(ones(length(subjects), 1)' + 1, squeeze(taus(1, 2, :)));

    for ii = 1:length(subjects)
        plot([1, 2], [taus(1, 1, ii), taus(1, 2, ii)], ...
             'Color', [0.2, 0.2, 0.2, 0.2], 'LineWidth', 0.1);
    end
    xlim([0.5, 2.5]);
    xticks([1, 2]);
    xticklabels({'left NAcc', 'right NAcc'});
    ylabel('\tau = (\alpha^+) / (\alpha^+ + \alpha^-)');
    title('Learning rate asymmetries (risk-sensitive TD)');

    subplot(2, 3, 5);

    mse_diff = sqrt(mse_rstd) - sqrt(mse_td);

    scatter(ones(length(subjects), 1)', squeeze(mse_diff(1, 1, :,  1)));
    hold on;
    scatter(ones(length(subjects), 1)' + 1, squeeze(mse_diff(1, 2, :,  1)));

    boxplot([squeeze(mse_diff(1, 1, :,  1)), squeeze(mse_diff(1, 2, :,  1))]);
    yline(0);
    xticklabels({'left NAcc', 'rightNAcc'});
    ylabel('rMSE_{RSTD} - rMSE_{TD}');
    title('MSE original models');
    subplot(2, 3, 6);

    mse_diff = sqrt(mse_fit(:, :, :, 2)) - sqrt(mse_fit(:, :, :, 1));

    scatter(ones(length(subjects), 1)', squeeze(mse_diff(1, 1, :,  1)));
    hold on;
    scatter(ones(length(subjects), 1)' + 1, squeeze(mse_diff(1, 2, :,  1)));

    boxplot([squeeze(mse_diff(1, 1, :,  1)), squeeze(mse_diff(1, 2, :,  1))]);
    xticklabels({'left NAcc', 'rightNAcc'});
    ylabel('rMSE_{RSTD} - rMSE_{TD}');
    title('MSE CPM');
    yline(0);

    sgtitle('Application to Niv, 2012');

    cpm_savefig(fig1, fullfile('realdata_niv2012', 'fig1_niv.png'));
end

function [cleanData, stimOns, rewardOns, rewards, ...
          chosen, rstdsig, tdsig, n_trials] = prep_data(sub_data, runs, dt, TR)

    cleanData = [];
    stimOns = [];
    rewardOns = [];
    rewards = [];
    chosen = [];
    rstdsig = [];
    tdsig = [];
    n_trials = 0;

    for run_idx = runs

        [n_times, n_voxels] = size(sub_data.TimeCourse{run_idx});
        nt = length(cleanData);
        % Very simple preprocessing: regressing out the motion parameter and a linear,
        % quadratic and an intercept term.
        regs = sub_data.Motion{run_idx};
        diff_motion = [zeros(1, size(regs, 2)); diff(regs)];
        intercept = ones(n_times, n_voxels);
        regs = [regs, [1:n_times]', [1:n_times]'.^2, diff_motion]; % , cs_onset', us_onset'];
        regs = [regs, intercept];
        [~, ~, cleanMRI, ~] = regress(sub_data.TimeCourse{run_idx}, regs);
        cleanData = [cleanData; cleanMRI];
        %% Create trial structrue
        % Onsets are already provided in micro time, we want them in seconds:
        stimOns = [stimOns, sub_data.CSonset{run_idx} .*  dt + nt * TR]; % Converting to seconds
        rewardOns = [rewardOns,  sub_data.USonset{run_idx} .* dt + nt * TR];
        rewards = [rewards, sub_data.Rewards{run_idx}];
        chosen = [chosen, sub_data.Chosen{run_idx}];
        n_trials = n_trials +  length(sub_data.Trials{run_idx});

        stick = zeros(1, n_times);
        stick(ceil((sub_data.TDtimes{run_idx} * dt) ./ TR)) = sub_data.RSTDerr{run_idx};
        canonical_hrf = spm_hrf(TR);
        orig_td_signal = conv(stick, canonical_hrf);
        orig_td_signal = orig_td_signal(1:n_times);
        rstdsig = [rstdsig, orig_td_signal];

        stick = zeros(1, n_times);
        stick(ceil((sub_data.TDtimes{run_idx} * dt) ./ TR)) = sub_data.TDerr{run_idx};
        canonical_hrf = spm_hrf(TR);
        orig_td_signal = conv(stick, canonical_hrf);
        orig_td_signal = orig_td_signal(1:n_times);
        tdsig = [tdsig, orig_td_signal];
    end

end

function save_wrapper(PRF, out)
    save(out, 'PRF');
end
