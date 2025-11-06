function simulations_td_decay(REDO, basedir)
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

    arguments
        REDO
        basedir = 'simulations_td_decay/'
    end
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%% SIMULATION STUDY %%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This function executes the simulation study published as (XXXXXX) by
    % Steinkamp, et al. 2024. For a more tutorial style and simpler example of CPM please see
    % example_rl.m in this folder.

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%% SETUP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % adding paths
    configs = simulation_configs();

    addpath(genpath('../toolbox'));
    addpath(configs.spm_path);
    addpath(genpath(configs.vba_path));

    simulationdir = fullfile(basedir, 'simulationfiles', filesep);
    mkdir(simulationdir)
    resultsdir = fullfile(basedir, 'results', filesep);
    mkdir(resultsdir)
    addpath('simulations_td/code');
    rng(23, 'twister'); % Set random seed
    % Auxiliary functions to transform tau parameters into alpha (learning rate) space
    logit_inv = @(x) 1 ./ (1 + exp(-x)); % Inverse logit
    logit = @(x)  log(x ./ (1 - x));

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ========================== SIMULATION SETTINGS ============================
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % First we need to simulate some data.
    % In this example we will simulate data using cpm for classic TD and for
    % risk-sensitive TD. Where classic TD has a single learning rate and
    % risk-sensitivie TD has a positive and a negative learning rate.
    % Setting for the experiment
    resolution = 41; % Resolution of the simulation grid
    TR = 0.592; % Repetition time of the fMRI data, i.e. the sampling rate
    nslices = 28; % Oversampling (often the number of slices to account for slice time correction
    test_centers = logit([0.25, 0.75]);
    test_widths = [abs(test_centers(1)) / 1.5, abs(test_centers(1)) / 2.5];
    grid_space = [-8.1259, 8.1259]; % Boundaries for numerical reasons
    reasonable_space = [-4.0, 4.0];  % [logit(0.0001), logit(1 - 0.0001)];
    % SNRs = [20, 10, 2, -2, -10, -20], signalvar ~ 0.0519 (estimate_sd_vor_snrs,
    % after VOI is estimated
    noise_levels = [0  0.0076    0.0240    0.0603    0.0956    0.2402    0.7594];
    fixedparams = struct('lambda', 1.0, 'gamma', 0.99);
    fixedparams_decay = struct();
    betap = 0.5;
    %%
    % CPM relies on it's core on a grid structure, that replaces the visual input in
    % PRF, which we can precalculate and use for both simulation and recovery, we
    % also need a stimulus structure, which we can also use for all our simulations.
    filename = readtable(fullfile(basedir, 'event_file.tsv'), 'FileType', 'text', ...
                         'TreatAsEmpty', 'n/a');
    events = cpm_events_to_trials(filename);
    trials = cpm_trials_to_csc(events, 5, 2:4, 1:3, 4);
    % The U structure uses  a data format, which includes the onsets, durations etc.
    % from each onset.
    data = {};
    data.trials = trials;
    data.ons = cat(2, trials.onsets);
    data.dur = zeros(size(data.ons)) + 0.1;
    data.dt =  zeros(size(data.ons)) +  TR / nslices;
    %%
    grid_crl = struct('tau', [grid_space(:)', resolution]);
    grid_drl = struct('taupos', [grid_space(:)', resolution], ...
                      'tauneg', [grid_space(:)', resolution]);
    grid_decay = struct('xi', [grid_space(:)', resolution]);
    %% Precompute
    rl_model = @cpm_td_learning;
    % Grid for classic RL
    U_crl = cpm_precompute(rl_model, grid_crl, fixedparams, data, ...
                           fullfile(simulationdir, 'U_crl.mat'), REDO);
    U_crl = cpm_set_constraints(U_crl, 'tau', reasonable_space);
    % Grid for distributional RL
    U_drl = cpm_precompute(rl_model, grid_drl, fixedparams, data, ...
                           fullfile(simulationdir, 'U_drl.mat'), REDO);
    U_drl = cpm_set_constraints(U_drl, 'taupos', reasonable_space);
    U_drl = cpm_set_constraints(U_drl, 'tauneg', reasonable_space);
    decay_model = @cpm_decay_learning;
    U_decay = cpm_precompute(decay_model, grid_decay, fixedparams_decay, data, ...
                           fullfile(simulationdir, 'U_decay.mat'), REDO);
    U_decay  = cpm_set_constraints(U_decay, 'xi', reasonable_space);
    %% Make dummy PRFs
    sim_crl_options = struct('model', 'spm_cpm_fcn_gaussian', ...
                             'outdir', simulationdir, ...
                             'name', 'simPRF_crl', 'params', {{'tau'}}, ...
                             'TE', 0.03, 'B0', 3);
    simPRF_crl = cpm_dummy_prf(U_crl, sim_crl_options, TR);

    sim_drl_options = struct('model', 'spm_cpm_fcn_gaussian', ...
                             'outdir', simulationdir, ...
                             'name', 'simPRF_drl', 'params', {{'taupos', 'tauneg'}}, ...
                             'TE', 0.03, 'B0', 3);

    simPRF_drl = cpm_dummy_prf(U_drl, sim_drl_options, TR);

    sim_decay_options = struct('model', 'spm_cpm_fcn_gaussian', ...
                             'outdir', simulationdir, ...
                             'name', 'simPRF_decay', 'params', {{'xi'}}, ...
                             'TE', 0.03, 'B0', 3);
    simPRF_decay = cpm_dummy_prf(U_decay, sim_decay_options, TR);

    %% Simulate VOI
    simulations_crl = struct();
    [simulations_crl.mu_tau, simulations_crl.width_tau] = ndgrid(test_centers, test_widths);

    simulations_drl = struct();
    [simulations_drl.mu_taupos, ...
     simulations_drl.mu_tauneg, ...
        simulations_drl.width_taupos, ...
        simulations_drl.width_tauneg] = ndgrid(test_centers, test_centers, ...
                                               test_widths, test_widths);

    simulations_decay = struct();
    [simulations_decay.mu_xi, simulations_decay.width_xi] = ndgrid(test_centers, test_widths);

    %% Unpacking
    Psim_crl =  unpack_cell_to_p_struct(simulations_crl, betap);
    Psim_drl =  unpack_cell_to_p_struct(simulations_drl, betap);
    Psim_decay = unpack_cell_to_p_struct(simulations_decay, betap);
    %%
    nnoise = length(noise_levels);
    nsims = length(Psim_crl)  + length(Psim_drl) + length(Psim_decay);
    nvoxels = nsims * nnoise;
    nsamps = simPRF_crl.M.ns;
    %%
    if ~isfile(fullfile(simulationdir, 'simVOI.mat')) || REDO
        VOI.Y = zeros(nsamps, 1);
        VOI.xY.y = nan(nsamps, nvoxels);
        VOI.xY.XYZmm = zeros(3, nvoxels);

        cc = 1;
        for nn = 1:nnoise
            for ii = 1:length(Psim_crl)
                ytmp = cpm_simulate(Psim_crl(ii), simPRF_crl, noise_levels(nn), false);
                VOI.xY.y(:, cc) = ytmp;
                VOI.xY.XYZmm(:, cc) = [1; nn; ii];
                cc = cc + 1;
            end
            for ii = 1:length(Psim_drl)
                ytmp = cpm_simulate(Psim_drl(ii), simPRF_drl, noise_levels(nn), false);
                VOI.xY.y(:, cc) = ytmp;
                VOI.xY.XYZmm(:, cc) = [2; nn; ii];
                cc = cc + 1;
            end
            for ii = 1:length(Psim_decay)
                ytmp = cpm_simulate(Psim_decay(ii), simPRF_decay, noise_levels(nn), false);
                VOI.xY.y(:, cc) = ytmp;
                VOI.xY.XYZmm(:, cc) = [4; nn; ii];
                cc = cc + 1;
            end
        end

    else
        VOI = load(fullfile(simulationdir, 'simVOI.mat'));
        VOI = VOI.VOI;
    end

    %% PRF
    SPM = {};
    SPM.xY.RT = TR;
    SPM.swd = simulationdir;
    %%
    crl_options = struct('model', 'spm_cpm_fcn_gaussian', 'name', 'crl', ...
                         'params', {{'tau'}}, ...
                         'TE', 0.03, 'B0', 3, 'voxel_wise', true, 'avg_sess', false);
    PRF_crl = spm_prf_analyse('specify', SPM, VOI, U_crl, crl_options);

    drl_options = struct('model', 'spm_cpm_fcn_gaussian', 'name', 'drl', ...
                         'params', {{'taupos', 'tauneg'}}, ...
                         'TE', 0.03, 'B0', 3, 'voxel_wise', true, 'avg_sess', false);
    PRF_drl = spm_prf_analyse('specify', SPM, VOI, U_drl, drl_options);

    decay_options = struct('model', 'spm_cpm_fcn_gaussian', 'name', 'decay', ...
                         'params', {{'xi'}}, ...
                         'TE', 0.03, 'B0', 3, 'voxel_wise', true, 'avg_sess', false);
    PRF_decay = spm_prf_analyse('specify', SPM, VOI, U_decay, decay_options);


    %%

    options = struct('use_parfor', true, 'init', 'NONE', 'nograph', true);

    % estimate_load_prf is a simple wrapping function to check if the PRFn file
    % exists, or if the file needs to be reestimated.
    PRFn_crl = estimate_load_prf(PRF_crl, 'PRFn_crl', simulationdir, options, REDO);
    PRFn_drl = estimate_load_prf(PRF_drl, 'PRFn_drl', simulationdir, options, REDO);
    PRFn_decay = estimate_load_prf(PRF_decay, 'PRFn_decay', simulationdir, options, REDO);

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%% PLOTTING STARTS HERE %%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    PRFn = {};
    PRFn{1} = PRFn_crl;
    PRFn{2} = PRFn_drl;
    PRFn{3} = PRFn_decay;
    %% Plotting setting up some helper variables
    simY = VOI.xY.y;
    nparams_crl = length(Psim_crl);
    nparams_drl = length(Psim_drl);
    nparams_decay = length(Psim_decay);

    noise_idx = VOI.xY.XYZmm(2, :);
    genp_idx = VOI.xY.XYZmm(1, :);
    param_idx = VOI.xY.XYZmm(3, :);
    nmodels = length(PRFn);
    base = 1; % Index, with 0 noise
    signal_var = var(simY(:, noise_idx == base)); % Voxel wise, signal variance
    snrs = zeros(size(simY, 2), 1); % pre allocate snrs
    mean_snr = zeros(length(noise_levels), 1); % pre - allocate

    for nidx = 1:nnoise
        noise_var = var(simY(:, noise_idx == nidx) - simY(:, noise_idx == base));
        snrs(noise_idx == nidx) = signal_var ./ (noise_var + eps);
        mean_snr(nidx) = mean(10 * log10(snrs(noise_idx == nidx)));
    end

    snr_label = round(mean_snr, 2);
    %% ===================== Extract F values ======================================
    modelF = zeros(nmodels, length(genp_idx(:)));
    for midx = 1:nmodels
        modelF(midx, :) = PRFn{midx}.F;
    end
    %%
    same_params = [zeros(1, nparams_crl), [Psim_drl(:).mu_taupos] == [Psim_drl(:).mu_tauneg], zeros(1, nparams_decay)] == 1;
    same_params = repmat(same_params, 1, nnoise);
    extended_model_idx = genp_idx;
    extended_model_idx(same_params) = 3;
    %%
    exceedanceEP = zeros(4, nmodels, nnoise);
    vba_options.DisplayWin = 0;

    for gp = 1 : 4
        for nn = 1:nnoise
            include_vec = extended_model_idx == gp & noise_idx == nn;

            [~, o] = VBA_groupBMC(modelF(1:3, include_vec), vba_options);
            exceedanceEP(gp, 1:3, nn) = o.ep;
        end
    end

    %% Plotting for model evidence:
    if true
        ep_titles = {'Classic TD', 'Risk-sensitive TD', '\tau^- = \tau^+', 'Decay'};
        fig1 = figure('Color', 'white', 'Units', 'pixels', 'Position', [0, 0, 800, 400]);
        % Plot RFX exceedence probabilities
        for rows = 1:4
            subplot(1, 4, rows);
            bar(1:nnoise, squeeze(exceedanceEP(rows, :, :)));
            legend({'Classic TD Model', 'Risk-sensitive TD Model', 'Decay model'});
            title(ep_titles{rows});
            xlabel('SNR');
            ylabel('Exceedance Probability');
            ylim([0, 1]);
            xticklabels(snr_label);
        end

        sgtitle('Model Recovery');
        cpm_savefig(fig1, fullfile(resultsdir, 'figSD_1_model_recovery.pdf'));
    end


    if true
        plot_noise = 4;
        pads = 40;
        ppd_samples = 500;
        plot_dimentions = 'posterior';
            sets = find(genp_idx == 4 & noise_idx == plot_noise);
            % make axes
            fig_x = 500;
            fig_y = 250;

        onedim_t  = linspace(PRFn{3}.U(1).grid.xi(1), ...
                            PRFn{3}.U(1).grid.xi(2), ...
                            PRFn{3}.U(1).grid.xi(3));

        fig2 = figure('Color', 'white', 'Units', 'pixels', 'Position', ...
                      [0, 0, fig_x + pads, fig_y]);

        cc = 1;

        for rows = 1:2
            for cols = 1:2
                subplot(2, 2, cc)
                hold on;
                gen_mu_xi = Psim_decay(cc).mu_xi;
                gen_width_xi = Psim_decay(cc).width_xi;

                plot_single_voxel(PRFn{3}, sets(cc), {'xi'}, {[]}, {[]},  ...
                                  ppd_samples, plot_dimentions);
                y = normpdf(onedim_t, gen_mu_xi, gen_width_xi);
                y = y ./ sum(y);

                plot(onedim_t, y, 'LineWidth', 1.5);
                xlim(PRFn{3}.U(1).grid.xi(1:2));
                ylim([0, 1]);
                tmp_title = sprintf('\\mu_\\xi= %4.2f, \\sigma_\\xi= %4.2f', ...
                                    gen_mu_xi, gen_width_xi);
                title(tmp_title);
                ylabel(['Probability']);

                % form [left bottom width height].
                cc = cc + 1;
            end
        end

        sgtitle({'Parameter Recovery: Decay', ['SNR:', num2str(snr_label(plot_noise))]});
        cpm_savefig(fig2, fullfile(resultsdir, 'figSD_2_parameter_recovery_decay.pdf'));
    end
end

function PRFn = estimate_load_prf(PRF, prfname, savedir, options, REDO)

    prfile = fullfile(savedir, [prfname, '.mat']);

    if ~isfile(prfile) || REDO
        PRFn = spm_prf_analyse('estimate', PRF, options);
        save(prfile, 'PRFn');
    else
        PRFn = load(prfile);
        PRFn = PRFn.PRFn;
    end

end

function P =  unpack_cell_to_p_struct(simparameters, betap)

    P = struct();
    for ii = fieldnames(simparameters)'
        ptmp = simparameters.(ii{1})(:);
        for jj = 1:length(ptmp)
            P(jj).(ii{1}) = ptmp(jj);
            P(jj).beta = betap;
        end
    end

end
