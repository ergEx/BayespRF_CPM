function simulations_td(REDO, basedir)
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
        basedir = 'simulations_td/'
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
    %% Simulate VOI
    simulations_crl = struct();
    [simulations_crl.mu_tau, simulations_crl.width_tau] = ndgrid(test_centers, test_widths);
    simulations_drl = struct();
    [simulations_drl.mu_taupos, ...
     simulations_drl.mu_tauneg, ...
        simulations_drl.width_taupos, ...
        simulations_drl.width_tauneg] = ndgrid(test_centers, test_centers, ...
                                               test_widths, test_widths);
    %% Unpacking
    Psim_crl =  unpack_cell_to_p_struct(simulations_crl, betap);
    Psim_drl =  unpack_cell_to_p_struct(simulations_drl, betap);
    %%
    nnoise = length(noise_levels);
    nsims = length(Psim_crl)  + length(Psim_drl);
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
    %%

    options = struct('use_parfor', true, 'init', 'NONE', 'nograph', true);

    % estimate_load_prf is a simple wrapping function to check if the PRFn file
    % exists, or if the file needs to be reestimated.
    PRFn_crl = estimate_load_prf(PRF_crl, 'PRFn_crl', simulationdir, options, REDO);
    PRFn_drl = estimate_load_prf(PRF_drl, 'PRFn_drl', simulationdir, options, REDO);
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%% PLOTTING STARTS HERE %%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    PRFn = {};
    PRFn{1} = PRFn_crl;
    PRFn{2} = PRFn_drl;
    %% Plotting setting up some helper variables
    simY = VOI.xY.y;
    nparams_crl = length(Psim_crl);
    nparams_drl = length(Psim_drl);
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
    same_params = [zeros(1, nparams_crl), [Psim_drl(:).mu_taupos] == [Psim_drl(:).mu_tauneg]] == 1;
    same_params = repmat(same_params, 1, nnoise);
    extended_model_idx = genp_idx;
    extended_model_idx(same_params) = 3;
    %%
    exceedanceEP = zeros(3, nmodels, nnoise);
    vba_options.DisplayWin = 0;

    for gp = 1:3
        for nn = 1:nnoise
            include_vec = extended_model_idx == gp & noise_idx == nn;
            [~, o] = VBA_groupBMC(modelF(:, include_vec), vba_options);
            exceedanceEP(gp, :, nn) = o.ep;
        end
    end

    %% Plotting for model evidence:
    if true
        ep_titles = {'Classic TD', 'Risk-sensitive TD', '\tau^- = \tau^+'};
        fig1 = figure('Color', 'white', 'Units', 'pixels', 'Position', [0, 0, 800, 400]);
        % Plot RFX exceedence probabilities
        for rows = 1:3
            subplot(1, 3, rows);
            bar(1:nnoise, squeeze(exceedanceEP(rows, :, :)));
            legend({'Classic TD Model', 'Risk-sensitive TD Model'});
            title(ep_titles{rows});
            xlabel('SNR');
            ylabel('Exceedance Probability');
            ylim([0, 1]);
            xticklabels(snr_label);
        end

        sgtitle('Model Recovery');
        cpm_savefig(fig1, fullfile(resultsdir, 'fig1_model_recovery.png'));
    end
    %% Plotting for shapes
    %% ==================== Recovery Plots options =================================
    ppd_samples = 250;
    plot_dimentions = 'posterior';
    plot_noise = 4;
    pads = 40;
    height_dims = [120, 225];
    row_height = sum(height_dims);

    onedim_t  = linspace(PRFn{1}.U(1).grid.tau(1), ...
                         PRFn{1}.U(1).grid.tau(2), ...
                         PRFn{1}.U(1).grid.tau(3));
    %% =================== Plot Recovery classical ==================================
    if true

        sets = find(genp_idx == 1 & noise_idx == plot_noise);
        % make axes
        fig_x = 500;
        fig_y = 700;
        fig2 = figure('Color', 'white', 'Units', 'pixels', 'Position', ...
                      [0, 0, fig_x + pads, fig_y + 2 * pads]);

        cc = 1;

        for rows = 1:2
            for cols = 1:2
                cl_axes{rows, cols} = axes('Units', 'pixels', 'Position', ...
                                           [0 + (fig_x / 2) * (cols - 1) + pads, ...
                                            fig_y - row_height * (rows - 1) - ...
                                            height_dims(1) - pads, ...
                                            fig_x / 2 - pads, ...
                                            height_dims(1) - 1.5 * pads]);
                hold on;
                gen_mu_tau = Psim_crl(cc).mu_tau;
                gen_width_tau = Psim_crl(cc).width_tau;

                plot_single_voxel(PRFn{1}, sets(cc), {'tau'}, {[]}, {[]},  ...
                                  ppd_samples, plot_dimentions);
                y = normpdf(onedim_t, gen_mu_tau, gen_width_tau);
                y = y ./ sum(y);

                plot(onedim_t, y, 'LineWidth', 1.5);
                xlim(PRFn{1}.U(1).grid.tau(1:2));
                ylim([0, 1]);
                tmp_title = sprintf('\\mu_\\tau= %4.2f, \\sigma_\\tau= %4.2f', ...
                                    gen_mu_tau, gen_width_tau);
                title(tmp_title);
                xticklabels([]);
                ylabel(['Probability']);

                dl_axes{rows, cols} = axes('Units', 'pixels',  'Position', ...
                                           [0 + (fig_x / 2) * (cols - 1) + pads, ...
                                            fig_y - rows * row_height - 0.25 * pads, ...
                                            fig_x / 2 - pads, ...
                                            height_dims(2) - pads]);
                hold on;
                plot_single_voxel(PRFn{2}, sets(cc), {'taupos', 'tauneg'}, ...
                                  { [], []}, {[], []}, ppd_samples, plot_dimentions);

                xlim(PRFn{2}.U(1).grid.tauneg(1:2));
                ylim(PRFn{2}.U(1).grid.tauneg(1:2));
                xlabel('\tau^-');
                ylabel('\tau^+');
                ci = [gen_mu_tau - 2 * gen_width_tau, gen_mu_tau + 2 * gen_width_tau];
                plot(ci, ci, 'color', 'yellow', 'LineWidth', 1.5);
                plot(gen_mu_tau, gen_mu_tau, 'o', 'color', 'white', 'MarkerFaceColor', 'white');
                % form [left bottom width height].
                cc = cc + 1;
                dl_tight =  dl_axes{rows, cols}.tightPosition();
                cl_axes{rows, cols}.InnerPosition([1, 3])  = dl_tight([1, 3]);
            end

        end

        for rows = 1:2
            for cols = 1:2
                cl_axes{rows, cols}.Position(2) = cl_axes{rows, cols}.Position(2) + 2 * pads;
                dl_axes{rows, cols}.Position(2) = dl_axes{rows, cols}.Position(2) + 2 * pads;
            end
        end

        sgt = sgtitle({'Parameter Recovery: Classic TD', ['SNR:', num2str(snr_label(plot_noise))]});
        cpm_savefig(fig2, fullfile(resultsdir, 'fig2_parameter_recovery_classic_rl.png'));
    end
    %% ============================ Plot recovery Distributional ===================
    if true
        fig_x = 1600;
        fig_y = 2.1 * row_height;
        pads = 35;
        fig3 = figure('Color', 'white', 'Units', 'pixels', 'Position', ...
                      [0, 0, fig_x + pads, fig_y + 6 * pads]);
        axis('off');
        normalize_vec = [1, 1, 1, 1];

        sets = find(genp_idx == 2 &  noise_idx == plot_noise);

        cl_axes = {};
        dl_axes = {};

        ncols = 8;
        nrows = 2;

        cc = 1;
        for cols = 1:ncols
            for rows = 1:nrows
                hold on;
                cl_axes{rows, cols} = axes('Units', 'pixels',  'Position', ...
                                           [0 + (fig_x / ncols) * (cols - 1) + pads, ...
                                            fig_y - row_height * (rows - 1) - ...
                                            height_dims(1) - 1.5 * pads, ...
                                            fig_x / ncols - pads, ...
                                            height_dims(1) - 1.75 * pads] ./ normalize_vec);

                plot_single_voxel(PRFn{1}, sets(cc), {'tau'}, { []}, {[]}, ...
                                  ppd_samples, plot_dimentions);

                xlim(PRFn{1}.U(1).grid.tau(1:2));
                ylim([0, 1]);

                gen_mu_taupos = Psim_drl(cc).mu_taupos;
                gen_mu_tauneg = Psim_drl(cc).mu_tauneg;
                gen_width_taupos = Psim_drl(cc).width_taupos;
                gen_width_tauneg = Psim_drl(cc).width_tauneg;

                tmp_title1 = sprintf('\\mu_{\\tau^+}= %4.2f, \\sigma_{\\tau^+}= %4.2f', ...
                                     gen_mu_taupos, gen_width_taupos);

                tmp_title2 = sprintf('\\mu_{\\tau^-}= %4.2f, \\sigma_{\\tau^-}= %4.2f', ...
                                     gen_mu_tauneg, gen_width_tauneg);
                tmpt =  title({tmp_title1; tmp_title2});
                tmpt.FontSize = 6;
                axis('off');

                dl_axes{rows, cols} = axes('Units', 'pixels',  'Position', ...
                                           [0 + (fig_x / ncols) * (cols - 1) + pads, ...
                                            fig_y - rows * row_height - 0.25 * pads, ...
                                            fig_x / ncols - pads, ...
                                            height_dims(2) - 0.25 * pads] ./ normalize_vec);

                plot_single_voxel(PRFn{2}, sets(cc), {'taupos', 'tauneg'}, ...
                                  { [], []}, {[], []}, ppd_samples, plot_dimentions);

                scatter(gen_mu_tauneg, gen_mu_taupos, ...
                        'filled', 'MarkerEdgeColor', [0.5 0.5 0.5], ...
                        'MarkerFaceColor', [1 1 1],  'LineWidth', 1.0);

                tp = -pi:0.01:pi;
                x_post = gen_mu_tauneg + 2 * gen_width_tauneg .* cos(tp);
                y_post = gen_mu_taupos + 2 * gen_width_taupos .* sin(tp);
                plot(x_post, y_post);
                xlabel('\tau^-');
                ylabel('\tau^+');

                xlim(PRFn{2}.U(1).grid.tauneg(1:2));
                ylim(PRFn{2}.U(1).grid.taupos(1:2));
                dl_tight =  dl_axes{rows, cols}.tightPosition();
                cl_axes{rows, cols}.InnerPosition([1, 3])  = dl_tight([1, 3]);
                cc = cc + 1;
            end
        end

        pads_move =  [40; 0];  % ; 40; 0] - 15;

        for rows = 1:nrows
            for cols = 1:ncols
                cl_axes{rows, cols}.Position(2) = (cl_axes{rows, cols}.Position(2) + ...
                                                   (40 + pads_move(rows, 1)) / normalize_vec(2));
                dl_axes{rows, cols}.Position(2) = (dl_axes{rows, cols}.Position(2) + ...
                                                   (pads_move(rows, 1)) / normalize_vec(2));

            end
        end

        sgtitle({'Parameter Recovery: Risk-sensitive TD', ...
                 ['SNR:', num2str(snr_label(plot_noise))]});

        cpm_savefig(fig3, fullfile(resultsdir, 'fig3_parameter_recovery_dist_rl.png'));
    end

    %%
    % Classic
    generators = {Psim_crl, Psim_drl};

    [trues, preds] = deal({zeros(6, nparams_crl * nnoise), zeros(8, nparams_drl * nnoise)}, ...
                          {zeros(6, nparams_crl * nnoise), zeros(8, nparams_drl * nnoise)});

    cc = ones(2, 1);

    for ii = 1:length(genp_idx)

        midx = genp_idx(ii);

        tmp_true = cpm_get_true_parameters(PRFn{midx}.M.pE{ii}, PRFn{midx}.M, PRFn{midx}.U);
        fn = fieldnames(generators{midx}(param_idx(ii)));

        for fi = 1:length(fn)
            tmp_true.(fn{fi}) = generators{midx}(param_idx(ii)).(fn{fi});
        end

        tmpcc = cc(midx);
        trues{midx}(:, tmpcc) = spm_vec(tmp_true);
        preds{midx}(:, tmpcc) = spm_vec(cpm_get_true_parameters(PRFn{midx}, ii));

        cc(midx) = cc(midx)  + 1;
    end

    error_fun =  @(true, pred) squeeze(sqrt(mean((true - pred).^2, 2)));
    mses{1} = error_fun(reshape(trues{1}, 6, [], nnoise), reshape(preds{1}, 6, [], nnoise));
    mses{2} = error_fun(reshape(trues{2}, 8, [], nnoise), reshape(preds{2}, 8, [], nnoise));
    %%
    fig4 = figure('Position', [0, 0, 1200, 600]);
    sub_titles = {'Classical TD', 'Risk-sensitive TD'};
    for nc = 1:2
        subplot(1, 2, nc);
        labels = fieldnames(cpm_get_true_parameters(PRFn{nc}, 1));
        labels = strrep(labels, '_', ' ');
        h = heatmap(round(mses{nc}, 4), 'YDisplayLabels', labels, 'XDisplayLabels', snr_label, ...
                    'XLabel', 'SNR', 'YLabel', 'Parameter');
        title(sub_titles{nc});
    end

    sgtitle('Parameter Recovery: RMSE');

    cpm_savefig(fig4, fullfile(resultsdir, 'fig4_rmse_parameter_recovery.png'));
    %%
    %% Recover Tau*

    cpm_logit_inv = @(x) 1 ./ (1 + exp(-x));

    alphas_pred = zeros(nsims * nnoise, 2);
    alphas_true = zeros(nsims * nnoise, 2);

    for  kk = 1:length(genp_idx)

        tmp_pred = cpm_get_true_parameters(PRFn{2}, kk);

        if genp_idx(kk) == 1
            tmp_true = [Psim_crl(param_idx(kk)).mu_tau,  Psim_crl(param_idx(kk)).mu_tau];
        elseif genp_idx(kk) == 2
            tmp_true = [Psim_drl(param_idx(kk)).mu_taupos, Psim_drl(param_idx(kk)).mu_tauneg];

        end

        tmp_pred = [tmp_pred.mu_taupos, tmp_pred.mu_tauneg];
        alphas_true(kk, :) = cpm_logit_inv(tmp_true);
        alphas_pred(kk, :) = cpm_logit_inv(tmp_pred);
    end
    %%
    laterality_true = (alphas_true(:, 1)) ./ sum(alphas_true, 2);
    laterality_pred = (alphas_pred(:, 1)) ./ sum(alphas_pred, 2);
    laterality_error = laterality_true - laterality_pred;
    %%
    laterality_true = [reshape(laterality_true(genp_idx == 1, :), [], nnoise); ...
                       reshape(laterality_true(genp_idx == 2, :),  [], nnoise)];
    laterality_pred = [reshape(laterality_pred(genp_idx == 1, :), [], nnoise); ...
                       reshape(laterality_pred(genp_idx == 2, :),  [], nnoise)];
    laterality_error = [reshape(laterality_error(genp_idx == 1, :), [], nnoise); ...
                        reshape(laterality_error(genp_idx == 2, :),  [], nnoise)];
    %
    noise_mat = ones(size(laterality_error)) .* [1:nnoise];
    model_index_res = [reshape(genp_idx(genp_idx == 1), [], nnoise); ...
                       reshape(genp_idx(genp_idx == 2), [], nnoise)];
    %%
    fig5 = figure('Position', [0, 0, 1200, 400]);
    subplot(1, 2, 1);
    for bl = unique(laterality_true)'
        scatter(reshape(noise_mat(laterality_true(:, 1) == bl, 1:end), [], 1), ...
                reshape(abs(laterality_error(laterality_true(:, 1) == bl, 1:end)), [], 1), ...
                30 * reshape(model_index_res(laterality_true(:, 1) == bl, 1:end), [], 1), 'filled');
        hold on;
    end

    xlabel('SNR');
    ylabel('Error');

    xticklabels(mean_snr(1:end));
    legend({'\tau^*=0.25', '\tau^*=0.50', '\tau^*=0.75'});

    title('Absolute error of learning rate asymmetry');
    subplot(1, 2, 2);
    plot(laterality_true(:, plot_noise), laterality_true(:, plot_noise), ...
         'Color',  [0.5, 0.5, 0.5]);
    hold on;

    markers = {'o', 'square'};
    for nc = 1:2
        scatter(laterality_true(model_index_res(:, plot_noise) == nc, plot_noise), ...
                laterality_pred(model_index_res(:, plot_noise) == nc, plot_noise), ...
                [], 'filled', 'Marker', markers{nc});
    end
    xlabel('True values');
    ylabel('Estimated values');
    xlim([0.0, 1.0]);
    ylim([0.0, 1.0]);
    title(sprintf('Estimation at SNR %4.2f', snr_label(plot_noise)));
    legend({'', 'Classic TD', 'Risk-sensitive TD'}, 'Location', 'northwest');

    cpm_savefig(fig5, fullfile(resultsdir, 'fig5_learning_assymetry_tau.png'));

    %% %% Classic BPA
    fig6 = figure('Position', [0, 0, 2000, 800]);
    prf_names = {'Classic TD', 'Risk-sensitive TD'};
    tiledlayout(2, 7, 'TileSpacing', 'tight', 'Padding', 'compact');

    for nc = 1:2
        for nn = 1:7
            nexttile();
            included =  find(genp_idx == nc & noise_idx == nn);
            nincluded = length(included);
            GCM = cell(nincluded, 1);
            i = 1;
            for v = included
                GCM{i}.Cp = PRFn{nc}.Cp{v};
                GCM{i}.Ep = PRFn{nc}.Ep{v};
                GCM{i}.M.pC = PRFn{nc}.M.pC{v};
                GCM{i}.M.pE = PRFn{nc}.M.pE{v};
                i = i + 1;
            end
            classic_BPA = spm_dcm_bpa(GCM);

            cl_labels = fieldnames(PRFn{nc}.M.pE{1});
            cl_labels = strrep(cl_labels, '_', ' ');
            tmp_mat = VBA_cov2corr(classic_BPA.Cp);
            idx = tril(tmp_mat);
            tmp_mat(~idx) = nan;
            t = heatmap(tmp_mat, 'MissingDataColor', 'w', 'GridVisible', ...
                        'off', 'MissingDataLabel', " ", 'ColorbarVisible', ...
                        'off', 'Colormap', colormap('parula'), 'XDisplayLabels', cl_labels, ...
                        'YDisplayLabels', cl_labels, 'FontSize', 10 - 3 *  nc, ...
                        'CellLabelFormat', '%0.2f');
            t.InnerPosition = [0, 0, 1, 1];

            if nc == 1
                title(sprintf('SNR %4.2f', snr_label(nn)));
            end
            if nn == 2
                ylabel(prf_names{nc});
            end

        end
    end

    sgtitle('Posterior Correlation after BPA');

    cpm_savefig(fig6, fullfile(resultsdir, 'fig6_posterior_correlation_bpa.png'));

    %%
    %% Predicted Y
    predYs = {};
    parfor nc =  1:nmodels
        predYs{nc} = zeros(size(PRFn{nc}.Y.y));
        for vx = 1:nnoise * nsims
            predYs{nc}(:, vx) = feval(PRFn{nc}.M.IS, PRFn{nc}.Ep{vx}, PRFn{nc}.M, PRFn{nc}.U);
        end
    end

    %% Simulated HRF accuracy RMSE / R2

    fits_r2 = zeros(2, nsims * nnoise);
    fits_mse = zeros(2, nsims * nnoise);

    for nc = 1:2
        for vx = 1:nnoise * nsims
            fits_r2(nc, vx) = VBA_r2(predYs{nc}(:, vx), simY(:, vx));
            fits_mse(nc, vx) = sqrt(mean((simY(:, vx) - predYs{nc}(:, vx)).^2));
        end
    end

    %% Plot classical RL

    fig7 = figure('Position', [0, 0, 1200, 500]);

    for nc = 1:2
        subplot(2, 2, nc);
        tmp_r2 = (fits_r2(:, genp_idx == nc));
        tmp_r2 =  [tmp_r2(1, :), tmp_r2(2, :)];
        model_g = [zeros(1, sum(genp_idx(:) == nc)), zeros(1, sum(genp_idx(:) == nc)) + 1];
        noise_g = [noise_idx(genp_idx(:) == nc), noise_idx(genp_idx(:) == nc)];

        boxchart(noise_g(noise_g > 1) - 2, tmp_r2(noise_g > 1), ...
                 'GroupByColor', model_g(noise_g > 1));
        legend({'Classic TD', 'Risk-sensitive TD'}, 'Location', 'SouthEast');
        xticks(0:5);
        xticklabels(snr_label(2:7));
        xlabel('SNR');
        ylabel('R^2');
        title('Generative Process:', prf_names{nc});
    end
    %
    subplot(2, 2, 3);

    % Fits of interest: At Model SNR largest difference
    tmp_idx =  find(genp_idx == 2 & noise_idx == plot_noise);
    [~, max_diff] = max(diff(fits_r2(:, tmp_idx))');
    max_diff = max_diff(1);
    hold on;
    t = (0:size(simY, 1) - 1) * 0.592;
    lh  = plot(t, simY(:, tmp_idx(max_diff)), 'Color', 'black');
    lh.Color(4) = 0.1;
    lh2 = plot(t, predYs{1}(:, tmp_idx(max_diff)), 'LineStyle', ':', 'LineWidth', 1.5);
    lh3 = plot(t, predYs{2}(:, tmp_idx(max_diff)), 'Color', 'black');

    lh2.Color  = [0, 0.4470, 0.7410, 1];
    lh3.Color = [0.8500    0.3250    0.0980, 1];

    xlabel('time');
    ylabel('signal, a.u.');
    title('Simulated and Recovered Signals', sprintf('SNR %4.2f', snr_label(plot_noise)));
    legend({'Simulated', 'Classic TD', 'Risk-sensitive TD'}, 'Location', 'SouthEast');

    subplot(2, 2, 4);
    hold on;

    sc1 = scatter(simY(:, tmp_idx(max_diff)), predYs{1}(:, tmp_idx(max_diff)), 0.8, 'filled');
    sc2 = scatter(simY(:, tmp_idx(max_diff)), predYs{2}(:, tmp_idx(max_diff)), 0.8, 'filled');

    l2  = lsline;

    l2(1).Color = sc1.CData;
    l2(2).Color = sc2.CData;

    xlim([min(simY(:, tmp_idx(max_diff))), max(simY(:, tmp_idx(max_diff)))]);
    ylim([min(predYs{2}(:, tmp_idx(max_diff))), max(predYs{2}(:, tmp_idx(max_diff)))]);

    ylabel('Simulated');
    xlabel('Recovered');
    rf = refline(0.5, 0);
    rf.Color = [0.2, 0.2, 0.2, 0.4];

    title_r2 = round(fits_r2(:, tmp_idx(max_diff)), 2);

    legend({['Classic TD (R^2 = ', num2str(round(title_r2(1), 2)) ')'], ...
            ['Risk-sensitive TD (R^2 = ', num2str(round(title_r2(2), 2)) ')']}, ...
           'Location', 'SouthEast');

    disp(generators{2}(max_diff));

    title('Simulated vs Recovered');

    sgtitle('Classical Model fit');

    cpm_savefig(fig7, fullfile(resultsdir, 'fig7_classic_model_fit.png'));

    %%
    cpm_prf_review(PRFn{2}, 110);
    %%
    sig = zeros(nmodels,  nsims * nnoise);

    for jj = 1:nsims * nnoise
        for nc = 1:nmodels
            sig(nc, jj) = feval(PRFn{nc}.M.IS, PRFn{nc}.Ep, PRFn{nc}.M, PRFn{nc}.U, ...
                                'is_above_threshold', PRFn{nc}.Cp, jj, 0.9);
        end
    end

    %% figure;
    %% mean_sig
    mean_sig = zeros(2, nmodels, nnoise);

    for ii = 1:nnoise
        for jj = 1:nmodels
            mean_sig(1, jj, ii) = mean(sig(jj, genp_idx == 1 & noise_idx == ii));
            mean_sig(2, jj, ii) = mean(sig(jj, genp_idx == 2 & noise_idx == ii));
        end
    end

    figure;
    subplot(1, 2, 1);
    bar(squeeze(mean_sig(1, :, :))');
    xticklabels(snr_label(1:7));

    subplot(1, 2, 2);
    bar(squeeze(mean_sig(2, :, :))');
    xticklabels(snr_label(1:7));

    %%
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
