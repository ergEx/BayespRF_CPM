function simulations_samsrf(REDO, basedir)
    % Script to perform simulation for retinotopy data.
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

    arguments
        REDO
        basedir = 'simulations_samsrf/'
    end
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%% SIMULATION STUDY %%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This function executes the simulation study published as (XXXXXX) by
    % Steinkamp, et al. 2023. For a more tutorial style and simpler example of CPM please see
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
    resultsdir = fullfile(basedir, 'results', filesep);

    addpath('simulations_samsrf/');
    rng(23, 'twister'); % Set random seed

    REDO = false; % REDO simulations or load previous files.
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ========================== SIMULATION SETTINGS ============================
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % First we need to simulate some data.
    % In this example we will simulate data using cpm for classical RL and for
    % distributional RL. Where classical RL has a single learning rate and
    % distributional RL has a positive and a negative learning rate.
    % Setting for the experiment
    resolution = 41; % Resolution of the simulation grid
    TR = 2.55; % Repetition time of the fMRI data, i.e. the sampling rate
    nslices = 22; % Oversampling (often the number of slices to account for slice time correction
    test_centers = [-5, 5];
    test_widths = [0.5, 2.2];
    grid_space = [-20, 20]; %
    reasonable_space = [-10, 10];  % [logit(0.0001), logit(1 - 0.0001)];
    % SNRs = [20, 10, 2, -2, -10, -20], signalvar ~ 0.0519 (estimate_sd_vor_snrs,
    % after VOI is estimated
    noise_levels = [0 0.1614    0.5103    1.2818    2.0314    5.1028   16.1363];
    fixedparams = struct('resolution', resolution);
    betap = 0.5;
    %%
    % CPM relies on it's core on a grid structure, that replaces the visual input in
    % PRF, which we can precalculate and use for both simulation and recovery, we
    % also need a stimulus structure, which we can also use for all our simulations.
    appertures = load(fullfile(basedir, 'example_data', 'aps_Bars.mat'));
    appertures = appertures.ApFrm;
    appertures_post = zeros(size(appertures, 3), resolution, resolution);
    % We will now resize and pad the appertures, so they can be better handled.
    % Because we use the default CPM model, we have to draw a inner box, which
    % limits the location. We further assume that the appertures a slightly larger
    % than the inner field.
    for ii = 1:size(appertures, 3)
        tmp_image = squeeze(appertures(:, :, ii));
        tmp_image = imresize(tmp_image, [resolution - 16, resolution - 16]);
        tmp_image = padarray(tmp_image, [8, 8], 0, 'both');
        appertures_post(ii, :, :) = tmp_image;
    end

    data = {};
    trial_n = size(appertures, 3);
    data.dur = ones(3 * trial_n, 1) .* TR;  % Duration of stimuli is as TR
    data.dt = ones(3 * trial_n, 1) .* (TR / nslices); % dt is the same across
    data.ons = [0:size(appertures, 3) - 1] .* TR;
    data.appertures = appertures_post;

    %%
    grid = struct('x', [grid_space(:)', resolution], 'y', [grid_space(:)', resolution]);
    %% Precompute
    rl_model = @cpm_grid_SAMSRF;
    % Grid for classic RL
    U = cpm_precompute(rl_model, grid, fixedparams, data, ...
                       fullfile(simulationdir, 'U_samsrf.mat'), REDO);
    U = cpm_set_constraints(U, 'x', reasonable_space);
    U = cpm_set_constraints(U, 'y', reasonable_space);
    %% Make dummy PRFs
    sim_options = struct('model', 'spm_cpm_fcn_gaussian', 'outdir', simulationdir, ...
                         'name', 'simPRF_drl', 'params', {{'x', 'y'}}, 'TE', 0.03, 'B0', 3);

    simPRF = cpm_dummy_prf(U, sim_options, TR);
    %% Simulate VOI
    simulations = struct();
    [simulations.mu_x, simulations.mu_y, ...
        simulations.width_x, simulations.width_y] = ndgrid(test_centers, test_centers, ...
                                                           test_widths, test_widths);
    %% Unpacking
    Psim =  unpack_cell_to_p_struct(simulations, betap);
    %%
    nnoise = length(noise_levels);
    nsims = length(Psim);
    nvoxels = nsims * nnoise;
    nsamps = simPRF.M.ns;
    %%
    if ~isfile(fullfile(simulationdir, 'simVOI.mat')) || REDO
        VOI.Y = zeros(nsamps, 1);
        VOI.xY.y = nan(nsamps, nvoxels);
        VOI.xY.XYZmm = zeros(3, nvoxels);

        cc = 1;
        for nn = 1:nnoise
            for ii = 1:length(Psim)
                ytmp = cpm_simulate(Psim(ii), simPRF, noise_levels(nn), false);
                VOI.xY.y(:, cc) = ytmp;
                VOI.xY.XYZmm(:, cc) = [1; nn; ii];
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
    prf_options = struct('model', 'spm_cpm_fcn_gaussian', 'name', 'samsrf', ...
                         'params', {{'x', 'y'}}, 'TE', 0.03, 'B0', 3, ...
                         'voxel_wise', true, 'avg_sess', false);
    PRF = spm_prf_analyse('specify', SPM, VOI, U, prf_options);
    %%
    options = struct('use_parfor', true, 'init', 'NONE', 'nograph', true);
    % estimate_load_prf is a simple wrapping function to check if the PRFn file
    % exists, or if the file needs to be reestimated.
    PRFn_sim = estimate_load_prf(PRF, 'PRFn', simulationdir, options, REDO);
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%% PLOTTING STARTS HERE %%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    PRFn = {};
    PRFn{1} = PRFn_sim;
    %% Plotting setting up some helper variables
    simY = VOI.xY.y;
    nparams = length(Psim);
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

    %% Plotting for shapes
    %% ==================== Recovery Plots options =================================
    ppd_samples = 250;
    plot_dimentions = 'response';
    plot_noise = 4;
    pads = 40;
    height_dims = [120, 225];
    row_height = sum(height_dims);

    %% ============================ Plot recovery Distributional ===================
    if true
        fig_x = 1600;
        fig_y = 2.0 * row_height;
        pads = 10;
        fig3 = figure('Color', 'white', 'Units', 'pixels', 'Position', ...
                      [0, 0, fig_x + pads, fig_y + 3 * pads]);
        axis('off');
        normalize_vec = [1, 1, 1, 1];

        sets = find(genp_idx == 1 &  noise_idx == plot_noise);

        dl_axes = {};

        ncols = 8;
        nrows = 2;

        cc = 1;
        for cols = 1:ncols
            for rows = 1:nrows
                hold on;

                gen_mu_x = Psim(cc).mu_x;
                gen_mu_y = Psim(cc).mu_y;
                gen_width_x = Psim(cc).width_x;
                gen_width_y = Psim(cc).width_y;

                tmp_title1 = sprintf('\\mu_{x}= %4.2f, \\sigma_{x}= %4.2f', ...
                                     gen_mu_x, gen_width_x);

                tmp_title2 = sprintf('\\mu_{y}= %4.2f, \\sigma_{y}= %4.2f', ...
                                     gen_mu_y, gen_width_y);
                %         tmpt.FontSize = 6;

                dl_axes{rows, cols} = axes('Units', 'pixels', ...
                                           'Position', [0 + (fig_x / ncols) * (cols - 1) + pads, ...
                                                        fig_y - rows * row_height - 0.25 * pads, ...
                                                        fig_x / ncols - pads, ...
                                                        height_dims(2) - 0.25 * pads] ./ ...
                                           normalize_vec);

                plot_single_voxel(PRFn{1}, sets(cc), {'y', 'x'}, ...
                                  { [], []}, {[], []}, ppd_samples, plot_dimentions);

                scatter(gen_mu_x, gen_mu_y, ...
                        'filled', 'MarkerEdgeColor', [0.5 0.5 0.5], 'MarkerFaceColor', [1 1 1],  ...
                        'LineWidth', 1.0);
                tmpt =  title({tmp_title1; tmp_title2});

                tp = -pi:0.01:pi;
                x_post = gen_mu_x + 2 * gen_width_x .* cos(tp);
                y_post = gen_mu_y + 2 * gen_width_y .* sin(tp);
                plot(x_post, y_post);
                xlabel('x');
                ylabel('y');

                xlim(PRFn{1}.U(1).grid.x(1:2));
                ylim(PRFn{1}.U(1).grid.y(1:2));

                axis('off');

                cc = cc + 1;
            end
        end

        pads_move =  [0; 0];

        for rows = 1:nrows
            for cols = 1:ncols
                dl_axes{rows, cols}.Position(2) = dl_axes{rows, cols}.Position(2) + ...
                 (pads_move(rows, 1)) / normalize_vec(2);
            end
        end

        sgtitle({'Parameter Recovery: Retinotopy', ['SNR:', num2str(snr_label(plot_noise))]});

        cpm_savefig(fig3, fullfile(resultsdir, 'fig3_parameter_recovery_retinotopy.png'));
    end

    %%
    % Classic
    generators = {Psim};

    [trues, preds] = deal({ zeros(8, nparams * nnoise)}, ...
                          {zeros(8, nparams * nnoise)});

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
    mses{1} = error_fun(reshape(trues{1}, 8, [], nnoise), reshape(preds{1}, 8, [], nnoise));
    %%
    fig4 = figure('Position', [0, 0, 1200, 600]);
    sub_titles = {'Field Recovery'};
    for nc = 1:1
        subplot(1, 1, nc);
        labels = fieldnames(cpm_get_true_parameters(PRFn{nc}, 1));
        labels = strrep(labels, '_', ' ');
        h = heatmap(round(mses{nc}, 4), 'YDisplayLabels', labels, 'XDisplayLabels', ...
                    snr_label, 'XLabel', 'SNR', 'YLabel', 'Parameter');
        title(sub_titles{nc});
    end

    sgtitle('Parameter Recovery: RMSE');

    cpm_savefig(fig4, fullfile(resultsdir, 'fig4_rmse_parameter_recovery_samsrf.png'));

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
