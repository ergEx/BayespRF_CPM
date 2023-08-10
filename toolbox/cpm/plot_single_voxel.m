function plot_single_voxel(PRF, idx, plot_names, transformers, ...
                           transform_names, num_samples, posterior)
    % Plotting function for a single PRF with many possibilities, needs to be
    % rewritten significantly in the near future.
    % TODO: Rewrite
    %
    % For details on why one might want to disable conditional dependencies,
    % see https://en.wikibooks.org/wiki/SPM/Bayesian_Parameter_Averaging_(BPA)
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
        PRF
        idx = 1
        plot_names = {'tau', 'eta'}
        transformers = {@cpm_logit_inv, []}
        transform_names = {'alpha', 'eta'}
        num_samples = 100
        posterior = 'posterior'
    end

    if islogical(posterior)
        if posterior
            posterior = 'posterior';
        else
            posterior = 'prior';
        end
    end

    for ii = 1:length(transformers)
        if isempty(transformers{ii})
            transformers{ii} = @(x) x;
            transform_names{ii} = plot_names{ii};
        end

    end

    if length(plot_names) > 2
        error("Not implemented!");
    end

    % Recover grid:
    grid = {};
    pmus = {};
    griddim = {};
    griddim_noname = {};
    resolution = {};

    cc = 1;

    param_names = fieldnames(PRF.U(1).grid);

    for pn = param_names(:)'
        grid.(pn{1}) = PRF.U(1).grid.(pn{1});
        pmus.(pn{1}) = [PRF.M.pmin.(['mu_' pn{1}]), PRF.M.pmax.(['mu_' pn{1}])];
        griddim.(pn{1}) = linspace(grid.(pn{1})(1), ...
                                   grid.(pn{1})(2), ...
                                   grid.(pn{1})(3));
        griddim_noname{cc} = griddim.(pn{1});
        resolution{cc} = grid.(pn{1})(3);
        cc = cc + 1;
    end

    xyz = {};
    [xyz{1:length(param_names)}] = meshgrid(griddim_noname{:});

    if length(griddim_noname) < 2
        xyz{1} =  xyz{1}(1, :);
    end

    XYZ = [];

    for ii = 1:length(param_names)
        XYZ = [XYZ, xyz{ii}(:)];
    end

    %%

    if strcmp(posterior, 'posterior')
        [~, ~, ~, z] = spm_cpm_get_ppd(PRF, XYZ, idx, num_samples);
    elseif strcmp(posterior, 'prior')
        [~, ~, z] = spm_cpm_get_ppd(PRF, XYZ, idx, num_samples);
    elseif strcmp(posterior, 'response')
        z = feval(PRF.options.model, PRF.Ep{idx}, PRF.M, PRF.U, 'get_response', XYZ);
    elseif strcmp(posterior, 'prior_response')
        z = feval(PRF.options.model, PRF.M.pE{idx}, PRF.M, PRF.U, 'get_response', XYZ);

    end

    if length(resolution) < 2
        resolution{end + 1} = [];
    end
    z = reshape(z, resolution{:});

    plt_idx = nan(1, length(plot_names));

    for cc = 1:length(plot_names)
        plt_idx(cc) = find(contains(param_names, plot_names{cc}));
    end

    zdims = 1:length(param_names);

    if length(zdims) > 1
        zp = permute(z, [plt_idx, zdims(~ismember(zdims, plt_idx))]);
    else
        zp = z;
    end

    param_names = param_names([plt_idx, zdims(~ismember(zdims, plt_idx))]);
    % transform_names = transform_names( [plt_idx, zdims(~ismember(zdims, plt_idx))]);
    % transformers = transformers( [plt_idx, zdims(~ismember(zdims, plt_idx))]);

    true_ep = feval(PRF.options.model, PRF.Ep{idx}, PRF.M, PRF.U, 'get_parameters');
    % true_pe = cpm_get_true_parameters(PRF.M.pE{idx}, PRF.M, PRF.U);
    % A is given, up to 8 dims here.

    cols = repmat({':'}, 1, ndims(zp));
    for ix = length(plot_names) + 1:length(param_names)
        % slice at posterior estimate of parameter
        [~, post] = min(abs(griddim.(param_names{ix}) - true_ep.(['mu_'  param_names{ix}])));
        cols{ix} = post;
    end

    zp_slc = zp(cols{:});

    if length(plot_names) == 1
        area(transformers{1}(griddim.(plot_names{1})), zp_slc);
        xline(pmus.(plot_names{1})(1), 'LineWidth', 1.5, 'LineStyle', '--');
        xline(pmus.(plot_names{1})(2), 'LineWidth', 1.5, 'LineStyle', '--');
    else

        hold on;

        if length(zdims) < 2
            zp_slc = zp_slc * zp_slc';
        end
        pmus.(plot_names{1}) = transformers{1}(pmus.(plot_names{1}));
        pmus.(plot_names{2}) = transformers{2}(pmus.(plot_names{2}));

        rec_pos = [sum(pmus.(plot_names{1})) / 2 - diff(pmus.(plot_names{1})) / 2, ...
                   sum(pmus.(plot_names{2})) / 2 - diff(pmus.(plot_names{2})) / 2, ...
                   diff(pmus.(plot_names{1})), ...
                   diff(pmus.(plot_names{2}))];

        pplot = pcolor(transformers{1}(griddim.(plot_names{1})), ...
                       transformers{2}(griddim.(plot_names{2})), ...
                       zp_slc);
        pplot.LineWidth = 1e-10;

        rectangle('Position', rec_pos, 'EdgeColor', 'white', 'LineWidth', 1, 'LineStyle', '--');
        axis square;
    end

end
