function U = cpm_precompute(model, grid, fixedparams, data, precomp_file, overwrite)
    % cpm_precompute Summary of this function goes here
    %   INPUTS:
    %       -model: name computational model of interest. The user would
    %               implement this function. See cpm_grid_template for details.
    %       -grid: a structure that defines the grid over the parameter space
    %       (grid.parami = [parami_min parami_max N_points]. The fieldnames
    %       should correspond with the parameters used in the computational model function
    %       -fixedparams: hyperparameters that are passed to the computational model
    %       -data: stimulus data to be used for computing neuronal signals. See
    %       cpm_grid_template for details
    %       -precomp_file: filename of the resulting population field U
    %       -overwrite: whether to overwrite existing files
    %
    %   OUTPUTS:
    %       -U: population response containing all possible responses.
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
    if (exist(precomp_file, 'file') == 2) && ~overwrite
        disp(['precomputation available: ' precomp_file]);
        load(precomp_file);

    else

        linspaces = {};
        p_names = fieldnames(grid);

        for i = 1:numel(p_names)
            dims(i) = grid.(p_names{i})(3);
            linspaces{i} = linspace(grid.(p_names{i})(1), ...
                                    grid.(p_names{i})(2), ...
                                    grid.(p_names{i})(3));
            idxs = linspace(1, ...
                            grid.(p_names{i})(3), ...
                            grid.(p_names{i})(3));
            idxcell{i} = idxs;

        end

        argmatrix = combvec(linspaces{:})';
        idxmatrix = combvec(idxcell{:})';

        p_values = num2cell(argmatrix(1, :), 1);
        freeparams = cell2struct(p_values, p_names, 2);
        tt = tic;
        signaltmp = model(freeparams, fixedparams, data);
        tt = toc(tt);

        disp(' PRECOMPUTING MESHGRID USING GRID: ');
        disp(grid);
        disp(' ETA :');
        disp(duration(0, 0, tt * size(idxmatrix, 1)));

        precomputed = zeros([length(signaltmp), dims]);

        for i = 1:size(idxmatrix, 1)

            p_values = num2cell(argmatrix(i, :), 1);
            freeparams = cell2struct(p_values, p_names, 2);

            signal = model(freeparams, fixedparams, data);
            idx = num2cell(idxmatrix(i, :));
            precomputed(:, idx{:}) = signal;
        end

        microtime = true;
        U(length(signal)) = struct();
        S = {};
        S.subs = repmat({':'}, 1, ndims(precomputed));
        S.type = '()';

        for nt = 1:length(signal)
            S.subs{1} = nt;
            signals = squeeze(subsref(precomputed, S));

            signals1D = zeros(length(idxmatrix), 1);
            for i = 1:size(idxmatrix, 1)
                idx = num2cell(idxmatrix(i, :));
                signals1D(i) = signals(idx{:});
            end

            U(nt).signals = signals;
            U(nt).signals1D = signals1D;

            try
                try
                    U(nt).ons = data.ons(nt);
                    U(nt).dur = data.dur(nt);
                    U(nt).dt = data.dt(nt);

                catch
                    U(nt).ons = data(nt).ons;
                    U(nt).dur = data(nt).dur;
                    U(nt).dt = data(nt).dt;
                end
            catch
                try
                    U(nt).ons = data.ons(1);
                    U(nt).dur = data.dur(1);
                    U(nt).dt = data.dt(1);
                    microtime  = false;
                catch
                    U(nt).ons = data(1).ons;
                    U(nt).dur = data(1).dur;
                    U(nt).dt = data(1).dt;
                    microtime  = false;
                end
            end
        end

        U(1).grid = grid;
        U(1).gridpoints = argmatrix;
        U(1).grididx = idxmatrix;
        U(1).microtime = microtime;
        U(1).model = model;
        if ~microtime
            disp('microtime bins are not correctly set');
        end

        % Default pmax / pmin for parameters

        for fn = 1:length(p_names)
            U(1).pmin.(['mu_' p_names{fn}]) = grid.(p_names{fn})(1);
            U(1).pmax.(['mu_' p_names{fn}]) = grid.(p_names{fn})(2);
            [wmin, wmax] = default_width_constraints(grid.(p_names{fn}));
            U(1).pmin.(['width_' p_names{fn}]) = wmin;
            U(1).pmax.(['width_' p_names{fn}]) = wmax;
        end

        U(1).rmin = 0; % For compatibility with PRF toolbox

        % Removing nan-entries from data.
        sum_sig = [U(:).ons];
        sum_sig = find(~isnan(sum_sig));
        if ~ismember(1, sum_sig)
            sum_sig = [1, sum_sig]; % First entry in U is importatn, so we keep it here.
        end
        U = U(sum_sig);

        save(precomp_file, "U", "grid", "precomp_file");

    end

end

function [width_min, width_max] = default_width_constraints(grid_mu)
    % Grid_mu = Grid over parameters (e.g. U(1).grid.eta = [-1, 1, 10]
    % Model_mu = Reasonable parameter space.

    % Find r the minimal radius given the resolution of the grid over P (e.g.
    % covering 95 % of the space around the space.

    % Rewriting variables for understanding:
    grid_min = grid_mu(1);
    grid_max = grid_mu(2);
    grid_n = grid_mu(3);

    grid_dist = grid_max - grid_min;

    r = grid_dist / grid_n;
    width_min = r / 2;

    width_max = grid_dist / 4;

    if width_min >= width_max
        error('Sigma_min is smaller than sigma_max, try to increase grid resolution.');
    end

end
