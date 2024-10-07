function auxiliary_td()
    % Creates auxiliary plots for TD learning
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

    % Create auxiliary figures for the td learning model, especially the CSC
    % representation we are using.

    addpath(genpath('simulations_td'));

    events = readtable('simulations_td/event_file.tsv',  'FileType', 'text', 'TreatAsEmpty', 'n/a');
    trials = cpm_events_to_trials(events);
    [trials_csc, stimid] = cpm_trials_to_csc(trials, 5, 2:4, 1:3, 4);

    fixedparams = struct('lambda', 1.0, 'gamma', 0.99);
    freeparams = struct('alpha', 0.2);

    [~, trials_td] = cpm_td_learning(freeparams, fixedparams, trials_csc);

    %%
    idx = [];
    cc = 1;
    for ii = 1:length(trials_td)
        if sum(diff(trials_td(ii).wealth)) == 428
            idx(cc) = ii;
            cc = cc + 1;
        end
    end
    %%
    fig1 =  figure('Position', [0, 0, 1200, 500]);

    BLUE = [0, 0.4470, 0.7410];
    ORANGE = [0.8500, 0.3250, 0.0980];

    endidx = idx(end);
    startidx = idx(2);

    ons = trials_td(endidx).onsets;
    ons(end) = ons(4) + 2;
    ons = ons - min(ons);
    ons(1) = -1;

    subplot(1, 2, 1);

    hold on;
    stimuli = trials_td(endidx).stimuli;

    sidx = find(sum(stimuli, 2) > 0);
    sid = unique(stimid(sidx));

    colors = {BLUE, ORANGE, 'black'};

    legend_names = {};

    for kk = 1:length(sidx)

        if ~ismember(legend_names, stimid(sidx(kk)))
            legend_names{kk} = stimid{sidx(kk)};
        else
            legend_names{kk} = '';
        end
        c = colors{ismember(sid, stimid(sidx(kk)))};
        tmpstim = stimuli(sidx(kk), :);
        tmpstim(tmpstim == 0) = nan;
        stem(ons + kk / 10, tmpstim, 'Color', c, 'LineWidth', 2);
    end

    for ii = 1:length(legend_names)
        legend_names{ii} = replace(legend_names{ii}, '_', ' ');
    end

    xlabel('Trial time in s');
    ylim([0, 1.1]);
    legend(legend_names, 'Location', 'southwest');
    title('CSC-representation');

    subplot(1, 2, 2);
    hold on;

    vs = trials_td(startidx).V;
    vs(1) = 0;
    ve = trials_td(endidx).V;
    ve(1) = 0;

    rs = trials_td(startidx).RPE;
    rs(1) = 0;
    re = trials_td(endidx).RPE;
    re(1) = 0;

    plot(ons, vs,  '--o', 'LineWidth', 2, 'Color', BLUE);
    plot(ons, ve, '-o',    'LineWidth', 2, 'Color', BLUE);
    plot(ons, rs, '--o', 'LineWidth', 2, 'Color', ORANGE);
    plot(ons, re,  '-o', 'LineWidth', 2, 'Color', ORANGE);

    we = [0, diff(trials_td(endidx).wealth)];
    we(we == 0) = nan;
    stem(ons, we, 'LineWidth', 2, 'Color', 'black');

    ylabel('Value');
    xlabel('Trial time in s');
    title('TD-error and Value');

    legend({'Value (early)', 'Value (late)', 'RPE (early)', ...
            'RPE (late)', 'Reward'}, 'Location', 'northwest');

    cpm_savefig(fig1, 'simulations_td/results/csc_visualization.pdf');
