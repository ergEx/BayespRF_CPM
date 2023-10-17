%% Simple RL example for CPM
% In this example script, I will show how to set up functions and data for the CPM and provide a guide through the results etc.
% The scenario
% Suppose you have behavioral and neuroimaging data from a simple multi-arm bandit experiment and you are interested in how reward prediction errors are encoded in the brain.
% The example here is losely based on XXX and reuses some code from
% https://github.com/AnneCollins/TenSimpleRulesModeling

%%% Simulation
% As in most cases for computational modeling we first want to simulate data to see if it our methods are actually valid for the task at hand.
% We are simulating data where the participant has to select between two bandits, with different probabilities, the learning rule and the model which is used to create behavior and choices is a Rescorla-Wagner model.
%
% For reuse later in CPM, we set the Rescorla-Wagner model up in a way that it can be used with some of the key functionality.
% First we do some set-up especially loading SPM 12
% Make sure, you have the correct paths set.

%%
addpath(genpath('toolbox'));
addpath('../../toolboxes/spm12/');
spm;
close all;

%%
% First we generate some data (building on the pymc example here: https://www.pymc.io/projects/examples/en/latest/case_studies/reinforcement_learning.html ):
% We can then generate a small experiment using extreme values at first (see function here generate data function):

alpha = 0.1;
beta = 1.0;
n = 100;
p_r = [0.2, 0.8];

[choices, rewards, rpe, Q] = generate_data(alpha, beta, n, p_r);

%%
% Visualize data:
figure;
t = 1:n;

title('Choices and Q trajectories');

plot(t(choices == 1 & rewards == 1), rewards(choices == 1 & rewards == 1), 'x', 'LineWidth', 2);
hold on;
plot(t(choices == 1 & rewards == 0), rewards(choices == 1 & rewards == 0) + 1, 'o', 'LineWidth', 2);

plot(t, Q(:, 1));
plot(t, Q(:, 2));

plot(t(choices == 2 & rewards == 1), rewards(choices == 2 & rewards == 1) - 1, 'x', 'LineWidth', 2);
hold on;
plot(t(choices == 2 & rewards == 0), rewards(choices == 2 & rewards == 0), 'o', 'LineWidth', 2);
snapnow;
%%
% Having some data to work with, we can now first test our model function (hyperlink to rescorla here).
% First we need to pack everything into the correct format:
freeparams = {};
freeparams.alpha = alpha;
fixedparams = {};
fixedparams.initQ = [0.5, 0.5];
data = {};
data.choices = choices;
data.reward = rewards;

rpe_rw = rescorla_wagner_model(freeparams, fixedparams, data);

% Do some basic sanity checks:
fprintf('Error between simulation and model %4.2f\n', mean((rpe_rw - rpe).^2));

%%
% Great! We have a simulation process and a model function.
% Now we go a step further and pretend it is a whole fMRI experiment. We will need this for CPM, so we need to add some information to data.
TR = 2.0; % Pretending that the TR is 2 s
nslices = 10; % Micro time sampling
onsets = [1:n] * 10; % 10 seconds between each trial
duration = zeros(size(onsets)) + 1; % Reward stimulus is presented for 2 s
deltaT = zeros(size(onsets)) + TR / nslices;

nscans = round((max(onsets) + 10) / TR);

data.ons = onsets;
data.dur = duration;
data.dt = deltaT;
%%
% Right now the CPM toolbox provides a few ways to setup models.
%
% # Using the full PRF inspired grid.
% # Directly infering the model from the data (slightly slower).
% # Doing something akin to a GLM model.
%
% To simulate data we will use the 3 option. There are no helper functions yet so we have to do some manual work to create the U structure.

U_sim = {};

for ii = 1:n
    U_sim(ii).ons = data.ons(ii);
    U_sim(ii).dur = data.dur(ii);
    U_sim(ii).dt = data.dt(ii);
    U_sim(ii).signals1D = rpe_rw(ii);
end

sim_options = {};
sim_options.model = 'spm_cpm_fcn_fullparametric';
sim_options.params = 1;
sim_options.TE = 0.3;
dPRF = cpm_dummy_prf(U_sim, sim_options, TR, nscans);

simParams = dPRF.M.pE{1};
simParams.beta1 = 0.1;

y = cpm_simulate(simParams, dPRF, 0, true);

figure;
plot(y);
snapnow;
%%
%
% Alternatively we could simulate data using the direct inference functions (TODO)
%
% Or even using the grid - although the data here starts to diverge now.
%
% Let us set up the real CPM now.

% For precomputation we set up a grid
grid = {};

grid.alpha = [0.0, 1.0, 21]; % Minium value, maximum value and resolution.

U_cpm = cpm_precompute(@rescorla_wagner_model, grid, fixedparams, data, 'ucpm', true);

cpm_options = {};
cpm_options.model = 'spm_cpm_fcn_gaussian';
cpm_options.params = fieldnames(grid);
cpm_options.TE = 0.3;

% We again need a dummy but now declare it explicitly:
% First, we need a few information from SPM (so using this might be the
% smoothest way);
SPM = {};
SPM.xY.RT = TR;
SPM.swd = ' ';
% Create a dummy VOI structure
xY = {};
xY.y = [y, y];
xY.XYZmm = zeros(2, 3);
VOI.Y = zeros(1, nscans)';
VOI.xY = xY;

PRF = spm_prf_analyse('specify', SPM, VOI, U_cpm, cpm_options);
%%
% That's it ... we are now ready to do model inversion

PRFn = spm_prf_analyse('estimate', PRF);
snapnow;
%%
% We now might to inspect our model:

cpm_prf_review(PRFn, 1);
snapnow;

%%
% If we for example used real data for this model and have the corresponding
% SPM.mat, we could for example extract and create maps over the voxels.

%% Functions

%%% Generate_data
% This function simulates data from a Q-learning process in a bandit task.

function [choices, rewards, rpes, Qs] = generate_data(alpha, beta, n, p_r)
    if nargin < 4
        p_r = [0.4, 0.6];
    end

    choices = zeros(1, n);
    rewards = zeros(1, n);
    Qs = zeros(n, 2);

    rpes = zeros(n, 1);
    % Initialize Q table
    Q = [0.5, 0.5];
    for i = 1:n
        % Apply the Softmax transformation
        exp_Q = exp(beta * Q);
        prob_a = exp_Q / sum(exp_Q);

        % Simulate choice and reward
        a = randsample([1, 2], 1, true, prob_a);
        r = rand() < p_r(a);

        rpes(i) = (r - Q(a));
        % Update Q table
        Q(a) = Q(a) + alpha * (r - Q(a));

        % Store values
        choices(i) = a;
        rewards(i) = r;
        Qs(i, :) = Q;
    end
end

%%% The Rescorla-Wagner model

function rpe = rescorla_wagner_model(freeparams, fixedparams, data)
    % unpacking parameters
    alpha = freeparams.alpha;
    initQ = fixedparams.initQ;

    % unpacking data
    reward = data.reward;
    choices = data.choices;

    ntrials = length(reward);

    rpe = zeros(ntrials, 1);
    Q = zeros(ntrials, length(initQ));
    Q(1, :) = initQ;

    % The learning process is a simple loop over the trials, updating the Q
    % value at each step.
    for t = 1:ntrials
        rpe(t) = reward(t) - Q(t, choices(t));
        Q(t + 1, :) = Q(t, :);
        Q(t + 1, choices(t)) = Q(t, choices(t)) + alpha * rpe(t);
    end
end
