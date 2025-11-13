function first_lvl_glm_VOIS(subidx, save_path, redo)

arguments
    subidx = 1;
    save_path = 'realdata_rewardMap/preprocessed_data/';
    redo = true;
end

addpath(genpath('/mnt/projects/CPM/BayespRF_CPM/SIMULATIONS/'))
cfgs = simulation_configs();

rewardmap_derivatives = [cfgs.derivativepath filesep 'fmriprep-25.0.0/'];
rewardmap_participants = cfgs.participants;

addpath(cfgs.spm_path);


sub = rewardmap_participants{subidx};

spm_folder = fullfile(save_path, ['sub-' char(sub)], filesep);

if ~exist(spm_folder, 'dir')
    mkdir(spm_folder)
end

raw_data = sprintf("sub-%s/ses-02/func/sub-%s_ses-02_task-risksensitive_acq-te14ipat2mb2me3pf78_dir-AP_run-04_part-mag_space-MNI152NLin2009cAsym_res-02_desc-preproc_bold.nii.gz", sub, sub);

confound_file = sprintf("sub-%s/ses-02/func/sub-%s_ses-02_task-risksensitive_acq-te14ipat2mb2me3pf78_dir-AP_run-04_part-mag_desc-confounds_timeseries.tsv", sub, sub);

in_file = fullfile(spm_folder, sprintf("sub-%s_bold.nii", sub))

if ~exist(in_file)
    outzip = gunzip(fullfile(rewardmap_derivatives, raw_data), spm_folder);
    movefile(outzip{1}, in_file)
end

%%
spmbatch = {};
    %Find which session was first

    niftis = spm_select('expand', in_file);

    confound_name = fullfile(rewardmap_derivatives, confound_file);

    [confnames, confmat] = create_confounds(confound_name);

    spmbatch{1}.spm.stats.fmri_spec.sess(1).scans = cellstr(niftis);
    spmbatch{1}.spm.stats.fmri_spec.sess(1).cond = struct('name', {}, 'onset', {}, ...
        'duration', {}, 'tmod', {}, ...
        'pmod',  ...
        struct('name', {}, 'param', {}, 'poly', {}), ...
        'orth', {});

    spmbatch{1}.spm.stats.fmri_spec.sess(1).multi = {[]};
    spmbatch{1}.spm.stats.fmri_spec.sess(1).regress = struct('name', {}, 'val', {});


    counter = 1;


    for kk = 1 : length(confnames)
        spmbatch{1}.spm.stats.fmri_spec.sess(1).regress(counter).name = confnames{kk};
        spmbatch{1}.spm.stats.fmri_spec.sess(1).regress(counter).val = confmat(:, kk);
        counter=counter+1;
    end

    spmbatch{1}.spm.stats.fmri_spec.sess(1).hpf = 128;


if exist([spm_folder, 'SPM.mat'], 'file') == 2 && redo
    delete([spm_folder 'SPM.mat'])
end
%So you always overwrite, you don't rename and keep last tries?

spmbatch{1}.spm.stats.fmri_spec.dir = {spm_folder};
spmbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
spmbatch{1}.spm.stats.fmri_spec.timing.RT = 1.8;
spmbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 44; % Number of slices for oversampling
spmbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 44 / 2; % middle slice;
spmbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
spmbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0, 0];
spmbatch{1}.spm.stats.fmri_spec.volt = 1;
spmbatch{1}.spm.stats.fmri_spec.global = 'None';
spmbatch{1}.spm.stats.fmri_spec.mthresh = 0.5;
spmbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';

spm_jobman('run',spmbatch)
%%
spmbatch = {};
spm_mat = [spm_folder 'SPM.mat'];
spmbatch{1}.spm.stats.fmri_est.spmmat = {spm_mat};
spmbatch{1}.spm.stats.fmri_est.write_residuals = 0;
spmbatch{1}.spm.stats.fmri_est.method.Classical = 1;
spm_jobman('run',spmbatch)


%%
matlabbatch = {};
matlabbatch{1}.spm.util.voi.spmmat = {fullfile(spm_folder, 'SPM.mat')};
matlabbatch{1}.spm.util.voi.adjust = NaN;
matlabbatch{1}.spm.util.voi.session = 1;
matlabbatch{1}.spm.util.voi.name = 'lnacc';
matlabbatch{1}.spm.util.voi.roi{1}.mask.image = {'realdata_rewardMap/masks/L_Accumbens.nii,1'};
matlabbatch{1}.spm.util.voi.roi{1}.mask.threshold = 0.01;
matlabbatch{1}.spm.util.voi.roi{2}.mask.image = {fullfile(spm_folder, 'mask.nii')};
matlabbatch{1}.spm.util.voi.roi{2}.mask.threshold = 0.01;
matlabbatch{1}.spm.util.voi.expression = 'i1 & i2';
spm_jobman('run',matlabbatch)

%%
matlabbatch = {};
matlabbatch{1}.spm.util.voi.spmmat = {fullfile(spm_folder, 'SPM.mat')};
matlabbatch{1}.spm.util.voi.adjust = NaN;
matlabbatch{1}.spm.util.voi.session = 1;
matlabbatch{1}.spm.util.voi.name = 'rnacc';
matlabbatch{1}.spm.util.voi.roi{1}.mask.image = {'realdata_rewardMap/masks/R_Accumbens.nii,1'};
matlabbatch{1}.spm.util.voi.roi{1}.mask.threshold = 0.01;
matlabbatch{1}.spm.util.voi.roi{2}.mask.image = {fullfile(spm_folder, 'mask.nii')};
matlabbatch{1}.spm.util.voi.roi{2}.mask.threshold = 0.01;
matlabbatch{1}.spm.util.voi.expression = 'i1 & i2';

spm_jobman('run',matlabbatch)

end



function [confound_names, confound_mat] = create_confounds(conf_file)

data = readtable(conf_file, 'FileType', 'text', 'Delimiter', '\t');
filteredData = data(:, ["framewise_displacement", "trans_x","trans_y","trans_z","rot_x","rot_y","rot_z"]);

filteredData{:,:}(ismissing(filteredData{:,:})) = 0;

confound_names = filteredData.Properties.VariableNames;
confound_mat = filteredData;
confound_mat = table2array(confound_mat);
end