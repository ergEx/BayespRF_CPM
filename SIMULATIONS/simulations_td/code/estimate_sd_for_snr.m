function sds = estimate_sd_for_snr(VOI, SNRdb)
    % estimate_sd_for_snr Estimates the SD for white noise, that added to the
    % pure signal in VOI should create the desired SNR in DB
    %   INPUTS:
    %       -VOI: a SPM VOI structure
    %       -SNRdb: a vector that defines for which SNRs (in DB) the standard
    %        deviation should be estimated for. Default=[20, 10, 2, -2, -10, -20]
    %   OUTPUTS:
    %       -sds: Standard deviations of noise that added to the signal should
    %        create the desired SNR
    % ---------------------------------------------------------------------
    % (c) Copyright 2023 Simon R. Steinkamp
    arguments
        VOI
        SNRdb = [20, 10, 2, -2, -10, -20]
    end

    signal_var = var(VOI.xY.y);
    mean_var = mean(signal_var);

    noise_var = mean_var .* 10.^(-SNRdb ./ 10);

    sds = sqrt(noise_var);
