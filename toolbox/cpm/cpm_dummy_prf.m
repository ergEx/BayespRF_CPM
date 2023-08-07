function dummyPRF = cpm_dummy_prf(U, options, TR, nscans)
% % cpm_dummy_prf, creates a PRF struct using only inputs, options and number of
% samples. Mostly a utility tool for simulations.
%
% U      - Experimental stimulus timing (see spm_prf_analyse.m for details).
% options - struct including options, most important is the output dir field.
% TR     - [optional] Echotime or sampling rate of the output signal, default=2.
% nscans - [optional] Number of samples to create (if not specified uses max(U.onset) / TR), 
%           default = nan
%
% Returns:
% dummyPRF    - pRF useful for simulations
%
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
    U 
    options
    TR = 2.0
    nscans = nan
end

% As we are using a simulation we can infer the number of scans from the onsets
% and the TR:
if isnan(nscans)
    nscans = ceil(max([U(:).ons]) / TR);
end

% Set the output directory (as save path)
try options.outdir; catch options.outdir = ''; end

% The BayesPRF requires an SPM struct, but only a few fields from there, which
% we generate here:
SPM = {};
SPM.xY.RT = TR;
SPM.swd = options.outdir; 
% Create a dummy VOI structure
xY = {};
xY.y = rand(2, nscans)';
xY.XYZmm = zeros(2, 3);
VOI.Y = zeros(1, nscans)';
VOI.xY = xY;

dummyPRF = spm_prf_analyse('specify', SPM, VOI, U, options);

end

