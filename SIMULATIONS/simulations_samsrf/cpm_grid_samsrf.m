function output = cpm_grid_samsrf(freeparams, fixedparams, data)
    % Creating a grid to precompute responses for retinotopy example. Essentially,
    % it is a look up function.
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

    appertures = data.appertures;
    resolution = fixedparams.resolution; % because 0, 0 is the middle, we need to add to
    offset = ceil(resolution / 2);

    signal = squeeze(appertures(:, freeparams.x + offset, freeparams.y + offset));

    abs_max = max(abs(signal(:))) + eps; % if everything is 0
    signal = signal(:) ./ abs_max;

    output = signal;

end
