function U = cpm_prepare_input_direct_model(data, fixedparams, model, constraints)
    % cpm_prepare_input_direct Creates an input structure for direct models of CPM,
    % Mostly re-allocating onsets etc. to the required U structure.
    %   INPUTS:
    %       -data: stimulus data to be used for computing neuronal signals. See
    %       cpm_grid_template for detail
    %       -fixedparams: hyperparameters that are passed to the computational model
    %       -model: name computational model of interest. The user would
    %               implement this function. See cpm_grid_template for details.
    %      - contraints: contraints on the parameterspace
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

    U = struct();
    for ii = 1:length(data.ons)
        U(ii).ons = data.ons(ii);
        U(ii).dur = data.dur(ii);
        U(ii).dt = data.dt(ii);
    end
    U(1).model = model;
    U(1).fixedparams = fixedparams;
    U(1).data = data;

    param_names = fieldnames(constraints);

    for fname = 1:length(param_names)
        U(1).pmin.(param_names{fname}) = constraints.(param_names{fname})(1);
        U(1).pmax.(param_names{fname}) = constraints.(param_names{fname})(2);
    end

    % For compatibility
    U(1).rmin = 1;
