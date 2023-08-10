function [U] = cpm_set_constraints(U, parameter, grid_constraints)
    % (Re)setting parameter constraints, using the initial constraints saved in
    % U.
    % U    - Stimulus estimates for the subject
    % parameter     - The parameter for which constraints have to be set.
    % grid_constraints - List of length 2 with the new upper and lower limits
    %                    of the parameter grid.
    %
    % Returns:
    %   U - Stimulus tructured
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
    grid = U(1).grid.(parameter);

    % Rewriting variables for understanding:
    grid_min = grid(1);
    grid_max = grid(2);
    grid_n = grid(3);

    cons_min = grid_constraints(1);
    cons_max = grid_constraints(2);

    % Check if symmetric:
    grid_dist = grid_max - grid_min;
    cons_dist = cons_max - cons_min;

    if abs((grid_max - grid_dist / 2) - (cons_max - cons_dist / 2)) > 1e-10
        warning(['Currently only defined under the assumption of', ...
                 'symmetric parameter spaces']);
    end

    r = grid_dist / grid_n;
    width_min = r / 2;

    % Optimal parameter space:  Sigma max is defined as the distance from mu to
    % P and is large enough to cover the mu space with 95 % probability. If not
    % using an sub-optimal definition of sigma max.
    % Sigma_max = P_max - Mu_max

    if (grid_max - cons_max) >= (cons_dist / 2)
        width_max = (grid_max - cons_max)  / 2;
    else
        width_max = cons_dist / 4;
        warning(["Distance between mu and P not large enough,", ...
                 "setting sigma_max to span 95 % of mu."]);
    end

    if width_min > width_max
        error(['Sigma_min is smaller than sigma_max,', ...
               ' try to increase grid resolution.']);
    end

    U(1).pmin.(['mu_' parameter]) = cons_min;
    U(1).pmax.(['mu_' parameter]) = cons_max;
    U(1).pmin.(['width_' parameter]) = width_min;
    U(1).pmax.(['width_' parameter]) = width_max;

end
