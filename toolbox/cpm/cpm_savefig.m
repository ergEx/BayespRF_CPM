function cpm_savefig(fig, fname)
    % % cpm_savefig, wrapper for better saving of functions. Using different
    % setting for pdf and png.
    %
    % fig - figure handle
    % fname - output name
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

    [~, ~, extension] = fileparts(fname);

    if strcmp(extension, '.png')

        print(fig, fname, '-dpng', '-r700', '-painters');

    elseif strcmp(extension, '.pdf')

        set(fig, 'Units', 'points');
        pos = get(fig, 'Position');
        set(fig, 'PaperPositionMode', 'Auto', 'PaperUnits', 'points', ...
            'PaperSize', [pos(3), pos(4)]);
        print(fig, fname, '-dpdf', '-r600',  '-bestfit');
    else
        error('%s - not yet implemented', extension);

    end

end
