function [v, outname] = cpm_write_to_nifti(image_file, coordinates, values, ...
                                           outname, set_dtype)
    % % cpm_write_to_nifti, utility function, to easily write arrays to niftis.
    %
    % http://andysbrainblog.blogspot.com/2013/06/creating-nifti-images-from-scratch.html
    %
    % image_file - reference nifti with the desired dimensions
    % coordinates - n x 3 array with the X Y Z coordinates (in cm/mni space)
    % values - n x 1 array, with corresponding values for coordinates
    % outname - name (and path) of the output file
    % set_dtype - [optional, default=64], datatype of nifti (see spm documentary)
    %
    % Returns:
    % v    - the nifti as SPM v struct
    % outname - the file name
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

    if nargin < 5
        set_dtype = 64;
    end

    % From:

    hdr = spm_vol(image_file);
    img = spm_read_vols(hdr);

    img = nan(size(img));
    hdr.fname = outname;
    hdr = rmfield(hdr, 'pinfo');
    nrows = size(coordinates, 2);

    hdr.dt(1) = set_dtype;

    mni = coordinates(:, :)';
    T = hdr.mat;

    coordinate = [mni(:, 1) mni(:, 2) mni(:, 3) ones(size(mni, 1), 1)] * (inv(T))';
    coordinate(:, 4) = [];
    coordinate = round(coordinate)';

    for i = 1:nrows

        img(coordinate(1, i), coordinate(2, i), coordinate(3, i)) = values(i);

    end

    v = spm_write_vol(hdr, img);
