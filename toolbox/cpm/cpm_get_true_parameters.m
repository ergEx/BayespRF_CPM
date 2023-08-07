function [true_P] = cpm_get_true_parameters(varargin)
% % cpm_get_true_parameters, transform parameters in PRF.Ep or arbitrary parameters
% from latent into real parameter space.
%
% true_P = cpm_get_true_parameters(PRF, idx) transforms posterior parameter 
%                                            estimates for voxel = idx
% true_P = cpm_get_true_parameters(P, M, U) transforms arbitrary parameters given
%                                           a PRF.M and PRF.U structures (for transformation)
%
% Returns:
% true_P    - transformed parameters in real space.
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

if nargin==2
    PRF=varargin{1};
    voxel=varargin{2};
    P = PRF.Ep{1,voxel};
    M = PRF.M;
    U = PRF.U;
elseif nargin==3
    P =varargin{1};
    M =varargin{2};
    U =varargin{3};     
end

true_P = feval(M.IS, P, M, U, 'get_parameters');

end
