function [J] = getJacobian(varargin)

% function [J]=getJacobian(X,Z,n,target_ind,leap)
% ------------------------------------------------------------------------

% X is an n dimensional cell array where each element is a vector
% indicating the parameter values in space

% Z is an n dimensional array of same size as X where each element is the
% objective function value for the set of parameters indexed by X.

% n: truncation error order O(h^n) (up to 8)

% target_ind: this is the position of the parameter point within the objective function space

% leap: distances between adjacent data points for calculation differences

%% Parse input

switch nargin
    case 2
        X=varargin{1};
        Z=varargin{2};
        target_ind=ceil(size(X)/2); %pick center point in space
        n=2; %default truncation error
        leap=1; %default leap
    case 3
        X=varargin{1};
        Z=varargin{2};
        n=varargin{3};
        target_ind=ceil(size(X)/2);
        leap=1;
    case 4
        X=varargin{1};
        Z=varargin{2};
        n=varargin{3};
        target_ind=varargin{4};
        leap=1;
    case 5
        X=varargin{1};
        Z=varargin{2};
        n=varargin{3};
        target_ind=varargin{4};
        leap=varargin{5};
end

%% Compute Hessian

switch n
    case 2
        W = [-1/2,0,1/2];
    case 4
        W = [1/12, -2/3, 0, 2/3, -1/12];
    case 6
        W = [-1/60, 3/20, -3/4, 0, 3/4, -3/20, 1/60];
    case 8
        a = [1/280 -4/105, 1/5, -4/5];
        W = [a, 0, -flip(a)];
end

r = n/2; %the amount of indexes from either side of the target point
J = zeros(1,length(target_ind)); %initialize Hessian matrix

for ii = 1:length(target_ind)
    try
        C=num2cell(target_ind);
        C{ii}=':';
        local_Z = Z(C{:}); %get 1D slice!
        local_X = X(C{:}); %get 1D slice!
        p_i = target_ind(ii);
        target = cell2mat(local_X(p_i));
        step = cell2mat(local_X(p_i+leap));
        h_pi = norm(step-target);
        drv_values=local_Z((p_i-r*leap):leap:(p_i+r*leap));
        dev=0;
        for q=1:length(drv_values)
            dev=dev+drv_values(q)*W(q);
        end
        J(ii) = dev/h_pi;

    catch ME % optimal parameter set is on the boundary of the parameters range
        J = [];
        warning('Optimal parameter set is on the boundary of the parameters range');
    end
end
end


%% 
% _*certainty footer text*_ 
%
% License: <https://github.com/SolavLab/certainty/LICENSE>
%
% Copyright (C) 2024 Amit Ashkenazi
% 
% This program is free software: you can redistribute it and/or modify
% it to fit your specifications.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.