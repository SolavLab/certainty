function [W] = secondDerWeights(n)

% function [W] = secondDerWeights(n)
% ------------------------------------------------------------------------
% This function calculates the weights for the second derivative finite 
% difference approximation based on the specified order of accuracy.
%
% n: Order of truncation error O(dt^n) (up to 8)
%
% W: Weights for the finite difference scheme from left (minus) to right (plus)
%
% The output weights can be used to approximate the second derivative as:
% df2dx2 = [f(t-n*dt), f(t-(n-1)*dt), ..., f(t), ..., f(t+(n-1)*dt), f(t+n*dt)] * W' / (dt^2)

    switch n
        case  2
            W = [1, -2, 1];
        case 4
            W = [-1/12, 4/3, -5/2, 4/3, -1/12];
        case 6
            W = [1/90, -3/20, 3/2, -49/18, 3/2, -3/20, 1/90];
        case 8
            W = [-1/560, 8/315, -1/5, 8/5, -205/72, 8/5, -1/5, 8/315, -1/560];
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