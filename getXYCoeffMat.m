function [N,D,B,K,W] = getXYCoeffMat(n_accuracy)

    syms x y hx hy f(x,y)
    n=n_accuracy;
    
    % Compute second-order mixed partial derivative
    df_dx_dy(x,y) = diff_t(diff_t(f,y,hy,n),x,hx,n);
    pos_arr = -(n/2):1:(n/2);
    
    for ii=1:length(pos_arr)
        for jj=1:length(pos_arr)
            % Compute weight matrix W
            W(ii,jj) = double(diff(df_dx_dy(x,y) * hx * hy, ...
                f(x+hx*pos_arr(ii),y+hy*pos_arr(jj))));
            [N(ii,jj), D(ii,jj)] = rat(W(ii,jj));
        end
    end
    
    K = max(D,[],'all');
    B =  N./D*K;
end

function [dgdt] = diff_t(g0,t,dt,n)
    g(t) = g0;
    % Differentiate by t based on accuracy level
    switch n
        case 2
            indices = -1:1;
            W = [-1/2,0,1/2];
        case 4
            indices = -2:2;
            W = [1/12, -2/3, 0, 2/3, -1/12];
        case 6
            indices = -3:3;
            W = [-1/60, 3/20, -3/4, 0, 3/4, -3/20, 1/60];
        case 8
            indices = -4:4;
            a = [1/280 -4/105, 1/5, -4/5];
            W = [a, 0, -flip(a)];
    end
    % Compute the weighted sum for differentiation
    for kk=1:length(indices)
        u(kk) = g(t+indices(kk)*dt);
    end
    % Calculate the derivative
    dgdt = u*W'/dt;
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