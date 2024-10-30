%% DEMO_confidence_bound_solution_conductivity_model
% Below is a demonstration for:
% 
% * Calculating objective function values at discrete parameter values
% * Finding the certainty level
% * Running the getHessian and getJacobian codes
% * Quantifying the certainty for the best-fit results
%
% In a 5 parameter model

%%
clear; close all; clc

% Plot settings
fontSize=7;
faceAlpha1=0.6;
faceAlpha2=0.2;
markerSize=7;

% Change default axes fonts
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize',fontSize)
set(0,'defaulttextinterpreter','latex');

%% Initialize model

k = @(P1, P2, P3, P4, n, m, T) (P1 * T + P2) * m^n * exp(-(P3 * m) / (T - P4));

T = 283.15:10:333.15;
m = [0.1,0.5,1,1.5,2,2.7,3.3,3.9,4.6];

P1 = linspace(2.8 - 1, 2.8 + 1, 9);
P2 = linspace(-663.3 - 250, -663.3 + 250, 5);
P3 = linspace(319.8 - 300, 319.8 + 300, 5);
P4 = linspace(-635.4 - 300, -635.4 + 300, 5);
n = linspace(0.977 - 4*0.03, 0.977 + 4*0.03, 5);

% Uncomment to create sliced plot
% P1 = linspace(2.8 - 1, 2.8 + 1, 3);
% P2 = linspace(-663.3 - 250, -663.3 + 250, 3);
% P3 = linspace(319.8 - 300, 319.8 + 50, 67);
% P4 = linspace(-635.4 - 200, -635.4 + 700, 67);
% n = linspace(0.977 - 4*0.03, 0.977 + 4*0.03, 3);

% Conductivity experimental results from Zhang, Chen, Wang, Wu, and Hu (2020)
% 9 rows for each m value, and 6 columns for each T value
expResults = [ 15.9288354	19.61808446	23.52433415	29.60070931	33.28995837	38.06425031
50.21704036	63.45487329	75.60766334	88.19445465	102.7343907	116.8402857
89.93057887	108.8107857	131.3802217	152.864595	176.0850727	198.8715292
115.3212019	140.2777838	166.5364689	192.795154	219.4878602	246.8316081
132.2482838	160.4600739	189.539946	214.7135484	244.8784832	275.9114602
137.2395963	166.3194484	201.6927161	228.1684019	259.8524404	292.6215317
133.3333466	164.149323	197.1354248	225.1302143	255.2951491	290.8854174
126.6059099	155.9027826	189.3229255	217.9687565	250.3038267	283.2899385
108.8107857	137.0225758	167.404531	197.569446	230.1215267	260.503482
];


% Plot theoretical model
% Create a figure with the specified width (190 mm) and font size (7 pt)
figure('Units', 'centimeters', 'Position', [0, 0, 19, 8]); % 190 mm width, 80 mm height;
% Adjust the position of the axes to shift the plot to the left
axes('Position', [0.1, 0.15, 0.65, 0.75]); % [left, bottom, width, height]

x = 0.1:0.1:4.6;

% loop for all T(i)
hold on;
for i = 2:length(T)
    k1 = arrayfun(@(x) k(2.8, -663.3, 319.8, -635.4, 0.977, x, T(i)), x);
    plot(x, k1, 'Color', [0.65    0.08    0.2], 'HandleVisibility','off'); 
    k2 = arrayfun(@(x) k(2.8, -663.3, 160, -160, 0.977, x, T(i)), x); 
    plot(x, k2, '--', 'Color', [0    0.5    0.75], 'HandleVisibility','off');
end

%theta3 = -0.33*theta4 + 110

k1 = arrayfun(@(x) k(2.8, -663.3, 319.8, -635.4, 0.977, x, T(1)), x); 
k2 = arrayfun(@(x) k(2.8, -663.3, 160, -160, 0.977, x, T(1)), x);

plot(x,k1, 'Color', [0.65    0.08    0.2], 'DisplayName', 'Best fit parameters $\hat{\theta}$');
plot(x,k2, '--', 'Color', [0    0.5    0.75], 'DisplayName', 'Arbitrary parameters');
hold off;

% Define marker types
markers = {'o', '<', '^', 'V', '>', 'd'};
hold on;
for i = 1:size(expResults, 2)
    plot(m, expResults(:, i), 'Marker', markers{i}, 'LineStyle', 'none', 'Markersize', markerSize, ...
        'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'DisplayName', [num2str(T(i)) 'K']);
end
hold off;

% Add labels and legend
xlabel('Salt Molality m (mol kg$^{-1}$)');
ylabel('Conductivity $\kappa$ (mS cm$^{-1}$)');
legend('show', 'Location', 'east', 'Interpreter', 'latex', 'FontSize', fontSize);
legend('boxoff')
ylim([0,300])
xlim([0,5])

pbaspect([2.5 1.2 1]);


%% Build objective function
% Store the variables in a cell array
variables = {P1, P2, P3, P4, n};

% Create a cell array with dimensions based on the lengths of the variables
dimensions = cellfun(@length, variables);
parameterCellArray = cell(dimensions);
conductivityCellArray = cell(dimensions);
Z = zeros(dimensions);
% Initialize the indices and start the recursion
indices = cell(1, numel(variables));
parameterCellArray = fillParameterArray(parameterCellArray, variables, indices, 1);

for i_space = 1:numel(parameterCellArray)
    resultMatrix = zeros(length(m), length(T));
    for i = 1:length(m)
        for j = 1:length(T)
            P1_temp = parameterCellArray{i_space}(1);
            P2_temp = parameterCellArray{i_space}(2);
            P3_temp = parameterCellArray{i_space}(3);
            P4_temp = parameterCellArray{i_space}(4);
            n_temp = parameterCellArray{i_space}(5);
     
            resultMatrix(i, j) = k(P1_temp, P2_temp, P3_temp, P4_temp, n_temp, m(i), T(j));
        end
    end
    conductivityCellArray{i_space} = resultMatrix;
    Z(i_space) = sum((resultMatrix-expResults).^2,'all');
end


%% Probability settings
alpha = 0.05;
p = length(variables);
N = numel(expResults);

Fstat = finv(1-alpha,p,N*3-p); % F-statistic for confidence ellipse
% Each reported measurement is the average of 3 experimental measurements!

%% Approximate confidence intervals

% Find minimum Z value and corresponding indices
[min_z, min_idx] = min(Z,[],'all');
[I1,I2,I3,I4,I5] = ind2sub(size(Z),min_idx);

% Calculate confidence interval
dF = (1 + p/(N*3-p)*Fstat)*min_z;
s = sqrt(dF-min_z);

% Normalize C values
mid_idx = ceil(size(parameterCellArray)/2);
mid_value = parameterCellArray{mid_idx(1),mid_idx(2),mid_idx(3),mid_idx(4),mid_idx(5)};

C_norm = parameterCellArray;
for i=1:numel(parameterCellArray)
    C_norm{i} = parameterCellArray{i}./mid_value;
end

% Calculate Hessian and errors
err = zeros(p,1);
err_J = zeros(p,1);

H=getHessian(C_norm,Z,4,[I1,I2,I3,I4,I5]);
A=inv(H);
for i=1:p
    err(i) = s*sqrt(2*A(i,i));
end

% Calculate Jacobian and errors
J = zeros(N,p);
f = zeros(dimensions);
for i = 1:N
    for i_space = 1:numel(parameterCellArray)
        f(i_space) = conductivityCellArray{i_space}(i);
    end
    J(i,:) = getJacobian(C_norm,f,4,[I1,I2,I3,I4,I5]);
end
J_matrix = (J'*J);
A=inv(J_matrix);
for i=1:p
    err_J(i) = s*sqrt(A(i,i));
end

%% Print results
ordinals = {'1st', '2nd', '3rd', '4th', '5th', '6th', '7th', '8th', '9th', '10th'};

for i = 1:p
    middle_value = variables{i}(ceil(length(variables{i})/2));
    certainty = err(i) * middle_value;
    if ~isreal(certainty)
        fprintf('Certainty of %s parameter is undefinable\n', ordinals{i});
    else
    fprintf('Certainty of %s parameter is %.4f ± %.4f\n', ordinals{i}, middle_value, abs(certainty));
    end
end

%% Plot slice of objective space H
% Fix dimensions 1, 2, and 5
dim1 = mid_idx(1);
dim2 = mid_idx(2);
dim5 = mid_idx(5);
% Extract the slice for dimensions 3 and 4
slice = squeeze(Z(dim1, dim2, :, :, dim5));
% Create a grid for dimensions 3 and 4
[X, Y] = meshgrid(P3, P4);

% Plot the contour map
figure('Units', 'centimeters', 'Position', [0, 0, 10.5, 9]); % 90 mm width, 90 mm height;
contour_lines = logspace(log10(10000), log10(max(Z,[],'all')), 20); % Adjusted contour levels
contourf(X, Y, slice', [0,contour_lines]); % Improved color map
% contourf(X, Y, slice);
colorbar;
h = gca;
h.ColorScale = 'log';
xlabel('$\theta_3$');
ylabel('$\theta_4$');

%% Functions
% Recursive function to fill the cell array
function cellArray = fillParameterArray(cellArray, variables, indices, currentDim)
    if currentDim > numel(variables)
        % Base case: all dimensions are filled
        positionVector = cellfun(@(var, idx) var(idx), variables, indices);
        cellArray{indices{:}} = positionVector;
    else
        % Recursive case: iterate over the current dimension
        for i = 1:length(variables{currentDim})
            indices{currentDim} = i;
            cellArray = fillParameterArray(cellArray, variables, indices, currentDim + 1);
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