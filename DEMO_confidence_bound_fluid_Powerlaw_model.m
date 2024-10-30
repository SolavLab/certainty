%% DEMO_confidence_bound_fluid_Powerlaw_model
% Below is a demonstration for:
% 
% * Calculating objective function values at discrete parameter values
% * Finding the certainty level
% * Running the getHessian and getJacobian codes
% * Quantifying the certainty for the best-fit results
%
% In a 2 parameter model

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
% Define the parameters
m = linspace(17.1250 - 4, 17.1250 + 4, 9); % [Pa s]
n = linspace(0.4022 - 0.1, 0.4022 + 0.1, 9); % power index

%Rheological experimental data from Laurencena & Williams (1974) for N liquid
experimentalData = [0.092321247	40.70189151
    1.206971993	14.09893698
    3.915131429	7.576043976
    12.48798749	3.770364694
    39.75474509	1.921678321
    125.3801763	0.945003122
    400.5953074	0.475582964
    1320.787205	0.217816736
    ];

shearRate = experimentalData(3:end,1);
expResults = experimentalData(3:end,2);

% ! Results are taken in a "range of confidence" over which it was judged that
%   the power-law is a plausible approximation.

%% Build objective function
% Store the variables in a cell array
variables = {m, n};

% Create a cell array with dimensions based on the lengths of the variables
dimensions = cellfun(@length, variables);
parameterCellArray = cell(dimensions);
EllisCellArray = cell(dimensions);
Z = zeros(dimensions);
% Initialize the indices and start the recursion
indices = cell(1, numel(variables));
parameterCellArray = fillParameterArray(parameterCellArray, variables, indices, 1);

for i_space = 1:numel(parameterCellArray)
    resultMatrix = zeros(length(expResults),1);
    for i = 1:length(shearRate)
        resultMatrix(i) = powerLawModel([parameterCellArray{i_space}(:)],shearRate(i));
    end
    EllisCellArray{i_space} = resultMatrix;
    Z(i_space) = sum((resultMatrix-expResults).^2,'all');
end


%% Probability settings
alpha = 0.05;
p = length(variables);
N = numel(expResults);

Fstat = finv(1-alpha,p,N-p); % F-statistic for confidence ellipse

%% Approximate confidence intervals

% Find minimum Z value and corresponding indices
[min_z, min_idx] = min(Z,[],'all');
[I1,I2] = ind2sub(size(Z),min_idx);

% Calculate confidence interval
dF = (1 + p/(N-p)*Fstat)*min_z;
s = sqrt(dF-min_z);

% Normalize C values
mid_idx = ceil(size(parameterCellArray)/2);
mid_value = parameterCellArray{mid_idx(1),mid_idx(2)};

C_norm = parameterCellArray;
for i=1:numel(parameterCellArray)
    C_norm{i} = parameterCellArray{i}./mid_value;
end

% Calculate Hessian and errors
err = zeros(p,1);
err_J = zeros(p,1);

H=getHessian(C_norm,Z,2,[I1,I2]);
A=inv(H);
for i=1:p
    err(i) = s*sqrt(2*A(i,i));
end

% Calculate Jacobian and errors
J = zeros(N,p);
f = zeros(dimensions);
for i = 1:N
    for i_space = 1:numel(parameterCellArray)
        f(i_space) = EllisCellArray{i_space}(i);
    end
    J(i,:) = getJacobian(C_norm,f,2,[I1,I2]);
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

function fittedVals = powerLawModel(x, shearRate)
    % x(1) corresponds to 'm' and x(2) corresponds to 'n'
    fittedVals = x(1) * shearRate.^(x(2) - 1);
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