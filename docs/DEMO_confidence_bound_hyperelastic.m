%% DEMO_confidence_bound_solution_conductivity_model
% Below is a demonstration for:
% 
% * Artifically adding noise to baseline simulated data
% * Calculating objective function values at discrete parameter values
% * Finding the certainty level
% * Running the getHessian and getJacobian codes
% * Plotting a comparison between the exact confidence region to the
%   approximations
%
% In a 2 parameter model

%%
clear; close all; clc

% Plot settings
fontSize=20;
faceAlpha1=0.6;
faceAlpha2=0.2;
markerSize=6;
markerSize2=3;
lineWidth=3;

% Change default axes fonts
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize',fontSize)
set(0,'defaulttextinterpreter','latex');

%%
load('hyperelastic_demo_data.mat'); %Oddes's simulation results generated using https://github.com/SolavLab/indentify

% Extract unique p1 and p2 values
p1_values = unique(cellfun(@(x) x.p1, amit_test));
p2_values = unique(cellfun(@(x) x.p2, amit_test));


%% Probability settings
alpha = 0.05;
p = 2;
N = numel(amit_test{1}.obj_fun_val.Ff(2:end));

Fstat = finv(1-alpha,p,N-p); % F-statistic for confidence ellipse

%% Initialize variables
n = length(p1_values);
m = length(p2_values);
C = cell(n, m);
f = zeros(n, m, N);
Z = zeros(n,m);
rand_err = 2*(rand([1,5])-0.5)
% rand_err = [-0.4619, 0.5310, -0.6227, -0.4250, -0.8178]; % Uncomment for fixed error
exp_force = amit_test{ceil(numel(C)/2)}.indenter_RB_out.Fz.data(end-N+1:end)-0.0015*rand_err;

% Populate the cell array C and calculate Z
for i = 1:n
    for j = 1:m
        idx = find(cellfun(@(x) x.p1 == p1_values(i) && x.p2 == p2_values(j), amit_test));
        if ~isempty(idx)
            % Fill C with [p1, p2]
            C{i, j} = [p1_values(i), p2_values(j)];

            % Find the normalized model indenter_RB_out.Fz.data(2:end) and fill f
            f(i,j,:) = amit_test{idx}.indenter_RB_out.Fz.data(end-N+1:end)./exp_force;

            % Calculate obj_fun_val_upd.Ff
            obj_fun_vals = (squeeze(f(i,j,:))'-1).^2;

            % Calculate the sum of obj_fun_val.Ff(2:end) and fill Z
            Z(i, j) = sum(obj_fun_vals);

        end
    end
end
%% Approximate confidence intervals

% Find minimum Z value and corresponding indices
[min_z, min_idx] = min(Z,[],'all');
[row,col] = ind2sub(size(Z),min_idx);

% Calculate confidence interval
dF = (1 + p/(N-p)*Fstat)*min_z;
s = sqrt(dF-min_z);

% Normalize C values
mid_idx = ceil(size(C)/2);
mid_value = C{mid_idx(1),mid_idx(2)};

C_norm = C;
for i=1:numel(C)
    C_norm{i} = C{i}./mid_value;
end

% Calculate Hessian and errors
H=getHessian(C_norm,Z,4,[row,col]);
A=inv(H);
for i=1:length(size(C))
    err(i) = s*sqrt(2*A(i,i));
end

% Calculate Jacobian and errors
J = zeros(N,2);
for i = 1:N
    J(i,:) = getJacobian(C_norm,f(:,:,i),4,[row,col]);
end
A=inv(J'*J);
for i=1:length(size(C))
    err_J(i) = s*sqrt(A(i,i));
end

%% Plotting

% Prepare data for plotting
xData=zeros(size(C));
yData=zeros(size(C));

for i=1:numel(xData) % loop over all the elements
    xData(i)=C_norm{i}(1); % X-coordinate
    yData(i)=C_norm{i}(2); % Y-coordinate
end

% Automatically adjust contour lines
numCurves = 10;
[Fx,Fy] = gradient(Z);
grad_temp = sqrt(Fx.^2+Fy.^2);

min_grad = min(abs(grad_temp),[],'all');
contour_lines = linspace(min(Z,[],'all'), 5*dF, numCurves);
[~, std_idx] = min(abs(contour_lines - dF)); % Find the index of the value closest to dF
[~, cont_idx] = min(abs(Z - 2*dF),[],'all'); % Find the index of the point closest to dF
contour_lines(std_idx) = dF;
contour_lines(std_idx:end) = dF:(5*abs(grad_temp(cont_idx))):(dF+(numCurves-std_idx)*5*abs(grad_temp(cont_idx)));
contour_lines(1:std_idx) = linspace(min(Z,[],'all'),dF,std_idx);
if std_idx == 1
    contour_lines = [0, contour_lines];
    std_idx = 2;
    numCurves = numCurves + 1;
end

% Plot contours
h = figure;
contourf(xData,yData,Z,contour_lines) %draw 2d contours with data points
hold on
axis manual
contour(xData, yData, Z, [dF dF], 'LineColor', 'w', 'LineWidth', 2);

% Draw estimated ellipse
theta_fun = @(x,y) ([x;y]-C_norm{min_idx}');

funcH = @(x,y) 0.5*dot((theta_fun(x,y)'*H)',theta_fun(x,y));
funcJ = @(x,y) dot((theta_fun(x,y)'*(J'*J))',theta_fun(x,y));
fcontour(funcH,'--m', 'LevelList', dF-min_z, 'Visible', 'on','tag','Contours','LineWidth', 1.5);
fcontour(funcJ,'--g', 'LevelList', dF, 'Visible', 'on','tag','Contours','LineWidth', 1.5);

% Draw global minimum
plot(C_norm{min_idx}(1),C_norm{min_idx}(2),'wo', 'MarkerSize',markerSize, 'MarkerFaceColor','w')
plot(C_norm{min_idx}(1),C_norm{min_idx}(2),'kx', 'MarkerSize',markerSize, 'LineWidth', markerSize/5)

% Draw true value
plot(1,1,'wo', 'MarkerSize',markerSize, 'MarkerFaceColor','w')
plot(1,1,'k.', 'MarkerSize',markerSize2)


hold off

xlabel('$\bar{\theta}_1$')
ylabel('$\bar{\theta}_2$')
axis equal
colorbar
clim([0, (dF+(numCurves-std_idx)*5*abs(grad_temp(cont_idx)))])

% Make plot be printable in PDF
set(h, 'Units', 'Inches');
pos = get(h, 'Position');
set(h, 'PaperPositionMode', 'Auto');
set(h, 'PaperUnits', 'Inches');
set(h, 'PaperSize', [pos(3), pos(4)]);


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
