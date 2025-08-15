%% Technical University of Catalonia (UPC)
%% Higher Technical School of Industrial Engineering of Barcelona (ETSEIB)
%% Centre of Technological Innovation in Static Converters and Drives (CITCEA)
%% Doctoral Program in Electrical Engineering
%% Developed by: Luis Angel Garcia Reyes, MSc
%% Eigenvalue analysis tool for impedance-based analysis

% GNU General Public License v3.0 (GPL-3.0)
% Copyright (C) 2025 Luis Angel Garcia Reyes, UPC-MSCA-ADOreD
% Email: luis.reyes@upc.edu
% This program is free software: you can redistribute it and/or modify it 
% under the terms of the GNU General Public License as published by the 
% Free Software Foundation, either version 3 of the License, or (at your 
% option) any later version.
% This program is distributed in the hope that it will be useful, but 
% WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
% or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License 
% for more details.
% You should have received a copy of the GNU General Public License along 
% with this program. If not, see <https://www.gnu.org/licenses/>.

% This work has received funding from the ADOreD project
% under the European Union’s Horizon Europe Research and 
% Innovation Programme under the Marie Skłodowska-Curie 
% Grant Agreement No. 101073554.


%% Eigenvalue descomposition

disp('---')
disp('SIaD Tool is applying the Impedance Modal Analysis.')
disp('The results are the next:')
disp('---')

% Scanner type
switch scanner_type
    case 1
%% Voltage perturbation
% Selection of reference frame:
switch scanner_selector
    case 1
%% ABC sequence
% Frequency vector of measurements: fd0
% System 1: Y_abc1
% System 2: Y_abc2
Z_sys_full = zeros(3,3,length(fd0)); % Initial zero L open loop matriz
for n=1:length(fd0)
    Y_system1 = Y_abc1(:,:,n);
    Y_system2 = Y_abc2(:,:,n);
    Y_sys_full(:,:,n) = Y_system1+Y_system2;
end
ZabcModalAnalysis(Y_sys_full,fd0,{'Sys1','Sys2'});
    case 2
%% dq0 sequence
% Set initial data and matrices:
% Frequency vector of measurements: fd0
% System 1: Yqq1, Yqd1, Ydq1 and Ydd1
% System 2: Yqq2, Yqd2, Ydq2 and Ydd2
Y_sys_full = zeros(2,2,length(fd0)); % Initial zero L open loop matriz
for n=1:length(fd0)
    Y_system1 = [Yqq1(n), Yqd1(n); 
                 Ydq1(n), Ydd1(n)];
    Y_system2 = [Yqq2(n), Yqd2(n); 
                 Ydq2(n), Ydd2(n)];
    Y_sys_full(:,:,n) = (Y_system1+Y_system2);
end
Zdq0ModalAnalysis(Y_sys_full,fd0,{'Sys1','Sys2'});
    case 3
%% 0pn sequence
% Frequency vector of measurements: fd0
% System 1: Y_0pn1
% System 2: Y_0pn2
Z_sys_full = zeros(3,3,length(fd0)); % Initial zero L open loop matriz
for n=1:length(fd0)
    Y_system1 = Y_0pn1(:,:,n);
    Y_system2 = Y_0pn2(:,:,n);
    Y_sys_full(:,:,n) = Y_system1+Y_system2;
end
ZabcModalAnalysis(Y_sys_full,fd0,{'Sys1','Sys2'});
end
    
    case 2
%% Current perturbation
% Selection of reference frame:
switch scanner_selector
    case 1
%% ABC sequence
% Frequency vector of measurements: fd0
% System 1: Z_abc1
% System 2: Z_abc2
Z_sys_full = zeros(3,3,length(fd0)); % Initial zero L open loop matriz
for n=1:length(fd0)
    Y_system1 = inv(Z_abc1(:,:,n));
    Y_system2 = inv(Z_abc2(:,:,n));
    Y_sys_full(:,:,n) = Y_system1+Y_system2;
end
ZabcModalAnalysis(Y_sys_full,fd0,{'Sys1','Sys2'});
    case 2
%% dq0 sequence
% Set initial data and matrices:
% Frequency vector of measurements: fd0
% System 1: Zqq1, Zqd1, Zdq1 and Zdd1
% System 2: Zqq2, Zqd2, Zdq2 and Zdd2
Z_sys_full = zeros(2,2,length(fd0)); % Initial zero L open loop matriz
for n=1:length(fd0)
    Y_system1 = inv([Zqq1(n), Zqd1(n); 
                     Zdq1(n), Zdd1(n)]);
    Y_system2 = inv([Zqq2(n), Zqd2(n); 
                     Zdq2(n), Zdd2(n)]);
    Y_sys_full(:,:,n) = Y_system1+Y_system2;
end
Zdq0ModalAnalysis(Y_sys_full,fd0,{'Sys1','Sys2'});
    case 3
%% 0pn sequence
% Frequency vector of measurements: fd0
% System 1: Z_0pn1
% System 2: Z_0pn2
Z_sys_full = zeros(3,3,length(fd0)); % Initial zero L open loop matriz
for n=1:length(fd0)
    Y_system1 = inv(Z_0pn1(:,:,n));
    Y_system2 = inv(Z_0pn2(:,:,n));
    Z_sys_full(:,:,n) = inv(Y_system1+Y_system2);
end
ZabcModalAnalysis(Z_sys_full,fd0,{'Sys1','Sys2'});
end
end

function Zdq0ModalAnalysis(A, frequencies, subsys_labels)
    % A: 2x2xN complex matrix, Y(f) for each frequency
    % frequencies: 1xN vector of frequencies (Hz or p.u.)
    % subsys_labels: labels for the inputs, e.g., {'Sys1','Sys2'}

    [n, m, N] = size(A);
    if n ~= 2 || m ~= 2
        error('Matrix must be of size 2x2xN');
    end
    if length(frequencies) ~= N
        error('Frequency vector must have length N');
    end
    if length(subsys_labels) ~= 2
        error('Two labels are required for the subsystems');
    end

    eigenvalues = zeros(2, N);
    participation_factors = zeros(2, 2, N);
    modal_impedance = zeros(2, N);
    critical_mode_idx = zeros(1, N);
    dominant_subsystem = strings(1, N);

    % === Eigenvalue plot ===
    figure('Name','Eigenvalues in the Complex Plane','NumberTitle','off');
    hold on; grid on;
    xlabel('Re(\lambda)'); ylabel('Im(\lambda)');
    title('Eigenvalues of Y(f) across Frequencies');
    xline(0, '--r', 'Stability Boundary (Re=0)', 'LineWidth', 2);

    for i = 1:N
        Yf = A(:, :, i);
        [V, D] = eig(Yf);
        W = inv(V);
        lambda = diag(D);
        eigenvalues(:, i) = lambda;

        % Modal impedance = 1 / lambda
        modal_impedance(:, i) = abs(1 ./ lambda);

        % Participation factors: L .* R
        participation_factors(:, :, i) = abs(W' .* V);

        % Identify critical mode (minimum |lambda|)
        [~, idx] = min(abs(lambda));
        critical_mode_idx(i) = idx;

        % Determine dominant subsystem
        pf = participation_factors(idx, :, i);
        if pf(1) > pf(2)
            dominant_subsystem(i) = subsys_labels{1};
        elseif pf(2) > pf(1)
            dominant_subsystem(i) = subsys_labels{2};
        else
            dominant_subsystem(i) = "Undetermined";
        end

        % Plot eigenvalues as circles
        plot(real(lambda), imag(lambda), 'o', 'MarkerSize', 8, ...
             'DisplayName', sprintf('f = %.2f', frequencies(i)));
    end
    legend show;
    hold off;

    % === Heatmap of participation factors of the critical mode ===
    heatmap_data = zeros(N, 2);
    row_labels = strings(N, 1);
    for i = 1:N
        pf = participation_factors(critical_mode_idx(i), :, i);
        heatmap_data(i, :) = pf;
        row_labels(i) = sprintf('f = %.2f', frequencies(i));
    end

    figure('Name','Critical Mode Participation Factors','NumberTitle','off');
    h = heatmap({'Input 1','Input 2'}, row_labels, heatmap_data, ...
                'Colormap', flipud(parula), ...
                'ColorbarVisible','on');
    h.Title = 'Participation Factors of Critical Mode';
    h.XLabel = 'Inputs';
    h.YLabel = 'Frequency';

    % === Summary table of metrics ===
    fprintf('\nModal Summary by Frequency:\n');
    fprintf('--------------------------------------------------------------------------------------------------\n');
    fprintf('| f (Hz/p.u.) | Critical Mode | Eigenvalue        | |Z_m|     | PF1     | PF2     | Dominant      |\n');
    fprintf('--------------------------------------------------------------------------------------------------\n');
    for i = 1:N
        idx = critical_mode_idx(i);
        lambda = eigenvalues(idx, i);
        Zm = modal_impedance(idx, i);
        pf = participation_factors(idx, :, i);
        fprintf('| %9.3f |      %d       | %8.4f%+8.4fj | %8.3f | %7.4f | %7.4f | %-13s |\n', ...
            frequencies(i), idx, real(lambda), imag(lambda), Zm, pf(1), pf(2), dominant_subsystem(i));
    end
    fprintf('--------------------------------------------------------------------------------------------------\n');
end

function ZabcModalAnalysis(A, frequencies, subsys_labels)
    % A: 3x3xN complex matrix, Y(f) for each frequency
    % frequencies: 1xN vector of frequencies (Hz or p.u.)
    % subsys_labels: labels for the two subsystems, e.g., {'Sys1','Sys2'}

    [n, m, N] = size(A);
    if n ~= 3 || m ~= 3
        error('Matrix must be of size 3x3xN');
    end
    if length(frequencies) ~= N
        error('Frequency vector must have length N');
    end
    if length(subsys_labels) ~= 2
        error('Two labels are required for the subsystems');
    end

    eigenvalues = zeros(3, N);
    participation_factors = zeros(3, 3, N);
    modal_impedance = zeros(3, N);
    critical_mode_idx = zeros(1, N);
    dominant_subsystem = strings(1, N);

    % === Eigenvalue plot ===
    figure('Name','Eigenvalues in the Complex Plane','NumberTitle','off');
    hold on; grid on;
    xlabel('Re(\lambda)'); ylabel('Im(\lambda)');
    title('Eigenvalues of Y(f) across Frequencies');
    xline(0, '--r', 'Stability Boundary (Re=0)', 'LineWidth', 2);

    for i = 1:N
        Yf = A(:, :, i);
        [V, D] = eig(Yf);
        W = inv(V);
        lambda = diag(D);
        eigenvalues(:, i) = lambda;

        % Modal impedance = 1 / lambda
        modal_impedance(:, i) = abs(1 ./ lambda);

        % Participation factors: L .* R
        participation_factors(:, :, i) = abs(W' .* V);

        % Identify critical mode (minimum |lambda|)
        [~, idx] = min(abs(lambda));
        critical_mode_idx(i) = idx;

        % Determine dominant subsystem (based on participation factors)
        pf = participation_factors(idx, :, i);
        pf_sys1 = sum(pf(1:3/2));  % Assuming Sys1 contributes to first half
        pf_sys2 = sum(pf(3/2+1:end));  % Sys2 contributes to second half

        if pf_sys1 > pf_sys2
            dominant_subsystem(i) = subsys_labels{1};
        elseif pf_sys2 > pf_sys1
            dominant_subsystem(i) = subsys_labels{2};
        else
            dominant_subsystem(i) = "Undetermined";
        end

        % Plot eigenvalues as circles
        plot(real(lambda), imag(lambda), 'o', 'MarkerSize', 8, ...
             'DisplayName', sprintf('f = %.2f', frequencies(i)));
    end
    legend show;
    hold off;

    % === Heatmap of participation factors of the critical mode ===
    heatmap_data = zeros(N, 3);
    row_labels = strings(N, 1);
    for i = 1:N
        pf = participation_factors(critical_mode_idx(i), :, i);
        heatmap_data(i, :) = pf;
        row_labels(i) = sprintf('f = %.2f', frequencies(i));
    end

    figure('Name','Critical Mode Participation Factors','NumberTitle','off');
    h = heatmap({'Input 1','Input 2','Input 3'}, row_labels, heatmap_data, ...
                'Colormap', flipud(parula), ...
                'ColorbarVisible','on');
    h.Title = 'Participation Factors of Critical Mode';
    h.XLabel = 'Inputs';
    h.YLabel = 'Frequency';

    % === Summary table of metrics ===
    fprintf('\nModal Summary by Frequency:\n');
    fprintf('----------------------------------------------------------------------------------------------------------------------\n');
    fprintf('| f (Hz/p.u.) | Critical Mode | Eigenvalue        | |Z_m|     | PF1     | PF2     | PF3     | Dominant System |\n');
    fprintf('----------------------------------------------------------------------------------------------------------------------\n');
    for i = 1:N
        idx = critical_mode_idx(i);
        lambda = eigenvalues(idx, i);
        Zm = modal_impedance(idx, i);
        pf = participation_factors(idx, :, i);
        fprintf('| %9.3f |      %d       | %8.4f%+8.4fj | %8.3f | %7.4f | %7.4f | %7.4f | %-15s |\n', ...
            frequencies(i), idx, real(lambda), imag(lambda), Zm, pf(1), pf(2), pf(3), dominant_subsystem(i));
    end
    fprintf('----------------------------------------------------------------------------------------------------------------------\n');
end

