%% Technical University of Catalonia (UPC)
%% Higher Technical School of Industrial Engineering of Barcelona (ETSEIB)
%% Centre of Technological Innovation in Static Converters and Drives (CITCEA)
%% Doctoral Program in Electrical Engineering
%% Developed by: Luis Angel Garcia Reyes, MSc
%% Multi-sequence Generalized Nyquist Criterion (GNC) Tool 

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


%% Stability by Generalized Nyquist Criterion (GNC)

disp('---')
disp('SIaD Tool is applying the GNC for stability assessment...')
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
L = zeros(3,3,length(fd0)); % Initial zero L open loop matriz
E = zeros(3,length(fd0)); % Initial zero E open loop eigenvalues
Y_sys1=zeros(3,3,length(fd0));
Y_sys2=zeros(3,3,length(fd0));
Z_sys1=zeros(3,3,length(fd0));
Z_sys2=zeros(3,3,length(fd0));
for n=1:length(fd0)
    % System 1
    Y_sys1(:,:,n) = Y_abc1(:,:,n);
    Z_sys1(:,:,n) = inv(Y_sys1(:,:,n));
    % System 2
    Y_sys2(:,:,n) = Y_abc2(:,:,n);
    Z_sys2(:,:,n) = inv(Y_sys2(:,:,n));
    % Open loop gain
    L(:,:,n) = Y_sys1(:,:,n)*Z_sys2(:,:,n);
    % Eigenvalues computation
    E(:,n) = eig(L(:,:,n));
end
MIMO_Nyquist3(fd0, E, 1);
    case 2
%% dq0 sequence

% Set initial data and matrices:
% Frequency vector of measurements: fd0
% System 1: Yqq1, Yqd1, Ydq1 and Ydd1
% System 2: Yqq2, Yqd2, Ydq2 and Ydd2
L = zeros(2,2,length(fd0)); % Initial zero L open loop matriz
E = zeros(2,length(fd0)); % Initial zero E open loop eigenvalues
for n=1:length(fd0)
    % System 1
    Y_sys1(:,:,n) = [Yqq1(n), Yqd1(n); 
                     Ydq1(n), Ydd1(n)];
    Z_sys1(:,:,n) = inv(Y_sys1(:,:,n));
    % System 2
    Y_sys2(:,:,n) = [Yqq2(n), Yqd2(n); 
                     Ydq2(n), Ydd2(n)];
    Z_sys2(:,:,n) = inv(Y_sys2(:,:,n));
    % Open loop gain
    L(:,:,n) = Y_sys1(:,:,n)*Z_sys2(:,:,n);
    % Eigenvalues computation
    E(:,n) = eig(L(:,:,n));
end
% MIMO Nyquist in dq0 reference
MIMO_Nyquist2(fd0, E, 1); % Put 0 1 to show -1 unstable point or 0 to hide it
    case 3
%% 0pn sequence
% Frequency vector of measurements: fd0
% System 1: Y_0pn1
% System 2: Y_0pn2
L = zeros(3,3,length(fd0)); % Initial zero L open loop matriz
E = zeros(3,length(fd0)); % Initial zero E open loop eigenvalues
Y_sys1=zeros(3,3,length(fd0));
Y_sys2=zeros(3,3,length(fd0));
Z_sys1=zeros(3,3,length(fd0));
Z_sys2=zeros(3,3,length(fd0));
for n=1:length(fd0)
    % System 1
    Y_sys1(:,:,n) = Y_0pn1(:,:,n);
    Z_sys1(:,:,n) = inv(Y_sys1(:,:,n));
    % System 2
    Y_sys2(:,:,n) = Y_0pn2(:,:,n);
    Z_sys2(:,:,n) = inv(Y_sys2(:,:,n));
    % Open loop gain
    L(:,:,n) = Y_sys1(:,:,n)*Z_sys2(:,:,n);
    % Eigenvalues computation
    E(:,n) = eig(L(:,:,n));
end
MIMO_Nyquist3(fd0, E, 1);
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
L = zeros(3,3,length(fd0)); % Initial zero L open loop matriz
E = zeros(3,length(fd0)); % Initial zero E open loop eigenvalues
Y_sys1=zeros(3,3,length(fd0));
Y_sys2=zeros(3,3,length(fd0));
Z_sys1=zeros(3,3,length(fd0));
Z_sys2=zeros(3,3,length(fd0));
for n=1:length(fd0)
    % System 1
    Z_sys1(:,:,n) = Z_abc1(:,:,n);
    Y_sys1(:,:,n) = inv(Z_sys1(:,:,n));
    % System 2
    Z_sys2(:,:,n) = Z_abc2(:,:,n);
    Y_sys2(:,:,n) = inv(Z_sys2(:,:,n));
    % Open loop gain
    L(:,:,n) = Y_sys1(:,:,n)*Z_sys2(:,:,n);
    % Eigenvalues computation
    E(:,n) = eig(L(:,:,n));
end
MIMO_Nyquist3(fd0, E, 1);
    case 2
%% dq0 sequence

% Set initial data and matrices:
% Frequency vector of measurements: fd0
% System 1: Zqq1, Zqd1, Zdq1 and Zdd1
% System 2: Zqq2, Zqd2, Zdq2 and Zdd2
L = zeros(2,2,length(fd0)); % Initial zero L open loop matriz
E = zeros(2,length(fd0)); % Initial zero E open loop eigenvalues
for n=1:length(fd0)
    % System 1
    Z_sys1(:,:,n) = [Zqq1(n), Zqd1(n); 
                     Zdq1(n), Zdd1(n)];
    Y_sys1(:,:,n) = inv(Z_sys1(:,:,n));
    % System 2
    Z_sys2(:,:,n) = [Zqq2(n), Zqd2(n); 
                     Zdq2(n), Zdd2(n)];
    Y_sys2(:,:,n) = inv(Z_sys2(:,:,n));
    % Open loop gain
    L(:,:,n) = Y_sys1(:,:,n)*Z_sys2(:,:,n);
    % Eigenvalues computation
    E(:,n) = eig(L(:,:,n));
end
% MIMO Nyquist in dq0 reference
MIMO_Nyquist2(fd0, E, 1); % Put 0 1 to show -1 unstable point or 0 to hide it
    case 3
%% 0pn sequence
% Frequency vector of measurements: fd0
% System 1: Z_0pn1
% System 2: Z_0pn2
L = zeros(3,3,length(fd0)); % Initial zero L open loop matriz
E = zeros(3,length(fd0)); % Initial zero E open loop eigenvalues
Y_sys1=zeros(3,3,length(fd0));
Y_sys2=zeros(3,3,length(fd0));
Z_sys1=zeros(3,3,length(fd0));
Z_sys2=zeros(3,3,length(fd0));
for n=1:length(fd0)
    % System 1
    Z_sys1(:,:,n) = Z_0pn1(:,:,n);
    Y_sys1(:,:,n) = inv(Z_sys1(:,:,n));
    % System 2
    Z_sys2(:,:,n) = Z_0pn2(:,:,n);
    Y_sys2(:,:,n) = inv(Z_sys2(:,:,n));
    % Open loop gain
    L(:,:,n) = Y_sys1(:,:,n)*Z_sys2(:,:,n);
    % Eigenvalues computation
    E(:,n) = eig(L(:,:,n));
end
MIMO_Nyquist3(fd0, E, 1);
end
end
disp('Graphical results from GNC are been displayed.')


function MIMO_Nyquist3(frequencies, eigenvalues, critical_visible)
    % PLOT_ORDERED_EIGENVALUES Plots ordered eigenvalues in the complex plane.
    % INPUT:
    %   frequencies       - Vector of frequencies in Hz (1xN).
    %   eigenvalues       - Matrix (3xN) containing the calculated eigenvalues.
    %   critical_visible  - 0/1 value to toggle the visibility of the (-1, 0) point.

    N = length(frequencies);
    ordered_eigenvalues = exchange3(eigenvalues);

    figure;
    hold on;
    grid on;
    xlabel('Re(\lambda)');
    ylabel('Im(\lambda)');

    colors = ["r", "b", "k"];
    mirror_colors = ["m", "c", "g"];
    
    for i = 1:3
        plot(real(ordered_eigenvalues(i, :)), imag(ordered_eigenvalues(i, :)), '-', ...
            'DisplayName', sprintf('Eigenvalue %d [0,+inf]', i), ...
            'Color', colors(i), 'LineWidth', 3);
        
        plot(real(ordered_eigenvalues(i, :)), -imag(ordered_eigenvalues(i, :)), '-', ...
            'DisplayName', sprintf('Eigenvalue %d [-inf,0]', i), ...
            'Color', mirror_colors(i), 'LineWidth', 3);
    end

    if critical_visible == 1
        plot(-1, 0, 'pentagram', 'MarkerSize', 8, "LineWidth", 3); 
    elseif critical_visible == 0
        disp("The (-1,0) point is hidden")
    else
        disp("Error in the input: 0- to hide (-1,0) point, 1- to show it")
    end

    legend show;
    hold off;
end

function E = exchange3(E)
    % EXCHANGE Reorders eigenvalues to maintain trajectory continuity.
    % INPUT:
    %   E - Matrix (3xN) of eigenvalues.
    
    [rows, cols] = size(E);
    
    for j = 2:cols
        prev = E(:, j-1);   % Eigenvalues at previous frequency
        curr = E(:, j);     % Eigenvalues at current frequency
        
        % Compute distance matrix: D(i,j) = |curr(i) - prev(j)|
        D = abs(curr - prev');  % 3x3 matrix
        
        % Greedy assignment: find closest match for each current eigenvalue
        assigned = zeros(rows, 1);
        used = false(rows, 1);
        
        for i = 1:rows
            [~, idx] = min(D(i, ~used));
            real_idx = find(~used);
            assigned(i) = real_idx(idx);
            used(real_idx(idx)) = true;
        end
        
        % Reorder current column based on assignment
        E(:, j) = curr(assigned);
    end
end



function MIMO_Nyquist2(frequencies, eigenvalues, critical_visible)
    % PLOT_ORDERED_EIGENVALUES Plots ordered eigenvalues in the complex plane.
    %
    % INPUT:
    %   frequencies       - Vector of frequencies in Hz (1xN).
    %   eigenvalues       - Matrix (2xN) containing the calculated eigenvalues.
    %   critical_visible  - 0/1 value to toggle the visibility of the (-1, 0) point.
    %
    % AUTHOR:
    %   Nicolae Darii (DTU and Siemens Gamesa)
    %
    % LICENSED UNDER:
    %   MIT License (see LICENSE file in the repository).
    %

    % Number of frequencies
    N = length(frequencies);

    %% Ordered case
    ordered_eigenvalues = exchange(eigenvalues);

    % Plot ordered eigenvalues in the complex plane
    figure;
    hold on;
    grid on;
    xlabel('Re(\lambda)');
    ylim([-250 250])
    xlim([-300 150])
    ylabel('Im(\lambda)');
    % title('Eigenvalues of MIMO open loop L(S)');

    % Plot eigenvalue trajectories
    h1 = plot(real(ordered_eigenvalues(1, :)), imag(ordered_eigenvalues(1, :)), '-', 'DisplayName', '\lambda_{1} [0,+\infty]', 'Color',"r","LineWidth",3);
    h2 = plot(real(ordered_eigenvalues(1, :)), -imag(ordered_eigenvalues(1, :)), '-','DisplayName', '\lambda_{1} [-\infty,0]','Color',"m","LineWidth",3);
    
    h3 = plot(real(ordered_eigenvalues(2, :)), imag(ordered_eigenvalues(2, :)), '-', 'DisplayName', '\lambda_{2} [0,+\infty]','Color',"b","LineWidth",3);
    h4 = plot(real(ordered_eigenvalues(2, :)), -imag(ordered_eigenvalues(2, :)), '-','DisplayName', '\lambda_{2} [-\infty,0]','Color',"k","LineWidth",3);
    
    if critical_visible == 1
        h5 = plot(-1, 0, 'pentagram', 'MarkerSize', 5, 'Color', "g","LineWidth",3); 
    elseif critical_visible == 0
        disp("The (-1,0) point is hidden")
    else
        disp("Error in the input: 0- to hide (-1,0) point, 1- to show it")
    end
    
       % --- Cuadro de zoom ---
    % Define los límites del zoom (modifica estos valores)
    X1 = -3; X2 = 1;
    Y1 = -0.5; Y2 = 0.5;

%     % Dibuja el rectángulo de zoom
%     rectangle('Position', [X1, Y1, X2 - X1, Y2 - Y1], ...
%               'EdgeColor', 'k', ...
%               'LineWidth', 1.5, ...
%               'LineStyle', '--');
    % --- Inset: lupa con zoom ---
    % Crea axes pequeños (lupa)
    axInset = axes('Position', [0.65, 0.65, 0.25, 0.25]); % [x y width height] en fracción de figura
    box on;
    hold on;
    grid on;

    % Traza los mismos datos en la lupa
    plot(real(ordered_eigenvalues(1, :)), imag(ordered_eigenvalues(1, :)), '-', 'Color',"r","LineWidth",1.5);
    plot(real(ordered_eigenvalues(1, :)), -imag(ordered_eigenvalues(1, :)), '-', 'Color',"m","LineWidth",1.5);
    plot(real(ordered_eigenvalues(2, :)), imag(ordered_eigenvalues(2, :)), '-', 'Color',"b","LineWidth",1.5);
    plot(real(ordered_eigenvalues(2, :)), -imag(ordered_eigenvalues(2, :)), '-', 'Color',"k","LineWidth",1.5);

    % Punto crítico si aplica
    if critical_visible == 1
        plot(-1, 0, 'p', 'MarkerSize', 5, 'Color', "g","LineWidth",1.5); 
    end

    % Ajusta límites del inset (zoom)
    xlim([X1, X2]);
    ylim([Y1, Y2]);
%     title('Zoom');

    legend([h2, h1, h4, h3]);
    hold off;
end

function E = exchange(E)
    % EXCHANGE Swaps components if necessary based on distance.
    %
    % INPUT:
    %   E   - Matrix containing the components to be swapped (2xN).
    %
    % OUTPUT:
    %   E   - Modified matrix after performing the swap operation.
    %
    % AUTHOR:
    %   Nicolae Darii (DTU and Siemens Gamesa)
    %
    % LICENSED UNDER:
    %   MIT License (see LICENSE file in the repository).
    %

    [rows, cols] = size(E);
    
    % Iterate over the columns of the matrix
    for j = 2:cols  % starts from the second column
        dist1 = distance(E(1,j), E(1,j-1));  % distance between E(1,j) and E(1,j-1)
        dist2 = distance(E(1,j), E(2,j-1));  % distance between E(1,j) and E(2,j-1)
        
        % If E(1,j) is farther from E(1,j-1) than from E(2,j-1), swap them
        if dist1 > dist2
            % Swap the values
            temp = E(1,j);
            E(1,j) = E(2,j);
            E(2,j) = temp;
        end
    end
end

function d = distance(z1, z2)
    % DISTANCE Calculates the distance between two complex numbers.
    %
    % INPUT:
    %   z1    - First complex number.
    %   z2    - Second complex number.
    %
    % OUTPUT:
    %   d     - Distance between the two complex numbers.
    %
    % AUTHOR:
    %   Nicolae Darii (DTU and Siemens Gamesa)
    %
    % LICENSED UNDER:
    %   MIT License (see LICENSE file in the repository).
    %

    d = abs(z1 - z2);
end