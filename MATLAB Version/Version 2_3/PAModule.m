%% Technical University of Catalonia (UPC)
%% Higher Technical School of Industrial Engineering of Barcelona (ETSEIB)
%% Centre of Technological Innovation in Static Converters and Drives (CITCEA)
%% Doctoral Program in Electrical Engineering
%% Developed by: Luis Angel Garcia Reyes, MSc
%% Passivity analysis tool

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

disp('---')
disp('SIaD Tool is applying the Passivity Assessment.')
disp('The results are the next:')
disp('---')

% Cálculo de valores propios mínimos
lambda_min_scan1 = zeros(1, length(fd0));
lambda_min_scan2 = zeros(1, length(fd0));

for fn = 1:length(fd0)
    Pmatrix_scan1 = Y_sys1(:,:,fn) + pagectranspose(Y_sys1(:,:,fn));
    Pmatrix_scan2 = Y_sys2(:,:,fn) + pagectranspose(Y_sys2(:,:,fn));
    lambda_min_scan1(fn) = min(eig(Pmatrix_scan1));
    lambda_min_scan2(fn) = min(eig(Pmatrix_scan2));
end

% === Tabla de violaciones de pasividad ===
fprintf('\nPassivity Violation Summary:\n');
fprintf('---------------------------------------------\n');
fprintf('| System   | Frequency Range [Hz] | Status   |\n');
fprintf('---------------------------------------------\n');

detect_passivity_ranges(fd0, lambda_min_scan1, 'System 1');
detect_passivity_ranges(fd0, lambda_min_scan2, 'System 2');

fprintf('---------------------------------------------\n');

% === Gráfica con formato mejorado ===
figure;
hold on;
grid on;
xlabel('Frequency (Hz)', 'FontSize', 12);
ylabel('Minimum Eigenvalue (\lambda_{min})', 'FontSize', 12);
title('Passivity Evaluation via Minimum Eigenvalues', 'FontSize', 14);

plot(fd0, lambda_min_scan1, '-', 'DisplayName', 'System 1', 'Color', 'b', 'LineWidth', 3);
plot(fd0, lambda_min_scan2, '-', 'DisplayName', 'System 2', 'Color', 'g', 'LineWidth', 3);
plot(fd0, zeros(1, length(fd0)), '--', 'Color', 'r', 'LineWidth', 2, 'DisplayName', 'Passive Limit');

legend('FontSize', 10);
hold off;

% === Función auxiliar para detectar rangos no pasivos ===
function detect_passivity_ranges(fd, eig_vals, system_name)
    is_non_passive = eig_vals < 0;
    start_idx = [];
    for i = 1:length(fd)
        if is_non_passive(i)
            if isempty(start_idx)
                start_idx = i;
            end
        elseif ~isempty(start_idx)
            fprintf('| %-8s | %5.1f to %5.1f     | Non-Passive |\n', system_name, fd(start_idx), fd(i-1));
            start_idx = [];
        end
    end
    if ~isempty(start_idx)
        fprintf('| %-8s | %5.1f to %5.1f     | Non-Passive |\n', system_name, fd(start_idx), fd(end));
    end
end
