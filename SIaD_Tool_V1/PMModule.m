%% Technical University of Catalonia (UPC)
%% Higher Technical School of Industrial Engineering of Barcelona (ETSEIB)
%% Centre of Technological Innovation in Static Converters and Drives (CITCEA)
%% Doctoral Program in Electrical Engineering
%% Developed by: Luis Angel Garcia Reyes, MSc
%% Phase Margin Analysis

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
disp('SIaD Tool is applying the Phase Margin Analysis.')
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

Phase_Margins_Table=Phase_Margins_ABC(fd0, Y_sys1, Y_sys2, 'abc');
    case 2
%% dq0 sequence

 Phase_Margins_Table=Phase_Margins(fd0, Y_sys1, Y_sys2, 'dq');
    case 3
%% 0pn sequence

Phase_Margins_Table=Phase_Margins_ABC(fd0, Y_sys1, Y_sys2, 'abc');
end
    
    case 2
%% Current perturbation
% Selection of reference frame:
switch scanner_selector
    case 1
%% ABC sequence

Phase_Margins_Table=Phase_Margins_ABC(fd0, Y_sys1, Y_sys2, 'abc');
    case 2
%% dq0 sequence

 Phase_Margins_Table=Phase_Margins(fd0, Y_sys1, Y_sys2, 'dq');
    case 3
%% 0pn sequence

Phase_Margins_Table=Phase_Margins_ABC(fd0, Y_sys1, Y_sys2, 'abc');
end
end

function Phase_Margins_Table = Phase_Margins_ABC(f_c, Y_c_abc, Y_g_abc, type)
    % PHASE_MARGINS_ABC Computes crossing frequency and phase margin (PM) for ABC systems.
    %
    % INPUT:
    %   f_c      - Frequency vector (Hz).
    %   Y_c_abc  - 3x3xN array of converter admittance in ABC domain.
    %   Y_g_abc  - 3x3xN array of grid admittance in ABC domain.
    %   type     - String, either 'abc' or 'pn' (for transformed cases).
    %
    % OUTPUT:
    %   Phase_Margins_Table - Table containing crossing frequency and phase margin for selected cases.

    % Compute crossing frequency and phase margin for each matrix element
    [crossing_frequency_aa, PM_aa] = PMs(f_c, squeeze(Y_c_abc(1,1,:)), squeeze(Y_g_abc(1,1,:)));
    [crossing_frequency_ab, PM_ab] = PMs(f_c, squeeze(Y_c_abc(1,2,:)), squeeze(Y_g_abc(1,2,:)));
    [crossing_frequency_ac, PM_ac] = PMs(f_c, squeeze(Y_c_abc(1,3,:)), squeeze(Y_g_abc(1,3,:)));

    [crossing_frequency_ba, PM_ba] = PMs(f_c, squeeze(Y_c_abc(2,1,:)), squeeze(Y_g_abc(2,1,:)));
    [crossing_frequency_bb, PM_bb] = PMs(f_c, squeeze(Y_c_abc(2,2,:)), squeeze(Y_g_abc(2,2,:)));
    [crossing_frequency_bc, PM_bc] = PMs(f_c, squeeze(Y_c_abc(2,3,:)), squeeze(Y_g_abc(2,3,:)));

    [crossing_frequency_ca, PM_ca] = PMs(f_c, squeeze(Y_c_abc(3,1,:)), squeeze(Y_g_abc(3,1,:)));
    [crossing_frequency_cb, PM_cb] = PMs(f_c, squeeze(Y_c_abc(3,2,:)), squeeze(Y_g_abc(3,2,:)));
    [crossing_frequency_cc, PM_cc] = PMs(f_c, squeeze(Y_c_abc(3,3,:)), squeeze(Y_g_abc(3,3,:)));

    % Define case labels
    switch type
        case "abc"
            cases = {'aa'; 'ab'; 'ac'; 'ba'; 'bb'; 'bc'; 'ca'; 'cb'; 'cc'};
        case "pn"
            cases = {'pp'; 'pn'; 'np'; 'nn'; 'zz'; 'pz'; 'nz'; 'zp'; 'zn'}; % Example for transformed domains
        otherwise
            error("Type error: enter 'abc' or 'pn'");
    end

    % Collect results
    crossing_frequencies = [crossing_frequency_aa; crossing_frequency_ab; crossing_frequency_ac;
                            crossing_frequency_ba; crossing_frequency_bb; crossing_frequency_bc;
                            crossing_frequency_ca; crossing_frequency_cb; crossing_frequency_cc];

    PM_values = [PM_aa; PM_ab; PM_ac;
                 PM_ba; PM_bb; PM_bc;
                 PM_ca; PM_cb; PM_cc];

    % Create result table
    Phase_Margins_Table = table(cases, crossing_frequencies, PM_values, ...
        'VariableNames', {'Case', 'Crossing_Frequency_Hz', 'PM_deg'});

    % Display table
    disp(Phase_Margins_Table);
end


function Phase_Margins_Table = Phase_Margins(f_c, Y_c_dq, Y_g_dq, type)
    % PHASE_MARGINS Computes crossing frequency and phase margin (PM) for given dq or pn cases.
    %
    % INPUT:
    %   f_c    - Frequency vector (Hz).
    %   Y_c_dq - 2x2xN array of converter admittance in dq domain.
    %   Y_g_dq - 2x2xN array of grid admittance in dq domain.
    %   type   - String, either 'dq' (for dq cases) or 'pn' (for pn cases).
    %
    % OUTPUT:
    %   Phase_Margins_Table - Table containing crossing frequency and phase margin for selected cases.
    %
    % AUTHOR:
    %   Nicolae Darii (DTU and Siemens Gamesa)
    %
    % LICENSED UNDER:
    %   MIT License (see LICENSE file in the repository).
    %

    % Calcolo della crossing frequency e del phase margin per ogni caso
    [crossing_frequency_dd, PM_dd] = PMs(f_c, squeeze(Y_c_dq(1,1,:)), squeeze(Y_g_dq(1,1,:)));
    [crossing_frequency_dq, PM_dq] = PMs(f_c, squeeze(Y_c_dq(1,2,:)), squeeze(Y_g_dq(1,2,:)));
    [crossing_frequency_qd, PM_qd] = PMs(f_c, squeeze(Y_c_dq(2,1,:)), squeeze(Y_g_dq(2,1,:)));
    [crossing_frequency_qq, PM_qq] = PMs(f_c, squeeze(Y_c_dq(2,2,:)), squeeze(Y_g_dq(2,2,:)));

    % Creazione della tabella in base al tipo di analisi ('dq' o 'pn')
    switch type
        case "dq"
            cases = {'dd'; 'dq'; 'qd'; 'qq'};
        case "pn"
            cases = {'pp'; 'pn'; 'np'; 'nn'};
        otherwise
            error("Type error: enter 'dq' or 'pn'");
    end

    % Assegnazione dei valori calcolati alla tabella
    crossing_frequencies = [crossing_frequency_dd; crossing_frequency_dq; crossing_frequency_qd; crossing_frequency_qq];
    PM_values = [PM_dd; PM_dq; PM_qd; PM_qq];

    % Creazione della tabella dei risultati
    Phase_Margins_Table = table(cases, crossing_frequencies, PM_values, ...
        'VariableNames', {'Case', 'Crossing_Frequency [Hz]', 'PM [deg]'});

    % Visualizzazione della tabella
    disp(Phase_Margins_Table);
end

function [crossing_frequency, PM] = PMs(f, Y_c, Y_g)
    % PMS Computes the crossing frequency and phase margin (PM) between converter and grid admittances.
    %
    % INPUT:
    %   f   - Frequency vector (Hz).
    %   Y_c - Converter admittance (complex values).
    %   Y_g - Grid admittance (complex values).
    %
    % OUTPUT:
    %   crossing_frequency - Frequency at which |Y_c / Y_g| = 1.
    %   PM                 - Phase margin (deg), difference in phase at crossing frequency.
    %
    % AUTHOR:
    %   Nicolae Darii (DTU and Siemens Gamesa)
    %
    % LICENSED UNDER:
    %   MIT License (see LICENSE file in the repository).
    %

    % Estrazione dei valori di ammettenza e calcolo del rapporto
    vec_c = squeeze(Y_c);
    vec_g = squeeze(Y_g);
    abs_ratio = abs(vec_c ./ vec_g);

    % Trova il punto di crossing: valore minimo di |1 - |Y_c / Y_g||    
    [min_value, min_index] = min(abs(1 - abs_ratio));
    crossing_frequency = f(min_index);  % Frequenza di crossing

    %% Calcolo del phase margin
    Phase_c = rad2deg(angle(vec_c(min_index))); % Fase dell'ammettenza del convertitore
    Phase_g = rad2deg(angle(vec_g(min_index))); % Fase dell'ammettenza della griglia

    % Calcolo del Phase Margin (PM)
    PM = 180 + (Phase_c - Phase_g);
end
