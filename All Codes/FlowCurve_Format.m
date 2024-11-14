clear; clc; clf; close all;

%% Test details -- EDIT DEPENDING ON YOUR SETTINGS %%

duration = 80; % Duration of every steady shear
spp = 1; % seconds/pt
p = duration/spp; % Calculated number of points
start = 5; % Default Don't change (data starts at row #5 in anton paar table)

%% Input Data file

for k = 1
    [file_t, path_t] = uigetfile('*.csv', 'Select AntonPaar Raw File');
    filename = fullfile(path_t, file_t);
    A = xlsread(filename);
    time_all = A(:,1);
    shear_rate_all = A(:,2);
    stress_all = A(:,3);
    viscosity_all = A(:,4);
end

%% File Rearrangement

j = 1;
i = 1;
while true
    idx_start = start+(i-1)*(p+1)+(start+1)*(i-1);
    idx_end = idx_start + (p-1);

    % Check if there's still data available in the current range
    if idx_start > length(time_all) || isempty(time_all(idx_start:idx_end))
        break; % Exit the loop if data is missing or beyond range
    end

    B(:,j) = time_all(idx_start:idx_end);
    B(:,j+1) = shear_rate_all(idx_start:idx_end);
    B(:,j+2) = stress_all(idx_start:idx_end);
    B(:,j+3) = viscosity_all(idx_start:idx_end);
    j = j+4;
    i = i+1;
end

%% Print

file_tx = erase(file_t, '.csv');
T_file = sprintf('Processed_%s',file_tx);
T_file_xlsx = [T_file, '.xlsx'];
T_name = fullfile(path_t, T_file_xlsx);
writematrix(B, T_name)
file_tx
disp('done')
