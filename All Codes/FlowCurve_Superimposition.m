clear; clc; clf; close all;

%% Initialize cell arrays to store data

chart_title = 'Input';
dims = [1 35];
prompt = "No. of Datasets to Plot";
datasetss = inputdlg(prompt,chart_title,dims);
datasets = str2double(datasetss);

time = cell(datasets, 1);
shearrate = cell(datasets, 1);
stress = cell(datasets, 1);
viscosity = cell(datasets, 1);
v = cell(datasets, 1);
v_std = cell(datasets, 1);
sr = cell(datasets, 1);
sigma = cell(datasets, 1);
sigma_std = cell(datasets, 1);
legendEntries = cell(datasets, 1);

%% Input Data Files

for d = 1:datasets
    [file_t, path_t] = uigetfile('*FlowCurve*', sprintf('Select FlowCurve File %d', d));
    file_t
    filename = fullfile(path_t, file_t);
    A = xlsread(filename);
    [numRows, numCols] = size(A);
    j = 1;
    for i = 1:4:numCols
        time{d}(:,j) = A(:,i);
        shearrate{d}(:,j) = A(:,i+1);
        stress{d}(:,j) = A(:,i+2);
        viscosity{d}(:,j) = A(:,i+3);
        j = j + 1;
    end
    legendEntries{d} = sprintf('%s', file_t); % Collect legend entries
end

%% Input Relaxation Times

prompt = "Implement Wi? [1] = yes, [0] = No";
Wiss = inputdlg(prompt,chart_title,dims);
Wis = str2double(Wiss);

if Wis == 1
    for d = 1:datasets
        prompt = sprintf('Relaxation time #%d',d);
        Lt = inputdlg(prompt,chart_title,dims);
        L(d) = str2double(Lt)
    end
else
    for d = 1:datasets
        L(d) = 1;
    end
end

%% Averaging Data

for d = 1:datasets
    [numRows, numCols] = size(viscosity{d});
    for i = 1:numCols
        v{d}(i) = mean(viscosity{d}(end*0.85:end, i));
        v_std{d}(i) = std(viscosity{d}(end*0.85:end, i));
        sr{d}(i) = mean(shearrate{d}(end*0.85:end, i));
        sigma{d}(i) = mean(stress{d}(end*0.85:end, i));
        sigma_std{d}(i) = std(stress{d}(end*0.85:end, i));
    end
end

%% Plotting Data

colors = {'#0072BD', '#D95319', '#EDB120', '#7E2F8E', '#77AC30', '#4DBEEE', '#A2142F', '#0072BD'};
f1 = figure;
a_handles = gobjects(1, datasets); % Preallocate array to store handles for legend
for d = 1:datasets
    a_handles(d) = semilogx(sr{d}*L(d), sigma{d}, 'o-', 'Color', colors{d});
    hold on
    eb = errorbar(sr{d}*L(d), sigma{d}, sigma_std{d}, 'vertical', 'LineStyle', 'none');
    set(eb, 'color', 'k', 'LineWidth', 1);
end
ylabel('\sigma [Pa]');
if L(1) == 1
    xlabel('Shear Rate [1/s]');
else
    xlabel('Wi');
end
legend(a_handles, legendEntries, 'Location', 'SouthEast');
grid on;

xlim([0.1 100])
ylim([0.1 40])

fullName = 'FlowCurve_sample';
saveas(f1,fullName,'png');

%% Export Selected

chart_title = 'Input';
dims = [1 35];
prompt = "Export?   [1] = Yes    [0] = Stop";
mode_exports = inputdlg(prompt,chart_title,dims);
mode_export = str2double(mode_exports);

if mode_export == 0
    disp('Run Section to Export Results.')
    return
else

    % Initialize the title row with headers for each dataset
    titleData = {};
    for d = 1:datasets
        titleData{1, end+1} = sprintf('Shear Rate [1/s] (Dataset %d)', d);
        titleData{1, end+1} = sprintf('Stress [Pa] (Dataset %d)', d);
        titleData{1, end+1} = sprintf('Viscosity [Pa s] (Dataset %d)', d);
    end

    % Determine the maximum length among all vectors
    maxLength = 0;
    for d = 1:datasets
        maxLength = max([maxLength, length(sr{d}), length(sigma{d}), length(v{d})]);
    end

    % Initialize exportData with NaN padding to handle different lengths
    exportData = NaN(maxLength, datasets * 3);

    % Fill exportData with the values, padding with NaNs where necessary
    for d = 1:datasets
        exportData(1:length(sr{d}), (d-1)*3 + 1) = sr{d};
        exportData(1:length(sigma{d}), (d-1)*3 + 2) = sigma{d};
        exportData(1:length(v{d}), (d-1)*3 + 3) = v{d};
    end

    % Combine title row with data matrix
    dataMatrixWithTitles = [titleData; num2cell(exportData)];

    % Define the file path and name
    name = sprintf('%sCombined_FlowC_%s.xlsx', path_t, file_t);

    % Write the cell array with titles and data to Excel
    writecell(dataMatrixWithTitles, name);

    winopen(name);
    disp('Results Exported');
end
