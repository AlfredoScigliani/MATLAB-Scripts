clear; clc; clf; close all
format shortg
warning off

%% Experiment Details

[file, path] = uigetfile('*.csv', 'Select totxtData File');
filename = fullfile(path, file);
data =  xlsread(filename);
file

chart_title = 'Input';
dims = [1 35];

prompt = "# of Rows";
numRowss = inputdlg(prompt,chart_title,dims);
numRows = str2double(numRowss)
if isempty(numRows)
    fprintf('Code terminated: Missing Input.')
    return
end

prompt = "# of Columns";
numColss = inputdlg(prompt,chart_title,dims);
numCols = str2double(numColss)
if isempty(numCols)
    fprintf('Code terminated: Missing Input.')
    return
end

prompt = "Left ppm";
Lppms = inputdlg(prompt,chart_title,dims);
Lppm = str2double(Lppms)
if isempty(Lppm)
    fprintf('Code terminated: Missing Input.')
    return
end

prompt = "Right ppm";
Rppms = inputdlg(prompt,chart_title,dims);
Rppm = str2double(Rppms)
if isempty(Rppm)
    fprintf('Code terminated: Missing Input.')
    return
end
x = linspace(Lppm, Rppm, numCols)';
P = data(:,1);
j = 1;
S = [];
for i = 1:numRows
    S{i}(:,1) = x;
    S{i}(:,2) = P(j:j+numCols-1); j = j+numCols+1;
end

prompt = "Data sets to plot";
dt = inputdlg(prompt,chart_title,dims);
d = str2double(dt)
if isempty(d)
    fprintf('Code terminated: Missing Input.')
    return
end

for i = 1:d
    prompt = sprintf('Slice No. of data set #%.0f (1 to %.0f) ', i, numRows);
    idxs = inputdlg(prompt,chart_title,dims);
    idx(i) = str2double(idxs)
    if isempty(idx(i))
        fprintf('Code terminated: Missing Input.')
        return
    end
end

for i = 1:d
    prompt = sprintf('B value for Slice #%.0f ', idx(i));
    Bs = inputdlg(prompt,chart_title,dims);
    B(i) = str2double(Bs)
    if isempty(idx(i))
        fprintf('Code terminated: Missing Input.')
        return
    end
end

prompt = "Number of Peaks (Vertical lines)";
numPeakss = inputdlg(prompt,chart_title,dims);
numPeaks = str2double(numPeakss)
if isempty(numPeaks)
    fprintf('Code terminated: Missing Input.')
    return
end

for i = 1:numPeaks
    prompt = sprintf('Chemical Shift #%.0f [ppm]', i);
    ppms = inputdlg(prompt,chart_title,dims);
    ppm(i,1) = str2double(ppms);
    if isempty(ppm(i,1))
        fprintf('Code terminated: Missing Input.')
        return
    end
end
ppm


%% Plotting ALL

f1 = figure;
for i = 1:numRows
    plot(S{i}(:,1),S{i}(:,2), '-')
    hold on
end
ylabel('Intensity [a.u.]')
xlabel('ppm')
xlim([0 10])
ylim([0 1.1*max(S{1}(:,2))])
set(gca, 'XDir','reverse')
title('All B Values')

[~, name, ~] = fileparts(file);
file_name = sprintf('Spectra_All_B_%s.fig', name);
fullName = fullfile(path, file_name);
saveas(f1, fullName);

%%  Plotting Selected

n=12;
colors = {'#0072BD','#D95319','#EDB120','#7E2F8E','#77AC30','#4DBEEE','#A2142F'}; %https://linkprotect.cudasvc.com/url?a=https%3a%2f%2fhtmlcolorcodes.com%2f&c=E,1,9tQNuXINuJsIRcPyGOJOZSUq4w_OCA6uMN4Q5_wxQ3QK3lHKETYaPH4_xQtCLLaFCLj0ASau4rV19hAm0ykD6HGOL43dlH_y55N5HJoSQA,,&typo=1
f2 = figure;
a = 1;
plot_handles = [];
plot_handles2 = [];
legend_labels = cell(1, d);
Bflip = flip(B);
idxflip = flip(idx);
shift = 0;

for i = d:-1:1
    h = plot(S{idx(i)}(:,1), S{idx(i)}(:,2) + shift*a, '-', 'Color', colors{a},'LineWidth',.8);
    plot_handles = [plot_handles; h];
    hold on
    a = a+1;
end
for i = 1:d
    legend_labels{i} = sprintf('B = %.0f', Bflip(i));
end

a = 1;
for i = d:-1:1
    h2 = plot(S{idx(i)}(:,1), S{idx(i)}(:,2) + shift*a, '-', 'Color', colors{a});
    plot_handles2 = [plot_handles2; h2];
    hold on
    a = a+1;
end

for i = 1:d
    legend_labels2{i} = sprintf('Slice %.0f', idxflip(i));
end

if shift == 0
    yline(0,'r--');
end
for i = 1:numPeaks
    xline(ppm(i), 'k--');
end
set(gca, 'XDir', 'reverse');

xlim([-1 10])
ylim([-0.5E7 1.1*max(S{1}(:,2))])

ylabel('Intensity [a.u.]', 'FontSize',14)
xlabel('ppm','fontsize',12)
title('Selected B Values', 'fontsize',14)

for k = 1 %legend
    plot_handles = flip(plot_handles);
    legend_labels = flip(legend_labels);
    legend1 = legend(plot_handles, legend_labels, 'Location', 'NorthWest','FontSize',n);

    plot_handles2 = flip(plot_handles2);
    legend_labels2 = flip(legend_labels2);

    ax2 = axes('Position', get(gca, 'Position'), 'Color', 'none', 'XTick', [], 'YTick', [], 'Box', 'off');
    set(ax2, 'color', 'none', 'XColor', 'none', 'YColor', 'none');
    legend2 = legend(ax2, plot_handles2, legend_labels2, 'Location', 'NorthEast','fontsize',n);
    set(legend2, 'Color', get(legend1, 'Color'));
end

file_name = sprintf('Spectra_Selected_B_%s.fig', name);
fullName = fullfile(path, file_name);
saveas(f2, fullName);

%% Noise stdv

std_S_noise = [];
SliceNo =[];
slice_start = 7;
for i = slice_start:numRows
    std_S_noise(i-slice_start+1,1) = std(S{i}(300:1000, 2));

end
SliceNo(:,1) = slice_start:1:numRows;
T = table(SliceNo,std_S_noise)
stdv = sqrt(sum(std_S_noise.^2)/length(std_S_noise));

%% Print

st = S{idx(1)};
for i = 2:d
    st = [st S{idx(i)}];
end

T_file = sprintf('Spectra_Selected_B_%s', name);
T_file_xlsx = [T_file, '.xlsx'];
T_name = fullfile(path, T_file_xlsx);
writematrix(st, T_name)

s = S{1};
for i = 2:numRows
    s = [s S{i}];
end

T_file = sprintf('Spectra_All_B_%s', name);
T_file_xlsx = [T_file, '.xlsx'];
T_name = fullfile(path, T_file_xlsx);
writematrix(s, T_name)

%% Save Workspace

Results_for = name
fileName = sprintf('work_Spectra_%s.mat', name);
fullFileName = fullfile(path, fileName);
save(fullFileName);
