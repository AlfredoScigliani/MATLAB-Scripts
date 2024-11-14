clear; clc; clf; close all;
format long

%% File Input

chart_title = 'Input';  %Title outside dialog box
dims = [1 35]; %Dimension of dialog box
[file_stress, path_1D] = uigetfile('*', 'Shear Stress File');
stress_file = fullfile(path_1D, file_stress);
shear_data = xlsread(stress_file);
file_stress

time(:,1) = shear_data(:,1);
shearrate(:,1) = shear_data(:,2);
viscosity(:,1) = shear_data(:,3);
stress(:,1) = shear_data(:,6);
shear_av = mean(shearrate(round(end*0.7):end));
strain = time*shear_av;

prompt = "Experiment Number";
nums = inputdlg(prompt,chart_title,dims);
num = str2double(nums);
if isempty(num)
    fprintf('Code terminated: Missing Input.')
    return
end

prompt = "No. of 2D Profiles (Videos)";
profiles2Ds = inputdlg(prompt,chart_title,dims);
profiles2D = str2double(profiles2Ds);
if isempty(profiles2D)
    fprintf('Code terminated: Missing Input.')
    return
end

for i = 1:profiles2D
    [file_range, path_range] = uigetfile('*range.csv', sprintf('2D Range File %.0f', i));
    range_file = fullfile(path_range, file_range);
    range_data2D = xlsread(range_file);
    range_data{i} = flipud(range_data2D);
    file_range
end

prompt = "Video Duration [s]";
vid_ts = inputdlg(prompt,chart_title,dims);
vid_t = str2double(vid_ts);
if isempty(vid_t)
    fprintf('Code terminated: Missing Input.')
    return
end

low_lim = [];
hi_lim = [];
for i = 1:profiles2D
    prompt = sprintf('Time for Video #%.0f [s]', i);
    time2D_tests = inputdlg(prompt,chart_title,dims);
    time2D_test(i) = str2double(time2D_tests);
    if isempty(time2D_test(i))
        fprintf('Code terminated: Missing Input.')
        return
    end
    low_lim(i) = time2D_test(i)*shear_av;
    hi_lim(i) = time2D_test(i)*shear_av + vid_t*shear_av;
end

[~, name, ~] = fileparts(file_stress);
fig_name = name;

%% Plot

close all
n = 16;
width = 1.6;
height = 1.7;
p = 256;
cmap = getPyPlot_cMap('seismic', p);
gr = 170;
nan_color = [gr/255, gr/255, gr/255]; % Dark gray color (#464646)
custom_cmap = [nan_color;cmap];

for i = 1:profiles2D
    f2 = figure('Units', 'inches', 'Position', [1, 4, width, height]);
    x_values = linspace(0, 1, size(range_data{i}, 2));
    y_values = linspace(low_lim(i), hi_lim(i), size(range_data{i}, 1));
    imagesc(x_values, y_values, imrotate(range_data{i}, -90));
    colormap(custom_cmap);
    clim([-1, 1]);
    set(gca, 'XAxisLocation', 'top');
    set(gca, 'FontSize', n)

    xlabh = get(gca, 'XLabel');
    ax = gca;
    ax.YAxis.Exponent = 0;
    set(gca, 'xtick', []);
    set(gca, 'xticklabel', []);

    if i>1
        yticks([round(y_values(round(length(y_values)*0.20))/10*10), round(y_values(round(length(y_values)*0.8))/10*10)]);
        set(gca, 'YTickLabelRotation', 90);
    else
        if round(y_values(round(length(y_values)*0.9))/10*10) < 50
            yticks([round(y_values(1)/10)*10, round(y_values(round(length(y_values)*0.9))/10*10)])
            set(gca, 'YTickLabelRotation', 90);
        else
            yticks([round(y_values(1)/10)*10, round(y_values(round(length(y_values)*0.8))/10)*10]);
            set(gca, 'YTickLabelRotation', 90);
        end
    end

    IM_file_axis = sprintf('%d_Final_Range_%d', i,num);
    IM_name_axis = fullfile(path_range, IM_file_axis);
    saveas(f2,IM_name_axis,'png')
end
