clear; clc; clf; close all;
format shortg

%% Delta File Input

chart_title = 'Input';  %Title outside dialog box
dims = [1 35]; %Dimension of dialog box
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
    [file_delta, path_delta] = uigetfile('*srd.csv', 'Delta File');
    delta_file = fullfile(path_delta, file_delta);
    delta_data{i} = xlsread(delta_file);
    file_delta

    strain{i}(:,1) = delta_data{i}(:,1);
    delta{i}(:,1) = delta_data{i}(:,2);
    piece{i}(:,1) = delta_data{i}(:,3);
    sig_ez{i}(:,1) = delta_data{i}(:,4);
    delta_data = [];
end


%% Delta Plot
close all
n = 14;
width = 4.5;
height = 2;
f1 = figure('Units', 'inches', 'Position', [1, 1, width, height]);
for i = 1:profiles2D
    semilogx(strain{i}, delta{i},'b-')
    hold on
end
ylim([-0.1 3])
%     xlim([10^-0 10^5])
xlabel('\gamma')
ylabel('\Delta')
set(gca, 'FontSize', n);


delta_name = sprintf('%d_Delta', num);
saveas(f1,delta_name,'png')
num

