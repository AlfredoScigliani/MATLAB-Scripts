clear; clc; clf; close all;
format shortg
warning off
width = 32; % Pixel width per frame

%% Video Input

disp('Select Video:')
[file, path] = uigetfile('*.avi', 'Select Video File');
if isequal(file, 0)
    fprintf('Code terminated: No video file selected.\n')
    return
end
filename = fullfile(path, file);
v = VideoReader(filename);
disp(['Selected File: ', file])

fps = 24;

prompt = "Shear Rate [1/s] (0 for ramp)";
srs = inputdlg(prompt, chart_title, dims);
sr = str2double(srs);
if isempty(sr) || isnan(sr)
    fprintf('Code terminated: Missing or invalid Shear Rate.\n')
    return
end

t = 0;
[file_stress, path_stress] = uigetfile('*.csv', 'Select Stress File');
if isequal(file_stress, 0)
    fprintf('Code terminated: No stress file selected.\n')
    return
end
stress_file = fullfile(path_stress, file_stress);
stress_data = xlsread(stress_file);
shear_rate = stress_data(:,2);

if sr > 0
    strain = stress_data(:,1)*sr;
else
    strain = stress_data(:,1).*shear_rate;
end

stress = stress_data(:,6);

%% Frame Extraction and Concatenation

num = 1:2*fps:v.NumFrames;
frames = cell(1, length(num));

% Extract and store frames at 1 frame per second
for i = 1:length(num)
    frames{1,i} = read(v, num(i));
end

% Concatenate the frames horizontally
con = frames{1,1};
for i = 2:length(frames)
    con = [con, frames{1,i}];
end


%% Plot

p = 14;
if sr > 0
    x_values = linspace(min(strain), max(strain), size(con, 2));
else
    x_values = linspace(min(shear_rate), max(shear_rate), size(con, 2));
end

y_min = 0;
y_max = 32;
num_rows = size(con, 1);
y_values = linspace(y_min, y_max, num_rows); % Increase from 0 mm to 32 mm

h = figure;
% subplot(2,1,1)
imagesc(x_values, y_values, con)
ylabel('x_3 [mm]','Fontsize', p)
if sr > 0
    xlabel('\gamma','Fontsize', p+2)
else
    xlabel('$\dot{\gamma}$', 'Interpreter', 'latex','Fontsize', p+2)

end
colormap gray
if sr > 0
    title(sprintf('Flow Visualization - sr = %0.1f', sr))
else
    title(sprintf('Flow Visualization - sr = %0.1f to %0.1f', min(shear_rate), max(shear_rate)))
end
set(gca, 'YDir', 'normal');


%% Save Workspace
save(sprintf('Mica_sr_%0.1f.mat', sr))