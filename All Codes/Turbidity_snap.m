clear; clc; clf; close all;
format shortg

%% Video Input

disp('Select Video:')
[file_t, path_t] = uigetfile('*.avi', 'Select Video File');
filename_t = fullfile(path_t, file_t);
v_t = VideoReader(filename_t);
file_t
fprintf('Video is %.2f seconds long \n', v_t.duration)
fprintf('with %.0f total frames. \n \n', v_t.NumFrames)

chart_title = 'Input';
dims = [1 35];
fps = 24

prompt = "Shear Rate [1/s]";
srs = inputdlg(prompt,chart_title,dims);
sr = str2double(srs)
if isempty(sr) == 1
    return
end

prompt = "No. of Snaps";
% snapss = inputdlg(prompt,chart_title,dims);
snaps = 5
if isempty(snaps) == 1
    return
end

prompt = "Seconds Per Snap [s]";
trs = inputdlg(prompt,chart_title,dims);
tr = str2double(trs)
fpsnap = ceil(tr*fps);
if isempty(fpsnap) == 1
    return
end

L = 3;

fprintf('\n')
fprintf('Frames per Snap: %.f \n \n', fpsnap)

disp('Select Shear Stress File:')
[file_stress, path_stress] = uigetfile('*.csv', 'Select Stress File');
stress_file = fullfile(path_stress, file_stress);
stress_data = xlsread(stress_file);
time = stress_data(:,1);
strain = time*sr;
stress = stress_data(:,6);
file_stress

fprintf('\nFind %.d desired times:', snaps)

u = figure;
loglog(time,stress, 'o-');
ylabel('\sigma [Pa]')
xlabel('Time [s]')
xlim([0,round(v_t.duration*0.9)])
ylim([min(stress)/5, max(stress)*10])
grid on
grid minor
datacursormode(u,'on');
fprintf('\n')
disp('Close figure to continue')
uiwait(u)

times(1) = 0;
fprintf('\nTime for profile #1 = %ds \n \n', times(1))
for i = 2:snaps
    prompt = sprintf('Time for profile #%.d [s]', i);
    slice = inputdlg(prompt,chart_title,dims);
    time_i = str2double(slice)
    times(i) = str2double(slice)*fps;
    if isempty(times(i))
        fprintf('Code terminated: Missing Input.')
        return
    end
end

colors = {'#000000', '#FF0000', '#0000FF', '#00FF00', '#FF00FF', '#00FFFF' '#000000', '#FF0000', '#0000FF', '#00FF00', '#FF00FF', '#00FFFF'};

[~, name, ~] = fileparts(filename_t);
parts = strsplit(name, '\');
file = parts{end};

%% Adjust Cylinder Location

k = 0;
bframe{1,1} = read(v_t, 1);
bsnap_frame{1,1} = uint8(bframe{1,1}); % Averages into 1 frame and then scales pixel intesity 1-256
bI{1,1} = mean(bframe{1,1},1);


disp('Find Cylinder Locations:')
while k == 0
    close all
    h = figure;
    img = bframe{1,1};
    imagesc(img);
    title('Find Inner and OuterCylinder location')
    colormap gray
    impixelinfo
    datacursormode(h,'on');
    xlabel('Pixel Position')
    grid on
    uiwait(h)

    prompt = "Inner Cylinder Position [px]";
    INs = inputdlg(prompt,chart_title,dims);
    in = str2double(INs);
    if isempty(in) == 1
        return
    end

    prompt = "Outer Cylinder Position [px]";
    OUTs = inputdlg(prompt,chart_title,dims);
    out = str2double(OUTs);
    if isempty(in) == 1
        return
    end

    if out > in
        dummy_i = in;
        dummy_o = out;
        in = dummy_o;
        out = dummy_i;
    end

    gap = abs(out-in) + 1;
    x_roi = linspace(0, 1, gap);
    left_gap = linspace(x_roi(1)-(out-1)*(1/gap), x_roi(1)-(1/gap), out-1);
    right_gap = linspace(x_roi(end)+(1/gap), x_roi(end) + (v_t.width-in)*(1/gap), v_t.width-in);
    x_values = [left_gap, x_roi, right_gap];
    x_3 = linspace(0,L,v_t.NumFrames);
    y_values = x_3;

    close all
    h = figure;
    img = bframe{1,1};
    subplot(2,2,1)
    imagesc(img);
    colormap gray
    axis image
    impixelinfo
    datacursormode(h,'on');
    xlabel('Pixel Position')
    grid on
    xline(in,'m-','LineWidth', 2);
    xline(out,'m-','LineWidth', 2);

    x =  linspace(1,length(bI{1,1}), length(bI{1,1}));
    subplot(2,2,3)
    plot(x, bI{1,1}, 'r-');
    ylim([-10 255])
    xlim([1, max(x)])
    yline(0,'k--');
    xline(in, 'm--');
    xline(out, 'm--');
    xlabel('Pixel Position')
    ylabel('Column Intensity')

    subplot(2,2,[2,4])
    imagesc(img);
    colormap gray
    impixelinfo
    xlabel('Pixel Position')
    grid on
    xline(in,'m-','LineWidth', 2);
    xline(out,'m-','LineWidth', 2);
    set(gcf, 'WindowState', 'Maximized')

    prompt = "0 = Re-adjust | 1 = Continue";
    is = inputdlg(prompt,chart_title,dims);
    k = str2double(is);
    if isempty(k) == 1
        return
    end
end
close all

%% Intensity Calculation

back = 1; %Background frame No.
frame_num{1,1}(:,:) = round(linspace(1, back, fpsnap));
check = fpsnap/2;

for i = 2:snaps % Determining number of frames centered at snap value
    if check > times(i)
        frame_num{1,i}(:,:) = round(linspace(2, times(i) + fpsnap/2, fpsnap));
    else
        frame_num{1,i}(:,:) = round(linspace(times(i) - fpsnap/2, times(i) + fpsnap/2, fpsnap));
    end
end

for j = 2:snaps %Converts frame num array into array of frames
    for i = 1:length(frame_num{1,j})
        frames{j,i} = read(v_t,frame_num{1,j}(i)); %#ok<*SAGROW>
        frames{j,i} = im2double(frames{j,i})*255;
    end
end

[numRows, numCols] = size(frames{1,1}); % Image size (ex: 1280px x 304px)

for j = 2:snaps
    sum_in{j,1} = frames{j,1}; % Initialize sums to average frames into snap for each snap
end

for j = 2:snaps %Sum all frames within temporal resolution
    for i = 1:length(frames)
        sum_in{j,1} = sum_in{j,1} + frames{j,i};
    end
end
snap_frame{1,1} = bsnap_frame{1,1};
for i = 2:snaps
    snap_frame{1,i} = uint8(sum_in{i,1}/fpsnap); % Averages into 1 frame and then scales pixel intesity 1-256
end

for i = 1:snaps
    I{1,i} = mean(snap_frame{1,i});
    % I{1,i} = mean(snap_frame{1,i})- bI{1,1};
end

for i = 1:snaps % Points in shear stress curve
    target_x = (times(i)/fps)*sr;
    valid_indices = find(strain >= target_x);

    if ~isempty(valid_indices)
        idx = valid_indices(1);
        closestX(i) = strain(idx);
        closestY(i) = stress(idx);

        if isnan(closestY(i))
            right = idx + 1;

            while right <= length(stress)
                if ~isnan(stress(right))
                    closestX(i) = strain(right);
                    closestY(i) = stress(right);
                    break;
                end
                right = right + 1;
            end
        end
    else
        closestX(i) = strain(end);
        closestY(i) = stress(end);
    end
end

%% Crop for Display

disp_snap_frame{1,1} = bsnap_frame{1,1};
for i = 2:snaps
    temp_snap_frame{1,i} = snap_frame{1,i};
    disp_snap_frame{1,i} = temp_snap_frame{1,i};
end

%% Plot

close all
p = 14;
h1 = figure;
hold on
for i = 1:snaps
    subplot(4,7,[i,i+7])
    imagesc(x_values, y_values, disp_snap_frame{1,i})
    colormap gray
    xlim([-inf 1.1])
    if i == 1
        ylabel('x_3 [mm]', 'FontSize', p)
        title(sprintf('Background \\gamma = %.f', 0));
    end
    if i > 1
        set(gca,'YTick',[])
        title(sprintf('\\gamma = %.f', closestX(i)));
    end
    xlabel('Radial Position')
end

for i = 1:snaps
    subplot(4,7, 14 + i)
    plot(x_values, I{1,i}, '-', 'Color', colors{i});
    if i == 1
        ylabel('Column Intensity')
    end
    if i >1
        % set(gca,'YTick',[])
    end
    xlabel('Radial Position')
    ylim([-255 255])
    xlim([-0.05 1.1])
    yline(0,'k--');
    xline(0, 'k--');
    xline(1, 'k--');
end

subplot(4,7,[6,7,13,14,20,21])
hold on
for i = 1:snaps
    plot(x_values, I{1,i}, '-', 'color', colors{i});
    title('Column Intensity Difference')
    xlabel('Radial Position')
    ylim([-50 255])
    xlim([-0.05 1.1])
end
yline(0, 'k--');
xline(0, 'k--');
xline(1, 'k--');
box on
grid on

subplot(4,7,[22,28])
semilogx(strain, stress, 'o-');
hold on
for i = 1:snaps
    semilogx(closestX(i), closestY(i), 'ko', 'MarkerFaceColor', sprintf('%s', colors{i}))
end
ylabel('\sigma [Pa]', 'FontSize', p)
xlabel('\gamma', 'FontSize', p+2)
ylim([0,50])
set(gcf, 'WindowState', 'maximized');

%% Wobble check

num_frames = 2400;

concatenated_frames = [];

for i = 1:num_frames
    frame = read(v_t, i); % Read the i-th frame
    concatenated_frames = [concatenated_frames; frame(100,:)]; % Concatenate frames vertically
end

h = figure;
imagesc(concatenated_frames);
title('Concatenated Frames');
xline(142, 'm')
xline(160, 'm')
colormap jet
xline(in, 'r-','LineWidth', 2)
xline(out, 'r-', 'LineWidth',2)

%% Save Final Workspace

sr
disp('Saving workspace...')
save(sprintf('Turb_Snap_sr%0.1f.mat', sr));
disp('Done.')
