clear; clc; clf; close all;
format shortg
warning off

%% Video Input

disp('Select Video:')
[file_t, path_t] = uigetfile('*', 'Select Video File');
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

prompt = "Initial Time[s]";
tis = inputdlg(prompt,chart_title,dims);
ti = str2double(tis)
if isempty(ti) == 1
    return
end

prompt = "Final Time[s]";
tfs = inputdlg(prompt,chart_title,dims);
tf = str2double(tfs)
if isempty(tf) == 1
    return
end

prompt = "Seconds Per Snap [s]";
trs = inputdlg(prompt,chart_title,dims);
tr = str2double(trs);
fpsnap = ceil(tr*fps);
if isempty(fpsnap) == 1
    return
end

snaps = floor((tf-ti)/tr)

fprintf('\n')
fprintf('Frames per Snap: %.f \n \n', fpsnap)

colors = {'#000000', '#FF0000', '#0000FF', '#00FF00', '#FF00FF', '#00FFFF' '#000000', '#FF0000', '#0000FF', '#00FF00', '#FF00FF', '#00FFFF'};

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

for j = 1:snaps %Converts frame num array into array of frames
    for i = 1:length(frame_num{1,j})
        frames{j,i} = read(v_t,frame_num{1,j}(i)); %#ok<*SAGROW>
        frames{j,i} = im2double(frames{j,i})*255;
    end
end

[numRows, numCols] = size(frames{1,1}); % Image size (ex: 1280px x 304px)

for j = 1:snaps
    sum_in{j,1} = frames{j,1}; % Initialize sums to average frames into snap for each snap
end

for j = 1:snaps %Sum all frames within temporal resolution
    for i = 1:size(frames,2)
        sum_in{j,1} = sum_in{j,1} + frames{j,i};
    end
end

for i = 1:snaps
    snap_frame{1,i} = uint8(sum_in{i,1}/fpsnap); % Averages into 1 frame and then scales pixel intesity 1-256
end

for i = 1:snaps
    I{1,i} = mean(snap_frame{1,i});
end

%% Plot

close all
figure
plot(1-x_values, bI{1,1}, 'b-'); % Reverse x_values here
hold on
for i = 1:snaps
    plot(1-x_values, I{1,i}, '-') % Reverse x_values here
    hold on
    ylabel('Column Intensity')
    ylim([-255 255])
    xlabel('Radial Position')
    xlim([-0.1 1.1])
    yline(0,'k--');
    xline(0, 'k--');
    xline(1, 'k--');
    set(gca, 'XDir','reverse')

    % legendEntries{i} = sprintf('t = %0.1fs', i*fpsnap/fps); % Create legend entry
end

% legend(legendEntries, 'Location', 'NorthWest'); % Add legend with all entries after loop


%% Wobble check

num_frames = 200;

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

