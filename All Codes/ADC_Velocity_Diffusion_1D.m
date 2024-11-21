clear; clc; clf; close all
format shortg
warning off
folder = 'C:\Users\alfre\OneDrive - Florida State University\Documents\FSU\Graduate\Research\PhD\RheoNMR\RheoNMR_Experiments\Results';
guessesPerDecade = 4;
options = optimset('Display', 'off');
pgs = 1; % 1 PGSE | 2 D-PGSE

%% Experiment Details

chart_title = 'Input';
dims = [1 35];
pgss = 1;

if pgs == 1
    disp('PGSE')
elseif pgs == 2
    disp('D-PGSE')
end
if isempty(pgss) == 1
    fprintf('Code terminated by user.')
    return
end
if pgs > 2
    fprintf('Invalid Input. Enter 1 for PGSE or 2 for D-PGSE')
    return
end
prompt = "Experiment Number";
exp_nums = inputdlg(prompt, chart_title, dims);
exp_num = str2double(exp_nums)
if isempty(exp_num) == 1
    fprintf('Code terminated by user.')
    return
end
prompt = "Number of Slices (B values)";
numSlicess = inputdlg(prompt, chart_title, dims);
numSlices = str2double(numSlicess)
if isempty(numSlices) == 1
    fprintf('Code terminated by user.')
    return
end
if pgs == 1
    prompt = "little delta [ms]";
    deltas = inputdlg(prompt, chart_title, dims);
    delta = str2double(deltas)/1000 %delta in s
    if isempty(delta) == 1
        fprintf('Code terminated by user.')
        return
    end
    prompt = "Big Delta [ms]";
    Deltas = inputdlg(prompt, chart_title, dims);
    Delta2 = str2double(Deltas)/1000 %Big Delta in s
    if isempty(Deltas) == 1
        fprintf('Code terminated by user.')
        return
    end
end
if pgs == 2
    prompt = "little delta [ms]";
    deltas = inputdlg(prompt, chart_title, dims);
    delta = str2double(deltas)/1000 %delta in s
    if isempty(delta) == 1
        fprintf('Code terminated by user.')
        return
    end
    prompt = "Big Delta 1 [ms]";
    Delta1s = inputdlg(prompt, chart_title, dims);
    Delta1 = str2double(Delta1s)/1000 %Big Delta1 in s
    if isempty(Delta1) == 1
        fprintf('Code terminated by user.')
        return
    end
    prompt = "Big Delta 2 [ms]";
    Delta2s = inputdlg(prompt, chart_title, dims);
    Delta2 = str2double(Delta2s)/1000 %Big Delta2 in s
    if isempty(Delta2) == 1
        fprintf('Code terminated by user.')
        return
    end
end
prompt = "Initial Gradient Strength [mT/m]";
g_ins = inputdlg(prompt, chart_title, dims);
g_in = str2double(g_ins)
if isempty(g_ins) == 1
    fprintf('Code terminated by user.')
    return
end
prompt = "Final Gradient Strength [mT/m]";
g_fs = inputdlg(prompt, chart_title, dims);
g_f = str2double(g_fs)
if isempty(g_fs) == 1
    fprintf('Code terminated by user.')
    return
end
prompt = "Shear Rate [1/s]";
srs = inputdlg(prompt, chart_title, dims);
sr = str2double(srs)
if isempty(sr) == 1
    fprintf('Code terminated by user.')
    return
end
d = 1;
V_imposed = sr*d
omega_rad = sr*d/(7.65);
omega_rps = sr*d/(7.65*2*pi);

if pgs == 2
    Delta = (Delta1-Delta2)/2;
else
    Delta = Delta2;
end

rotate = 1; %Assume scanning top of TC gap

bin = 0; %No bin for now
resultFolder = exp_nums{1};

mkdir(resultFolder);
g = linspace(g_in, g_f, numSlices)*0.001; % Gradient strengths in T/m
gZeroMax = max(g); % Gradient strength of last zero filled q
gamma = 2.675e8; % Gyromagnetic ratio in rad/(s*T)
slices_lin = linspace(1, 2*numSlices, 2*numSlices); %DICOMs x2 (Real and Im)
real_qSpaceData = cell(1, numSlices);
im_qSpaceData = cell(1, numSlices);

%% Read DICOMs

for i = 1:numSlices*2
    dicomData = dicomread(sprintf('MRIm%d_%d.dcm', exp_num, slices_lin(i)));
    dicomInfo = dicominfo(sprintf('MRIm%d_%d.dcm', exp_num, slices_lin(i)));

    scalingFactor(i) = 511.98;
    intercept(i) = dicomInfo.RescaleIntercept;

    if i <= numSlices
        real_qSpaceData{i} = double(dicomData)*scalingFactor(i) + intercept(i); % First half correspond to real data
    else
        im_qSpaceData{i-numSlices} = double(dicomData)*scalingFactor(i) + intercept(i); % Last half correspond to imaginary data
    end
end

for i = 1:numSlices
    if rotate == 1
        real_qSpaceData{i} = rot90(real_qSpaceData{i},-1);
        im_qSpaceData{i} = rot90(im_qSpaceData{i},-1);
    end
end
for i = 1:numSlices
    real_qSpaceData{i} = real_qSpaceData{i}(1,1:end);
    im_qSpaceData{i} = im_qSpaceData{i}(1,1:end);
end

for i = 1:numSlices
    qSpaceData{i} = real_qSpaceData{i} + 1i*im_qSpaceData{i};
end

test_real_slice1 = real_qSpaceData{1}; %Confirm with PV that scaling is correct
test_im_slice1 = im_qSpaceData{1};

for i = 1:numSlices
    dicomData = dicomread(sprintf('MRIm%d_%d_Mag.dcm', exp_num, slices_lin(i)));
    dicomInfo = dicominfo(sprintf('MRIm%d_%d_Mag.dcm', exp_num, slices_lin(i)));

    intercept(i) = dicomInfo.RescaleIntercept;
    mag_qSpaceData{i} = double(dicomData)*scalingFactor(i) + intercept(i);
end

for i = 1:numSlices
    if rotate == 1
        mag_qSpaceData{i} = rot90(mag_qSpaceData{i}, -1);
    end
end
for i = 1:numSlices
    mag_qSpaceData{i} = mag_qSpaceData{i}(1,1:end);
end
test_mag_slice1 = mag_qSpaceData{1}; %Confirm with PV that scaling is correct

%% Partial Save

pre_saved = 1;
disp('Saving Raw Workspace...')
if bin == 1
    save(sprintf('Raw_work_%d_bin.mat', exp_num));
else
    save(sprintf('Raw_work_%d.mat', exp_num));
end
disp('Done (Raw).')

%% Adjust Cylinder Location

i = 0;
while i == 0
    close all
    u = figure;
    imagesc(test_mag_slice1);
    title('Find Inner and OuterCylinder location')
    colormap gray
    xlabel('Pixel Position')
    grid on
    axis image
    set(gcf, 'WindowState', 'maximized');
    center = 1;
    prompt = "Right Wall Position [px]";
    INs = inputdlg(prompt,chart_title,dims);
    in = str2double(INs)
    if isempty(in) == 1
        return
    end

    prompt = "Left Wall Position [px]";
    OUTs = inputdlg(prompt,chart_title,dims);
    out = str2double(OUTs)
    if isempty(in) == 1
        return
    end

    h = 1;

    close all
    u = figure;
    imagesc(test_mag_slice1);
    colormap gray
    axis image
    impixelinfo
    datacursormode(u,'on');
    xlabel('Pixel Position')
    grid on
    set(gcf, 'WindowState', 'maximized');

    if out > in
        dummy_i = in;
        dummy_o = out;
        in = dummy_o;
        out = dummy_i;
    end

    line([in, in],[center + round(h/2), center - round(h/2)], 'Color', 'm', 'LineWidth', 2);
    line([out, out],[center + round(h/2), center - round(h/2)], 'Color', 'm', 'LineWidth', 2);
    yline(center,'m-','LineWidth', 2);

    prompt = "0 = Re-adjust | 1 = Continue";
    is = inputdlg(prompt,chart_title,dims);
    i = str2double(is);
    if isempty(i) == 1
        return
    end
end
close all

%% Bin Gap Area

if bin == 1
    for k = 1:numSlices
        result = [];
        a = real_qSpaceData{k};
        a = flip(a);
        bin_start = 64-in;
        bin_end = 64-out;

        bin_elements = zeros(1, 10);
        if bin_end > bin_start

            bin_size = floor((bin_end - bin_start + 1) / num_bins);
            remainder = rem((bin_end - bin_start + 1), num_bins);
            idx = bin_start;
            for i = 1:num_bins
                if i <= remainder
                    next_idx = idx + bin_size;
                    bin_elements(i) = bin_size + 1;
                else
                    next_idx = idx + bin_size - 1;
                    bin_elements(i) = bin_size;
                end
                result = [result, mean(a(idx:next_idx))];
                idx = next_idx + 1;
            end
        end
        result = flip([a(1:bin_start-1), result, a(bin_end+1:end)]);
        real_qSpaceData{k} = [];
        real_qSpaceData{k} = result;

        result = [];
        a = []; %#ok<NASGU>
        a = im_qSpaceData{k};
        a = flip(a);
        bin_start = 64-in;
        bin_end = 64-out;

        bin_elements = zeros(1, 10);
        if bin_end > bin_start

            bin_size = floor((bin_end - bin_start + 1) / num_bins);
            remainder = rem((bin_end - bin_start + 1), num_bins);
            idx = bin_start;
            for i = 1:num_bins
                if i <= remainder
                    next_idx = idx + bin_size;
                    bin_elements(i) = bin_size + 1;
                else
                    next_idx = idx + bin_size - 1;
                    bin_elements(i) = bin_size;
                end
                result = [result, mean(a(idx:next_idx))];
                idx = next_idx + 1;
            end
        end
        result = flip([a(1:bin_start-1), result, a(bin_end+1:end)]);
        im_qSpaceData{k} = [];
        im_qSpaceData{k} = result;

        result = [];
        a = []; %#ok<NASGU>
        a = mag_qSpaceData{k};
        a = flip(a);
        bin_start = 64-in;
        bin_end = 64-out;

        bin_elements = zeros(1, 10);
        if bin_end > bin_start

            bin_size = floor((bin_end - bin_start + 1) / num_bins);
            remainder = rem((bin_end - bin_start + 1), num_bins);
            idx = bin_start;
            for i = 1:num_bins
                if i <= remainder
                    next_idx = idx + bin_size;
                    bin_elements(i) = bin_size + 1;
                else
                    next_idx = idx + bin_size - 1;
                    bin_elements(i) = bin_size;
                end
                result = [result, mean(a(idx:next_idx))];
                idx = next_idx + 1;
            end
        end
        result = flip([a(1:bin_start-1), result, a(bin_end+1:end)]);
        mag_qSpaceData{k} = [];
        mag_qSpaceData{k} = result;
    end
end

%% Define ROI

if h == 1
    topRow = center;
    bottomRow = center;
else
    topRow = center - round(h/2);
    bottomRow = center + round(h/2);
end
leftCol = out;
if bin == 1
    rightCol = out+num_bins-1;
else
    rightCol = in;
end

%% Pixel Selection Display

selectedPixels = zeros((bottomRow - topRow + 1)*(rightCol - leftCol + 1), 2);
index = 1;

for col = leftCol:1:rightCol
    for row = bottomRow:-1:topRow
        selectedPixels(index,:) = [row, col];
        index = index + 1;
    end
end

numPixels = size(selectedPixels, 1);
height = bottomRow - topRow + 1;
gap = rightCol - leftCol + 1;
loc_row = linspace(topRow, bottomRow, bottomRow - topRow + 1);
loc_col = linspace(rightCol, leftCol, (rightCol - leftCol));
a_row = round(linspace(loc_row(1), loc_row(end), height + 1));
a_col = round(linspace(loc_col(1), loc_col(1), height + 1));

f1 = figure;
imagesc(mag_qSpaceData{1});
colormap gray
axis image
xlabel('Pixel Position')
grid on
set(gcf, 'WindowState', 'maximized');

line([out, out],[center + round(h/2), center - round(h/2)], 'Color', 'm', 'LineWidth', 2);
if bin == 1
    line([out+num_bins-1, out+num_bins-1],[center + round(h/2), center - round(h/2)], 'Color', 'm', 'LineWidth', 2);
else
    line([in, in],[center + round(h/2), center - round(h/2)], 'Color', 'm', 'LineWidth', 2);
end

yline(center,'m-','LineWidth', 2);
title('Area Selcted', 'FontSize', 20)
set(gcf, 'WindowState', 'maximized');

file_name = sprintf('Region_selection_%d', exp_num);
fullName = fullfile(resultFolder, file_name);
saveas(f1,fullName,'fig');
saveas(f1,fullName,'png');

%% q-Space Data Images

f2 = figure;
for slice = 1:round(numSlices/2)
    subplot(round(numSlices/2), 2, slice*2 - 1)
    imagesc(real_qSpaceData{slice})
    colormap gray
    title(['Real Slice ', num2str(slice)], 'FontSize', 14)
    axis on

    subplot(round(numSlices/2), 2, slice*2)
    imagesc(im_qSpaceData{slice})
    colormap gray
    title(['Imaginary Slice ', num2str(slice)], 'FontSize', 14)
    axis on
end
sgtitle('Pixel Intensity in q-Space', 'FontSize', 16)
set(gcf, 'WindowState', 'maximized');

file_name = sprintf('q-images_1_%d', exp_num);
fullName = fullfile(resultFolder, file_name);
saveas(f2,fullName,'fig');
saveas(f2,fullName,'png');

f3 = figure;
for slice = round(numSlices/2) + 1:numSlices
    win = 1:2:numSlices;
    win2 = 2:2:numSlices;
    subplot(round(numSlices/2), 2, win(slice-round(numSlices/2)))
    imagesc(real_qSpaceData{slice})
    colormap gray
    title(['Real Slice ', num2str(slice)], 'FontSize', 14)
    axis on

    subplot(round(numSlices/2), 2, win2(slice-round(numSlices/2)))
    imagesc(im_qSpaceData{slice})
    colormap gray
    title(['Imaginary Slice ', num2str(slice)], 'FontSize', 14)
    axis on
end
sgtitle('Pixel Intensity in q-Space')
set(gcf, 'WindowState', 'maximized');
file_name = sprintf('q-images_2_%d', exp_num);
fullName = fullfile(resultFolder, file_name);
saveas(f3,fullName,'fig');
saveas(f3,fullName,'png');


%% Build FID

N = 4048; % Total points in FID and DDP

for pixel = 1:numPixels
    row = selectedPixels(pixel, 1);
    col = selectedPixels(pixel, 2);
    for slice = 1:N
        if slice <= numSlices
            FID_Real(slice,:) = real_qSpaceData{slice}(row, col);
            FID_IM(slice,:) = im_qSpaceData{slice}(row, col);
            pixelFID(slice,:) = FID_Real(slice,:) + 1i*FID_IM(slice,:);
        else
            pixelFID(slice,:) = 0; %Zero fill
        end
        FID{pixel} = pixelFID;
    end
end

qValues = zeros(1, numSlices);

for slice = 1:numSlices
    if pgs == 1
        q = gamma*g(slice)*delta;
        B = q^2*(Delta-delta/3)*10^-6;
    end
    if pgs == 2
        q = gamma*g(slice)*delta;
        B = q^2*(2*Delta-delta/3)*10^-6
    end
    qValues(slice) = q;
    BValues(slice) = B/1.102378; %Adjust for actual PV output
end
BValues = BValues';

qValuesZeroFilled = zeros(1, N);
g_step = g(2)-g(1);
g_0 = g(end);

for i = 1:N-length(g)
    g_extra(i) = g_0 + g_step;
    g_0 = g_extra(i);
end

g_ZeroFilled = [g, g_extra];

for slice = 1:N
    if pgs == 1
        qZero = gamma*g_ZeroFilled(slice)*delta;
    end
    if pgs == 2
        qZero = gamma*g_ZeroFilled(slice)*delta;
    end
    qValuesZeroFilled(slice) = qZero;
end

%% Display FID

pix = linspace(1,numPixels,numPixels);
window = 1:2:numPixels*2;
f4 = figure;
if numPixels < 16
    j = numPixels;
else
    j = 16;
end

for i = 1:j
    subplot(4,8, window(i))
    semilogx(qValuesZeroFilled, real(FID{pix(i)}), 'o - ')
    xlabel('q [m^{-1}]')
    if bin == 1
        title(sprintf('Real - Bin %d', pix(i)))
    else
        title(sprintf('Real - Pixel %d', pix(i)))
    end

    subplot(4,8, window(i)+1)
    semilogx(qValuesZeroFilled, imag(FID{pix(i)}), 'o - ')
    xlabel('q [m^{-1}]')
    if bin == 1
        title(sprintf('maginary - Bin %d', pix(i)))
    else
        title(sprintf('Imaginary - Pixel %d', pix(i)))
    end
    sgtitle('FID of pixels of q-Space Images')
end
set(gcf, 'WindowState', 'maximized');

file_name = sprintf('FIDs_%d_1', exp_num);
fullName = fullfile(resultFolder, file_name);
saveas(f4,fullName,'fig');
saveas(f4,fullName,'png');

if numPixels < 32
    j = numPixels;
else
    j = 32;
end

if numPixels > 16
    f4_2 = figure;
    for i = 17:j
        subplot(4,8, window(i)-16*2)
        semilogx(qValuesZeroFilled, real(FID{pix(i)}), 'o - ')
        xlabel('q [m^{-1}]')
        if bin == 1
            title(sprintf('Real - Bin %d', pix(i)))
        else
            title(sprintf('Real - Pixel %d', pix(i)))
        end

        subplot(4,8, window(i)+1-16*2)
        semilogx(qValuesZeroFilled, imag(FID{pix(i)}), 'o - ')
        xlabel('q [m^{-1}]')
        if bin == 1
            title(sprintf('maginary - Bin %d', pix(i)))
        else
            title(sprintf('Imaginary - Pixel %d', pix(i)))
        end
        sgtitle('FID of pixels of q-Space Images')
    end
    set(gcf, 'WindowState', 'maximized');

    file_name = sprintf('FIDs_%d_2', exp_num);
    fullName = fullfile(resultFolder, file_name);
    saveas(f4_2,fullName,'fig');
    saveas(f4_2,fullName,'png');
end

if numPixels < 48
    j = numPixels;
else
    j = 48;
end
if numPixels > 32
    f4_3 = figure;
    for i = 33:j
        subplot(4,8, window(i)-32*2)
        semilogx(qValuesZeroFilled, real(FID{pix(i)}), 'o - ')
        xlabel('q [m^{-1}]')
        if bin == 1
            title(sprintf('Real - Bin %d', pix(i)))
        else
            title(sprintf('Real - Pixel %d', pix(i)))
        end

        subplot(4,8, window(i)+1-32*2)
        semilogx(qValuesZeroFilled, imag(FID{pix(i)}), 'o - ')
        xlabel('q [m^{-1}]')
        if bin == 1
            title(sprintf('maginary - Bin %d', pix(i)))
        else
            title(sprintf('Imaginary - Pixel %d', pix(i)))
        end
        sgtitle('FID of pixels of q-Space Images')
    end
    set(gcf, 'WindowState', 'maximized');

    file_name = sprintf('FIDs_%d_3', exp_num);
    fullName = fullfile(resultFolder, file_name);
    saveas(f4_3,fullName,'fig');
    saveas(f4_3,fullName,'png');
end

if numPixels < 64
    j = numPixels;
else
    j = 64;
end

if numPixels > 48
    f4_4 = figure;
    for i = 49:j
        subplot(4,8, window(i)-48*2)
        semilogx(qValuesZeroFilled, real(FID{pix(i)}), 'o - ')
        xlabel('q [m^{-1}]')
        if bin == 1
            title(sprintf('Real - Bin %d', pix(i)))
        else
            title(sprintf('Real - Pixel %d', pix(i)))
        end

        subplot(4,8, window(i)+1-48*2)
        semilogx(qValuesZeroFilled, imag(FID{pix(i)}), 'o - ')
        xlabel('q [m^{-1}]')
        if bin == 1
            title(sprintf('maginary - Bin %d', pix(i)))
        else
            title(sprintf('Imaginary - Pixel %d', pix(i)))
        end
        sgtitle('FID of pixels of q-Space Images')
    end
    set(gcf, 'WindowState', 'maximized');

    file_name = sprintf('FIDs_%d_3', exp_num);
    fullName = fullfile(resultFolder, file_name);
    saveas(f4_4,fullName,'fig');
    saveas(f4_4,fullName,'png');
end

%% Dynamic displacement profile in Z-space

Naxis = linspace(-N/2, N/2-1, N);
npts = 2*N;
x_continuous = linspace(min(Naxis),max(Naxis),npts);

v = [];
D = [];
k_fwhm = [];
k_v = [];

for pixel = 1:numPixels
    ddp = fft(FID{pixel});
    ddp = fftshift(ddp);
    DDP{pixel} = ddp;
    [Maxc, idx1c] = max(abs(DDP{pixel}));

    y_continuous = pchip(Naxis, abs(DDP{pixel}), x_continuous);
    x1 = x_continuous(1:idx1c);
    x2 = x_continuous(idx1c:end);
    y1 = y_continuous(1:idx1c);
    y2 = y_continuous(idx1c:end);
    [Max, idx1] = max(abs(DDP{pixel}));
    [Maxref, idx1c] = max(y_continuous);
    halfMax = Maxref/2;

    intercepts_x = [];
    intercepts_y = [];

    for i = 2:length(x_continuous)
        if (y_continuous(i-1) > halfMax && y_continuous(i) <= halfMax) || (y_continuous(i-1) < halfMax && y_continuous(i) >= halfMax)
            m = (y_continuous(i) - y_continuous(i-1)) / (x_continuous(i) - x_continuous(i-1));
            c = y_continuous(i) - m*x_continuous(i);
            intercept_x = (halfMax-c)/ m;
            intercepts_x = [intercepts_x intercept_x];
            intercepts_y = [intercepts_y halfMax];
        end
    end

    reference_x = x_continuous(idx1c);
    reference_y = halfMax;
    distances_left = abs(intercepts_x(intercepts_x < reference_x) - reference_x);
    distances_right = abs(intercepts_x(intercepts_x > reference_x) - reference_x);
    w1 = intercepts_x(intercepts_x < reference_x);
    w2 = intercepts_x(intercepts_x > reference_x);
    [~, idx_left] = min(abs(w1 - reference_x));
    w1 = w1(idx_left);
    [~, idx_right] = min(abs(w2 - reference_x));
    w2 = w2(idx_right);

    if isempty(w1) == 1 || isempty(w2) == 1
        w1 = NaN;
        w2 = NaN;
        fwhm = NaN;
    else
        fwhm = abs(w2 - w1);
    end

    halfMaxValues{pixel} = halfMax;
    k_fwhm(pixel) = fwhm;
    k_v(pixel) = Naxis(idx1);

    %         if isnan(k_fwhm(pixel))
    %             k_v(pixel) = NaN;
    %         else
    %             k_v(pixel) = Naxis(idx1);
    %         end

    if pixel == 1
        f5 = figure;
    end

    if pixel >= 1 && pixel <= 24
        subplot(4,6, pixel)
        plot(Naxis, abs(DDP{pixel}), 'o')
        hold on
        yline(halfMax, 'k--');
        xline(Naxis(idx1), 'k--');
        plot(x_continuous, y_continuous, '-')
        sgtitle('Dynamic Displacement Profile','FontSize', 20)
        if ~isnan(w1)
            plot(w1, halfMax, 'kx')
            plot(w2, halfMax, 'kx')
            plot([w1, w2], [halfMax, halfMax], 'r', 'LineWidth', 2)
            text(Naxis(idx1), halfMax/1.5, sprintf('FWHM: %.1f', fwhm), 'HorizontalAlignment', 'center')
        end
        xlabel('Z')
        ylabel('Magnitude')
        if bin == 1
            title(sprintf('Bin #%d',pixel),'FontSize', 14)
        else
            title(sprintf('Pixel #%d',pixel),'FontSize', 14)
        end
        ylim([0 max(abs(DDP{pixel})*1.3)])
    end

    if pixel == 24 || numPixels < 24 && pixel == numPixels
        set(gcf, 'WindowState', 'maximized');
        file_name = sprintf('DDPs_%d_1', exp_num);
        fullName = fullfile(resultFolder, file_name);
        saveas(f5,fullName,'fig');
        saveas(f5,fullName,'png');
    end

    if pixel == 25
        f6 = figure;
    end

    if pixel >= 25 && pixel <= 48
        subplot(4,6, pixel-24)
        plot(Naxis, abs(DDP{pixel}), 'o')
        hold on
        yline(halfMax, 'k--');
        xline(Naxis(idx1), 'k--');
        plot(x_continuous, y_continuous, '-')
        sgtitle('Dynamic Displacement Profile','FontSize', 20)
        if ~isnan(w1)
            plot(w1, halfMax, 'kx')
            plot(w2, halfMax, 'kx')
            plot([w1, w2], [halfMax, halfMax], 'r', 'LineWidth', 2)
            text(Naxis(idx1), halfMax/1.5, sprintf('FWHM: %.1f', fwhm), 'HorizontalAlignment', 'center')
        end
        xlabel('Z')
        ylabel('Magnitude')
        if bin == 1
            title(sprintf('Bin #%d',pixel),'FontSize', 14)
        else
            title(sprintf('Pixel #%d',pixel),'FontSize', 14)
        end

        ylim([0 max(abs(DDP{pixel})*1.3)])
        set(gcf, 'WindowState', 'maximized');
    end

    if pixel == 48 || numPixels > 24 && numPixels < 48 && pixel == numPixels
        set(gcf, 'WindowState', 'maximized');
        file_name = sprintf('DDPs_%d_2', exp_num);
        fullName = fullfile(resultFolder, file_name);
        saveas(f6,fullName,'fig');
        saveas(f6,fullName,'png');
    end
    if pixel == 49
        f6_2 = figure;
    end

    if pixel >= 49 && pixel <= numPixels
        subplot(4,6, pixel-48)
        plot(Naxis, abs(DDP{pixel}), 'o')
        hold on
        yline(halfMax, 'k--');
        xline(Naxis(idx1), 'k--');
        plot(x_continuous, y_continuous, '-')
        sgtitle('Dynamic Displacement Profile','FontSize', 20)
        if ~isnan(w1)
            plot(w1, halfMax, 'kx')
            plot(w2, halfMax, 'kx')
            plot([w1, w2], [halfMax, halfMax], 'r', 'LineWidth', 2)
            text(Naxis(idx1), halfMax/1.5, sprintf('FWHM: %.1f', fwhm), 'HorizontalAlignment', 'center')
        end
        xlabel('Z')
        ylabel('Magnitude')
        if bin == 1
            title(sprintf('Bin #%d',pixel),'FontSize', 14)
        else
            title(sprintf('Pixel #%d',pixel),'FontSize', 14)
        end

        ylim([0 max(abs(DDP{pixel})*1.3)])
        set(gcf, 'WindowState', 'maximized');
    end

    if pixel == 64 || numPixels > 48 && numPixels < 64 && pixel == numPixels
        set(gcf, 'WindowState', 'maximized');
        file_name = sprintf('DDPs_%d_3', exp_num);
        fullName = fullfile(resultFolder, file_name);
        saveas(f6,fullName,'fig');
        saveas(f6,fullName,'png');
    end
    if pgs == 1
        v(pixel) = (2*pi*numSlices*k_v(pixel))/(N*gamma*delta*Delta*max(g));
        D(pixel) = (3.56*(numSlices*k_fwhm(pixel))^2)/(gamma^2*delta^2*(max(g))^2*N^2*Delta);
    end

    if pgs == 2
        v(pixel) = (2*pi*numSlices*k_v(pixel))/(N*gamma*delta*Delta*2*max(g));
        D(pixel) = (3.56*(numSlices*k_fwhm(pixel))^2)/(gamma^2*delta^2*(max(g))^2*N^2*Delta*2);
    end
end

pixel_number = linspace(1, numPixels, numPixels)';
v = v';
D = D';
k_fwhm = k_fwhm';
k_v = k_v';
if pixel == numPixels
    T1 = table(pixel_number, k_fwhm, k_v, v, D)
else
    T1 = table(k_fwhm, k_v, v, D)
end
T_file = sprintf('DDP_table_%d',exp_num);
T_file_csv = [T_file, '.csv'];
T_name = fullfile(resultFolder, T_file_csv);
delete(T_name)
writetable(T1, T_name)

%% Analytical Newtonian Velocity

Ri = 7.65; %mm
Ro = 8.65; %mm
d = Ro-Ri;
V_imposed = sr*d;
omega_rps = sr*d/(Ri*2*pi);
w = sr*d/(Ri);
r_gap = linspace(Ri,Ro,gap);
C1 = 2*w*Ri^2/(Ri^2-Ro^2);
C2 = -w*Ri^2*Ro^2/(Ri^2-Ro^2);
v_fun = @(x) (C1/2).*x + C2./x;
rgap_norm = (r_gap-min(r_gap))/(max(r_gap)-min(r_gap));

%% Lever Rule Prediction

LR_gamma_app = sr;
LR_d = 1;
LR_v_wall = LR_gamma_app / LR_d;
LR_gamma_L = 1.2;
LR_gamma_H = 20;
LR_Ri = 7.65;
LR_Ro = 8.65;

LR_a_H = (LR_gamma_app - LR_gamma_L) / (LR_gamma_H - LR_gamma_L);
LR_r_b = LR_Ri + LR_a_H * (LR_Ro - LR_Ri);
LR_r = linspace(LR_Ri, LR_Ro, 1000);
LR_v = zeros(size(LR_r));
LR_v(LR_r <= LR_r_b) = LR_v_wall - LR_gamma_H * (LR_r(LR_r <= LR_r_b) - LR_Ri);
LR_v(LR_r > LR_r_b) = (LR_v_wall - LR_gamma_H * (LR_r_b - LR_Ri)) - LR_gamma_L * (LR_r(LR_r > LR_r_b) - LR_r_b);
LR_r_norm = (LR_r - LR_Ri) / (LR_Ro - LR_Ri);

%% Averaging Columns (resets shift)

numAverages = floor(length(v)/height);
for i = 1:numAverages
    startIndex = (i-1)*height + 1;
    endIndex = i*height;
    averaged_v(i) = mean(v(startIndex:endIndex), 'omitnan')*1000;
    std_v(i) = std(v(startIndex:endIndex),'omitnan')*1000;
end

numAverages = floor(length(D)/height);
for i = 1:numAverages
    startIndex = (i-1)*height + 1;
    endIndex = i*height;
    averaged_D(i) = mean(D(startIndex:endIndex), 'omitnan');
    std_D(i) = std(D(startIndex:endIndex), 'omitnan');
end
if V_imposed == 0
    a = mean(averaged_v);
else
    a = averaged_v(end);
end
if V_imposed == 0
    V_norm = 1;
else
    V_norm = V_imposed;
end

v_an = v_fun(r_gap);
r = linspace(0,1,gap);
rgap_norm = (r-min(r))/(max(r)-min(r));
v_plot = averaged_v-a;

%% Velocity and Diffusion - Plot

close all
v_plot = averaged_v-a;
norm = 0; % 1 = normalized, 0 = absolute
f7 = figure;
subplot(1,2,1)
if norm == 0
    plot(r, v_plot, 'o')
    hold on
    plot(rgap_norm, v_an, 'k--')
    if V_imposed > LR_gamma_L
        plot(LR_r_norm, LR_v, 'm--');
    end
else
    plot(r, averaged_v/V_norm, 'o')
    hold on
    % plot(rgap_norm, v_an/V_norm, 'k--')
    if V_imposed > LR_gamma_L
        plot(LR_r_norm, LR_v, 'm--');
    end
end
yline(0, 'k--', 'LineWidth', 1.5);
xlabel('Radial Position')

if norm == 0
    ylabel('Absolute Velocity [mm/s]')
else
    ylabel('Normalized Velocity')
end

title('Velocity')

xlim([-0.01 1.01])
if norm ==1
    ylim([-0.01 1.01])
else
    % ylim([-0.05 inf])
end
if V_imposed == 0
    ylabel('Absolute Velocity [mm/s]')
    ylim([-0.1 1])
end
% ---------------------------------------------------------------------------------
subplot(1,2,2)
semilogy(r, averaged_D, 'o')
hold on
xlabel('Radial Position')
ylabel('ADC [m^2/s]')
title('Diffusion')
yline(2.29*10^-9, 'k--', 'LineWidth', 2);
yline(pi*2.29*10^-9, 'r--', 'LineWidth', 1);
text(0.5, 1.5*averaged_D(round(end/2)), sprintf('D_{r=0.5} = %.2e m^2/s', averaged_D(round(end/2))), 'HorizontalAlignment', 'center', 'FontSize', 12)

xlim([-0.01 1.01])
ylim([10^-11, 10^-8])
% ylim([10^-10, 10^-7])
grid on

if bin == 1
    if pgs == 1
        sgtitle(sprintf('Experiment #%.d - PGSE     V = %.1f mm/s     Shear Rate = %.1f/s      \\Delta = %d ms       Bins = %d ', exp_num,V_imposed, sr, Delta*1000, num_bins))
    elseif pgs == 2
        sgtitle(sprintf('Experiment #%.d - D-PGSE    V = %.1f mm/s     Shear Rate = %.1f/s      \\Delta = %d ms      Bins = %d', exp_num,V_imposed, sr, Delta*1000, num_bins))
    end
else
    if pgs == 1
        sgtitle(sprintf('Experiment #%.d - PGSE     V = %.1f mm/s     Shear Rate = %.1f/s      \\Delta = %d ms', exp_num,V_imposed, sr, Delta*1000))
    elseif pgs == 2
        sgtitle(sprintf('Experiment #%.d - D-PGSE    V = %.1f mm/s     Shear Rate = %.1f/s      \\Delta = %d ms', exp_num,V_imposed, sr, 2*Delta*1000))
    end
end

set(gcf, 'WindowState', 'maximized');
file_name = sprintf('Velocity_Diffusion_%d', exp_num);
fullName = fullfile(resultFolder, file_name);
saveas(f7,fullName,'fig');
saveas(f7,fullName,'png');

radial_position = r';
velocity = v_plot';
velocity_norm = averaged_v./V_norm; velocity_norm = velocity_norm';
velocity_norm_std =  std_v./V_norm; velocity_norm_std = velocity_norm_std';
Diffusion = averaged_D; Diffusion = Diffusion';
Diffusion_std = std_D';
r_an = rgap_norm';
v_analy = v_an';

T2 = table(radial_position, velocity, velocity_norm, velocity_norm_std, Diffusion, Diffusion_std,r_an,v_analy)
T_file = sprintf('Velocity_Diffusion_table_%d',exp_num);
T_file_csv = [T_file, '.csv'];
T_name = fullfile(resultFolder, T_file_csv);
delete(T_name)
writetable(T2, T_name)
exp_num

%% Magnitude - Build FID

for pixel = 1:numPixels
    row = selectedPixels(pixel, 1);
    col = selectedPixels(pixel, 2);
    for slice = 1:numSlices
        pixelFID_Mag(slice,:) = mag_qSpaceData{slice}(row, col);
        FID_Mag{pixel} = pixelFID_Mag;
        FID_Mag{pixel} = FID_Mag{pixel}/FID_Mag{pixel}(1);
    end
end

%% Magnitude - Fit Monoexponential Plot

f8 = figure;
funm = @(D,x) exp(-(x-BValues(1)).*D(1));
D0 = 0.001;
options = optimoptions('lsqcurvefit', 'Display', 'off');
for pixel = 1:numPixels
    [fit, rsmm] = lsqcurvefit(funm, D0, BValues, FID_Mag{pixel},[],[], options);
    fitm(pixel) = fit*10^-6;
    rsm(pixel) = rsmm;
    if pixel >= 1 && pixel <=24
        subplot(4,6, pixel), semilogy(BValues, FID_Mag{pixel}, 'ko', BValues, funm(fitm(pixel)/10^-6, BValues), 'b-')
        mag_data(:,pixel) =  FID_Mag{pixel};
        mag_B(:,pixel) = BValues;
        % ylim([0.1,1])
        xlabel('B [s/mm^2]')
        ylabel('Relative Intensity')
        if bin == 1
            title(sprintf('Bin #%d',pixel))
        else
            title(sprintf('Pixel #%d',pixel))
        end
        sgtitle('Magnitude - Signal Attenuation')
    end
    if pixel == 24 || numPixels < 24 && pixel == numPixels
        set(gcf, 'WindowState', 'maximized');
        file_name = sprintf('Mag_fit_1_%d', exp_num);
        fullName = fullfile(resultFolder, file_name);
        saveas(f8,fullName,'fig');
        saveas(f8,fullName,'png');
    end

    if pixel == 25
        f9 = figure;
    end
    if pixel >= 25 && pixel <=48
        subplot(4,6, pixel-24), semilogy(BValues, FID_Mag{pixel}, 'ko', BValues, funm(fitm(pixel)/10^-6, BValues), 'b-')
        ylim([0.001,1])
        xlabel('B [s/mm^2]')
        ylabel('Relative Intensity')
        if bin == 1
            title(sprintf('bin #%d',pixel))
        else
            title(sprintf('Pixel #%d',pixel))
        end
        sgtitle('Magnitude - Signal Attenuation')
    end


    if pixel == 48 ||numPixels > 24 && numPixels < 48 && pixel == numPixels
        set(gcf, 'WindowState', 'maximized');
        file_name = sprintf('Mag_fit_2_%d', exp_num);
        fullName = fullfile(resultFolder, file_name);
        saveas(f9,fullName,'fig');
        saveas(f9,fullName,'png');
    end

    if pixel == 49
        f9_2 = figure;
    end
    if pixel >= 49 && pixel <=numPixels
        subplot(4,6, pixel-48), semilogy(BValues, FID_Mag{pixel}, 'ko', BValues, funm(fitm(pixel)/10^-6, BValues), 'b-')
        ylim([0.001,1])
        xlabel('B [s/mm^2]')
        ylabel('Relative Intensity')
        if bin == 1
            title(sprintf('bin #%d',pixel))
        else
            title(sprintf('Pixel #%d',pixel))
        end
        sgtitle('Magnitude - Signal Attenuation')
    end

    if pixel == 64 ||numPixels > 48 &&  numPixels < 64 && pixel == numPixels
        set(gcf, 'WindowState', 'maximized');
        file_name = sprintf('Mag_fit_2_%d', exp_num);
        fullName = fullfile(resultFolder, file_name);
        saveas(f9_2,fullName,'fig');
        saveas(f9_2,fullName,'png');
    end
end

numCols = size(mag_B, 2);
new_matrix = zeros(size(mag_B, 1), numCols * 2);
new_matrix(:, 1:2:end) = mag_B;
new_matrix(:, 2:2:end) = mag_data;
column_titles = cell(1, numCols * 2);
for i = 1:numCols
    column_titles{2*i-1} = sprintf('B%d', i);
    column_titles{2*i} = sprintf('S%d', i);
end
T9 = array2table(new_matrix, 'VariableNames', column_titles);
T_file = sprintf('MagFit_%d', exp_num);
T_file_csv = [T_file, '.csv'];
T_name = fullfile(resultFolder, T_file_csv);
delete(T_name)
writetable(T9, T_name)
winopen(T_name)

%% Magnitude - Diffusion Plot

if pixel == numPixels
    for i = 1:numAverages
        startIndex = (i-1)*height + 1;
        endIndex = i*height;
        averaged_Dmag(i) = mean(fitm(startIndex:endIndex), 'omitnan');
        std_Dmag(i) = std(fitm(startIndex:endIndex),'omitnan');
    end

    f9 = figure;
    subplot(1,1,1)
    semilogy(r, averaged_Dmag, 'o')
    hold on
    %         eb = errorbar(r, averaged_Dmag, std_Dmag,'vertical', 'LineStyle', 'none');
    %         set(eb, 'color', 'k', 'LineWidth', 1)
    yline(2.29*10^-9, 'k--', 'LineWidth', 2);

    xlabel('Radial Position')
    xlim([-0.01 1.01])

    ylabel('ADC [m^2/s]')
    %         set(gca, 'XDir','reverse')
    % ylim([10^-10 10^-7])
    ylim([10^-12, 10^-9])
    title('Magnitude - Diffusion')
end

if pgs == 1
    sgtitle(sprintf('Experiment #%.d - PGSE     V = %.1f mm/s     Shear Rate = %.1f/s', exp_num,V_imposed, sr))
elseif pgs == 2
    sgtitle(sprintf('Experiment #%.d - D-PGSE    V = %.1f mm/s     Shear Rate = %.1f/s', exp_num,V_imposed, sr))
end

file_name = sprintf('Diffusion_Mag_%d', exp_num);
fullName = fullfile(resultFolder, file_name);
saveas(f9,fullName,'fig');
saveas(f9,fullName,'png');

%% Visualize DICOMs

freq = linspace(0,100, length(mag_qSpaceData{1}));
f10 = figure;
for i = 1:numSlices
    subplot(1,3,1)
    plot(freq, mag_qSpaceData{i},'-')
    hold on
    axis square
    title('Magnitude')
    xlabel('Frequency [ppm]')
    set(gca, 'XDir','reverse')
    grid

    subplot(1,3,2)
    plot(freq, real_qSpaceData{i})
    hold on
    axis square
    title('Real')
    xlabel('Frequency [ppm]')
    set(gca, 'XDir','reverse')

    subplot(1,3,3)
    plot(freq, im_qSpaceData{i})
    hold on
    axis square
    xlabel('Frequency [ppm]')
    set(gca, 'XDir','reverse')
    title('Imaginary')

    sgtitle('DICOM Spectra')
end

file_name = sprintf('DICOMs_Check_%d', exp_num);
fullName = fullfile(resultFolder, file_name);
saveas(f10,fullName,'fig');
saveas(f10,fullName,'png');

%% Save Final Workspace

close all
Results_for = sprintf('Exp #%d, V = %.1f mm/s, Shear Rate = %.1f/s', exp_num,V_imposed, sr)
disp('Saving workspace...')
if bin == 1
    save(sprintf('work_%d_bin.mat', exp_num));
else
    save(sprintf('work_%d.mat', exp_num));
end
disp('Done.')
