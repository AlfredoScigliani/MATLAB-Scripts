clear; clc; clf; close all;
format shortg
resultFolder = 'C:\Users\alfre\OneDrive - Florida State University\Documents\FSU\Graduate\Research\PhD\RheoPTV\Results Summary\Paper_Figures';
path_stress = 'C:\Users\alfre\OneDrive - Florida State University\Documents\FSU\Graduate\Research\PhD\RheoPTV\Results Summary\All Data - E\';

%% Data Input

chart_title = 'Input';  %Title outside dialog box
dims = [1 35]; %Dimension of dialog box
prompt = "Enter Wi Number";
Wis = inputdlg(prompt,chart_title,dims);
Wi = str2double(Wis)
if isempty(Wi)
    fprintf('Code terminated: Missing Input.')
    return
end

num = 2;
for j = 1:num
    if j == 1
        s = 'CTAB';
    else
        s = 'CPCl';
    end
    prompt = sprintf('Exp. No. for %s Wi = %0.1f', s, Wi);
    exp_nums = inputdlg(prompt,chart_title,dims);
    exp_num(j) = str2double(exp_nums)
    if isempty(exp_num)
        fprintf('Code terminated: Missing Input.')
        return
    end
end

for j = 1:num
    if j == 1
        s = 'CTAB';
    else
        s = 'CPCl';
    end
    gss = '4';
    gs(j) = str2double(gss);
end
for j = 1
    if j == 1
        s = 'CTAB';
    else
        s = 'CPCl';
    end
    for k = 1
        prompt = sprintf('Gap Size [mm] #%d for %s',k,s);
        d{1} = [4.11; 2.23; 1.18; 0.7];
        d{2} = [4.11; 2.23; 1.18; 0.7];
        if isempty(d)
            fprintf('Code terminated: Missing Input.')
            return
        end
    end
end

for j = 1:num
    if j == 1
        s = 'CTAB';
    else
        s = 'CPCl';
    end
    for i = 1:gs(j)
        stress_file = fullfile(path_stress, sprintf('C%d-E_%d.xlsx',i,exp_num(j)));
        stress_filename = sprintf('C%d-E_%d.xlsx',i,exp_num(j))

        shear_data = xlsread(stress_file);
        time{j,i}(:,1) = shear_data(:,1);
        shearrate{j,i}(:,1) = shear_data(:,2);
        viscosity{j,i}(:,1) = shear_data(:,3);
        stress{j,i}(:,1) = shear_data(:,6);
        shear_av{j,i} = mean(shearrate{j,i}(round(end*0.7):end));
        strain{j,i} = time{j,i}*shear_av{j,i};
    end
end


if exp_num(1) == 1
    Dss = {'2'};
elseif exp_num(1) == 2
    Dss = {'2'};
elseif exp_num(1) == 3
    Dss = {'3'};
elseif exp_num(1) == 4
    Dss = {'4'};
elseif exp_num(1) == 5
    Dss = {'4'};
elseif exp_num(1) == 6
    Dss = {'4'};
end

for i = 1:sum(gs)
    Ds(i) = str2double(Dss);
    if isempty(Ds)
        fprintf('Code terminated: Missing Input.')
        return
    end
end

N{1} = 'T';
N{2} = 'S1';
if Ds(1) == 2
    N{2} = 'S';
end
N{3} = 'S2';
N{4} = 'S3';

a = 0;
for j = 1:num
    for k = 1:gs(j)
        a = a+1;
        for i = 1:Ds(a)
            delta_file = fullfile(path_stress,sprintf('C%d-E_%d%s_delta_srd.csv',k,exp_num(j),N{i}));
            delta_filename = sprintf('C%d-E_%d%s_delta_srd.csv',k,exp_num(j),N{i})
            delta_data{i} = xlsread(delta_file);
            delta{k,i,j}(:,1) = delta_data{i}(:,1);
            delta{k,i,j}(:,2) = delta_data{i}(:,2);
        end
    end
end
colors_CTAB = {'#210087', '#1f00ff', '#0382f5', '#03e2f5','#b9d3f6'};
colors_CPCl= {'#5C3317', '#CF6112', '#CD853F', '#FFA54F', '#fbd8b6'};
colors = {'#0072BD', '#D95319', '#7E2F8E', '#77AC30', '#4DBEEE', '#A2142F', '#EDB120'};
colorsr = {'#8C0F28', '#BD2B4C', '#EE4770', '#FF71B1', '#FFA8D0', '#FFD4E6'};
colorsb = {'#00517F', '#487BAC', '#8CA3D8', '#CFCCF3', '#E3EBFA', '#F5F8FF'};

%% PCHIP
for j = 1:num
    for i = 1:gs(j)
        st = []; %#ok<NASGU>
        str = []; %#ok<NASGU>
        st = strain{j,i};
        str = stress{j,i};
        [st, sorted_indices] = sort(st);
        str = str(sorted_indices);
        interpolated_func = pchip(st, str);
        strain_c{j,i} = logspace(log10(min(strain{j,i})), log10(max(strain{j,i})), 100);
        stress_cc = interp1(strain{j,i}, stress{j,i}, strain_c{j,i}, 'pchip');
        stress_c{j,i} = stress_cc;

    end
end

%% Plot

close all
for b = 1
    x1_L = 10^-2;
    x1_H = 10^2;
    y1_L = 10^-2;
    y1_H = 2*10^2;

    x2_L = x1_L;
    x2_H = x1_H;
    y2_L = y1_L;
    y2_H = y1_H;

    x3_L = x1_L;
    x3_H = x1_H;
    y3_L = 0;
    y3_H = 3;

    x4_L = x1_L;
    x4_H = x1_H;
    y4_L = y3_L;
    y4_H = y3_H;
end  % Plot Limits

n = 17; % Axis label
p = 17; % Legends
l = 2;

a = 0;
for j = 1
    f1 = figure;
    % Stress CTAB --------------------------------------------------------------------------------------------------------------
    subplot(2,2,1)
    for i = 1:gs(j)
        a = a+1;
        loglog(strain_c{j,i}, stress_c{j,i}, '-', 'color', sprintf('%s', colors_CTAB{i}), 'LineWidth', l)
        hold on
        legend_entries{i} = sprintf('d = %0.2f mm', d{j}(i));
    end
    xlabel('\gamma')
    ylabel('\sigma [Pa]')
    set(gca, 'FontSize', n);
    box on
    xlim([x1_L x1_H])
    ylim([y1_L y1_H])
    legend(legend_entries, 'Location', 'NorthEast', 'FontSize', p);
end

% Stress CPCl--------------------------------------------------------------------------------------------------------------

for j = 2
    subplot(2,2,2)
    for i = 1:gs(j)
        a = a+1;
        loglog(strain_c{j,i}, stress_c{j,i}, '-', 'color', sprintf('%s', colors_CPCl{i}), 'LineWidth', l)
        hold on
        legend_entries{i} = sprintf('d = %0.2f mm', d{j}(i));

    end
    xlabel('\gamma')
    set(gca, 'FontSize', n);
    box on
    xlim([x2_L x2_H])
    ylim([y2_L y2_H])
    legend(legend_entries, 'Location', 'NorthEast', 'FontSize', p);
    yticks([]);
end

% Delta CTAB --------------------------------------------------------------------------------------------------------------

for j = 1
    subplot(2,2,3)
    a = 0;
    for k = 1:gs(j)
        a = a+1;
        for i = 1:Ds(a)
            semilogx(delta{k,i,j}(:,1),  delta{k,i,j}(:,2), '-', 'color', sprintf('%s', colors_CTAB{k}), 'LineWidth', l);
            hold on
        end
        xlabel('\gamma')
        ylabel('\Delta')
        set(gca, 'FontSize', n);
        box on
        xlim([x3_L x3_H])
        ylim([y3_L y3_H])
    end
end
D_N1 = 2*4.11/(92.96/2);
yline(D_N1,'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off')

% Delta CPCl--------------------------------------------------------------------------------------------------------------

for j = 2
    subplot(2,2,4)
    for k = 1:gs(j)
        a = a+1;
        for i = 1:Ds(a)
            semilogx(delta{k,i,j}(:,1),  delta{k,i,j}(:,2), '-', 'color', sprintf('%s', colors_CPCl{k}),'LineWidth', l);
            hold on
        end
        xlabel('\gamma')
        set(gca, 'FontSize', n);
        box on
        xlim([x4_L x4_H])
        ylim([y4_L y4_H])
    end
end
yticks([]);
D_N1 = 2*4.11/(92.96/2);
yline(D_N1,'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off')

set(gcf, 'WindowState', 'maximized');

file_name = sprintf('E-Wi%0.1f.fig', Wi);
fig_folder = fullfile(resultFolder, 'Figs');
if ~exist(fig_folder, 'dir')
    mkdir(fig_folder);
end
fullName = fullfile(fig_folder, file_name);
saveas(f1,fullName)

file_name = sprintf('E-Wi%0.1f.png', Wi);
fullName = fullfile(resultFolder, file_name);
fig = gcf;
saveas(f1,fullName)

%% Save

folder_path = fullfile(resultFolder, 'Workspaces');
if ~exist(folder_path, 'dir')
    mkdir(folder_path);
end
if Wi < 1
    save(fullfile(folder_path, sprintf('E-Wi%0.1f.mat', Wi)));
else
    save(fullfile(folder_path, sprintf('E-Wi%d.mat', Wi)));
end

disp('Done.')
Wi
