clear; clc; clf; close all
format shortg
initial_params = [1.43 41.70 0.11 1.92 0.44 16.9 0.12 0.34 3.84 0.005 0.01 0.02]; % INITIAL GUESSES [G1, L1, G2, L2, ...]
max_modes = 6;
start = 0.1;% FREQ TO START
threshold = 100; % FREQ TO STOP FITTING

%% File Input

[file, path] = uigetfile('*SAOS.csv', 'Select SAOS Data File');
stress_file = fullfile(path, file);
data = xlsread(stress_file);
file

t_all = data(:,1); % Time [s]
w_all = data(:,2); % Frequency [rad/s]
storage_all = data(:,3);
loss_all = data(:,4);
torque_all = data(:,5);

differences = abs(w_all-start);
[~, index_start] = min(differences);

differences = abs(w_all-threshold);
[~, index] = min(differences);

last = length(w_all) - index_start;
first = index - 1;

t = t_all(index:index_start);
w = w_all(index:index_start);
storage = storage_all(index:index_start);
loss = loss_all(index:index_start);
torque = torque_all(index:index_start);

modes_list = 1:max_modes;
num_modes = length(modes_list);
residuals_all = zeros(num_modes, 1);
S_model_opt_all = cell(num_modes, 1);
L_model_opt_all = cell(num_modes, 1);
S_model_plot_all = cell(num_modes, 1);
L_model_plot_all = cell(num_modes, 1);
G_values_all = cell(num_modes, 1);
L_values_all = cell(num_modes, 1);


%% Bounds of Iteration

G_low = 0;
L_low = 0;

G_up = 1;
L_up = 500;
% To activate upper bound, add it in lsqcurvefit line ub
fprintf('Finding optimal parameters... \n \n')


%% Iterate over different modes

for idx = 1:num_modes
    modes = modes_list(idx);

    if modes == 1
        model_fun = @(params, w) [params(1)*((w*params(2)).^2)./(1+(w*params(2)).^2); params(1)*(w*params(2))./(1+(w*params(2)).^2)];

    elseif modes == 2
        model_fun = @(params, w) [params(1)*((w*params(2)).^2)./(1+(w*params(2)).^2) + params(3)*((w*params(4)).^2)./(1+(w*params(4)).^2); params(1)*(w*params(2))./(1+(w*params(2)).^2) + params(3)*(w*params(4))./(1+(w*params(4)).^2)];

    elseif modes == 3
        model_fun = @(params, w) [params(1)*((w*params(2)).^2)./(1+(w*params(2)).^2) + params(3)*((w*params(4)).^2)./(1+(w*params(4)).^2) + params(5)*((w*params(6)).^2)./(1+(w*params(6)).^2); params(1)*(w*params(2))./(1+(w*params(2)).^2) + params(3)*(w*params(4))./(1+(w*params(4)).^2) + params(5)*(w*params(6))./(1+(w*params(6)).^2)];

    elseif modes == 4
        model_fun = @(params, w) [params(1)*((w*params(2)).^2)./(1+(w*params(2)).^2) + params(3)*((w*params(4)).^2)./(1+(w*params(4)).^2) + params(5)*((w*params(6)).^2)./(1+(w*params(6)).^2) + params(7)*((w*params(8)).^2)./(1+(w*params(8)).^2); params(1)*(w*params(2))./(1+(w*params(2)).^2) + params(3)*(w*params(4))./(1+(w*params(4)).^2) + params(5)*(w*params(6))./(1+(w*params(6)).^2) + params(7)*(w*params(8))./(1+(w*params(8)).^2)];

    elseif modes == 5
        model_fun = @(params, w) [params(1)*((w*params(2)).^2)./(1+(w*params(2)).^2) + params(3)*((w*params(4)).^2)./(1+(w*params(4)).^2) + params(5)*((w*params(6)).^2)./(1+(w*params(6)).^2) + params(7)*((w*params(8)).^2)./(1+(w*params(8)).^2) + params(9)*((w*params(10)).^2)./(1+(w*params(10)).^2); params(1)*(w*params(2))./(1+(w*params(2)).^2) + params(3)*(w*params(4))./(1+(w*params(4)).^2) + params(5)*(w*params(6))./(1+(w*params(6)).^2) + params(7)*(w*params(8))./(1+(w*params(8)).^2) + params(9)*(w*params(10))./(1+(w*params(10)).^2)];

    elseif modes == 6
        model_fun = @(params, w) [params(1)*((w*params(2)).^2)./(1+(w*params(2)).^2) + params(3)*((w*params(4)).^2)./(1+(w*params(4)).^2) + params(5)*((w*params(6)).^2)./(1+(w*params(6)).^2) + params(7)*((w*params(8)).^2)./(1+(w*params(8)).^2) + params(9)*((w*params(10)).^2)./(1+(w*params(10)).^2) + params(11)*((w*params(12)).^2)./(1+(w*params(12)).^2); params(1)*(w*params(2))./(1+(w*params(2)).^2) + params(3)*(w*params(4))./(1+(w*params(4)).^2) + params(5)*(w*params(6))./(1+(w*params(6)).^2) + params(7)*(w*params(8))./(1+(w*params(8)).^2) + params(9)*(w*params(10))./(1+(w*params(10)).^2) + params(11)*(w*params(12))./(1+(w*params(12)).^2)];
    end

    q = sprintf('Current mode = %d of %d', idx, num_modes);
    qs = sprintf('done.');
    disp(q)

    options = optimoptions('lsqcurvefit', 'Display', 'final', 'TolFun', 1e-10, 'TolX', 1e-10, 'MaxIterations', 50000, 'MaxFunctionEvaluations', 50000);
    if idx == num_modes
        options = optimoptions('lsqcurvefit', 'Display', 'final', 'TolFun', 1e-10, 'TolX', 1e-10, 'MaxIterations', 50000, 'MaxFunctionEvaluations', 50000);
    end

    lb = zeros(size(initial_params));
    for i = 1:2:length(initial_params)
        lb(i) = lb(i) + G_low; %G LOWER
    end
    for i = 2:2:length(initial_params)
        lb(i) = lb(i) + L_low; %L LOWER
    end

    ub = zeros(size(initial_params));
    for i = 1:2:length(initial_params)
        ub(i) = ub(i) + G_up; %G UPPER
    end
    for i = 2:2:length(initial_params)
        ub(i) = ub(i) + L_up; %L UPPER
    end
    [optimal_params, ~, residual, ~, ~, ~, jacobian] = lsqcurvefit(@(params, w) model_fun(params, w), initial_params, w, [storage; loss], lb, [], options);

    if modes == 1
        G_values = optimal_params(1);
        L_values = optimal_params(2);

    elseif modes == 2
        G_values = [optimal_params(1), optimal_params(3)];
        L_values = [optimal_params(2), optimal_params(4)];

    elseif modes == 3
        G_values = [optimal_params(1), optimal_params(3), optimal_params(5)];
        L_values = [optimal_params(2), optimal_params(4), optimal_params(6)];

    elseif modes == 4
        G_values = [optimal_params(1), optimal_params(3), optimal_params(5), optimal_params(7)];
        L_values = [optimal_params(2), optimal_params(4), optimal_params(6), optimal_params(8)];

    elseif modes == 5
        G_values = [optimal_params(1), optimal_params(3), optimal_params(5), optimal_params(7), optimal_params(9)];
        L_values = [optimal_params(2), optimal_params(4), optimal_params(6), optimal_params(8), optimal_params(10)];

    elseif modes == 6
        G_values = [optimal_params(1), optimal_params(3), optimal_params(5), optimal_params(7), optimal_params(9), optimal_params(11)];
        L_values = [optimal_params(2), optimal_params(4), optimal_params(6), optimal_params(8), optimal_params(10), optimal_params(12)];
    end
    w_mod = logspace(log10(min(w)), log10(max(w)), 10000)';
    S_model_opt = model_fun(optimal_params, w);
    L_model_opt = S_model_opt((end/2)+1:end);
    S_model_opt = S_model_opt(1:end/2);

    S_model_plot = model_fun(optimal_params, w_mod);
    L_model_plot = S_model_plot((end/2)+1:end);
    S_model_plot = S_model_plot(1:end/2);

    residuals_all(idx) = sum((log10(storage)-log10(S_model_opt)).^2 + (log10(loss)-log10(L_model_opt)).^2);
    S_model_opt_all{idx} = S_model_opt;
    L_model_opt_all{idx} = L_model_opt;
    S_model_plot_all{idx} = S_model_plot;
    L_model_plot_all{idx} = L_model_plot;
    G_values_all{idx} = G_values;
    L_values_all{idx} = L_values;
    disp(qs)
end

modes_table = table(modes_list', residuals_all, 'VariableNames', {'Modes', 'Sum of Log Squared Errors'});
disp(modes_table);

for idx = 1:num_modes
    fprintf('\n');
    disp(['Parameters for No. Modes = ', num2str(modes_list(idx)), ':']);
    results_table = table(G_values_all{idx}', L_values_all{idx}', 'VariableNames', {'G [Pa]', 'L [s]'});
    disp(results_table);
end


%% Plot Results

close all
f1 = figure;
for idx = 1:num_modes
    subplot(2,3,idx);
    loglog(w_all, storage_all, 'ko', 'MarkerFaceColor', 'b')
    hold on
    loglog(w_all, loss_all, 'k^', 'MarkerFaceColor', 'r')
    loglog(w_mod, S_model_plot_all{idx}, 'k.', 'LineWidth', 1)
    loglog(w_mod, L_model_plot_all{idx}, 'k.', 'LineWidth', 1)
    title(['SAOS Data Fitting - Modes = ', num2str(modes_list(idx))])
    xlabel('Angular Frequency (rad/s)')
    ylabel('G'', G'''' (Pa)')
    grid on
    ylim([0.1 100])
    set(gcf, 'WindowState', 'maximized');
end
file

fullName = 'SAOS_sample';
saveas(f1,fullName,'png');


%% Export Selected

chart_title = 'Input';
dims = [1 35];
prompt = "Enter Mode to Export, [0] = Stop";
mode_exports = inputdlg(prompt,chart_title,dims);
mode_export = str2double(mode_exports);
if mode_export == 0
    disp('Run Section to Export Results.')
    return
elseif isempty(mode_export) == 1
    disp('Run Section to Export Results.')
    return
else
    S_model = S_model_opt_all{mode_export};
    L_model = L_model_opt_all{mode_export};

    maxLength = max([length(w_all), length(storage_all), length(loss_all), length(L_model), length(S_model)]);

    w_all(length(w_all)+1:maxLength) = 0;
    storage_all(length(storage_all)+1:maxLength) = 0;
    loss_all(length(loss_all)+1:maxLength) = 0;
    S_model = [zeros(first, 1); S_model; zeros(last, 1)];
    L_model = [zeros(first, 1); L_model; zeros(last, 1)];
    titles = {'Angular Frequency [rad/s]', 'G''', 'G''''', 'G'' model', 'G'''' model'};
    dataMatrix = [w_all, storage_all, loss_all, S_model, L_model];
    dataMatrixWithTitles = [titles; num2cell(dataMatrix)];
    name = sprintf('%sSAOS_results_%s.xlsx', path, file);

    if exist(name, 'file') == 2
        delete(name)
    end

    writecell(dataMatrixWithTitles, name);

    G_ex = G_values_all{mode_export}';
    L_ex = L_values_all{mode_export}';
    res_ex = residuals_all(mode_export);

    writecell({'G0 [Pa]'}, name, 'Range', 'G1');
    writematrix(G_ex, name, 'Range', 'G2');

    writecell({'L [s]'}, name, 'Range', 'H1');
    writematrix(L_ex, name, 'Range', 'H2');

    writecell({'Sum log error^2'}, name, 'Range', 'I1');
    writematrix(res_ex, name, 'Range', 'I2');

    winopen(name)
    disp('Results Exported')
end


