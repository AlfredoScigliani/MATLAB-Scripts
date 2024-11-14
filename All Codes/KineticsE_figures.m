clear; clc; clf; close all
format shortg

%% Input Data

[file, path] = uigetfile('E_data.xlsx', 'Select File');
filename = fullfile(path, file);
data = xlsread(filename);
file
i = 1;

% Cell(:,1) = data(:,i); i = i+1;
q = data(:,i); i = i+1;
System = data(:,i); i = i+1;
Wi = data(:,i); i = i+1;
shear_rate = data(:,i); i = i+1;
t_smax = data(:,i); i = i+1;
s_max = data(:,i); i = i+1;
y_max = data(:,i); i = i+1;
s_ss = data(:,i); i = i+1;
TFR = data(:,i); i = i+1;
SB = data(:,i); i = i+1;
Vi_2k = data(:,i); i = i+1;
Vo_2k = data(:,i); i = i+1;
Vi_5k = data(:,i); i = i+1;
Vo_5k = data(:,i); i = i+1;
Vi_10k = data(:,i); i = i+1;
Vo_10k = data(:,i); i = i+1;
D_max = data(:,i); i = i+1;
D_2k = data(:,i); i = i+1;
D_2k_std = data(:,i); i = i+1;
D_5k = data(:,i); i = i+1;
D_5k_std = data(:,i); i = i+1;
D_10k = data(:,i); i = i+1;
D_10k_stq = data(:,i); i = i+1;
D_SS = data(:,i); i = i+1;
D_SS_std = data(:,i); i = i+1;
Vi_SS = data(:,i); i = i+1;
Vo_SS = data(:,i); i = i+1;
LR = data(:,i); i = i+1;


%% TFR and D

close all; clc;
colors_hex = {'#210087','#5C3317', '#1f00ff','#CF6112', '#0382f5','#CD853F', '#03e2f5','#FFA54F','#b9d3f6', '#fbd8b6'};

%olors_hex = {'#0072BD','#D95319','#EDB120','#7E2F8E','#77AC30','#4DBEEE','#A2142F', '#9933ff', '#5A9BD4', '#FF6600'};
colors = cellfun(@(x) sscanf(x(2:end)','%2x%2x%2x',[1 3])/255, colors_hex, 'UniformOutput', false);
shape = {'o', 's', 'o', 's', 'o', 's', 'o', 's', 'o', 's'};

p = 20; % labels
t = 20; %ticks
n = 14; % Font size for legends
m = 10;  % Marker size
b = 1.5;  % Thickness of plot edges
L = 0.92; % x position for letter a) b)...
lw = 3;

figure
j = 1;
subplot(2,2,1); % TFR CTAB
hold on;
for i = 1:2:8
    a1(j) = semilogx(Wi((i-1)*6+1:(i-1)*6+1+5)', TFR((i-1)*6+1:(i-1)*6+1+5)', sprintf('k%s-', shape{i}), 'Color', colors{i}, 'MarkerFaceColor', colors{i}, 'MarkerEdgeColor', 'k', 'MarkerSize', m);
    j = j+1;
end
xlabel('Wi', 'FontSize', p, 'FontName', 'Times');
legend('d = 4.11 mm', 'd = 2.23 mm', 'd = 1.18 mm', 'd = 0.70 mm', 'Location', 'SouthWest', 'FontSize', n, 'FontName', 'Times');
ylim([-0.4 0.3]);
yline(0, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
text(L, 0.99, '(a)', 'Units', 'normalized', 'VerticalAlignment', 'top', 'FontWeight', 'bold', 'FontSize', t, 'FontName', 'Times');
ax = gca;
ax.FontSize = t;
ax.FontName = 'Times';
ax.LineWidth = b;
set(ax, 'Box', 'on', 'TickLength', [0.02, 0.02]);
ylabel('v/U_i', 'FontName', 'Times', 'FontSize', p);
legendHandle = legend;
set(legendHandle, 'LineWidth', 0.5);
set(ax, 'XScale', 'log');
xlim([10^-1 3*10^3])
set(gca, 'XTick', [10^-1 10^0 10^1 10^2 10^3]);
set(gca, 'XTickLabel', {'10^{-1}','10^{0}', '10^{1}','10^{2}','10^{3}'});
% -----------------------------------------------------------------------------------------------------------------------------------------------
subplot(2,2,2); % TFR CPCL
hold on;
for i = 2:2:8
    a1(j) = semilogx(Wi((i-1)*6+1:(i-1)*6+1+5)', TFR((i-1)*6+1:(i-1)*6+1+5)', sprintf('k%s-', shape{i}), 'Color', colors{i}, 'MarkerFaceColor', colors{i}, 'MarkerEdgeColor', 'k', 'MarkerSize', m);
    j = j+1;
end
xlabel('Wi', 'FontSize', p, 'FontName', 'Times');
legend('d = 4.11 mm', 'd = 2.23 mm', 'd = 1.18 mm', 'd = 0.70 mm', 'Location', 'SouthWest', 'FontSize', n, 'FontName', 'Times');
ylim([-0.4 0.3]);
yline(0, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
text(L, 0.99, '(b)', 'Units', 'normalized', 'VerticalAlignment', 'top', 'FontWeight', 'bold', 'FontSize', t, 'FontName', 'Times');
ax = gca;
ax.FontSize = t;
ax.FontName = 'Times';
ax.LineWidth = b;
set(ax, 'Box', 'on', 'TickLength', [0.02, 0.02]);
ylabel('v/U_i', 'FontName', 'Times', 'FontSize', p);
legendHandle = legend;
set(legendHandle, 'LineWidth', 0.5);
set(ax, 'XScale', 'log');
xlim([10^-1 3*10^3])
set(gca, 'XTick', [10^-1 10^0 10^1 10^2 10^3]);
set(gca, 'XTickLabel', {'10^{-1}','10^{0}', '10^{1}','10^{2}','10^{3}'});

% -----------------------------------------------------------------------------------------------------------------------------------------------
subplot(2,2,3); % Delta SS CTAB
yyaxis left
for i = 1:2:8
    eb = errorbar(Wi((i-1)*6+1:(i-1)*6+1+5)', D_SS((i-1)*6+1:(i-1)*6+1+5)', D_SS_std((i-1)*6+1:(i-1)*6+1+5)', 'vertical', 'LineStyle', 'none');
    set(eb, 'color', colors{i}, 'LineWidth', 1);
    a(i) = loglog(Wi((i-1)*6+1:(i-1)*6+1+5)', D_SS((i-1)*6+1:(i-1)*6+1+5)', sprintf('%s-', shape{i})', 'Color', colors{i}, 'MarkerFaceColor', colors{i}, 'MarkerSize', m, 'LineWidth', lw);
    hold on;
end
set(gca, 'FontSize', t, 'FontName', 'Times', 'LineWidth', b, 'YColor', 'k');
ylim([0.1 10]);
xlim([10^-1 3*10^3]);
xlabel('Wi', 'FontSize', p, 'FontName', 'Times');
ylabel('\Delta_S_S', 'FontSize', t, 'FontName', 'Times', 'Color', 'k');
% set(gca, 'YTick', [0.1 1 10]);
% set(gca, 'YTickLabel', {'0.1','1','10'});
D_N1 = 2*4.11/(92.96/2);
D_N2 = 2*2.23/(50.55/2);
D_N3 = 2*1.18/(26.72/2);
D_N4 = 2*0.70/(15.85/2);
yline(D_N1,'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off')
yline(D_N2,'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off')
yline(D_N3,'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off')
yline(D_N4,'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off')

yyaxis right % Delta max CTAB
for i = 1:2:8
    semilogx(Wi((i-1)*6+1:(i-1)*6+1+5)', D_max((i-1)*6+1:(i-1)*6+1+5)', sprintf('%s-', shape{i}), 'Color', colors{i}, 'MarkerEdgeColor', colors{i}, 'MarkerSize', m, 'LineWidth', b, 'LineWidth', lw);
    % semilogx(Wi((i-1)*6+1:(i-1)*6+1+5)', D_max((i-1)*6+1:(i-1)*6+1+5)','-', 'Color', colors{i});
end

yyaxis right
set(gca, 'FontSize', t, 'FontName', 'Times', 'LineWidth', b, 'YColor', 'k');
ylim([-2 30]);
ylabel('\Delta_m_a_x', 'FontSize', t, 'FontName', 'Times', 'Color', 'k');
xlim([10^-1 3*10^3])
% set(gca, 'XScale', 'log', 'YScale', 'log');
set(gca, 'XTick', [10^-1 10^0 10^1 10^2 10^3]);
set(gca, 'XTickLabel', {'10^{-1}','10^{0}', '10^{1}','10^{2}','10^{3}'});
% set(gca, 'YTick', [0 10 20]);
% set(gca, 'YTickLabel', {'0','10', '20'});
text(L-0.01, 0.98, '(c)', 'Units', 'normalized', 'VerticalAlignment', 'top', 'FontWeight', 'bold', 'FontSize', t, 'FontName', 'Times');
annotation('arrow', [0.39 0.41], [0.24 0.24], 'Units', 'normalized', 'LineWidth', 2);
annotation('arrow', [0.16 0.14], [0.31 0.31], 'Units', 'normalized', 'LineWidth', 2);
% -----------------------------------------------------------------------------------------------------------------------------------------------
subplot(2,2,4); % Delta SS CPCL
yyaxis left
for i = 2:2:8
    eb = errorbar(Wi((i-1)*6+1:(i-1)*6+1+5)', D_SS((i-1)*6+1:(i-1)*6+1+5)', D_SS_std((i-1)*6+1:(i-1)*6+1+5)', 'vertical', 'LineStyle', 'none');
    set(eb, 'color', colors{i}, 'LineWidth', 1);
    a(i) = loglog(Wi((i-1)*6+1:(i-1)*6+1+5)', D_SS((i-1)*6+1:(i-1)*6+1+5)', sprintf('%s-', shape{i})', 'Color', colors{i}, 'MarkerFaceColor', colors{i}, 'MarkerSize', m, 'LineWidth', lw);
    hold on;
end
set(gca, 'FontSize', t, 'FontName', 'Times', 'LineWidth', b, 'YColor', 'k');
ylim([0.1 10]);
xlim([10^-1 3*10^3]);
xlabel('Wi', 'FontSize', p, 'FontName', 'Times');
ylabel('\Delta_S_S', 'FontSize', t, 'FontName', 'Times', 'Color', 'k');
% set(gca, 'YTick', [0.1 1 10]);
% set(gca, 'YTickLabel', {'0.1','1','10'});
D_N1 = 2*4.11/(92.96/2);
D_N2 = 2*2.23/(50.55/2);
D_N3 = 2*1.18/(26.72/2);
D_N4 = 2*0.70/(15.85/2);
yline(D_N1,'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off')
yline(D_N2,'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off')
yline(D_N3,'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off')
yline(D_N4,'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off')

yyaxis right % Delta max CPCL
for i = 2:2:8
    semilogx(Wi((i-1)*6+1:(i-1)*6+1+5)', D_max((i-1)*6+1:(i-1)*6+1+5)', sprintf('%s-', shape{i}), 'Color', colors{i}, 'MarkerEdgeColor', colors{i}, 'MarkerSize', m, 'LineWidth', b,'LineWidth', lw);
    % semilogx(Wi((i-1)*6+1:(i-1)*6+1+5)', D_max((i-1)*6+1:(i-1)*6+1+5)','-', 'Color', colors{i});
end

yyaxis right
set(gca, 'FontSize', t, 'FontName', 'Times', 'LineWidth', b, 'YColor', 'k');
ylim([-2 30]);
ylabel('\Delta_m_a_x', 'FontSize', t, 'FontName', 'Times', 'Color', 'k');
xlim([10^-1 3*10^3])
% set(gca, 'XScale', 'log', 'YScale', 'log');
set(gca, 'XTick', [10^-1 10^0 10^1 10^2 10^3]);
set(gca, 'XTickLabel', {'10^{-1}','10^{0}', '10^{1}','10^{2}','10^{3}'});
% set(gca, 'YTick', [0 10 20]);
% set(gca, 'YTickLabel', {'0','10', '20'});
text(L-0.01, 0.98, '(d)', 'Units', 'normalized', 'VerticalAlignment', 'top', 'FontWeight', 'bold', 'FontSize', t, 'FontName', 'Times');
annotation('arrow', [0.60 0.58], [0.32 0.32], 'Units', 'normalized', 'LineWidth', 2);
annotation('arrow', [0.83 0.85], [0.23 0.23], 'Units', 'normalized', 'LineWidth', 2);
set(gcf, 'WindowState', 'maximized');

% print(gcf, 'KineticsE1.png', '-dpng', '-r600');

%% SB and V

figure
subplot(2,2,1); % SB CTAB
for i = 1:2:8
    plot(Wi((i-1)*6+2:(i-1)*6+1+5)', SB((i-1)*6+2:(i-1)*6+1+5)', sprintf('k%s', shape{i})', 'Color', colors{i}, 'MarkerFaceColor', colors{i}, 'MarkerEdgeColor', 'k', 'MarkerSize', m);
    hold on
end
plot(Wi((i-1)*6+1:(i-1)*6+1+5)', LR((i-1)*6+1:(i-1)*6+1+5)', 'k--', 'LineWidth', 1.5);

ylim([-0.1 1]);
text(L, 0.98, '(e)', 'Units', 'normalized', 'VerticalAlignment', 'top', 'FontWeight', 'bold', 'FontSize', t, 'FontName', 'Times');
set(gca, 'FontSize', t, 'FontName', 'Times', 'LineWidth', b);
%     legend('q = 0.024 (CTAB)', 'q = 0.044 (CTAB)', 'q = 0.085 (CTAB)', 'q = 0.147 (CTAB)', 'q = 0.171 (CTAB)','q = 0.024 (CPCl)', 'q = 0.044 (CPCl)', 'q = 0.085 (CPCl)', 'q = 0.147 (CPCl)', 'q = 0.171 (CPCl)', 'Location', 'NorthWest', 'FontSize', n, 'FontName', 'Times', 'NumColumns', 2);
%     legendHandle = legend;
%     set(legendHandle, 'LineWidth', 0.5);
xlabel('Wi', 'FontName', 'Times', 'FontSize', t);
ylabel('$\frac{\alpha_H}{d}$', 'Interpreter', 'latex', 'FontSize', p+5, 'Rotation', 0);
xlim([-20 450])
% set(gca, 'XTick', [10^-1 10^0 10^1 10^2 10^3]);
% set(gca, 'XTickLabel', {'10^{-1}','10^{0}', '10^{1}','10^{2}','10^{3}'});
% -----------------------------------------------------------------------------------------------------------------------------------------------
subplot(2,2,2); % SB CPCL
for i = 2:2:8
    plot(Wi((i-1)*6+2:(i-1)*6+1+5)', SB((i-1)*6+2:(i-1)*6+1+5)', sprintf('k%s', shape{i})', 'Color', colors{i}, 'MarkerFaceColor', colors{i}, 'MarkerEdgeColor', 'k', 'MarkerSize', m);
    hold on
end
plot(Wi((i-1)*6+1:(i-1)*6+1+5)', LR((i-1)*6+1:(i-1)*6+1+5)', 'k--', 'LineWidth', 1.5);

ylim([-0.1 1]);
text(L, 0.98, '(f)', 'Units', 'normalized', 'VerticalAlignment', 'top', 'FontWeight', 'bold', 'FontSize', t, 'FontName', 'Times');
set(gca, 'FontSize', t, 'FontName', 'Times', 'LineWidth', b);
%     legend('q = 0.024 (CTAB)', 'q = 0.044 (CTAB)', 'q = 0.085 (CTAB)', 'q = 0.147 (CTAB)', 'q = 0.171 (CTAB)','q = 0.024 (CPCl)', 'q = 0.044 (CPCl)', 'q = 0.085 (CPCl)', 'q = 0.147 (CPCl)', 'q = 0.171 (CPCl)', 'Location', 'NorthWest', 'FontSize', n, 'FontName', 'Times', 'NumColumns', 2);
%     legendHandle = legend;
%     set(legendHandle, 'LineWidth', 0.5);
xlabel('Wi', 'FontName', 'Times', 'FontSize', t);
ylabel('$\frac{\alpha_H}{d}$', 'Interpreter', 'latex', 'FontSize', p+5, 'Rotation', 0);
xlim([-20 450])
% set(gca, 'XTick', [10^-1 10^0 10^1 10^2 10^3]);
% set(gca, 'XTickLabel', {'10^{-1}','10^{0}', '10^{1}','10^{2}','10^{3}'});
% -----------------------------------------------------------------------------------------------------------------------------------------------
subplot(2,2,3) % Vi CTAB
yyaxis left
for i = 1:2:8
    semilogx(Wi((i-1)*6+1:(i-1)*6+1+5)', Vi_SS((i-1)*6+1:(i-1)*6+1+5)', sprintf('k%s-', shape{i})', 'Color', colors{i}, 'MarkerFaceColor', colors{i}, 'MarkerEdgeColor', 'k', 'MarkerSize', m);
    hold on;
end
set(gca, 'FontSize', t, 'FontName', 'Times', 'LineWidth', b, 'YColor', 'k');
ylabel('v_i/U_i', 'FontName', 'Times', 'Color', 'k', 'FontSize', p);
xlabel('Wi', 'FontSize', t, 'FontName', 'Times');

ax = gca;
ax.YColor = 'k';
ylim([-1 1.05])
set(gca, 'YTick', [ 0 0.5 1]);
set(gca, 'YTickLabel', {'0','0.5','1'});

yyaxis right % Vo CTAB
for i = 1:2:8
    semilogx(Wi((i-1)*6+1:(i-1)*6+1+5)', Vo_SS((i-1)*6+1:(i-1)*6+1+5)', sprintf('k%s-', shape{i})', 'Color', colors{i}, 'MarkerFaceColor', colors{i}, 'MarkerEdgeColor', 'k', 'MarkerSize', m);
end
% ylabel('v_o_u_t_e_r/U_i', 'FontSize', p, 'FontName', 'Times', 'Color', 'k');
ax.YColor = 'k';
ylim([-0.05 2])
xlim([10^-1 3*10^3])
set(gca, 'XTick', [10^-1 10^0 10^1 10^2 10^3]);
set(gca, 'XTickLabel', {'10^{-1}','10^{0}', '10^{1}','10^{2}','10^{3}'});
text(L, 0.98, '(g)', 'Units', 'normalized', 'VerticalAlignment', 'top', 'FontWeight', 'bold', 'FontSize', t, 'FontName', 'Times');
annotation('arrow', [0.16 0.14], [0.43 0.43], 'Units', 'normalized', 'LineWidth', 2);
annotation('arrow', [0.4 0.42], [0.18 0.18], 'Units', 'normalized', 'LineWidth', 2);
ylabel('v_o/U_i', 'FontName', 'Times', 'Color', 'k', 'FontSize', p);
set(gca, 'YTick', [0 0.5 1]);
set(gca, 'YTickLabel', {'0','0.5', '1'});
% -----------------------------------------------------------------------------------------------------------------------------------------------
subplot(2,2,4) % Vi CPCL
yyaxis left
for i = 2:2:8
    semilogx(Wi((i-1)*6+1:(i-1)*6+1+5)', Vi_SS((i-1)*6+1:(i-1)*6+1+5)', sprintf('k%s-', shape{i})', 'Color', colors{i}, 'MarkerFaceColor', colors{i}, 'MarkerEdgeColor', 'k', 'MarkerSize', m);
    hold on;
end
set(gca, 'FontSize', t, 'FontName', 'Times', 'LineWidth', b, 'YColor', 'k');
ylabel('v_i/U_i', 'FontName', 'Times', 'Color', 'k', 'FontSize', p);
xlabel('Wi', 'FontSize', t, 'FontName', 'Times');

ax = gca;
ax.YColor = 'k';
ylim([-1 1.05])
set(gca, 'YTick', [ 0 0.5 1]);
set(gca, 'YTickLabel', {'0','0.5','1'});

yyaxis right % Vo CPCL
for i = 2:2:8
    semilogx(Wi((i-1)*6+1:(i-1)*6+1+5)', Vo_SS((i-1)*6+1:(i-1)*6+1+5)', sprintf('k%s-', shape{i})', 'Color', colors{i}, 'MarkerFaceColor', colors{i}, 'MarkerEdgeColor', 'k', 'MarkerSize', m);
end
% ylabel('v_o_u_t_e_r/U_i', 'FontSize', p, 'FontName', 'Times', 'Color', 'k');
ax.YColor = 'k';
ylim([-0.05 2])
xlim([10^-1 3*10^3])
set(gca, 'XTick', [10^-1 10^0 10^1 10^2 10^3]);
set(gca, 'XTickLabel', {'10^{-1}','10^{0}', '10^{1}','10^{2}','10^{3}'});
text(L, 0.98, '(h)', 'Units', 'normalized', 'VerticalAlignment', 'top', 'FontWeight', 'bold', 'FontSize', t, 'FontName', 'Times');
annotation('arrow', [0.16 0.14], [0.43 0.43], 'Units', 'normalized', 'LineWidth', 2);
annotation('arrow', [0.4 0.42], [0.18 0.18], 'Units', 'normalized', 'LineWidth', 2);
ylabel('v_o/U_i', 'FontName', 'Times', 'Color', 'k', 'FontSize', p);

set(gca, 'YTick', [0 0.5 1]);
set(gca, 'YTickLabel', {'0','0.5', '1'});
set(gcf, 'WindowState', 'maximized');

set(gcf, 'WindowState', 'maximized');

% print(gcf, 'KineticsE2.png', '-dpng', '-r600');

