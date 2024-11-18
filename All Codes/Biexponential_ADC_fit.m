clear; clc; clf; close all;
format shortg
warning off
options = optimset('Display', 'off');
guessesPerDecade = 4; %Modify for more or less initial guess random combinations

%% Experiment Details
[file, path] = uigetfile('*.xlsx*', 'Select Data File');
filename = fullfile(path, file);
og =  xlsread(filename);

chart_title = 'Input';  %Title outside dialog box
dims = [1 35]; %Dimension of dialog box
numDelta = 1;

prompt = "# of Peaks per BigDelta";
numPeakss = inputdlg(prompt,chart_title,dims);
numPeaks = str2double(numPeakss);
if isempty(numPeaks)
    fprintf('Code terminated: Missing Input.')
    return
end

prompt = "# of Points per Scan";
numPointss = inputdlg(prompt,chart_title,dims);
numPoints = str2double(numPointss);
if isempty(numPoints)
    fprintf('Code terminated: Missing Input.')
    return
end
prompt = sprintf('Arrange %d plots. Enter No. ROWS:',numPeaks);
pnumRowss = inputdlg(prompt,chart_title,dims);
pnumRows = str2double(pnumRowss);
if isempty(pnumRows)
    fprintf('Code terminated: Missing Input.')
    return
end

prompt = sprintf('Arrange %d plots. Enter No. Cols:',numPeaks);
pnumColss = inputdlg(prompt,chart_title,dims);
pnumCols = str2double(pnumColss);
if isempty(pnumCols)
    fprintf('Code terminated: Missing Input.')
    return
end

for i = 1:numPeaks
    prompt = sprintf('ppm for Peak #%.0f', i);
    ppms = inputdlg(prompt,chart_title,dims);
    ppm(i) = str2double(ppms);
end
peakNo = 1:1:numPeaks;
peakNo = peakNo';
ppm = ppm';


%% Organizing Data
data = og(1:numPoints+1, numDelta+2:numDelta+2+(2*numPeaks*numDelta)-1);
[numRows,numCols] = size(data);
for i = 1:numDelta
    Delta(i) = og(1,i);
end
j = 1;
for i = 1:2:numCols
    B{j} = data(:,i);
    B{j}(isnan(B{j})) = [];
    S{j} = data(:,i+1);
    S{j}(isnan(S{j})) = [];
    j = j+1;
end

%% Limit for Guesses

lb = [0.01, 1*10^-9, 1*10^-9]; %(10^-9 mm^2/s = 10^-15 m^2/s)
ub = [1, 0.01, 0.01]; %(10^-3 mm^2/s = 10^-9 m^2/s)
lb_m = 1e-9;  %(10^-9 mm^2/s = 10^-15 m^2/s)
ub_m = 0.1;   %(10^-1 mm^2/s = 10^-7 m^2/s)
lb_b = lb(2);
ub_b = ub(2);

%% Iterating for best fit Mono

k = 1;
for j = 1:numDelta
    for i = 1:numPeaks
        bestResidual_b = Inf;
        bestResidual_m = Inf;
        D_guess_m = [];
        lowerBound_m = lb_m;
        upperBound_m = ub_m;

        logRange_m = log10(upperBound_m) - log10(lowerBound_m);
        numIntervals_m = ceil(logRange_m)*guessesPerDecade;

        for g_m = 1:numIntervals_m
            intervalLowerBound_m = 10^(log10(lowerBound_m) + (g_m -1)*logRange_m/numIntervals_m);

            if g_m == numIntervals_m
                intervalUpperBound_m = min(10^(log10(lowerBound_m) + g_m*logRange_m/numIntervals_m), upperBound_m);
            else
                intervalUpperBound_m = 10^(log10(lowerBound_m) + g_m*logRange_m/numIntervals_m);
            end

            D_guess_m(end+1) = intervalLowerBound_m + (intervalUpperBound_m - intervalLowerBound_m)*rand();

            funm = @(D, x) exp(-(B{k}-B{k}(1)).*D(1));
            [fit_m{k}, residual_m(k,1)] = lsqcurvefit(funm, D_guess_m(g_m), B{k}, S{k}, lb_m, ub_m, options);

            if residual_m(k, 1) < bestResidual_m
                bestResidual_m(k) = residual_m(k, 1);
                bestGuess_m(k) = D_guess_m(g_m);
            end
        end
        k = k + 1;
    end
end

%% Iterating for best fit Biexp

k = 1;
for j = 1:numDelta
    for i = 1:numPeaks
        bestResidual_b = Inf;
        bestResidual_m = Inf;

        D_guess_b2 = [];
        D_guess_b4 = [];
        D_guess_b = [];

        lowerBound_b = lb_b;
        upperBound_b = ub_b;

        logRange_b = log10(upperBound_b) - log10(lowerBound_b);
        numIntervals_b = ceil(logRange_b)*guessesPerDecade;

        for g_b = 1:numIntervals_b
            intervalLowerBound_b = 10^(log10(lowerBound_b) + (g_b - 1)*logRange_b/numIntervals_b);

            if g_b == numIntervals_b
                intervalUpperBound_b = min(10^(log10(lowerBound_b) + g_b*logRange_b/numIntervals_b), upperBound_b);
            else
                intervalUpperBound_b = 10^(log10(lowerBound_b) + g_b*logRange_b/numIntervals_b);
            end

            D_guess_b2(end+1) = intervalLowerBound_b + (intervalUpperBound_b - intervalLowerBound_b)*rand();
            D_guess_b4(end+1) = intervalLowerBound_b + (intervalUpperBound_b - intervalLowerBound_b)*rand();

            D_guess_b{g_b}(1) = lb(1) + (ub(1) - lb(1))*rand();
            D_guess_b{g_b}(2) = D_guess_b2(g_b);
            D_guess_b{g_b}(3) = D_guess_b4(g_b);

            funb = @(D, x) D(1)*exp(-(B{k}-B{k}(1)).*D(2))+(1-D(1))*(exp(-(B{k}-B{k}(1)).*D(3)));
            [fit_b{k}, residual_b(k,1)] = lsqcurvefit(funb, D_guess_b{g_b}, B{k}, S{k}, lb, ub, options);

            if residual_b(k, 1) < bestResidual_b
                bestResidual_b(k) = residual_b(k, 1);
                bestGuess_b{k} = D_guess_b{g_b};
            end
        end
        k = k+1;

    end
end


%% Plotting Best fit & Save Work

close all
p=14;
ADC_m = [];  residual_m = [];  A_b = [];  ADC1_b = [];  B_b = [];  ADC2_b = [];  residual_b = [];
k = 1;
for j = 1:numDelta
    f1 = figure;
    for i = 1:numPeaks
        funm = @(D,x) exp(-(B{k}-B{k}(1))*D(1));
        [best_fit_m{k}, best_residual_m(k,1)] = lsqcurvefit(funm, bestGuess_m(k), B{k}, S{k}, lb_m, ub_m, options);

        funb = @(D,x) D(1)*exp(-(B{k}-B{k}(1)).*D(2))+(1-D(1))*(exp(-(B{k}-B{k}(1)).*D(3)));
        [best_fit_b{k}, best_residual_b(k,1)] = lsqcurvefit(funb, bestGuess_b{k}, B{k}, S{k}, lb, ub, options);

        if numDelta == 1 && numPeaks == 1
        else
            subplot(pnumRows, pnumCols, i)
        end

        hold on
        a1 = semilogy(B{k}, S{k}, 'ko');
        a3 = semilogy(B{k}, funb(best_fit_b{k}), 'r-'); %BIEXPONENTIAL
        a2 = semilogy(B{k}, funm(best_fit_m{k}), 'b-'); y = 1;
        set(gca, 'Yscale', 'log')
        box on
        xlabel('B [s/mm^2]','FontSize',p)
        ylabel('Relative Intensity','FontSize',p)
        legend([a1,a2,a3],{'Data', 'MonoExp.','BiExp.'},'Box','on', 'Location','SouthWest','FontSize',p-3)
        sgtitle(sprintf('\\Delta = %d ms', Delta(j)),'fontsize',p);
        ylim([min(S{k})-0.1 1])
        xlim([min(B{k}),max(B{k})+1000])
        title(sprintf('%0.1f ppm', ppm(i)),'FontSize',p);

        ADC_m(k,1) = best_fit_m{k}/10^6;
        residual_m(k,1) = best_residual_m(k,1);

        if  best_fit_b{k}(2)> best_fit_b{k}(3)
            A_b(k,1) = best_fit_b{k}(1);
            ADC1_b(k,1) = best_fit_b{k}(2)/10^6;
            B_b(k,1) = 1 - best_fit_b{k}(1);
            ADC2_b(k,1) = best_fit_b{k}(3)/10^6;
            residual_b(k,1) = best_residual_m(k,1);
            k = k+1;
        elseif best_fit_b{k}(2)< best_fit_b{k}(3)
            B_b(k,1) = best_fit_b{k}(1);
            ADC2_b(k,1) = best_fit_b{k}(2)/10^6;
            A_b(k,1) = 1 - best_fit_b{k}(1);
            ADC1_b(k,1) = best_fit_b{k}(3)/10^6;
            residual_b(k,1) = best_residual_m(k,1);
            k = k+1;
        end
    end
end

set(gcf, 'WindowState', 'maximized');
% set(gca,'fontsize',20);

[~, name, ~] = fileparts(file);
file_name = sprintf('Signal_Attenuation_%s', name);
f_file_fig = [file_name, '.fig'];
f_name = fullfile(path, f_file_fig);
savefig(f1, f_name);

fileName = sprintf('work_ADC_fit_%s.mat', name);
fullFileName = fullfile(path, fileName);
save(fullFileName);


%% Print Table/Export section

clc; BigDelta = [];
n = numPeaks*numDelta;
list = cell(n, 1);
for i = 1:n
    BigDelta(i, 1) = Delta(ceil(i/numPeaks));
end

fprintf('\n \n')
disp('                          A*exp(-b*D)')
T = table(ppm, BigDelta, ADC_m, residual_m, A_b, ADC1_b, B_b, ADC2_b, residual_b);
fprintf('\n \n')
disp(T)

T_file = sprintf('ADC_Table_%s',name);
T_file_csv = [T_file, '.csv'];
T_name = fullfile(path, T_file_csv);
delete(T_name)
writetable(T, T_name)
winopen(T_name) % Opens Processed Table with ADCs
winopen(filename)% Re-opens Raw data file
file
