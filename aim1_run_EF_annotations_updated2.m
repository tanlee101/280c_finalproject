% run_EF_from_annotations_csv.m
clear; clc; close all;


T = readtable('annotations.csv'); %load csv

T = T(strcmp(T.view, 'SAX'), :); % Keep only SAX rows

%%Setup Variables
retention_percent = [100 75 50 25]; % define slice retention levels
patientIDs = unique(T.patient_id);

optStrategyNames = {'Uniform','MiddleWeighted','EdgePreserving','BasalWeighted'};
regionalNames = {'Uniform_75','Apical_75','Middle_75','Basal_75'};

EF_all = nan(length(patientIDs), length(retention_percent));
EF_opt = nan(length(patientIDs), length(retention_percent), length(optStrategyNames));
EF_regional = nan(length(patientIDs), length(regionalNames));

EDV_full_all = nan(length(patientIDs), 1); % save full-volume EDV, ESV, disease label, and # of usable slices
ESV_full_all = nan(length(patientIDs), 1);
disease_all = strings(length(patientIDs), 1);
n_valid_slices_all = nan(length(patientIDs), 1);

%%Loop through patients
for p = 1:length(patientIDs) %process one patient at a time

    patientID = patientIDs(p);
    patientData = T(T.patient_id == patientID, :); % extract single patient's rows and sort slices in anatomical order
    patientData = sortrows(patientData, 'z_slice');

    ED_area = patientData.lv_area_ED_mm2;
    ES_area = patientData.lv_area_ES_mm2;

    valid_idx = ~isnan(ED_area) & ~isnan(ES_area);  %remove slices where ED or ES area is missing
    patientData = patientData(valid_idx, :);
    ED_area = ED_area(valid_idx);
    ES_area = ES_area(valid_idx);

    dz = mean(patientData.z_spacing_mm, 'omitnan');  %mean slice spacing, convert slice areas into volumes
    disease_all(p) = string(patientData.disease{1});

    nSlices = length(ED_area);
    n_valid_slices_all(p) = nSlices;

    %% Percent-based reduced-slice EF: uniform baseline
    for k = 1:length(retention_percent) %loop through different slice retentions

        percent = retention_percent(k);
        nKeep = max(round(nSlices * percent / 100), 1);

        slice_idx = unique(round(linspace(1, nSlices, nKeep))); %choose evenly spaced slices
        spacing_factor = nSlices / length(slice_idx);

        EDV = sum(ED_area(slice_idx), 'omitnan') * dz * spacing_factor / 1000; %estimate EDV and ESV and convert to mL
        ESV = sum(ES_area(slice_idx), 'omitnan') * dz * spacing_factor / 1000;

        EF_all(p,k) = (EDV - ESV) / EDV * 100; %EF equation

        if percent == 100
            EDV_full_all(p) = EDV;
            ESV_full_all(p) = ESV;
        end
    end

    %% Optimization experiment: clinically feasible sampling protocols
    for k = 1:length(retention_percent)

        percent = retention_percent(k);
        nKeep = max(round(nSlices * percent / 100), 1);
        mid = round(nSlices/2);

        % Uniform
        idx_uniform = unique(round(linspace(1, nSlices, nKeep)));

        % Middle-weighted
        [~, priority_mid] = sort(abs((1:nSlices) - mid), 'ascend');
        idx_middle_weighted = sort(priority_mid(1:nKeep));

        % Edge-preserving
        if nKeep == 1
            idx_edge = mid;
        elseif nKeep == 2
            idx_edge = [1 nSlices];
        else
            idx_edge = unique(round(linspace(1, nSlices, nKeep)));
            idx_edge(1) = 1;
            idx_edge(end) = nSlices;
        end

        % Basal-weighted
        basal_center = round(0.65 * nSlices);
        [~, priority_basal] = sort(abs((1:nSlices) - basal_center), 'ascend');
        idx_basal_weighted = sort(priority_basal(1:nKeep));

        all_strategy_idx = {idx_uniform, idx_middle_weighted, idx_edge, idx_basal_weighted};

        for s = 1:length(optStrategyNames)
            idx = all_strategy_idx{s};
            spacing_factor = nSlices / length(idx);

            EDV = sum(ED_area(idx), 'omitnan') * dz * spacing_factor / 1000;
            ESV = sum(ES_area(idx), 'omitnan') * dz * spacing_factor / 1000;

            EF_opt(p,k,s) = (EDV - ESV) / EDV * 100;
        end
    end

    %% Fair regional comparison at 75% slice retention
    nKeep_region = max(round(nSlices * 0.75), 1); %estimate EF at 75% slice retention at different locations

    idx_uniform_region = unique(round(linspace(1, nSlices, nKeep_region)));
    idx_apical = 1:nKeep_region;

    mid_center = round(nSlices/2);
    mid_start = max(round(mid_center - nKeep_region/2), 1);
    mid_end = min(mid_start + nKeep_region - 1, nSlices);
    idx_middle = mid_start:mid_end;

    idx_basal = (nSlices - nKeep_region + 1):nSlices;

    regional_idx = {idx_uniform_region, idx_apical, idx_middle, idx_basal};

    for r = 1:length(regionalNames)
        idx = regional_idx{r};
        spacing_factor = nSlices / length(idx);

        EDV = sum(ED_area(idx), 'omitnan') * dz * spacing_factor / 1000;
        ESV = sum(ES_area(idx), 'omitnan') * dz * spacing_factor / 1000;

        EF_regional(p,r) = (EDV - ESV) / EDV * 100;
    end

    fprintf('Finished patient %03d\n', patientID);
end

%% Errors relative to full SAX
EF_full = EF_all(:,1);
EF_error = EF_all - EF_full;
abs_EF_error = abs(EF_error);

%% Disease-group analysis
diseaseGroups = unique(disease_all);

figure; hold on;
for d = 1:length(diseaseGroups)
    group = diseaseGroups(d);  %select current disease group
    group_idx = disease_all == group; %find patients in this group
    mean_group_error = mean(abs_EF_error(group_idx,:), 1, 'omitnan'); %calculate mean absolute EF error for this group

    plot(retention_percent, mean_group_error, '-o', 'LineWidth', 2, ...
        'DisplayName', char(group));
end

set(gca, 'XDir', 'reverse');
xlabel('Percent of SAX Slices Retained (%)');
ylabel('Mean Absolute EF Error (percentage points)');
title('Disease Group Comparison: EF Error vs SAX Undersampling');
legend('Location','northwest');
grid on;

%% Results table
resultsTable = array2table(EF_all, ...
    'VariableNames', {'Full_SAX','Retain_75pct','Retain_50pct','Retain_25pct'});

resultsTable.PatientID = patientIDs;
resultsTable.Disease = disease_all;
resultsTable.ValidSlices = n_valid_slices_all;
resultsTable.EDV_Full_mL = EDV_full_all;
resultsTable.ESV_Full_mL = ESV_full_all;

resultsTable = movevars(resultsTable, {'PatientID','Disease','ValidSlices'}, 'Before', 1);

disp(resultsTable);

%% Plot 1: Mean absolute EF error
mean_abs_error = mean(abs_EF_error, 1, 'omitnan');
std_abs_error = std(abs_EF_error, 0, 1, 'omitnan');

figure;
errorbar(retention_percent, mean_abs_error, std_abs_error, '-o', 'LineWidth', 2);
set(gca, 'XDir', 'reverse');
xlabel('Percent of SAX Slices Retained (%)');
ylabel('Mean Absolute EF Error (percentage points)');
title('EF Error vs Percent SAX Slices Retained');
grid on;

%% Plot 2: Mean EF estimate
mean_EF = mean(EF_all, 1, 'omitnan');
std_EF = std(EF_all, 0, 1, 'omitnan');

figure;
errorbar(retention_percent, mean_EF, std_EF, '-o', 'LineWidth', 2);
set(gca, 'XDir', 'reverse');
xlabel('Percent of SAX Slices Retained (%)');
ylabel('Estimated EF (%)');
title('Mean EF vs Percent SAX Slices Retained');
grid on;

%% Fair regional comparison plot
EF_regional_error = abs(EF_regional - EF_all(:,1));

mean_regional_error = mean(EF_regional_error, 1, 'omitnan');
std_regional_error = std(EF_regional_error, 0, 1, 'omitnan');

figure;
bar(mean_regional_error);
hold on;
errorbar(1:length(regionalNames), mean_regional_error, std_regional_error, ...
    'k.', 'LineWidth', 1.5);

set(gca, 'XTick', 1:length(regionalNames));
set(gca, 'XTickLabel', regionalNames);
xtickangle(30);

ylabel('Mean Absolute EF Error (percentage points)');
title('Fair Regional Comparison at 75% SAX Slice Retention');
grid on;


%% Optimization experiment plot
EF_opt_error = abs(EF_opt - EF_all(:,1));

figure; hold on;

lineStyles = {'--o','-s',':^','-.d'};

for s = 1:length(optStrategyNames)
    strategy_error = squeeze(EF_opt_error(:,:,s));

    mean_error = mean(strategy_error, 1, 'omitnan');
    std_error = std(strategy_error, 0, 1, 'omitnan');

    errorbar(retention_percent, mean_error, std_error, lineStyles{s}, ...
        'LineWidth', 2, 'MarkerSize', 8, ...
        'DisplayName', optStrategyNames{s});
end

set(gca, 'XDir', 'reverse');
xlabel('Percent of SAX Slices Retained (%)');
ylabel('Mean Absolute EF Error (percentage points)');
title('Optimization of Reduced SAX Sampling');
legend('Location','northwest');
grid on;

%% Save results
resultsDir = 'results';
if ~exist(resultsDir, 'dir')
    mkdir(resultsDir);
end

writetable(resultsTable, fullfile(resultsDir, 'EF_results_from_annotations.csv'));

fprintf('\nSaved results to results/EF_results_from_annotations.csv\n');
