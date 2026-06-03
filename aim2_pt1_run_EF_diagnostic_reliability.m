% run_EF_diagnostic_reliability.m
clear; clc; close all;

T = readtable('annotations.csv'); % Load CSV

T = T(strcmp(T.view, 'SAX'), :); %keep only SAX rows

%%Setup
retention_percent = [100 75 50 25];
patientIDs = unique(T.patient_id);

EF_all = nan(length(patientIDs), length(retention_percent));
disease_all = strings(length(patientIDs), 1);

%%Loop through patients
for p = 1:length(patientIDs)

    patientID = patientIDs(p);
    patientData = T(T.patient_id == patientID, :);
    patientData = sortrows(patientData, 'z_slice');

    ED_area = patientData.lv_area_ED_mm2;
    ES_area = patientData.lv_area_ES_mm2;


    valid_idx = ~isnan(ED_area) & ~isnan(ES_area) & (ED_area > 0) & (ES_area > 0);  % remove slices missing either ED or ES area
    patientData = patientData(valid_idx, :);
    ED_area = ED_area(valid_idx);
    ES_area = ES_area(valid_idx);

    dz = mean(patientData.z_spacing_mm, 'omitnan');
    disease_all(p) = string(patientData.disease{1});

    nSlices = length(ED_area);

    %%Calculate EF at each retention level
    for k = 1:length(retention_percent)

        percent = retention_percent(k);
        nKeep = max(round(nSlices * percent / 100), 1);

        % Uniform reduced sampling
        slice_idx = unique(round(linspace(1, nSlices, nKeep)));
        spacing_factor = nSlices / length(slice_idx);

        EDV = sum(ED_area(slice_idx), 'omitnan') * dz * spacing_factor / 1000;
        ESV = sum(ES_area(slice_idx), 'omitnan') * dz * spacing_factor / 1000;

        EF_all(p,k) = (EDV - ESV) / EDV * 100;
    end

   fprintf('Finished patient %03d\n', patientID);
end

%% Clinical EF category reclassification
% Categories: >70 (hyperdynamic), 50-70 (normal), 41-49 (mildly reduced), <=40 (reduced)
category_labels = {'Hyperdynamic (>70)','Normal (50-70)','Mildly Reduced (41-49)','Reduced (<=40)'};

EF_full = EF_all(:,1);
classify = @(EF) (EF > 70)*1 + (EF >= 50 & EF <= 70)*2 + (EF > 40 & EF < 50)*3 + (EF <= 40)*4;
n_patients = length(patientIDs);
reclass_rate = nan(1, length(retention_percent));

class_full = classify(EF_full);

for k = 1:length(retention_percent)
    class_k = classify(EF_all(:,k));
    reclass_rate(k) = mean(class_full ~= class_k, 'omitnan') * 100;
end

% Stacked bar: for each retention level, count patients per category
figure;
category_counts = nan(length(retention_percent), 4);
for k = 1:length(retention_percent)
    class_k = classify(EF_all(:,k));
    for cat = 1:4
        category_counts(k,cat) = sum(class_k == cat, 'omitnan');
    end
end

b = bar(category_counts, 'stacked');
for k = 2:length(retention_percent)
    text(k, n_patients + 4, sprintf('%.1f%% reclass.', reclass_rate(k)), ...
        'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold');
end
set(gca, 'XTickLabel', arrayfun(@(x) sprintf('%d%%',x), retention_percent, 'UniformOutput', false));
xlabel('Percent of SAX Slices Retained (%)', 'FontSize', 15);
ylabel('Number of Patients', 'FontSize', 15);
title('EF Category Distribution vs SAX Undersampling', 'FontSize', 15);
legend(category_labels, 'Location', 'eastoutside', 'FontSize', 15);
grid on;

% Print reclassification rate
fprintf('\nCategory Reclassification Rate:\n');
for k = 1:length(retention_percent)
    fprintf('Retain %d%%: %.1f%% patients reclassified\n', retention_percent(k), reclass_rate(k));
end

%% Reclassification direction bar chart
n_upgraded = nan(1, length(retention_percent)-1);
n_downgraded = nan(1, length(retention_percent)-1);
n_to_normal = nan(1, length(retention_percent)-1);

for k = 2:length(retention_percent)
    class_k = classify(EF_all(:,k));
    changed = class_full ~= class_k & class_full > 0 & class_k > 0;
    n_upgraded(k-1) = sum(class_k(changed) < class_full(changed));
    n_downgraded(k-1) = sum(class_k(changed) > class_full(changed));
    n_to_normal(k-1) = sum(class_k(changed) == 2 & class_full(changed) ~= 2);
end

figure;
b = bar([n_upgraded; n_downgraded]', 'grouped');
b(1).DisplayName = 'Upgraded';
b(2).DisplayName = 'Downgraded';
hold on;

x = 1:length(retention_percent)-1;
scatter(x - 0.15, n_to_normal, 80, 'k^', 'filled', 'DisplayName', 'Reclassified to Normal');
xticklabels(arrayfun(@(x) sprintf('%d%%', x), retention_percent(2:end), 'UniformOutput', false));
xlabel('Percent of SAX Slices Retained (%)', 'FontSize', 15);
ylabel('Number of Patients', 'FontSize', 15);
title('EF Category Reclassification Direction', 'FontSize', 15);
legend({'Upgraded','Downgraded','Reclassified to Normal'}, 'Location', 'best', 'FontSize', 15);
grid on;


%%Define healthy vs diseased groups
healthy_idx = disease_all == "NOR";
diseased_idx = ~healthy_idx;

%%Boxplot: EF distributions by retention level and diagnosis group
%convert EF results so MATLAB can group by retention level and diagnosis
%category
figure;

groupLabels = strings(length(patientIDs) * length(retention_percent), 1);
EF_values = nan(length(patientIDs) * length(retention_percent), 1);
retentionLabels = strings(length(patientIDs) * length(retention_percent), 1);

row = 1;
%populate plotting arrays with EF values, diagnosis labels and slice
%retention labels
for k = 1:length(retention_percent)

    for p = 1:length(patientIDs)

        EF_values(row) = EF_all(p,k);

        if healthy_idx(p)
            groupLabels(row) = "Healthy";
        else
            groupLabels(row) = "Diseased";
        end

        retentionLabels(row) = string(retention_percent(k)) + "%";
        row = row + 1;
    end
end
%plot
boxchart(categorical(retentionLabels), EF_values, 'GroupByColor', categorical(groupLabels));
xlabel('Percent of SAX Slices Retained');
ylabel('Estimated EF (%)');
title('Healthy vs Diseased EF Distributions Across SAX Retention Levels');
legend;
grid on;

%%Statistical comparison: healthy vs diseased at each retention level
%test significance difference in EF between healthy and diseased patients
%at each retention level 
p_values = nan(length(retention_percent), 1);
mean_difference = nan(length(retention_percent), 1);

fprintf('\nHealthy vs Diseased EF Comparison:\n');

%independent two-sample t-tests
for k = 1:length(retention_percent)

    EF_healthy = EF_all(healthy_idx,k);
    EF_diseased = EF_all(diseased_idx,k);

    EF_healthy = EF_healthy(~isnan(EF_healthy));
    EF_diseased = EF_diseased(~isnan(EF_diseased));

    [~, p_val] = ttest2(EF_healthy, EF_diseased);

    p_values(k) = p_val;
    mean_difference(k) = mean(EF_healthy, 'omitnan') - mean(EF_diseased, 'omitnan');

    fprintf('%d%% retained: mean difference = %.2f EF points, p = %.4f\n', ...
        retention_percent(k), mean_difference(k), p_val);
end
%%Mean EF by group at each retention level
%solve for mean EF and variability for healthy and diseased patients
healthy_mean = mean(EF_all(healthy_idx,:), 1, 'omitnan');
healthy_std = std(EF_all(healthy_idx,:), 0, 1, 'omitnan');

diseased_mean = mean(EF_all(diseased_idx,:), 1, 'omitnan');
diseased_std = std(EF_all(diseased_idx,:), 0, 1, 'omitnan');
%plot
figure;
errorbar(retention_percent, healthy_mean, healthy_std, '-o', 'LineWidth', 2, ...
    'DisplayName', 'Healthy');
hold on;
errorbar(retention_percent, diseased_mean, diseased_std, '-o', 'LineWidth', 2, ...
    'DisplayName', 'Diseased');

set(gca, 'XDir', 'reverse');
xlabel('Percent of SAX Slices Retained (%)');
ylabel('Estimated EF (%)');
title('Mean EF: Healthy vs Diseased Across Retention Levels');
legend('Location','best');
grid on;
% Add significance markers

ymax = max([healthy_mean + healthy_std, ...
            diseased_mean + diseased_std]);

for k = 1:length(retention_percent)

    if p_values(k) < 0.001
        sigText = '***';
    elseif p_values(k) < 0.01
        sigText = '**';
    elseif p_values(k) < 0.05
        sigText = '*';
    else
        sigText = 'ns';
    end

    text(retention_percent(k), ...
         68, ...
         sigText, ...
         'HorizontalAlignment','center', ...
         'FontSize',12, ...
         'FontWeight','bold');
end
% significance key
annotation('textbox',[0.15 0.15 0.15 0.10], ...
    'String',{'* p<0.05','** p<0.01','*** p<0.001'}, ...
    'FitBoxToText','on', ...
    'EdgeColor','none');

%%Save results table
resultsTable = array2table(EF_all, ...
    'VariableNames', {'Full_SAX','Retain_75pct','Retain_50pct','Retain_25pct'});

resultsTable.PatientID = patientIDs;
resultsTable.Disease = disease_all;

if any(healthy_idx)
    diagnosis = strings(length(patientIDs), 1);
    diagnosis(healthy_idx) = "Healthy";
    diagnosis(diseased_idx) = "Diseased";
    resultsTable.DiagnosisGroup = diagnosis;
end

resultsTable = movevars(resultsTable, {'PatientID','Disease','DiagnosisGroup'}, 'Before', 1);

resultsDir = 'results';
if ~exist(resultsDir, 'dir')
    mkdir(resultsDir);
end

writetable(resultsTable, fullfile(resultsDir, 'EF_diagnostic_reliability_results.csv'));

fprintf('\nSaved results to results/EF_diagnostic_reliability_results.csv\n');
