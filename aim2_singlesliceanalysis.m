%% run_EF_selected_figures.m
clear; clc; close all;

%% Load CSV
T = readtable('annotations.csv');
T_orig = T;

%% Keep only SAX rows
T = T(strcmp(T.view, 'SAX'), :);

%% Setup
retention_percent = [100 90 75 50 25];
patientIDs = unique(T.patient_id);

EF_all = nan(length(patientIDs), length(retention_percent));
disease_all = strings(length(patientIDs), 1);

%% Loop through patients
for p = 1:length(patientIDs)

    patientID = patientIDs(p);
    patientData = T(T.patient_id == patientID, :);
    patientData = sortrows(patientData, 'z_slice');

    ED_area = patientData.lv_area_ED_mm2;
    ES_area = patientData.lv_area_ES_mm2;

    valid_idx = ~isnan(ED_area) & ~isnan(ES_area) & (ED_area > 0) & (ES_area > 0);
    patientData = patientData(valid_idx, :);
    ED_area = ED_area(valid_idx);
    ES_area = ES_area(valid_idx);

    dz = mean(patientData.z_spacing_mm, 'omitnan');
    disease_all(p) = string(patientData.disease{1});

    nSlices = length(ED_area);

    for k = 1:length(retention_percent)
        percent = retention_percent(k);
        nKeep = max(round(nSlices * percent / 100), 1);
        slice_idx = unique(round(linspace(1, nSlices, nKeep)));
        spacing_factor = nSlices / length(slice_idx);
        EDV = sum(ED_area(slice_idx), 'omitnan') * dz * spacing_factor / 1000;
        ESV = sum(ES_area(slice_idx), 'omitnan') * dz * spacing_factor / 1000;
        EF_all(p,k) = (EDV - ESV) / EDV * 100;
    end

    fprintf('Finished patient %03d\n', patientID);
end

%% Shared derived variables
EF_full = EF_all(:,1);
diseaseGroups = unique(disease_all);
diseaseGroups = ["NOR"; diseaseGroups(diseaseGroups ~= "NOR")];
group_list = cellstr(diseaseGroups);
colors = lines(length(diseaseGroups));
nor_pos = find(diseaseGroups == "NOR");
n_non_nor = sum(diseaseGroups ~= "NOR");

%% LAX and mid-SAX reduction
T_lax = T_orig(strcmp(T_orig.view, 'LAX'), :);

lax_reduction = nan(length(patientIDs), 1);
sax_mid_reduction = nan(length(patientIDs), 1);

for p = 1:length(patientIDs)
    patientID = patientIDs(p);

    lax = T_lax(T_lax.patient_id == patientID, :);
    if ~isempty(lax)
        lax_reduction(p) = (lax.lv_area_ED_mm2 - lax.lv_area_ES_mm2) / lax.lv_area_ED_mm2 * 100;
    end

    sax = T(T.patient_id == patientID, :);
    sax = sortrows(sax, 'z_slice');
    sax_ed = sax.lv_area_ED_mm2;
    sax_es = sax.lv_area_ES_mm2;
    valid = sax_ed > 0 & sax_es > 0 & ~isnan(sax_ed) & ~isnan(sax_es);
    sax_ed = sax_ed(valid);
    sax_es = sax_es(valid);
    if ~isempty(sax_ed)
        mid = round(length(sax_ed) / 2);
        sax_mid_reduction(p) = (sax_ed(mid) - sax_es(mid)) / sax_ed(mid) * 100;
    end
end

%% Middle SAX Slice LV Reduction Boxplot
valid = ~isnan(sax_mid_reduction);
plot_boxplot(sax_mid_reduction(valid), disease_all(valid), diseaseGroups, group_list, colors, nor_pos, n_non_nor, ...
    'LV Area Reduction (%)', 'Middle SAX Slice LV Reduction by Disease Group');

%% LAX LV Reduction Boxplot
valid = ~isnan(lax_reduction);
plot_boxplot(lax_reduction(valid), disease_all(valid), diseaseGroups, group_list, colors, nor_pos, n_non_nor, ...
    'LV Area Reduction (%)', 'LAX LV Reduction by Disease Group');

%% Full SAX EF Boxplot
valid = ~isnan(EF_full);
plot_boxplot(EF_full(valid), disease_all(valid), diseaseGroups, group_list, colors, nor_pos, n_non_nor, ...
    'Ejection Fraction (%)', 'Full SAX EF by Disease Group');

%% ROC curves and AUC comparison
auc_table_groups = {};
auc_table_sax_mid = [];
auc_table_lax = [];
auc_table_full_ef = [];

figure; hold on;
for d = 1:length(diseaseGroups)
    group = diseaseGroups(d);
    if group == "NOR"; continue; end

    pair_idx = disease_all == group | disease_all == "NOR";
    labels = double(disease_all(pair_idx) == group);

    scores_sax = sax_mid_reduction(pair_idx);
    valid = ~isnan(scores_sax);
    [~, ~, ~, auc_sax] = perfcurve(labels(valid), scores_sax(valid), 1);

    scores_lax = lax_reduction(pair_idx);
    valid = ~isnan(scores_lax);
    [~, ~, ~, auc_lax] = perfcurve(labels(valid), scores_lax(valid), 1);

    scores_ef = EF_all(pair_idx, 1);
    valid = ~isnan(scores_ef);
    [X, Y, ~, auc_ef] = perfcurve(labels(valid), scores_ef(valid), 1);

    auc_table_groups{end+1} = char(group);
    auc_table_sax_mid(end+1) = max(auc_sax, 1-auc_sax);
    auc_table_lax(end+1) = max(auc_lax, 1-auc_lax);
    auc_table_full_ef(end+1) = max(auc_ef, 1-auc_ef);

    auc_ef_corrected = max(auc_ef, 1-auc_ef);
    if auc_ef < 0.5
        plot(1-X, 1-Y, '-', 'LineWidth', 2, 'Color', colors(d,:), ...
            'DisplayName', sprintf('%s (AUC=%.2f)', char(group), auc_ef_corrected));
    else
        plot(X, Y, '-', 'LineWidth', 2, 'Color', colors(d,:), ...
            'DisplayName', sprintf('%s (AUC=%.2f)', char(group), auc_ef_corrected));
    end
end
plot([0 1], [0 1], 'k--', 'HandleVisibility', 'off');
xlabel('False Positive Rate');
ylabel('True Positive Rate');
title('ROC Curves: Full SAX EF vs NOR');
legend('Location', 'southeast');
grid on;

figure;
auc_matrix = [auc_table_full_ef; auc_table_sax_mid; auc_table_lax]';
b = bar(auc_matrix, 'grouped');
colors_bar = lines(3);
for k = 1:3; b(k).FaceColor = colors_bar(k,:); end
set(gca, 'XTick', 1:length(auc_table_groups), 'XTickLabel', auc_table_groups);
ylabel('AUC');
title('AUC Comparison: Full SAX EF vs Single Slice');
legend({'Full SAX EF','Mid SAX Slice','LAX Slice'}, 'Location', 'southwest');
yline(0.5, 'k--', 'Chance', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off');
yline(0.7, 'k:', 'Moderate', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off');
yline(0.8, 'k-.', 'Good', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off');
ylim([0 1]); grid on;

%% Figure 6: SAX vs LAX Correlation
fprintf('\n=== SAX vs LAX Correlation ===\n');
figure; hold on;

for d = 1:length(diseaseGroups)
    group = diseaseGroups(d);
    idx = disease_all == group;
    x = sax_mid_reduction(idx);
    y = lax_reduction(idx);
    valid = ~isnan(x) & ~isnan(y);
    x = x(valid); y = y(valid);
    if length(x) >= 3
        [r, p_corr] = corr(x, y, 'Type', 'Spearman');
        fprintf('%s: r=%.3f, p=%.4f (n=%d)\n', char(group), r, p_corr, length(x));
    end
    scatter(x, y, 60, colors(d,:), 'filled', 'DisplayName', char(group));
end

valid = ~isnan(sax_mid_reduction) & ~isnan(lax_reduction);
[r_all, p_all] = corr(sax_mid_reduction(valid), lax_reduction(valid), 'Type', 'Spearman');
fprintf('Overall: r=%.3f, p=%.4f (n=%d)\n', r_all, p_all, sum(valid));

x_all = sax_mid_reduction(valid);
y_all = lax_reduction(valid);
p_fit = polyfit(x_all, y_all, 1);
x_fit = linspace(min(x_all), max(x_all), 100);
plot(x_fit, polyval(p_fit, x_fit), 'k-', 'LineWidth', 2, 'HandleVisibility', 'off');
xl = xlim; plot(xl, xl, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');

xlabel('Middle SAX Slice LV Reduction (%)');
ylabel('LAX LV Reduction (%)');
title(sprintf('SAX vs LAX Correlation (Overall r=%.3f)', r_all));
legend('Location', 'northwest');
grid on;