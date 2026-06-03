% run_EF_hybrid_SAX_LAX_kfold.m
clear; clc; close all;

T = readtable('annotations.csv'); % load CSV

SAX = T(strcmp(T.view, 'SAX'), :); % separate SAX and LAX rows
LAX = T(strcmp(T.view, 'LAX'), :);

%% Setup
retention_percent = [100 75 50 25];
patientIDs = intersect(unique(SAX.patient_id), unique(LAX.patient_id)); %keep only patients with both SAX and LAX images

SAX_EF = nan(length(patientIDs), length(retention_percent));
Hybrid_EF_CV = nan(length(patientIDs), length(retention_percent));
LAX_FAC = nan(length(patientIDs), 1);
disease_all = strings(length(patientIDs), 1);

%% First pass: calculate SAX EF and LAX FAC
for p = 1:length(patientIDs)

    patientID = patientIDs(p);

    %% SAX data
    %load SAX slices and remove incomplete annotations
    saxData = SAX(SAX.patient_id == patientID, :);
    saxData = sortrows(saxData, 'z_slice');

    ED_area = saxData.lv_area_ED_mm2;
    ES_area = saxData.lv_area_ES_mm2;

   valid_idx = ~isnan(ED_area) & ...
            ~isnan(ES_area) & ...
            (ED_area > 0) & ...
            (ES_area > 0);
    saxData = saxData(valid_idx, :);
    ED_area = ED_area(valid_idx);
    ES_area = ES_area(valid_idx);

    dz = mean(saxData.z_spacing_mm, 'omitnan');
    nSlices = length(ED_area);

    disease_all(p) = string(saxData.disease{1});

    %% SAX EF at each retention level
    %estimate EF using progressively fewer SAX slices
    for k = 1:length(retention_percent)

        percent = retention_percent(k);
        nKeep = max(round(nSlices * percent / 100), 1);

        slice_idx = unique(round(linspace(1, nSlices, nKeep)));
        spacing_factor = nSlices / length(slice_idx);

        EDV = sum(ED_area(slice_idx), 'omitnan') * dz * spacing_factor / 1000;
        ESV = sum(ES_area(slice_idx), 'omitnan') * dz * spacing_factor / 1000;

        SAX_EF(p,k) = (EDV - ESV) / EDV * 100;
    end

    %% LAX data
    %calculate LAX fractional area change
    laxData = LAX(LAX.patient_id == patientID, :);

    lax_ED = mean(laxData.lv_area_ED_mm2, 'omitnan');
    lax_ES = mean(laxData.lv_area_ES_mm2, 'omitnan');

    LAX_FAC(p) = (lax_ED - lax_ES) / lax_ED * 100;

    fprintf('Finished patient %03d\n', patientID);
end

%% 5-fold cross-validation hybrid model
%use LAX FAC and reduced-SAX EF together to predict the full SAX EF
%reference value
Full_SAX_EF = SAX_EF(:,1);

validPatients = ~isnan(Full_SAX_EF) & ~isnan(LAX_FAC);
validIdx = find(validPatients);
%randomly split patients into training and testing groups using 5-fold
%cross validation partitions
rng(1);
K = 5;
cv = cvpartition(length(validIdx), 'KFold', K);

fprintf('\nRunning %d-fold cross-validation on %d patients.\n', K, length(validIdx));
%train and test hybrid model
for k = 1:length(retention_percent)

    reduced_EF = SAX_EF(:,k);

    for fold = 1:K
%select training and testing patients for current fold
        train_local = training(cv, fold);
        test_local = test(cv, fold);

        train_idx = validIdx(train_local);
        test_idx = validIdx(test_local);

        train_idx = train_idx(~isnan(reduced_EF(train_idx)));
        test_idx = test_idx(~isnan(reduced_EF(test_idx)));
%construct predictor matrix using reduced SAX Ef and LAX FAC
        X_train = [reduced_EF(train_idx), LAX_FAC(train_idx)];
        y_train = Full_SAX_EF(train_idx);

        X_test = [reduced_EF(test_idx), LAX_FAC(test_idx)];
%train linear regression model and predict test patients
        mdl = fitlm(X_train, y_train);

        Hybrid_EF_CV(test_idx,k) = predict(mdl, X_test);
    end

    fprintf('Completed %d%% retention hybrid CV model.\n', retention_percent(k));
end

%% Error comparison across all cross-validated test predictions
%calculate prediction errors for SAX-only and hybrid approach
SAX_error = abs(SAX_EF(validIdx,:) - repmat(Full_SAX_EF(validIdx), 1, length(retention_percent)));
Hybrid_error = abs(Hybrid_EF_CV(validIdx,:) - repmat(Full_SAX_EF(validIdx), 1, length(retention_percent)));

mean_SAX_error = mean(SAX_error, 1, 'omitnan');
std_SAX_error = std(SAX_error, 0, 1, 'omitnan');

mean_Hybrid_error = mean(Hybrid_error, 1, 'omitnan');
std_Hybrid_error = std(Hybrid_error, 0, 1, 'omitnan');

disp('Mean SAX-only error:')
disp(mean_SAX_error)

disp('Mean Hybrid CV error:')
disp(mean_Hybrid_error)

%% Percent improvement
improvement = 100 * (mean_SAX_error - mean_Hybrid_error) ./ mean_SAX_error;

fprintf('\nHybrid improvement with %d-fold cross-validation:\n', K);
for k = 1:length(retention_percent)
    fprintf('%d%% retained: %.1f%% error reduction\n', ...
        retention_percent(k), improvement(k));
end

%% Plot: SAX-only vs Hybrid CV error

figure;
errorbar(retention_percent, mean_SAX_error, std_SAX_error, '-o', ...
    'LineWidth', 2, 'DisplayName', 'Reduced SAX only');
hold on;

errorbar(retention_percent, mean_Hybrid_error, std_Hybrid_error, '-s', ...
    'LineWidth', 2, 'DisplayName', 'Hybrid SAX + LAX');

set(gca, 'XDir', 'reverse');
xlabel('Percent of SAX Slices Retained (%)');
ylabel('Mean Absolute EF Error (percentage points)');
title('5-Fold Cross-Validated Hybrid SAX + LAX EF Estimation');
legend('Location','northwest');
grid on;

%% Plot: LAX FAC vs Full SAX EF

figure;
scatter(LAX_FAC(validIdx), Full_SAX_EF(validIdx), 'filled');
xlabel('LAX Fractional Area Change (%)');
ylabel('Full SAX EF (%)');
title('Relationship Between LAX FAC and Full SAX EF');
grid on;

%% Save results

resultsTable = array2table([SAX_EF Hybrid_EF_CV], ...
    'VariableNames', {'SAX_100','SAX_75','SAX_50','SAX_25', ...
                      'HybridCV_100','HybridCV_75','HybridCV_50','HybridCV_25'});

resultsTable.PatientID = patientIDs;
resultsTable.Disease = disease_all;
resultsTable.LAX_FAC = LAX_FAC;

resultsTable = movevars(resultsTable, {'PatientID','Disease','LAX_FAC'}, 'Before', 1);

resultsDir = 'results';
if ~exist(resultsDir, 'dir')
    mkdir(resultsDir);
end

writetable(resultsTable, fullfile(resultsDir, 'EF_hybrid_SAX_LAX_kfold_results.csv'));

fprintf('\nSaved results to results/EF_hybrid_SAX_LAX_kfold_results.csv\n');