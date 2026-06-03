%% Helper function to plot boxplot with significance markers
function plot_boxplot(feat_clean, dis_clean, diseaseGroups, group_list, colors, nor_pos, n_non_nor, y_label, fig_title)
    nor_vals = feat_clean(strcmp(cellstr(dis_clean), 'NOR'));
    figure; hold on;
    for d = 1:length(diseaseGroups)
        group = diseaseGroups(d);
        idx = strcmp(cellstr(dis_clean), char(group));
        vals = feat_clean(idx);
        boxplot(vals, repmat(group_list(d), sum(idx), 1), ...
            'Positions', d, 'Widths', 0.5, 'Colors', colors(d,:), 'Symbol', '');
        jitter = (rand(sum(idx), 1) - 0.5) * 0.3;
        scatter(d + jitter, vals, 40, colors(d,:), 'filled', 'MarkerFaceAlpha', 0.6);
    end
    y_min = min(feat_clean); y_max = max(feat_clean) * 1.05;
    y_step = (max(feat_clean) - min(feat_clean)) * 0.07;
    for d = 1:length(diseaseGroups)
        group = diseaseGroups(d);
        if group == "NOR"; continue; end
        grp_vals = feat_clean(strcmp(cellstr(dis_clean), char(group)));
        [~, p_val] = ttest2(nor_vals, grp_vals);
        if p_val < 0.001; sig = '***'; elseif p_val < 0.01; sig = '**'; elseif p_val < 0.05; sig = '*'; else; sig = 'ns'; end
        y_line = y_max + y_step * (d - 1);
        plot([nor_pos d], [y_line y_line], 'k-', 'LineWidth', 1);
        text(mean([nor_pos d]), y_line + y_step * 0.3, sig, 'HorizontalAlignment', 'center', 'FontSize', 11);
    end
    ylim([y_min * 0.95, y_max + y_step * (n_non_nor + 1)]);
    set(gca, 'XTick', 1:length(diseaseGroups), 'XTickLabel', group_list);
    ylabel(y_label);
    title(fig_title);
    grid on;
end