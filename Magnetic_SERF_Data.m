clc; 
clear; 
close all;

file_list = {
    '20250624_105949_synced_cut_aligned.mat'; 
    '20250624_114455_synced_cut_aligned.mat';
    '20250624_143311_synced_cut_aligned.mat'; 
    '20250624_151055_synced_cut_aligned.mat';
    '20250624_154038_synced_change_cut_aligned.mat'; 
    '20250624_160333_synced_cut_aligned.mat';
    '20250624_162816_synced_cut_aligned.mat'; 
    '20250624_165356_synced_cut_aligned.mat';
};
[H_all_axes, dist_vec, num_freqs, target_freqs, comp_names, numPoints] = process_all_core(file_list);
for ax = 1:3
    H_all_axes(:,:,ax) = fillmissing(H_all_axes(:,:,ax), 'linear', 2);
end

pseudo_Z = 1 ./ sqrt(target_freqs);

[z_sort, sort_idx] = sort(pseudo_Z, 'ascend'); 
freq_sort = target_freqs(sort_idx); 
[X_grid, Z_grid] = meshgrid(dist_vec, z_sort);
xq = linspace(0, dist_vec(end), 300); 
zq = linspace(min(z_sort), max(z_sort), 300);
[Xq, Zq] = meshgrid(xq, zq);

figure( 'Color', 'w', 'Position',[50, 50, 1600, 1050]);
fontN = 'Times New Roman'; fontS = 12; 

anom_x_start = 50;
anom_x_width = 30; 

x_ticks_vec = 0:14:84; 
x_max = 84; 

dummy_y_labels = {'16', '8', '4', '2', '1'};

h_ax = zeros(3, 3); 

for ax_idx = 1:3 
    H_current = H_all_axes(:, :, ax_idx);
    Imags = imag(H_current);
    Phases_raw = angle(H_current);
    
    Decoupled_Phases = zeros(num_freqs, numPoints);
    for f_i = 1:num_freqs
        ref_phi = Phases_raw(f_i, 1);
        Decoupled_Phases(f_i, :) = wrapToPi(Phases_raw(f_i, :) - ref_phi) * 1000;
    end
    
    Spatial_PD = Decoupled_Phases(3, :) - Decoupled_Phases(1, :);
    
    Imags_sorted = Imags(sort_idx, :); 
    Phase_sorted = Decoupled_Phases(sort_idx, :);
    
    Imag_interp = interp2(X_grid, Z_grid, Imags_sorted, Xq, Zq, 'spline'); 
    Phase_interp = interp2(X_grid, Z_grid, Phase_sorted, Xq, Zq, 'spline');
    
    h_ax(1, ax_idx) = subplot(3, 3, ax_idx); hold on;
    plot(dist_vec, Spatial_PD, '-o', 'LineWidth', 2, 'Color', '#D95319', 'MarkerFaceColor', '#D95319');
    yline(0, '--k', 'LineWidth', 1.2); 
    grid on; box on; 
    
    xlim([0, x_max]); xticks(x_ticks_vec); xticklabels(string(x_ticks_vec)); 
    
    ylabel([comp_names{ax_idx}, ' Spatial Phase Diff (mrad)'], 'FontName', fontN, 'FontSize', fontS+1, 'FontWeight', 'bold');
    
    if ax_idx == 2
        xlabel({'Distance (m)', '(a)'}, 'FontName', fontN, 'FontSize', fontS+2, 'FontWeight', 'bold');
    else
        xlabel({'Distance (m)', ' '}, 'FontName', fontN, 'FontSize', fontS+2, 'FontWeight', 'bold'); 
    end
    
    yl1 = ylim; ylim([yl1(1)-0.5, yl1(2)+0.5]); yl1 = ylim;
    rectangle('Position',[anom_x_start, yl1(1)+0.05*(yl1(2)-yl1(1)), anom_x_width, 0.9*(yl1(2)-yl1(1))], ...
              'EdgeColor', 'm', 'LineStyle', '--', 'LineWidth', 2.5);
              
    if ax_idx == 1
        text(anom_x_start + anom_x_width/4, yl1(2) - 0.8*(yl1(2)-yl1(1)), 'Known Anomaly Zone', ...
            'Color', 'r', 'FontName', fontN, 'FontSize', fontS, 'FontWeight', 'bold', ...
            'HorizontalAlignment', 'center', 'BackgroundColor',[1 1 1 0.8], 'Margin', 1);
    end
    
    h_ax(2, ax_idx) = subplot(3, 3, ax_idx + 3);
    pcolor(Xq, Zq, Imag_interp); shading interp; hold on; 
    contour(Xq, Zq, Imag_interp, 10, 'LineColor',[0.3 0.3 0.3], 'LineWidth', 0.5); 
    
    colormap(h_ax(2, ax_idx), 'jet'); c2 = colorbar;
    c2.Label.String =[comp_names{ax_idx}, '/I Imaginary (nT/A)']; 
    c2.Label.FontName = fontN; c2.Label.FontSize = fontS; c2.Label.FontWeight = 'bold';
    
    set(gca, 'YDir', 'reverse', 'YTick', z_sort, 'YTickLabel', dummy_y_labels);
    xlim([0, x_max]); xticks(x_ticks_vec); xticklabels(string(x_ticks_vec)); 
    
    if ax_idx == 1, ylabel('Frequency (Hz)', 'FontName', fontN, 'FontSize', fontS+1, 'FontWeight', 'bold'); end
    box on; set(gca, 'Layer', 'top');
    
    if ax_idx == 2
        xlabel({'Distance (m)', '(b)'}, 'FontName', fontN, 'FontSize', fontS+2, 'FontWeight', 'bold');
    else
        xlabel({'Distance (m)', ' '}, 'FontName', fontN, 'FontSize', fontS+2, 'FontWeight', 'bold'); 
    end
    
    yl2 = ylim;
    rectangle('Position',[anom_x_start, yl2(1), anom_x_width, yl2(2)-yl2(1)], 'EdgeColor', 'w', 'LineStyle', '-', 'LineWidth', 3.5);
    rectangle('Position',[anom_x_start, yl2(1), anom_x_width, yl2(2)-yl2(1)], 'EdgeColor', 'm', 'LineStyle', '--', 'LineWidth', 2.5);
    
    h_ax(3, ax_idx) = subplot(3, 3, ax_idx + 6);
    pcolor(Xq, Zq, Phase_interp); shading interp; hold on;
    contour(Xq, Zq, Phase_interp, 10, 'LineColor',[0.3 0.3 0.3], 'LineWidth', 0.5);
    
    colormap(h_ax(3, ax_idx), 'turbo'); c3 = colorbar;
    c3.Label.String =[comp_names{ax_idx}, '/I Absolute Phase (mrad)'];
    c3.Label.FontName = fontN; c3.Label.FontSize = fontS; c3.Label.FontWeight = 'bold';
    
    set(gca, 'YDir', 'reverse', 'YTick', z_sort, 'YTickLabel', dummy_y_labels);
    xlim([0, x_max]); xticks(x_ticks_vec); xticklabels(string(x_ticks_vec)); 
    
    if ax_idx == 1, ylabel('Frequency (Hz)', 'FontName', fontN, 'FontSize', fontS+1, 'FontWeight', 'bold'); end
    box on; set(gca, 'Layer', 'top');
    
    if ax_idx == 2
        xlabel({'Distance (m)', '(c)'}, 'FontName', fontN, 'FontSize', fontS+2, 'FontWeight', 'bold');
    else
        xlabel({'Distance (m)', ' '}, 'FontName', fontN, 'FontSize', fontS+2, 'FontWeight', 'bold'); 
    end
    
    yl3 = ylim;
    rectangle('Position',[anom_x_start, yl3(1), anom_x_width, yl3(2)-yl3(1)], 'EdgeColor', 'w', 'LineStyle', '-', 'LineWidth', 3.5);
    rectangle('Position',[anom_x_start, yl3(1), anom_x_width, yl3(2)-yl3(1)], 'EdgeColor', 'm', 'LineStyle', '--', 'LineWidth', 2.5);
end

all_axes = findobj(gcf, 'Type', 'axes');
set(all_axes, 'FontName', fontN, 'FontSize', fontS);

drawnow; 
for c = 1:3 
    pos2 = get(h_ax(2, c), 'Position'); 
    
    pos1 = get(h_ax(1, c), 'Position'); 
    pos1(1) = pos2(1); pos1(3) = pos2(3); 
    set(h_ax(1, c), 'Position', pos1);
    
    pos3 = get(h_ax(3, c), 'Position'); 
    pos3(1) = pos2(1); pos3(3) = pos2(3); 
    set(h_ax(3, c), 'Position', pos3);
end
