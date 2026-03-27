
clc;
clear;
close all;

fs = 96000; 
target_freqs = [16, 8, 4, 2, 1]; 
num_freqs = length(target_freqs);

fprintf('Please choose .mat \n');[fileNames, pathName] = uigetfile('*.mat',  'MultiSelect', 'on');

if isequal(fileNames, 0), disp('取消操作'); return; end
if ischar(fileNames), fileNames = {fileNames}; end

fileNames = sort(fileNames); 
numFiles = length(fileNames);
numPoints = numFiles * 2 - 1; 



[All_Phases_mrad, All_Amps_V, x_labels, global_ref_sig, valid_point_idx] ...
    = process_core(fileNames, pathName, fs, target_freqs);


Ref_Phase = All_Phases_mrad(1, :);
Rel_Phases_mrad = All_Phases_mrad - repmat(Ref_Phase, numPoints, 1);

idx_4Hz = find(target_freqs == 4);
idx_1Hz = find(target_freqs == 1);
DFIP_1D = Rel_Phases_mrad(:, idx_4Hz) - Rel_Phases_mrad(:, idx_1Hz);

dx = 7; 
dist_vec = 0:dx:(numPoints-1)*dx; 

pseudo_Z_vec = 1 ./ sqrt(target_freqs); 
[X_grid, Z_grid] = meshgrid(dist_vec, pseudo_Z_vec); % X轴改用距离
Amp_matrix = All_Amps_V';
Phase_matrix = Rel_Phases_mrad';

xq = linspace(0, dist_vec(end), 300); 
zq_main = linspace(min(pseudo_Z_vec), max(pseudo_Z_vec), 300);[Xq_main, Zq_main] = meshgrid(xq, zq_main);

Amp_interp   = interp2(X_grid, Z_grid, Amp_matrix, Xq_main, Zq_main, 'spline');
Phase_interp = interp2(X_grid, Z_grid, Phase_matrix*5, Xq_main, Zq_main, 'spline');

figure('Color', 'w', 'Position',[50, 50, 1000, 1050]);

fontName = 'Times New Roman';
fontSize = 12;

anom_x_start = dist_vec(8);
anom_x_width = dist_vec(13) - dist_vec(8);

ax1 = subplot(3, 1, 1); hold on;
plot(dist_vec, DFIP_1D*10, '-o', 'Color', '#D95319', 'LineWidth', 2, ...
    'MarkerFaceColor', '#D95319', 'MarkerSize', 6);
yline(0, '--k', 'LineWidth', 1.2); 
set(gca, 'FontName', fontName, 'FontSize', fontSize, 'LineWidth', 1);

xlim([0, dist_vec(end)]); 
xticks(dist_vec); 
xticklabels(string(dist_vec));


ylabel('Elect D-Freq Phase Dif (mrad)', 'FontName', fontName, 'FontSize', fontSize+1, 'FontWeight', 'bold'); 
xlabel({'Distance (m)', '(a)'}, 'FontName', fontName, 'FontSize', fontSize+2, 'FontWeight', 'bold'); 
yl1 = ylim; ylim([yl1(1)-1, yl1(2)+1]); yl1 = ylim; 
grid on; box on; 

rectangle('Position',[anom_x_start, yl1(1)+0.02, anom_x_width, yl1(2)-yl1(1)-0.04], ...
          'EdgeColor', 'm', 'LineStyle', '--', 'LineWidth', 2);
text(anom_x_start + anom_x_width/1.4, yl1(2) - 0.15*(yl1(2)-yl1(1)), 'Known Anomaly Zone', ...
    'Color', 'r', 'FontName', fontName, 'FontSize', fontSize, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center', 'BackgroundColor', [1 1 1 0.8], 'Margin', 2);

ax2 = subplot(3, 1, 2); hold on;
pcolor(Xq_main, Zq_main, Amp_interp); shading interp; 
contour(Xq_main, Zq_main, Amp_interp, 15, 'LineColor',[0.3 0.3 0.3], 'LineWidth', 0.5); 
colormap(ax2, 'jet'); c2 = colorbar; 
c2.Label.String = 'Voltage Amplitude (V)'; 
c2.Label.FontName = fontName; c2.Label.FontSize = fontSize;

set(gca, 'YDir', 'reverse', 'FontName', fontName, 'FontSize', fontSize, 'LineWidth', 1); 
xlim([0, dist_vec(end)+1]); 
yticks(pseudo_Z_vec); yticklabels(string(target_freqs)); 
xticks(dist_vec); 
xticklabels(string(dist_vec)); 

ylabel('Frequency (Hz)', 'FontName', fontName, 'FontSize', fontSize+1, 'FontWeight', 'bold'); 
xlabel({'Distance (m)', '(b)'}, 'FontName', fontName, 'FontSize', fontSize+2, 'FontWeight', 'bold'); 
axis tight; box on; set(gca, 'Layer', 'top');

yl2 = ylim;
rectangle('Position',[anom_x_start, yl2(1), anom_x_width, yl2(2)-yl2(1)], ...
          'EdgeColor', 'w', 'LineStyle', '-', 'LineWidth', 4);
rectangle('Position',[anom_x_start, yl2(1), anom_x_width, yl2(2)-yl2(1)], ...
          'EdgeColor', 'm', 'LineStyle', '--', 'LineWidth', 2);

ax3 = subplot(3, 1, 3); hold on;
pcolor(Xq_main, Zq_main, Phase_interp); shading interp; 
contour(Xq_main, Zq_main, Phase_interp, 15, 'LineColor',[0.3 0.3 0.3], 'LineWidth', 0.5); 
colormap(ax3, 'turbo'); c3 = colorbar; 
c3.Label.String = 'Electrical Absolute Phase (mrad)'; 
c3.Label.FontName = fontName; c3.Label.FontSize = fontSize;

set(gca, 'YDir', 'reverse', 'FontName', fontName, 'FontSize', fontSize, 'LineWidth', 1);
xlim([0, dist_vec(end)+1]); 
yticks(pseudo_Z_vec); yticklabels(string(target_freqs)); 
xticks(dist_vec); 
xticklabels(string(dist_vec)); 

ylabel('Frequency (Hz)', 'FontName', fontName, 'FontSize', fontSize+1, 'FontWeight', 'bold'); 
xlabel({'Distance (m)', '(c)'}, 'FontName', fontName, 'FontSize', fontSize+2, 'FontWeight', 'bold'); 
axis tight; box on; set(gca, 'Layer', 'top');

yl3 = ylim;
rectangle('Position',[anom_x_start, yl3(1), anom_x_width, yl3(2)-yl3(1)], ...
          'EdgeColor', 'w', 'LineStyle', '-', 'LineWidth', 4);
rectangle('Position',[anom_x_start, yl3(1), anom_x_width, yl3(2)-yl3(1)], ...
          'EdgeColor', 'm', 'LineStyle', '--', 'LineWidth', 2);

drawnow; 

pos2 = get(ax2, 'Position'); 
pos1 = get(ax1, 'Position'); 
pos1(1) = pos2(1); pos1(3) = pos2(3); 
set(ax1, 'Position', pos1);

pos3 = get(ax3, 'Position');
pos3(1)