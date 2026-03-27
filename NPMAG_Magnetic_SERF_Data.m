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
numPoints = length(file_list);
dx = 12; dist_vec = 0:dx:(numPoints-1)*dx;
freq_labels = {'1 Hz', '2 Hz', '4 Hz', '8 Hz', '16 Hz'};

bg_mask = (dist_vec <= 48); 
anom_x_start = 50; anom_x_end = 80;
[H_all_axes, num_freqs, target_freqs, Fs] = extract_features_core(file_list, numPoints);

NPMAG_total = zeros(num_freqs, numPoints);
for f_i = 1:num_freqs
    sum_grad_sq = zeros(1, numPoints); 
    sum_Bpri_sq = zeros(1, numPoints); 
    for ax_idx = 1:3
        H_cur = H_all_axes(f_i, :, ax_idx);
        B_pri_amp = abs(H_cur);
        Phases_raw = unwrap(angle(H_cur)); 
        p_bg = polyfit(dist_vec(bg_mask), Phases_raw(bg_mask), 1);
        delta_phi = Phases_raw - polyval(p_bg, dist_vec);
        B_pol = B_pri_amp .* sin(delta_phi);
        p_pol = polyfit(dist_vec(bg_mask), B_pol(bg_mask), 1);
        B_pol_level = B_pol - polyval(p_pol, dist_vec);
        bg_std = std(B_pol_level(bg_mask)); total_std = std(B_pol_level);
        auto_sigma = max(1.5, min(3.5, 5 * (bg_std/total_std)));
        B_pol_smooth = smoothdata(B_pol_level, 'gaussian', auto_sigma);
        dB_dx = gradient(B_pol_smooth, dx);
        dB_dz = calc_vertical_derivative_1D(B_pol_smooth, dx);
        sum_grad_sq = sum_grad_sq + dB_dx.^2 + dB_dz.^2;
        sum_Bpri_sq = sum_Bpri_sq + B_pri_amp.^2;
    end
    AS_vector = sqrt(sum_grad_sq);
    B_pri_vector = sqrt(sum_Bpri_sq);
    epsilon_base = 2 * mean(AS_vector(bg_mask));
    epsilon = max(epsilon_base, 0.05 * max(B_pri_vector)); 
        NPMAG_total(f_i, :) = (AS_vector ./ (B_pri_vector + epsilon)) * 1e6;

end

%% 4. 绘制折线图
figure('Name', 'Refined NPMAG Profile', 'Color', 'w', 'Position', [150, 150, 1100, 550]);
fontN =  'Times New Roman'; fontS = 16; 
colors = lines(num_freqs); 

hold on; box on; grid on;

max_val = max(NPMAG_total(:, :), [], 'all'); 
yl = [0, max_val * 1.25]; 

p1 = patch([14 28 28 14], [yl(1) yl(1) yl(2) yl(2)], [0.93 0.93 0.93], 'FaceAlpha', 0.6, 'EdgeColor', 'none');
p2 = patch([50 84 84 50], [yl(1) yl(1) yl(2) yl(2)], [1 0.88 0.88], 'FaceAlpha', 0.6, 'EdgeColor', 'none'); 

text(20, yl(2)*0.94, ['Goaf'], 'HorizontalAlignment', 'center', ...
     'FontName', fontN, 'FontSize', fontS, 'FontWeight', 'bold', 'Color', [0.3 0.3 0.3]);
text(65, yl(2)*0.94, 'Known Anomaly Zone', 'HorizontalAlignment', 'center', ...
     'FontName', fontN, 'FontSize', fontS, 'FontWeight', 'bold', 'Color', [0.7 0 0]);

h_lines = zeros(1, num_freqs);
for f_i = 1:num_freqs
if  f_i == 5
    NPMAG_total(f_i, :)=NPMAG_total(f_i, :);
end
    h_lines(f_i) = plot(dist_vec, NPMAG_total(f_i, :), '-o', 'LineWidth', 2.8, 'MarkerSize', 9, ...
                        'Color', colors(f_i,:), 'MarkerFaceColor', colors(f_i,:));
end

set(gca, 'FontName', fontN, 'FontSize', fontS, 'LineWidth', 1.8, 'TickDir', 'in');
xlim([0, dist_vec(end)]); xticks(0:dx:dist_vec(end)); ylim(yl);

xlabel('Distance (m)', 'FontName', fontN, 'FontSize', fontS+2, 'FontWeight', 'bold');
ylabel('Bx NPMAG Operator (ppm/m)', 'FontName', fontN, 'FontSize', fontS+2, 'FontWeight', 'bold');

legend(h_lines, freq_labels, 'Location', 'northwest', 'FontName', fontN, 'FontSize', fontS, ...
       'Box', 'on', 'LineWidth', 1.2, 'EdgeColor', [0.5 0.5 0.5]);


function dz = calc_vertical_derivative_1D(data, dx)
    N = length(data); if N < 3, dz = zeros(1, N); return; end
    data_detrend = data - mean(data);
    data_pad = [fliplr(data_detrend), data_detrend, fliplr(data_detrend)];
    N_pad = length(data_pad);
    win = ones(1, N_pad); t_len = floor(N_pad * 0.1);
    taper = sin(linspace(0, pi/2, t_len)).^2;
    win(1:t_len) = taper; win(end-t_len+1:end) = fliplr(taper);
    data_taper = data_pad .* win;
    F_data = fft(data_taper);
    k = (2*pi / (N_pad * dx)) * [0:(ceil(N_pad/2)-1), -floor(N_pad/2):-1];
    F_dz = F_data .* abs(k);
    dz_pad = real(ifft(F_dz));
    dz = dz_pad(N+1 : 2*N);
end