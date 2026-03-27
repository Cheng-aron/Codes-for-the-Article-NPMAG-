
clc; 
clear; 
close all;


data_filename = 'SimData_Selectivity_2D_FFT.mat';
load(data_filename);



idx_low = find(freq_list == 1, 1);
idx_high = find(freq_list == 4, 1);
dx = x_map_vec(2) - x_map_vec(1);
dy = y_map_vec(2) - y_map_vec(1);

B_pri_mod_2D = sqrt(abs(MapData.Pri.Bx{idx_low}).^2 + abs(MapData.Pri.By{idx_low}).^2 + abs(MapData.Pri.Bz{idx_low}).^2);
Abs_By_surf_2D = abs(MapData.Tot.By{idx_low} - MapData.Pri.By{idx_low});

sigma_smooth = 2; 
S_By_pol = imgaussfilt(Abs_By_surf_2D, sigma_smooth);

[Grad_X_2D, Grad_Y_2D] = gradient(S_By_pol, dx, dy);
[Ny, Nx] = size(S_By_pol);
kx = 2*pi * (0:Nx-1) / (Nx*dx);  kx(kx > pi/dx) = kx(kx > pi/dx) - 2*pi/dx;
ky = 2*pi * (0:Ny-1) / (Ny*dy);  ky(ky > pi/dy) = ky(ky > pi/dy) - 2*pi/dy;

[KX, KY] = meshgrid(kx, ky);
K_radial = sqrt(KX.^2 + KY.^2); 
F_dz = fft2(S_By_pol) .* K_radial;         
Grad_Z_2D = real(ifft2(F_dz));       

eps0 = 0.05 * max(B_pri_mod_2D(:)); 
NPMAG_2D = sqrt(Grad_X_2D.^2 + Grad_Y_2D.^2 + Grad_Z_2D.^2) ./ (B_pri_mod_2D + eps0);


trim_margin_x = round(30 / dx); 
trim_margin_y = round(30 / dy); 
NPMAG_2D(1:trim_margin_y, :) = 0;              
NPMAG_2D(end-trim_margin_y+1:end, :) = 0;      
NPMAG_2D(:, 1:trim_margin_x) = 0;              
NPMAG_2D(:, end-trim_margin_x+1:end) = 0;      

Diff_Phase_Ex_2D = angle(MapData.Tot.Ex{idx_high} ./ MapData.Tot.Ex{idx_low});

[~, ref_x_idx] = min(abs(x_map_vec - (-250)));
[~, ref_y_idx] = min(abs(y_map_vec - 0));
Ref_By_Value = MapData.Tot.By{idx_low}(ref_y_idx, ref_x_idx);
dPhi_H_space_2D  = angle(MapData.Tot.By{idx_low} ./ Ref_By_Value);

Norm_Diff_Phase_Ex_2D = Diff_Phase_Ex_2D ./ max(abs(Diff_Phase_Ex_2D(:)));
Norm_dPhi_H_space_2D  = dPhi_H_space_2D  ./ max(abs(dPhi_H_space_2D(:)));
Norm_NPMAG_2D         = (NPMAG_2D ./ max(abs(NPMAG_2D(:)))).^2; 

[~, y0_idx] = min(abs(y_map_vec - 0));
x_line = x_map_vec;

Raw1D_Ex_Phase = Diff_Phase_Ex_2D(y0_idx, :);
Raw1D_H_Space  = dPhi_H_space_2D(y0_idx, :);
Raw1D_NPMAG    = NPMAG_2D(y0_idx, :);

Norm_1D_Ex_Phase = Raw1D_Ex_Phase ./ max(abs(Raw1D_Ex_Phase));
Norm_1D_H_Space  = Raw1D_H_Space  ./ max(abs(Raw1D_H_Space));
max_np = max(abs(Raw1D_NPMAG)); if max_np == 0; max_np = 1; end 
Norm_1D_NPMAG = (Raw1D_NPMAG ./ max_np).^2;


fwhm_vals = zeros(1,3); 
true_sci_vals = zeros(1,3); 
false_fpr_vals = zeros(1,3);

data_1d = {abs(Norm_1D_Ex_Phase), abs(Norm_1D_H_Space), abs(Norm_1D_NPMAG)};
param_names = {'Elect. Phase Diff', 'Mag. Spatial Phase', 'Proposed NPMAG'};


for k = 1:3
    y_tmp = data_1d{k} / max(data_1d{k});
    
    idx_left = find(x_line < 0);
    x_lp = x_line(idx_left); y_lp = y_tmp(idx_left);
    y_lp = y_lp / max(y_lp); 
    
    idx_half = find(y_lp >= 0.5);
    if ~isempty(idx_half)
        fwhm_vals(k) = x_lp(idx_half(end)) - x_lp(idx_half(1));
    else
        fwhm_vals(k) = NaN;
    end
    
    mask_true = (x_line >= -100) & (x_line <= -50);
    area_true = trapz(x_line(mask_true), y_tmp(mask_true));
    
    mask_false = (x_line >= 50) & (x_line <= 100);
    area_false = trapz(x_line(mask_false), y_tmp(mask_false));
    
    area_total = trapz(x_line, y_tmp);
    
    true_sci_vals(k) = (area_true / area_total) * 100;
    false_fpr_vals(k) = (area_false / area_total) * 100; 
    
    fprintf('     %-20s | Left FWHM = %5.1f m | True SCI = %5.1f %% | False FPR = %5.1f %%\n', ...
        param_names{k}, fwhm_vals(k), true_sci_vals(k), false_fpr_vals(k));
end


fs_tick = 13;  
fs_label = 15; 
font_name = 'Times New Roman'; 


figure('Name', '2D Normalized Selectivity Map', 'Color', 'w', 'Position',[100, 50, 1400, 450]); 
t2 = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
colormap(jet);


bx_true =[-100, -50, -50, -100, -100]; by_true =[-50, -50, 50, 50, -50];
bx_false =[50, 100, 100, 50, 50];      by_false =[-50, -50, 50, 50, -50];

data_2d = {Norm_Diff_Phase_Ex_2D, Norm_dPhi_H_space_2D, Norm_NPMAG_2D};
titles_2d = {'(a)', '(b)', '(c)'};

for k = 1:3
    nexttile;
    pcolor(X_grid, Y_grid, data_2d{k});
    shading interp; caxis([-1, 1]); hold on;


    plot3(bx_true, by_true, ones(size(bx_true)), 'w-', 'LineWidth', 4);   
    plot3(bx_true, by_true, ones(size(bx_true)), 'r--', 'LineWidth', 2);   

    plot3(bx_false, by_false, ones(size(bx_true)), 'w-', 'LineWidth', 4); 
    plot3(bx_false, by_false, ones(size(bx_true)), 'b--', 'LineWidth', 2); 
    

    axis equal;     
    axis tight;     
    box on;         
    
    set(gca, 'LineWidth', 1, 'FontSize', fs_tick, 'FontName', font_name, 'TickDir', 'in'); 
    set(gca, 'XTick',[-200 -100 0 100 200]); 
    set(gca, 'YTick',[-200 -100 0 100 200]);
    
    if k == 1
        ylabel('y (m)', 'FontWeight', 'bold', 'FontSize', fs_label, 'FontName', font_name);
    else
        yticklabels({});
    end

    xlabel({'x (m)', titles_2d{k}}, 'FontWeight', 'bold', 'FontSize', fs_label, 'FontName', font_name);
end

cb = colorbar;
cb.Layout.Tile = 'east'; 
cb.FontSize = fs_tick;
cb.FontName = font_name;

ylabel(cb, 'Normalized Response', 'FontWeight', 'bold', 'FontSize', fs_label, 'FontName', font_name);
