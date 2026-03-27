clc; 
clear; 
close all;

sigma_list =[0.00001, 0.0001, 0.001, 0.01];

freq_list = [1, 4];

x_range = 250; 
grid_step = 5;  
final_data_filename = 'SimData_Sigma_Sweep_2D_FFT.mat';

load(final_data_filename);  

Results = struct();
idx_low = 1; idx_high = length(freq_list);
dx = x_map_vec(2) - x_map_vec(1);
dy = y_map_vec(2) - y_map_vec(1);

[~, y0_idx] = min(abs(y_map_vec - 0));
x_line = x_map_vec;[~, ref_x_idx] = min(abs(x_map_vec - (-x_range)));

for k = 1:length(Data_Store)
    curr_sigma = Data_Store(k).Sigma;
    PriMap = Data_Store(k).PriMap;
    TotMap = Data_Store(k).TotMap;
    
  
    B_pri_mod_2D = sqrt(abs(PriMap.Bx{idx_low}).^2 + abs(PriMap.By{idx_low}).^2 + abs(PriMap.Bz{idx_low}).^2);
    Abs_By_surf_2D = abs(TotMap.By{idx_low} - PriMap.By{idx_low});
    

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
    
    Diff_Phase_Ex_2D = angle(TotMap.Ex{idx_high} ./ TotMap.Ex{idx_low});
    Ref_By_Value = TotMap.By{idx_low}(y0_idx, ref_x_idx);
    dPhi_H_space_2D  = angle(TotMap.By{idx_low} ./ Ref_By_Value);
    
    Norm2D_Ex = Diff_Phase_Ex_2D ./ max(abs(Diff_Phase_Ex_2D(:)));
    Norm2D_H_Space = dPhi_H_space_2D ./ max(abs(dPhi_H_space_2D(:)));
    Norm2D_NPMAG = (NPMAG_2D ./ max(abs(NPMAG_2D(:)))).^2; % 二阶锐化
    
    Raw1D_Ex = Diff_Phase_Ex_2D(y0_idx, :);
    Raw1D_H  = dPhi_H_space_2D(y0_idx, :);
    Raw1D_NP = NPMAG_2D(y0_idx, :);
    
    Norm1D_Ex = Raw1D_Ex ./ max(abs(Raw1D_Ex));
    Norm1D_H  = Raw1D_H  ./ max(abs(Raw1D_H));
    Norm1D_NP = (Raw1D_NP ./ max(abs(Raw1D_NP))).^2;

    fwhm_vals = zeros(1,3); sci_vals = zeros(1,3);
    data_1d_abs = {abs(Norm1D_Ex), abs(Norm1D_H), abs(Norm1D_NP)};
    param_names = {'Elect. Phase Diff', 'Mag. Spatial Phase', 'Proposed NPMAG'};
    
    
    for m = 1:3
        y_tmp = data_1d_abs{m} / max(data_1d_abs{m});
        
    
        idx_half = find(y_tmp >= 0.5);
        if ~isempty(idx_half)
            fwhm_vals(m) = x_line(idx_half(end)) - x_line(idx_half(1)); 
        else
            fwhm_vals(m) = NaN; 
        end
        

        mask_in = (x_line >= -50) & (x_line <= 50);
        sci_vals(m) = (trapz(x_line(mask_in), y_tmp(mask_in)) / trapz(x_line, y_tmp)) * 100;
        
        fprintf('     %-20s | FWHM = %5.1f m | SCI = %5.1f %%\n', param_names{m}, fwhm_vals(m), sci_vals(m));
    end
    
    Results(k).Norm2D = {Norm2D_Ex, Norm2D_H_Space, Norm2D_NPMAG};
    Results(k).Norm1D = {Norm1D_Ex, Norm1D_H, Norm1D_NP};
    Results(k).Metrics = struct('FWHM', fwhm_vals, 'SCI', sci_vals);
    Results(k).Raw2D = {PriMap, TotMap}; 
end
disp('======================================================');




fs_tick = 13;  
fs_label = 15; 
font_name = 'Times New Roman'; 

titles_2d = {'(a)', '(b)', '(c)'};


N_sigmas = length(sigma_list);

box_x =[-50, 50, 50, -50, -50];
box_y =[-50, -50, 50, 50, -50];


fig_height = min(180 * N_sigmas, 1600);
figure('Name', 'All Sigma 2D Maps (Merged View)', 'Color', 'w', 'Position',[100, 50, 800, fig_height]);


t2 = tiledlayout(N_sigmas, 3, 'TileSpacing', 'compact', 'Padding', 'compact'); 
colormap(jet);

for row = 1:N_sigmas
    curr_sig = sigma_list(row);
    P2D = Results(row).Norm2D;
    
    for col = 1:3
        nexttile;
        
      
        pcolor(X_grid, Y_grid, P2D{col});
        shading interp; caxis([-1, 1]); hold on;
        
        plot(box_x, box_y, 'w-', 'LineWidth', 3); 
        plot(box_x, box_y, 'r--', 'LineWidth', 1.5); 
    
        axis equal; axis tight; box on;
        set(gca, 'LineWidth', 1, 'FontSize', fs_tick, 'FontName', font_name, 'TickDir', 'in'); 
        set(gca, 'XTick',[-200 -100 0 100 200]); 
        set(gca, 'YTick',[-200 -100 0 100 200]);
        
    
        if row ~= N_sigmas 
            xticklabels({}); 
        end
        if col ~= 1
            yticklabels({}); 
        end
        

        if col == 1
             ylabel({['\sigma = ', num2str(curr_sig)],  'y (m)'}, ...
                    'FontWeight', 'bold', 'FontSize', fs_label, 'FontName', font_name);
        end
        
        if row == N_sigmas

            xlabel({'x (m)',  titles_2d{col}}, 'FontWeight', 'bold', 'FontSize', fs_label, 'FontName', font_name);
        end
    end
end

cb = colorbar;
cb.Layout.Tile = 'east'; 
cb.FontSize = fs_tick;
cb.FontName = font_name;

ylabel(cb, 'Normalized Response', 'FontWeight', 'bold', 'FontSize', fs_label, 'FontName', font_name);




