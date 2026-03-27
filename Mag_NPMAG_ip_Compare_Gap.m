clc; 
clear; 
close all;

Gap_List =[20, 50, 100, 150, 200]; 

W_BLK = 50;      
x_range = 250;   
grid_step = 2.5; 

data_filename = 'SimData_Exp2_Gap_Scan_2D_FFT.mat';
load(data_filename);


Results = struct();
idx_low = 1; idx_high = length(freq_list);
dx = x_map_vec(2) - x_map_vec(1);
dy = y_map_vec(2) - y_map_vec(1);

[~, y0_idx] = min(abs(y_map_vec - 0));
x_line = x_map_vec;[~, ref_x_idx] = min(abs(x_map_vec - (-x_range)));

for k = 1:length(Data_Store)
    curr_gap = Data_Store(k).Gap;
    PriMap = Data_Store(k).PriMap;
    TotMap = Data_Store(k).TotMap;
    
    
    B_pri_mod_2D = sqrt(abs(PriMap.Bx{idx_low}).^2 + abs(PriMap.By{idx_low}).^2 + abs(PriMap.Bz{idx_low}).^2);
    Abs_By_surf_2D = abs(TotMap.By{idx_low} - PriMap.By{idx_low});
    S_By_pol = imgaussfilt(Abs_By_surf_2D, 2);
    [Grad_X_2D, Grad_Y_2D] = gradient(S_By_pol, dx, dy);
    
    [Ny, Nx] = size(S_By_pol);
    kx = 2*pi * (0:Nx-1) / (Nx*dx);  kx(kx > pi/dx) = kx(kx > pi/dx) - 2*pi/dx;
    ky = 2*pi * (0:Ny-1) / (Ny*dy);  ky(ky > pi/dy) = ky(ky > pi/dy) - 2*pi/dy;
    [KX, KY] = meshgrid(kx, ky);
    F_dz = fft2(S_By_pol) .* sqrt(KX.^2 + KY.^2);         
    Grad_Z_2D = real(ifft2(F_dz));       
    
  
    eps0 = 0.05 * max(B_pri_mod_2D(:)); 
    NPMAG_2D = sqrt(Grad_X_2D.^2 + Grad_Y_2D.^2 + Grad_Z_2D.^2) ./ (B_pri_mod_2D + eps0);
    Diff_Phase_Ex_2D = angle(TotMap.Ex{idx_high} ./ TotMap.Ex{idx_low});
    dPhi_H_space_2D  = angle(TotMap.By{idx_low} ./ TotMap.By{idx_low}(y0_idx, ref_x_idx));
    

    Norm2D_Ex = Diff_Phase_Ex_2D ./ max(abs(Diff_Phase_Ex_2D(:)));
    Norm2D_H_Space = dPhi_H_space_2D ./ max(abs(dPhi_H_space_2D(:)));
    Norm2D_NPMAG = (NPMAG_2D ./ max(abs(NPMAG_2D(:)))).^2; 
    

    Norm1D_Ex = Diff_Phase_Ex_2D(y0_idx, :) ./ max(abs(Diff_Phase_Ex_2D(y0_idx, :)));
    Norm1D_H  = dPhi_H_space_2D(y0_idx, :)  ./ max(abs(dPhi_H_space_2D(y0_idx, :)));
    Norm1D_NP = (NPMAG_2D(y0_idx, :) ./ max(abs(NPMAG_2D(y0_idx, :)))).^2;
 
    offset = (curr_gap + W_BLK)/2;
    left_b1 = -offset - W_BLK/2;  right_b1 = -offset + W_BLK/2;
    left_b2 =  offset - W_BLK/2;  right_b2 =  offset + W_BLK/2;
    
    fwhm_vals = zeros(1,3); sci_vals = zeros(1,3);
    data_1d_abs = {abs(Norm1D_Ex), abs(Norm1D_H), abs(Norm1D_NP)};
    param_names = {'Elect. Phase Diff', 'Mag. Spatial Phase', 'Proposed NPMAG'};
    
    for m = 1:3
        y_tmp = data_1d_abs{m} / max(data_1d_abs{m});
        
        idx_right = find(x_line > 0);
        x_rp = x_line(idx_right); y_rp = y_tmp(idx_right);
        y_rp = y_rp / max(y_rp); 
        idx_half = find(y_rp >= 0.5);
        if ~isempty(idx_half)
            fwhm_vals(m) = x_rp(idx_half(end)) - x_rp(idx_half(1));
        else
            fwhm_vals(m) = NaN;
        end
        
        mask_body1 = (x_line >= left_b1) & (x_line <= right_b1);
        mask_body2 = (x_line >= left_b2) & (x_line <= right_b2);
        area_in = trapz(x_line(mask_body1), y_tmp(mask_body1)) + trapz(x_line(mask_body2), y_tmp(mask_body2));
        sci_vals(m) = (area_in / trapz(x_line, y_tmp)) * 100;
        
        fprintf('     %-20s | FWHM = %5.1f m | SCI = %5.1f %%\n', param_names{m}, fwhm_vals(m), sci_vals(m));
    end
    
    Results(k).Norm2D = {Norm2D_Ex, Norm2D_H_Space, Norm2D_NPMAG};
    Results(k).Norm1D = {Norm1D_Ex, Norm1D_H, Norm1D_NP};
    Results(k).Metrics = struct('FWHM', fwhm_vals, 'SCI', sci_vals);
    Results(k).Bounds =[ left_b1, right_b1, left_b2, right_b2 ];
end

leg_str = arrayfun(@(x) sprintf('Gap = %d m', x), Gap_List, 'UniformOutput', false);


fs_tick = 13;  
fs_label = 15; 
font_name = 'Times New Roman'; 
titles_2d = {'(a)', '(b)', '(c)'};

N_gaps = length(Gap_List);
W_Y = 50; 


figure('Name', 'All Gaps 2D Maps (Merged View)', 'Color', 'w', 'Position', [100, 100, 1000, 1200]);
t2 = tiledlayout(N_gaps, 3, 'TileSpacing', 'compact', 'Padding', 'compact'); 
colormap(jet);

for row = 1:N_gaps
    curr_gap = Gap_List(row);
    P2D = Results(row).Norm2D;
    bnds = Results(row).Bounds; 
    
    bx1 =[bnds(1), bnds(2), bnds(2), bnds(1), bnds(1)]; by1 =[-W_Y, -W_Y, W_Y, W_Y, -W_Y];
    bx2 =[bnds(3), bnds(4), bnds(4), bnds(3), bnds(3)]; by2 =[-W_Y, -W_Y, W_Y, W_Y, -W_Y];
    
    for col = 1:3
        nexttile;
        pcolor(X_grid, Y_grid, P2D{col});
        shading interp; caxis([-1, 1]); hold on;
        
 
        plot(bx1, by1, 'w-', 'LineWidth', 3); 
        plot(bx2, by2,  'w-', 'LineWidth', 3);
        plot(bx1, by1, 'r--', 'LineWidth', 1.5); 
        plot(bx2, by2, 'r--', 'LineWidth', 1.5);

        axis equal; axis tight; box on;
        set(gca, 'LineWidth', 1, 'FontSize', fs_tick, 'FontName', font_name, 'TickDir', 'in'); 
        set(gca, 'XTick',[-200 -100 0 100 200]);
        set(gca, 'YTick',[-200 -100 0 100 200]);
        
       
        if row ~= N_gaps, xticklabels({}); end
        if col ~= 1, yticklabels({}); end
        
        if col == 1

             ylabel({sprintf('Gap = %d m', curr_gap), 'y (m)'}, ...
                    'FontWeight', 'bold', 'FontSize', fs_label, 'FontName', font_name);
        end
        
        if row == N_gaps
            xlabel({'x (m)', titles_2d{col}}, 'FontWeight', 'bold', 'FontSize', fs_label, 'FontName', font_name);
        end
    end
end


cb = colorbar;
cb.Layout.Tile = 'east'; 
cb.FontSize = fs_tick;
cb.FontName = font_name;
ylabel(cb, 'Normalized Response', 'FontWeight', 'bold', 'FontSize', fs_label, 'FontName', font_name);