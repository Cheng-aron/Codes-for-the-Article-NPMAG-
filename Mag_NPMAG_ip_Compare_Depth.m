clc; 
clear; 
close all;

target_depth_list =[50, 100, 150, 200, 250,300];
freq_list =[1, 4]; 
x_range = 250; 
grid_step = 5;
data_filename = 'SimData_Depth_Sweep_2D_FFT.mat';
load(data_filename, 'Data_Store', 'x_map_vec', 'y_map_vec', 'X_grid', 'Y_grid', 'freq_list');        
depths_in_data = [Data_Store.Depth];
keep_idx = ismember(depths_in_data, target_depth_list);
Data_Store = Data_Store(keep_idx);

depth_list = target_depth_list; 
        
Results = struct();
idx_low = 1; idx_high = length(freq_list);
dx = x_map_vec(2) - x_map_vec(1);
dy = y_map_vec(2) - y_map_vec(1);

[~, y0_idx] = min(abs(y_map_vec - 0));
x_line = x_map_vec;
[~, ref_x_idx] = min(abs(x_map_vec - (-x_range)));

for k = 1:length(Data_Store)
    PriMap = Data_Store(k).Pri.Map;
    TotMap = Data_Store(k).Tot.Map;
    
    B_pri_mod_2D = sqrt(abs(PriMap.Bx{idx_low}).^2 + abs(PriMap.By{idx_low}).^2 + abs(PriMap.Bz{idx_low}).^2);
    Abs_By_surf_2D = abs(TotMap.By{idx_low} - PriMap.By{idx_low});
  
    S_By_pol = imgaussfilt(Abs_By_surf_2D, 2);
    
 
    [Grad_X_2D, Grad_Y_2D] = gradient(S_By_pol, dx, dy);
    [Ny, Nx] = size(S_By_pol);
    kx = 2*pi * (0:Nx-1) / (Nx*dx);  kx(kx > pi/dx) = kx(kx > pi/dx) - 2*pi/dx;
    ky = 2*pi * (0:Ny-1) / (Ny*dy);  ky(ky > pi/dy) = ky(ky > pi/dy) - 2*pi/dy;
    
    [KX, KY] = meshgrid(kx, ky);
    K_radial = sqrt(KX.^2 + KY.^2); 
    F_dz = fft2(S_By_pol) .* K_radial;         
    Grad_Z_2D = real(ifft2(F_dz));       

    eps0 = 0.03 * max(B_pri_mod_2D(:)); 
    NPMAG_2D = sqrt(Grad_X_2D.^2 + Grad_Y_2D.^2 + Grad_Z_2D.^2) ./ (B_pri_mod_2D + eps0);
    
    Diff_Phase_Ex_2D = angle(TotMap.Ex{idx_high} ./ TotMap.Ex{idx_low});
    Ref_By_Value = TotMap.By{idx_low}(y0_idx, ref_x_idx);
    dPhi_H_space_2D  = angle(TotMap.By{idx_low} ./ Ref_By_Value);
    

    Norm2D_Ex = Diff_Phase_Ex_2D ./ max(abs(Diff_Phase_Ex_2D(:)));
    Norm2D_H_Space = dPhi_H_space_2D ./ max(abs(dPhi_H_space_2D(:)));
    Norm2D_NPMAG = (NPMAG_2D ./ max(abs(NPMAG_2D(:)))).^2; 

    Raw1D_Ex = Diff_Phase_Ex_2D(y0_idx, :);
    Raw1D_H  = dPhi_H_space_2D(y0_idx, :);
    Raw1D_NP = NPMAG_2D(y0_idx, :);
    
    Norm1D_Ex = Raw1D_Ex ./ max(abs(Raw1D_Ex));
    Norm1D_H  = Raw1D_H  ./ max(abs(Raw1D_H));
    Norm1D_NP = (Raw1D_NP ./ max(abs(Raw1D_NP))).^2;

    data_1d_abs = {abs(Norm1D_Ex), abs(Norm1D_H), abs(Norm1D_NP)};
    for m = 1:3
        y_tmp = data_1d_abs{m} / max(data_1d_abs{m});
        idx_half = find(y_tmp >= 0.5);
        mask_in = (x_line >= -50) & (x_line <= 50);
        sci_vals(m) = (trapz(x_line(mask_in), y_tmp(mask_in)) / trapz(x_line, y_tmp)) * 100;
    end
    
    Results(k).Norm2D = {Norm2D_Ex, Norm2D_H_Space, Norm2D_NPMAG};
    Results(k).Norm1D = {Norm1D_Ex, Norm1D_H, Norm1D_NP};
end



titles_1d = {'(a)', '(b)', '(c)'};
x_b1 = -50; x_b2 = 50; 

show_indices_m2 = [1, 2,3]; 
show_indices_m3 = [1, 2,3, 4,5];    

figure('Name', 'Type 2: 1D Layout (1x3)', 'Color', 'w', 'Position',[50, 450, 1600, 400]);
colors = jet(length(depth_list));


for m = 1:3
    subplot(1, 3, m); hold on; 
    
 
    for k = 1:length(depth_list)
     
        should_plot = false;
        if m == 1
            should_plot = true; 
        elseif m == 2 && ismember(k, show_indices_m2)
            should_plot = true;
        elseif m == 3 && ismember(k, show_indices_m3)
            should_plot = true;
        end
        
        if should_plot
            P1D = Results(k).Norm1D;
            c = colors(k, :);
            y_data = P1D{m};
            
         
            y_data = smoothdata(y_data, 'gaussian', 5); 
            
       
            plot(x_line, y_data, 'Color', c, 'LineWidth', 1.5, ...
                 'DisplayName', sprintf('Depth = %d m', depth_list(k)));
        end
    end
    

    grid on; box on; 
    set(gca, 'FontSize', 13, 'FontName', 'Times New Roman');
    ylabel('Normalized Response', 'FontWeight', 'bold'); 
    xlabel({'x (m)', titles_1d{m}}, 'FontWeight', 'bold'); 
    xlim([-200, 200]); ylim([-1.1, 1.1]);
    
    xline(x_b1, '--r', 'LineWidth', 2, 'HandleVisibility', 'off'); 
    xline(x_b2, '--r', 'LineWidth', 2, 'HandleVisibility', 'off');
    
    legend('show', 'Location', 'Best', 'FontSize', 11);
end
box_x =[-50, 50, 50, -50, -50];
box_y =[-50, -50, 50, 50, -50];
num_depths = length(depth_list);
