clc; 
clear; 
close all;

data_filename = 'SimData_2D_Surface_FFT.mat';
load(data_filename);


idx_low = find(freq_list == 1, 1);
idx_high = find(freq_list == 4, 1);
dx = x_map_vec(2) - x_map_vec(1);
dy = y_map_vec(2) - y_map_vec(1);


B_pri_mod_2D = sqrt(abs(MapData.Pri.Bx{idx_low}).^2 + abs(MapData.Pri.By{idx_low}).^2 + abs(MapData.Pri.Bz{idx_low}).^2);

By_pol_surf_2D = MapData.Tot.By{idx_low} - MapData.Pri.By{idx_low};
Abs_By_surf_2D = abs(By_pol_surf_2D);

sigma_smooth = 3; 
S_By_pol = imgaussfilt(Abs_By_surf_2D, sigma_smooth);

[Grad_X_2D, Grad_Y_2D] = gradient(S_By_pol, dx, dy);

[Ny, Nx] = size(S_By_pol);

kx = 2*pi * (0:Nx-1) / (Nx*dx);  
kx(kx > pi/dx) = kx(kx > pi/dx) - 2*pi/dx;

ky = 2*pi * (0:Ny-1) / (Ny*dy);  
ky(ky > pi/dy) = ky(ky > pi/dy) - 2*pi/dy;

[KX, KY] = meshgrid(kx, ky);


K_radial = sqrt(KX.^2 + KY.^2); 

F_Bpol = fft2(S_By_pol);            
F_dz   = F_Bpol .* K_radial;        
Grad_Z_2D = real(ifft2(F_dz));      

eps0 = 1e-18; 
NPMAG_2D = sqrt(Grad_X_2D.^2 + Grad_Y_2D.^2 + Grad_Z_2D.^2) ./ (B_pri_mod_2D + eps0);

Diff_Phase_Ex_2D = angle(MapData.Tot.Ex{idx_high} ./ MapData.Tot.Ex{idx_low});
Diff_Phase_By_2D = angle(MapData.Tot.By{idx_high} ./ MapData.Tot.By{idx_low});


[~, ref_x_idx] = min(abs(x_map_vec - (-250)));
[~, ref_y_idx] = min(abs(y_map_vec - (-250)));
Ref_By_Value = MapData.Tot.By{idx_low}(ref_y_idx, ref_x_idx);
dPhi_H_space_2D  = angle(MapData.Tot.By{idx_low} ./ Ref_By_Value) * 1000;



[~, y0_idx] = min(abs(y_map_vec - 0));
x_line = x_map_vec;
NPMAG_1D         = NPMAG_2D(y0_idx, :);
Diff_Phase_Ex_1D = Diff_Phase_Ex_2D(y0_idx, :);
Diff_Phase_By_1D = Diff_Phase_By_2D(y0_idx, :);
dPhi_H_space_1D  = dPhi_H_space_2D(y0_idx, :);


figure('Name', 'Profile Analysis (Line y=0)', 'Color', 'w', 'Position',[50, 50, 1200, 800]);
x_b1 = -50; x_b2 = 50; 

subplot(2, 2, 1);
plot(x_line, Diff_Phase_Ex_1D * 1000, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 4); grid on;
xline(x_b1,'--r'); xline(x_b2,'--r'); title('(a) Electrical Phase Diff (Ex)');
xlabel('x (m)'); ylabel('Phase (mrad)', 'FontWeight', 'bold');xlim([-200, 200]);

subplot(2, 2, 2);
plot(x_line, Diff_Phase_By_1D * 1000, 'k-s', 'LineWidth', 1.5, 'MarkerSize', 4); grid on;
xline(x_b1,'--r'); xline(x_b2,'--r'); title('(b) Magnetic Phase Diff (By)');
xlabel('x (m)'); ylabel('Phase (mrad)', 'FontWeight', 'bold');xlim([-200, 200]);

subplot(2, 2, 3);
plot(x_line, dPhi_H_space_1D, 'm-^', 'LineWidth', 1.5, 'MarkerFaceColor', 'm', 'MarkerSize', 4); grid on;
xline(x_b1,'--r'); xline(x_b2,'--r'); title('(c) Magnetic Spatial Phase (By)');
xlabel('x (m)'); ylabel('Phase (mrad)', 'FontWeight', 'bold');xlim([-200, 200]);

subplot(2, 2, 4);
NPMAG_scale = 1e6;  
area(x_line, NPMAG_1D * NPMAG_scale, 'FaceColor',[1 0.2 0.2], 'FaceAlpha', 0.5, 'EdgeColor', 'r', 'LineWidth', 2); hold on;
xline(x_b1,'--k'); xline(x_b2,'--k'); title('(d) NPMAG (FFT Reconstructed)');
xlabel('x (m)'); ylabel('NPMAG (ppm/m)', 'FontWeight', 'bold'); xlim([-200, 200]);

sgtitle('1D Profile Analysis Extract from 2D Plane (y=0)', 'FontSize', 16, 'FontWeight', 'bold');


figure('Name', '2D Normalized Core Parameters Map (1x4 Layout)', 'Color', 'w', 'Position',[50, 50, 1200, 800]);

fs_tick = 13;  
fs_label = 15; 
font_name = 'Times New Roman'; 

box_x =[-50, 50, 50, -50, -50];
box_y =[-50, -50, 50, 50, -50];

ax1 = subplot(2, 2, 1);
pcolor(X_grid, Y_grid, Diff_Phase_Ex_2D * 1000); 
shading interp; 
axis equal; 
xlim([-200, 200]); ylim([-200, 200]); 
colormap(ax1, jet); hold on;


set(ax1, 'FontSize', fs_tick, 'FontName', font_name); 

cb1 = colorbar; 
cb1.FontSize = fs_tick; 
cb1.FontName = font_name; 
ylabel(cb1, 'EPD mrad', 'FontWeight', 'bold', 'FontSize', fs_label, 'FontName', font_name); 

plot3(box_x, box_y, ones(size(box_x)), 'w-', 'LineWidth', 4); 
plot3(box_x, box_y, ones(size(box_x)), 'r--', 'LineWidth', 2);
xlabel({'x (m)', '(a)'}, 'FontWeight', 'bold', 'FontSize', fs_label, 'FontName', font_name); 
ylabel('y (m)', 'FontWeight', 'bold', 'FontSize', fs_label, 'FontName', font_name);



ax2 = subplot(2, 2, 2);
pcolor(X_grid, Y_grid, Diff_Phase_By_2D * 1000); 
shading interp; 
axis equal; 
xlim([-200, 200]); ylim([-200, 200]); 
colormap(ax2, jet); hold on;


set(ax2, 'FontSize', fs_tick, 'FontName', font_name);

cb2 = colorbar; 
cb2.FontSize = fs_tick;
cb2.FontName = font_name;
ylabel(cb2, 'MPD mrad', 'FontWeight', 'bold', 'FontSize', fs_label, 'FontName', font_name);

plot3(box_x, box_y, ones(size(box_x)), 'w-', 'LineWidth', 4);
plot3(box_x, box_y, ones(size(box_x)), 'r--', 'LineWidth', 2); 
xlabel({'x (m)', '(b)'}, 'FontWeight', 'bold', 'FontSize', fs_label, 'FontName', font_name); 
ylabel('y (m)', 'FontWeight', 'bold', 'FontSize', fs_label, 'FontName', font_name);

drawnow; 
cb2_pos = cb2.Position; 
rect_x = cb2_pos(1) + cb2_pos(3) + 0.002; 
rect_y = cb2_pos(2)+ 0.01;
rect_w = 0.03; 
rect_h = cb2_pos(4);
annotation('rectangle',[rect_x, rect_y, rect_w, rect_h], ...
    'EdgeColor',[0.9, 0.1, 0.5], 'LineWidth', 2.5);


ax3 = subplot(2, 2, 3);
pcolor(X_grid, Y_grid, dPhi_H_space_2D);
shading interp; 
axis equal; 
xlim([-200, 200]); ylim([-200, 200]); 
colormap(ax3, jet); hold on;

set(ax3, 'FontSize', fs_tick, 'FontName', font_name);

cb3 = colorbar; 
cb3.FontSize = fs_tick;
cb3.FontName = font_name;
ylabel(cb3, 'MSPD mrad', 'FontWeight', 'bold', 'FontSize', fs_label, 'FontName', font_name);

plot3(box_x, box_y, ones(size(box_x)), 'w-', 'LineWidth', 4);
plot3(box_x, box_y, ones(size(box_x)), 'r--', 'LineWidth', 2);
xlabel({'x (m)', '(c)'}, 'FontWeight', 'bold', 'FontSize', fs_label, 'FontName', font_name); 
ylabel('y (m)', 'FontWeight', 'bold', 'FontSize', fs_label, 'FontName', font_name);


ax4 = subplot(2, 2, 4);
pcolor(X_grid, Y_grid, NPMAG_2D * NPMAG_scale);
shading interp; 
axis equal; 
xlim([-200, 200]); ylim([-200, 200]); 
colormap(ax4, jet); hold on;

set(ax4, 'FontSize', fs_tick, 'FontName', font_name);

cb4 = colorbar; 
cb4.FontSize = fs_tick;
cb4.FontName = font_name;
ylabel(cb4, 'NPMAG ppm/m', 'FontWeight', 'bold', 'FontSize', fs_label, 'FontName', font_name); 

plot3(box_x, box_y, ones(size(box_x)), 'w-', 'LineWidth', 4);
plot3(box_x, box_y, ones(size(box_x)), 'r--', 'LineWidth', 2); 
xlabel({'x (m)', '(d)'}, 'FontWeight', 'bold', 'FontSize', fs_label, 'FontName', font_name); 
ylabel('y (m)', 'FontWeight', 'bold', 'FontSize', fs_label, 'FontName', font_name);

figure('Name', '2D Map: Pure Secondary Fields', 'Color', 'w', 'Position',[150, 150, 1400, 900]);

unit_E = 1e6;  % uV/m
unit_B = 1e12; % pT
map_idx = idx_low;

ax11 = subplot(2, 3, 1);
pcolor(X_grid, Y_grid, abs(MapData.Tot.Ex{map_idx} - MapData.Pri.Ex{map_idx}) * unit_E);
shading interp; axis equal tight; colorbar; colormap(ax11, jet); hold on;
plot(box_x, box_y, 'r--', 'LineWidth', 2); 
title('Secondary Ex (uV/m)'); xlabel('x (m)'); ylabel('y (m)');

ax12 = subplot(2, 3, 2);
pcolor(X_grid, Y_grid, abs(MapData.Tot.Ey{map_idx} - MapData.Pri.Ey{map_idx}) * unit_E);
shading interp; axis equal tight; colorbar; colormap(ax12, jet); hold on;
plot(box_x, box_y, 'r--', 'LineWidth', 2); 
title('Secondary Ey (uV/m)'); xlabel('x (m)'); ylabel('y (m)');

subplot(2, 3, 3); axis off;
text(0.05, 0.6, {'Secondary Fields Info:', ...['Freq: ', num2str(freq_list(map_idx)), ' Hz'], ...
    'Plane: z = 0 (Ground)'}, 'FontSize', 12, 'BackgroundColor',[0.95 0.95 0.95]);

ax14 = subplot(2, 3, 4);
pcolor(X_grid, Y_grid, abs(MapData.Tot.Bx{map_idx} - MapData.Pri.Bx{map_idx}) * unit_B);
shading interp; axis equal tight; colorbar; colormap(ax14, jet); hold on;
plot(box_x, box_y, 'r--', 'LineWidth', 2); 
title('Secondary Bx (pT)'); xlabel('x (m)'); ylabel('y (m)');

ax15 = subplot(2, 3, 5);
pcolor(X_grid, Y_grid, abs(MapData.Tot.By{map_idx} - MapData.Pri.By{map_idx}) * unit_B);
shading interp; axis equal tight; colorbar; colormap(ax15, jet); hold on;
plot(box_x, box_y, 'r--', 'LineWidth', 2); 
title('Secondary By (pT)'); xlabel('x (m)'); ylabel('y (m)');

ax16 = subplot(2, 3, 6);
pcolor(X_grid, Y_grid, abs(MapData.Tot.Bz{map_idx} - MapData.Pri.Bz{map_idx}) * unit_B);
shading interp; axis equal tight; colorbar; colormap(ax16, jet); hold on;
plot(box_x, box_y, 'r--', 'LineWidth', 2); 
title('Secondary Bz (pT)'); xlabel('x (m)'); ylabel('y (m)');

sgtitle('Ground Surface Secondary Field Maps', 'FontSize', 16, 'FontWeight', 'bold');