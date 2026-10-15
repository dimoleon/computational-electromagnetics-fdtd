clear;
clc;
close all;

addpath("./pml");
common;

n_pml = 8;
R_pml = 1e-6;
order_pml = 2;
courant_coeff = 0.9;
num_periods = 60;
freq_sim = 3.75e9;
omega = 2 * pi * freq_sim;
cells_per_lambda = 10;
beta_x = (2 * pi * freq_sim / c0) * sqrt(1 - (cutoff_freq / freq_sim)^2);
v_phase = omega / beta_x;

% A3
% grid to big, wave not reflected in time to cause interference
[ez0, times] = a34_sim(3.75e9, 30, 0.9, 150, 1, 2, 10, false);
% % proper parameters
[ez1, ~] = a34_sim(3.75e9, 30, 0.9, 8, 1e-6, 2, 10, false);
[ez1a, ~] = a34_sim(3.75e9, 30, 0.9, 8, 1e-6, 1, 10, false);
[ez1b, ~] = a34_sim(3.75e9, 30, 0.9, 8, 1e-6, 3, 10, false);
[ez2, ~] = a34_sim(3.75e9, 30, 0.9, 16, 1e-6, 2, 10, false);
[ez3, ~] = a34_sim(3.75e9, 30, 0.9, 16, 1e-9, 2, 10, false);
[ez4, ~] = a34_sim(3.75e9, 30, 0.9, 32, 1e-9, 2, 10, false);
[ez1_error, ez1_err_max] = pml_error(ez1, ez0)
[ez1a_error, ez1a_err_max] = pml_error(ez1a, ez0)
[ez1b_error, ez1b_err_max] = pml_error(ez1b, ez0)
[ez2_error, ez2_err_max] = pml_error(ez2, ez0)
[ez3_error, ez3_err_max] = pml_error(ez3, ez0)
[ez4_error, ez4_err_max] = pml_error(ez4, ez0)


% A4
tic
[ez5, times5] = a34_sim(3.75e9, 30, 0.9, 8, 1e-9, 2, 5, false);
toc
v_phase_num5 = num_phase_velocity(4, 4.4, times5, freq_sim, 5, ez5, 10);
err5 = abs(v_phase_num5 - v_phase)/v_phase
tic
[ez6, times6] = a34_sim(3.75e9, 30, 0.9, 16, 1e-9, 2, 10, false);
toc
v_phase_num6 = num_phase_velocity(4, 4.4, times6, freq_sim, 10, ez6, 10);
err6 = abs(v_phase_num6 - v_phase)/v_phase
tic
[ez7, times7] = a34_sim(3.75e9, 30, 0.9, 32, 1e-9, 2, 20, false);
toc
v_phase_num7 = num_phase_velocity(4, 4.4, times7, freq_sim, 20, ez7, 10);
err7 = abs(v_phase_num7 - v_phase)/v_phase


function [ez_physical, times] = a34_sim(freq_sim, num_periods, courant_coeff,...
    n_pml, R_pml, order_pml, cells_per_lambda, video)
    common;
    dx = lambda0 / cells_per_lambda;
    nx = ceil(a / dx);
    dx = a / nx;
    ny = ceil(b / dx);
    dy = b / ny;

    y_int = linspace(0.0, b, ny + 1);
    Ezx = zeros(nx + 1 + 2*n_pml, ny + 1);
    Ezy = zeros(nx + 1 + 2*n_pml, ny + 1);
    Hx = zeros(nx + 1 + 2*n_pml, ny);
    Hy = zeros(nx + 2*n_pml, ny + 1);
    source_x_idx = nx / 2 + 1;                      % f matlab
    source_profile = sin(pi * y_int / b);
    dt_courant = 1.0 / (c0 * sqrt(1/dx^2 + 1/dy^2));

    x_int = linspace(-n_pml * dx, a + n_pml * dx, nx + 1 + 2*n_pml).';
    % real_idxs = (x_int > -1e-12) & (x_int < a + 1e-12);
    real_idxs = (n_pml + 1):(n_pml + nx + 1);
    x_half = (x_int(1:end-1) + x_int(2:end)) * 0.5;

    dt = dt_courant * courant_coeff;
    omega = 2 * pi * freq_sim;
    time_steps = ceil(num_periods / (freq_sim * dt));
    times = (0:time_steps)*dt;
    ez_tracker = zeros(time_steps + 1, nx + 1 + 2*n_pml, ny + 1);

    Ezx = hard_source_split(Ezx, Ezy, source_x_idx + n_pml, source_profile, omega, times(1));
    ez_tracker(1, :, :) = Ezx + Ezy;

    rate_const_e = rate_constant(x_int, a, n_pml * dx, c0, order_pml, R_pml);
    rate_const_h = rate_constant(x_half, a, n_pml * dx, c0, order_pml, R_pml);
    [self_e, curl_e] = update_coeffs(rate_const_e, dt, eps0);
    [self_h, curl_h] = update_coeffs(rate_const_h, dt, mu0);

    for n = 1:time_steps
        % first H, then E
        Hx = hx_update_pml(Ezx, Ezy, Hx, dt, dy, mu0);
        Hy = hy_update_pml(Ezx, Ezy, Hy, dx, self_h, curl_h);
        Ezx = ezx_update_pml(Ezx, Hy, dx, self_e, curl_e);
        Ezy = ezy_update_pml(Ezy, Hx, dt, dy, eps0);
        Ezx = hard_source_split(Ezx, Ezy, source_x_idx + n_pml, source_profile, omega, times(n + 1));
        ez_tracker(n + 1, :, :) = Ezx + Ezy;
    end

    ez_physical = ez_tracker(:, real_idxs, :);
    stamps = (1:4)*floor(time_steps/4) + 1;
    stamp_times = times(stamps);
    stamp_periods = stamp_times * freq_sim;

    figure;
    t = tiledlayout(4, 1, "TileSpacing", "compact");
    for i=1:4
        ax = nexttile;
        bo = zeros(nx + 1 + 2*n_pml, ny + 1);
        bo(:, :) = ez_tracker(stamps(i), :, :);
        imagesc(x_int, y_int, bo');
        xline(x_int(source_x_idx + n_pml), 'k:', 'LineWidth', 1.5);
        xline(x_int(n_pml + 1), 'k:', 'LineWidth', 1.5);
        xline(x_int(n_pml + nx + 1), 'k:', 'LineWidth', 1.5);
        xlabel("x [m]")
        ylabel("y [m]")
        title(sprintf("E_z, t=%.2f [ns], periods=%.2f", stamp_times(i)*1e9, stamp_periods(i)));
    end
    cb = colorbar;
    cb.Label.String = 'E_z [V/m]';
    cb.Layout.Tile = 'east';

    % max_ez = max(abs(ez_tracker(:)))
    if video
        writerObj = VideoWriter('ElectricField.mp4', 'MPEG-4');
        writerObj.FrameRate = 10;
        open(writerObj);

        figure;
        for k = 1:time_steps
            bo = zeros(nx + 1 + 2*n_pml, ny + 1);
            bo(:, :) = ez_tracker(k, :, :);
            imagesc(x_int, y_int, bo(:, :)');
            xline(x_int(source_x_idx + n_pml), 'k:', 'LineWidth', 1.5);
            xline(x_int(n_pml + 1), 'k:', 'LineWidth', 1.5);
            xline(x_int(n_pml + nx + 1), 'k:', 'LineWidth', 1.5);
            xlabel("x [m]")
            ylabel("y [m]")
            title(sprintf("E_z, t=%.2f [ns], periods=%.2f", times(k)*1e9, times(k)*freq_sim));
            cb = colorbar;
            cb.Label.String = 'E_z [V/m]';

            currFrame = getframe(gcf);
            writeVideo(writerObj, currFrame);
        end

        close(writerObj);
    end
end
