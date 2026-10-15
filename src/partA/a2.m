clear;
clc;
close all;

addpath("./pec");

% courant_coeff = 1.1;
% num_periods = 12;
% freq_sim = 3.75e9;

% (a)
a2_sim(3.75e9, 12, 1.1);
% (b)
a2_sim(3.75e9, 12, 0.9);
% (c)
a2_sim(2.5e9, 12, 0.9);

function a2_sim(freq_sim, num_periods, courant_coeff)
    common;

    Ez = zeros(nx + 1, ny + 1);
    Hx = zeros(nx + 1, ny);
    Hy = zeros(nx, ny + 1);

    dt = dt_courant * courant_coeff;
    omega = 2 * pi * freq_sim;
    time_steps = ceil(num_periods / (freq_sim * dt));
    times = (0:time_steps) * dt;
    ez_tracker = zeros(time_steps + 1, nx + 1, ny + 1);

    % starting point
    Ez = hard_source(Ez, source_x_idx, source_profile, omega, times(1));
    ez_tracker(1, :, :) = Ez;

    for n = 1:time_steps
        % first H, then E
        Hy = hy_update(Ez, Hy, dt, dx, mu0);
        Hx = hx_update(Ez, Hx, dt, dy, mu0);
        Ez = ez_update(Ez, Hx, Hy, dt, dx, dy, eps0);
        Ez = hard_source(Ez, source_x_idx, source_profile, omega, times(n + 1));
        ez_tracker(n + 1, :, :) = Ez;
    end

    stamps = (1:4)*floor(time_steps/4) + 1;
    stamp_times = times(stamps);
    stamp_periods = stamp_times * freq_sim;

    figure;
    t = tiledlayout(4, 1, "TileSpacing", "compact");
    for i=1:4
        ax = nexttile;
        bo = zeros(nx + 1, ny + 1);
        bo(:, :) = ez_tracker(stamps(i), :, :);
        imagesc(x_int, y_int, bo');
        xline(x_int(source_x_idx), 'k:', 'LineWidth', 1.5);
        xlabel("x [m]")
        ylabel("y [m]")
        title(sprintf("E_z, t=%.2f [ns], periods=%.2f", stamp_times(i)*1e9, stamp_periods(i)));
    end
    cb = colorbar;
    cb.Label.String = 'E_z [V/m]';
    cb.Layout.Tile = 'east';
    writerObj = VideoWriter('a2.mp4', 'MPEG-4');
    writerObj.FrameRate = 10;
    open(writerObj);

    n_pml = 0;
    figure;
    for k = 1:time_steps
        bo = zeros(nx + 1 + 2*n_pml, ny + 1);
        bo(:, :) = ez_tracker(k, :, :);
        imagesc(x_int, y_int, bo');
        xline(x_int(source_x_idx + n_pml), 'k:', 'LineWidth', 1.5);
        % xline(x_int(n_pml + 1), 'k:', 'LineWidth', 1.5);
        % xline(x_int(n_pml + nx + 1), 'k:', 'LineWidth', 1.5);
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
