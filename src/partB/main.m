clear;
clc;
close all;

ez0 = b_sim(, 1, )
ez1 = b_sim(0, 1e-6, 2);
ez2 = b_sim(8, 1e-6, 2);
ez3 = b_sim(16, 1e-6, 2);
ez4 = b_sim(100, 1e-6, 2);

function ez_physical = b_sim(n_pml, R_pml, order_pml)
    c0 = 299792458;
    eps0 = 8.85418781762e-12;
    mu0 = 1.256637061436e-06;

    freq = 1.5e9;
    lambda0 = c0 / freq;

    L = 1;
    Ly = 0.1;

    eps_rho = 5.25;
    delta_eps = 0.85;
    omega_eps = 2 * pi * 1e9;
    eps_relative = @(t) eps_rho + delta_eps * sin(omega_eps * t);

    lambda_ref = lambda0 / sqrt(eps_rho + delta_eps);
    n_dielectric = ceil(10 * Ly / lambda_ref);
    dy = Ly / n_dielectric;
    ny = round(2 * L / dy);
    y_int = linspace(-L - n_pml*dy, L + n_pml*dy, ny + 1 + 2*n_pml)';
    y_half = (y_int(2:end) + y_int(1:end-1)) * 0.5;
    real_idxs = (n_pml + 1):(n_pml + ny + 1);

    dielectric_mask = zeros(size(y_int));
    % take average eps_r in dielectric material edges
    dielectric_mask(ny/2 + n_pml + 1) = 0.5;
    dielectric_mask(ny/2 + n_pml + 1 + n_dielectric) = 0.5;
    % properly inside the dielectric material
    dielectric_mask(ny/2 + n_pml + 2: ny/2 + n_pml + n_dielectric) = 1.0;

    duration_source = 3 / freq;
    y_source = -L + 2 * lambda0;
    [~, source_idx] = min(abs(y_int - y_source));

    dt_courant = dy / c0;
    courant_coeff = 0.9;

    Ez = zeros(size(y_int));
    Hx = zeros(size(y_half));

    dt = dt_courant * courant_coeff;
    num_periods = 10;
    time_steps = ceil(num_periods / (freq * dt));
    times = (0:time_steps)*dt;
    ez_tracker = zeros(time_steps + 1, size(Ez, 1));

    Ez = hard_source(Ez, source_idx, freq, times(1), duration_source);
    ez_tracker(1, :) = Ez;

    rate_const_e = rate_constant(y_int + L, 2*L, n_pml*dy, c0, order_pml, R_pml);
    rate_const_h = rate_constant(y_half + L, 2*L, n_pml*dy, c0, order_pml, R_pml);
    [self_e, curl_e] = update_coeffs(rate_const_e, dt);
    [self_h, curl_h] = update_coeffs(rate_const_h, dt);
    eps_old = eps0 * dielectric_relative(dielectric_mask, eps_relative, times(1));

    for n = 1:time_steps
        Hx = hx_update(Ez, Hx, dy, mu0, self_h, curl_h);
        eps_new = eps0 * dielectric_relative(dielectric_mask, eps_relative, times(n + 1));
        Ez = ez_update(Ez, Hx, dy, eps_new, eps_old, self_e, curl_e);
        Ez = hard_source(Ez, source_idx, freq, times(n + 1), duration_source);
        eps_old = eps_new;
        ez_tracker(n + 1, :) = Ez;
    end

    ez_physical = ez_tracker(:, real_idxs);
    figure;
    imagesc(times*1e9, y_int(real_idxs), ez_tracker(:, real_idxs)');
    xline(y_int(source_idx), 'k:', 'LineWidth', 1.5);
    % xline(y_int(n_pml + 1), 'k:', 'LineWidth', 1.5);
    % xline(y_int(n_pml + ny + 1), 'k:', 'LineWidth', 1.5);
    hold on;
    t_limits = xlim;
    patch([t_limits(1) t_limits(1) t_limits(2) t_limits(2)], ...
           [0 Ly Ly 0], ...
          [0.4 0.4 0.4], 'FaceAlpha', 0.15, 'EdgeColor', 'none');

    yline(0, 'k--');
    yline(Ly, 'k--');

    xlabel("t [ns]")
    ylabel("y [m]")
    cb = colorbar;
    cb.Label.String = 'E_z [V/m]';
    hold off;
end
