function v_phase = num_phase_velocity(x1_w, x2_w, times, freq, cells_per_l, ez, period_begin)
    common;
    omega = 2 * pi * freq;
    lambda = c0 / freq;
    nx = size(ez, 2) - 1;
    ny = size(ez, 3) - 1;
    x_mid_idx = round(nx/2) + 1;
    y_mid_idx = round(ny/2) + 1;
    x1_idx = round((x1_w - 3) * cells_per_l) + x_mid_idx;
    x2_idx = round((x2_w - 3) * cells_per_l) + x_mid_idx;
    times_passed = (times > (period_begin / freq));
    t = times(times_passed);
    t = t(:);
    ez1 = ez(times_passed, x1_idx, y_mid_idx);
    % size(ez1)
    ez2 = ez(times_passed, x2_idx, y_mid_idx);
    A = [cos(omega*t), sin(omega*t)];
    % size(A)
    trig_phi1 = A \ ez1;
    r1 = trig_phi1(1);
    s1 = trig_phi1(2);
    trig_phi2 = A \ ez2;
    r2 = trig_phi2(1);
    s2 = trig_phi2(2);
    delta_phi = atan2(-r1*s2 + r2*s1, r1*r2 + s1*s2);
    v_phase = -omega * (x2_w - x1_w) * lambda / delta_phi;
end
