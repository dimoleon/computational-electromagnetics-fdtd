function hx = hx_update(ez, hx, dt, dy, mu)
    % ez : nx + 1, ny + 1
    % hx : nx + 1, ny
    % hy : nx, ny + 1
    dEz = ez(:, 2:end) - ez(:, 1:end-1);
    hx = hx - (dt / (mu * dy)) * dEz;
end
