function hy = hy_update(ez, hy, dt, dx, mu)
    % ez : nx + 1, ny + 1
    % hx : nx + 1, ny
    % hy : nx, ny + 1
    dEz = ez(2:end, :) - ez(1:end-1, :);
    hy = hy + (dt / (mu * dx)) * dEz;
end
