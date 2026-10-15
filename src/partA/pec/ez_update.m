function ez = ez_update(ez, hx, hy, dt, dx, dy, eps)
    % ez : nx + 1, ny + 1
    % hx : nx + 1, ny
    % hy : nx, ny + 1
    dHx = hx(2:end-1, 2:end) - hx(2:end-1, 1:end-1);
    dHy = hy(2:end, 2:end-1) - hy(1:end-1, 2:end-1);
    ez(2:end-1, 2:end-1) = ez(2:end-1, 2:end-1) + (dt/eps) * (dHy/dx - dHx/dy);
end
