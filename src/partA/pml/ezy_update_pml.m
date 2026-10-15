function ezy = ezy_update_pml(ezy, hx, dt, dy, eps)
    % ezx, ezy: nx + 1 + 2*pad, ny + 1
    % hx: nx + 2*pad + 1, ny
    dHx = hx(2:end-1, 2:end) - hx(2:end-1, 1:end-1);
    ezy(2:end-1, 2:end-1) = ezy(2:end-1, 2:end-1) - (dt / eps) * (dHx / dy);
end
