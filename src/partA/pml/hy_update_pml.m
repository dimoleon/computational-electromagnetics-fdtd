function hy = hy_update_pml(ezx, ezy, hy, dx, self_h, curl_h)
    % ezx, ezy: nx + 1 + 2*pad, ny + 1
    % hy: nx + 2*pad, ny + 1
    % self_h, curl_h: nx + 2*pad, 1
    ez = ezx + ezy;
    dEz = ez(2:end, :) - ez(1:end-1, :);
    hy = hy .* self_h + dEz .* (curl_h / dx);
end
