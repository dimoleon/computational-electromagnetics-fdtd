function ezx = ezx_update_pml(ezx, hy, dx, self_e, curl_e)
    % ezx, ezy: nx + 1 + 2*pad, ny + 1
    % hy: nx + 2*pad, ny + 1
    % self_e, curl_e: nx + 1 + 2*pad, 1
    dHy = hy(2:end, 2:end-1) - hy(1:end-1, 2:end-1);
    ezx(2:end-1, 2:end-1) = ezx(2:end-1, 2:end-1) .* self_e(2:end-1)...
        + dHy .* (curl_e(2:end-1) / dx);
end
