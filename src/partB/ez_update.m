function ez = ez_update(ez, hx, dy, eps_new, eps_old, self_e, curl_e)
    dHx = hx(2:end) - hx(1:end-1);
    ez(2:end-1) = (self_e(2:end-1) .* eps_old(2:end-1) ./ eps_new(2:end-1)) .* ez(2:end-1)...
                    -(curl_e(2:end-1) ./ eps_new(2:end-1)) .* dHx / dy;
end
