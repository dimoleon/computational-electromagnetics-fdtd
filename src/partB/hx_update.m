function hx = hx_update(ez, hx, dy, mu, self_h, curl_h)
    dEz = ez(2:end) - ez(1:end-1);
    hx = self_h .* hx - curl_h .* dEz / (mu * dy);
end
