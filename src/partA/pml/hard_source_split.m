function ezx = hard_source_split(ezx, ezy, source_idx_padded, profile, omega, t)
    amplitude = 1;
    source_sum = amplitude * profile(2:end-1) * cos(omega * t);
    ezx(source_idx_padded, 2:end-1) = source_sum - ezy(source_idx_padded, 2:end-1);
end
