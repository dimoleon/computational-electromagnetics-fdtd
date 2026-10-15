function ez = hard_source(ez, source_x_idx, profile, omega, t)
    amplitude = 1;
    ez(source_x_idx, 2:end-1) = amplitude * profile(2:end-1) * cos(omega * t);
end
