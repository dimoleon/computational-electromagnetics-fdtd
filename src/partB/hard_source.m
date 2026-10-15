function ez = hard_source(ez, source_idx, freq, time, duration)
    if time <= duration
        ez(source_idx) = sin(2 * pi * freq * time);
    end
end
