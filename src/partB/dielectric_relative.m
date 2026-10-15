function eps_r = dielectric_relative(mask, eps_func, time)
    eps_r = (1 - mask)*1 + mask*eps_func(time);
end
