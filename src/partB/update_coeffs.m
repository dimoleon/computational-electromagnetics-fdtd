function [self, curl] = update_coeffs(rate_constant, dt)
    expon = -rate_constant * dt;
    self = exp(expon);
    curl = ones(size(rate_constant)) * dt;
    nonzero = (rate_constant ~= 0);
    curl(nonzero) = -expm1(expon(nonzero)) ./ rate_constant(nonzero);
end
