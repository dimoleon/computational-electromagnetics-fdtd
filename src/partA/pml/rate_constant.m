function alpha = rate_constant(x_pos, real_length, pml_length, c_in, order, R_goal)
    pml_depth = zeros(size(x_pos));
    pml_depth(x_pos < 0) = -x_pos(x_pos < 0);
    pml_depth(real_length - x_pos < 0) = x_pos(real_length - x_pos < 0) - real_length;
    alpha_max = -(order + 1) * c_in * log(R_goal) / (2 * pml_length);
    alpha = zeros(size(x_pos));
    alpha(pml_depth > 0) = alpha_max * (pml_depth(pml_depth > 0) / pml_length).^order;
end
