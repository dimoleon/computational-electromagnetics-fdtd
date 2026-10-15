function [err_sq, err_max] = pml_error(ezA, ezR)
    err_sq = sqrt(sum(ezA(round(0.8*end):end)-ezR(round(0.8*end):end))^2) / sqrt(sum(ezR(round(0.8*end):end))^2);
    err_max = max(abs(ezA - ezR), [], "all");
end
