function [c, ceq] = constraint(x, sig, r)
c  = [];
ceq = x(1) * sig + x(2) * sig - r;
end
