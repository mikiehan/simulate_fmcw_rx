function [c,ceq] = nlcon(x)
  c = []; %sum(abs(exp(1i*x))) - 1.0;
  ceq = x(1) - 0;