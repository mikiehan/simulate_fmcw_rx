function [c,ceq] = nlcon(x)
  c = sum(abs(x)) - 1.0;
  ceq = [];