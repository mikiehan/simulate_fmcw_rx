function x = opt_beam(Sr, Nr) % for now num mics are hard-coded to 8


objective = @(x) func_obj(x, Sr);

x0 = ones(1, Nr)/Nr;

nonlincon = @nlcon;
x = fmincon(objective,x0,[],[],[],[],[],[],nonlincon);

for i=1:Nr
    disp(['x' num2str(i) ' = ' num2str(x(i))])
end


% Maximize RSS (sum of each sample's absolute value)
function obj = func_obj(x,r)
    obj = sum(abs(x * r)) * -1;
end


end

