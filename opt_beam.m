function w = opt_beam(Sr, Nr) 

lb = zeros(1, Nr-1);
ub = 2 * pi * ones (1, Nr-1);
objective = @(x) func_obj(x, Sr);

x0 = zeros(1, Nr-1);

%nonlincon = @nlcon;
%x = fmincon(objective,x0,[],[],[],[],lb,ub,nonlincon);
x = fmincon(objective,x0,[],[],[],[],lb,ub,[]);

for i=1:Nr-1
    disp(['x' num2str(i) ' = ' num2str(x(i))])
end

w = [1 exp(1i*x)];

% Maximize RSS (sum of each sample's absolute value)
function obj = func_obj(x,r)
    obj = sum(abs([1 exp(1i*x)] * r)) * -1;
end

end

