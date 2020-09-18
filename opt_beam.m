function w = opt_beam(Sr) % for now num mics are hard-coded to 8

lb = zeros(1, 8);
ub = 2 * pi * ones (1, 8);

objective = @(x) func_obj(x, Sr);

x0 = zeros(1, 8);

x = fmincon(objective,x0,[],[],[],[],lb,ub,[]);

disp(['x1 = ' num2str(x(1))])
disp(['x2 = ' num2str(x(2))])
disp(['x3 = ' num2str(x(3))])
disp(['x4 = ' num2str(x(4))])
disp(['x5 = ' num2str(x(5))])
disp(['x6 = ' num2str(x(6))])
disp(['x7 = ' num2str(x(7))])
disp(['x8 = ' num2str(x(8))])

w = exp(1i*x); %weights

% Maximize RSS (sum of each sample's absolute value)
function obj = func_obj(x,r)
   obj = sum(abs(exp(1i*x(1))*r(1, :)+ exp(1i*x(2))*r(2, :) + exp(1i*x(3))*r(3, :) ...
         + exp(1i*x(4))*r(4, :)+ exp(1i*x(5))*r(5, :) + exp(1i*x(6))*r(6, :) ...
         + exp(1i*x(7))*r(7, :)+ exp(1i*x(8))*r(8, :)))*-1;
end


end

