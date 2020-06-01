function [q_res, sigma_revs] = optimize_combined(unit_direction,weighting_factor)
global kc;
global w_tau;
global q_initial;

A=[];
b=[];
Aeq=[];
Beq=[];
lb= [ -pi+pi/60, -pi+pi/30, -pi+pi/30, -pi+pi/30];
ub=[   pi-pi/60,    pi-pi/30,    pi-pi/30,   pi-pi/30];

% weighting_factor  = [ 0.6  0.4];
cost_func_gravity =@(q)  ((kc\Gq(q))' * (kc\Gq(q)));
cost_func_sfe       =@(q)  (unit_direction' * kc * get_jacob(q)* w_tau^2 *  get_jacob(q)' * kc * unit_direction)^0.5;
cost_func  = @(q)   weighting_factor(1)*cost_func_sfe(q) + weighting_factor(2)*cost_func_gravity(q);

options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp');%, 'OutputFcn',@outfunc);
[q_res, sigma_revs] = fmincon(cost_func, q_initial, A, b, Aeq, Beq, lb, ub, @nonlinear_constr, options);
end


% function stop = outfunc(x,optimValues,state)
% 
% [~, tau_g] = Gq(x);
% global norm_list;
% norm_list = [norm_list, norm(tau_g)];
% 
% stop= 0;
% end
