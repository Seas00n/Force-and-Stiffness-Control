function [q_res, sigma_revs] = optimize_gravity(~)
global kc;
global q_initial;


cost_func =@(q)  ((kc\Gq(q))' * (kc\Gq(q)));
A=[];
b=[];
Aeq=[];
Beq=[];
lb= [ -pi+pi/60, -pi+pi/30, -pi+pi/30, -pi+pi/30];
ub=[   pi-pi/60,    pi-pi/30,    pi-pi/30,   pi-pi/30];
options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp', 'OutputFcn',@outfunc);
disp('optimizing gravity!!!!!!!!!!!!!!!!!!!!!!!!')
[q_res, sigma_revs] = fmincon(cost_func, q_initial, A, b, Aeq, Beq, lb, ub, @nonlinear_constr,options);
%Gq_out = Gq(q_res);
end


function stop = outfunc(x,optimValues,state)

% [~, tau_g] = Gq(x);
% global norm_list;
% norm_list = [norm_list, norm(tau_g)];
global dispacement_bygravity;
dispacement_bygravity = [dispacement_bygravity, optimValues.fval];
stop= 0;
end



