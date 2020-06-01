clc
clear all
close all
%% four links robot arm parameters
    global l1;
    global l2;
    global l3;
    global l4;
%     l1 = 0.14; 
%     l2 = 0.12;
%     l3 = 0.12;
%     l4 = 0.14;
    l1 = 0.2; 
    l2 = 0.2;
    l3 = 0.2;
    l4 = 0.2;

    global kc;
    global w_tau;
%     kc = [1000,200;
%              200,1000];
        kc = [100,20;
                20,100];
    w_tau =    [1/5, 0,    0,   0;
                       0,    1/4, 0,   0;
                       0,    0,    1/3,0;
                       0,    0,    0,   1/2];
    global mass_each_link;
    global q_initial;
    mass_each_link = 0.25;
    q_initial = [30, 30, -70, 45]/180*pi; 

%% initial position and jacobian
global ee_initial_position;
global points_initial;
[ee_initial_position, points_initial] = forward_kine(q_initial);

%%  sfe optimize
% unit_direction = [-1/sqrt(2); 1/sqrt(2)];
unit_direction = [-1; 0];
[q_res_sfe, sigma_revs] = optimize_sfe(unit_direction);
optimized_sigma_SFE = 1 / sigma_revs;
[Hend, points_res_SFE] = forward_kine(q_res_sfe);

%% gravity optimize
% global norm_list;
% norm_list = [];
global dispacement_bygravity;
dispacement_bygravity =  [];
[q_res_gravity, sigma_revs] = optimize_gravity();

[Hend, points_res_gravity] = forward_kine(q_res_gravity);
gravity_direction = mapminmax(Gq(q_res_gravity)', 0, 1);
%optimized_sigma_gravity = 1 / sigma_revs;
%% combined optimize
% unit_direction = [-1; 0];
weighting_factor  = [ 0.3  0.7];
[q_res_combined, sigma_revs] = optimize_combined(unit_direction, weighting_factor);

[Hend, points_res_combined] = forward_kine(q_res_combined);
gravity_direction2 = mapminmax(Gq(q_res_combined)', 0, 1);
%optimized_sigma_combined = 1 / sigma_revs;
%% visualization sfe
f1 = figure;
subplot(1,2,1);
plot_data(f1, 'optimized SFE', q_res_sfe, points_res_SFE, unit_direction);
subplot(1,2,2);
plot_data2(f1, 'optimized SFE', q_res_sfe)
title('optimized SFE and SFR', 'FontSize',15);
%% visualization gravity
f2 = figure;
subplot(1,2,1);
plot_data(f2, 'optimized gravity', q_res_gravity, points_res_gravity, gravity_direction');
subplot(1,2,2);
plot(dispacement_bygravity);
title('e-e displacement by gravity', 'FontSize',15);

%% visualizaion combined
f3 = figure;
% subplot(1,2,1);
plot_data(f3, 'optimized combined', q_res_combined, points_res_combined, [unit_direction, gravity_direction2']);
% subplot(1,2,2);
% plot_data2(f3, 'optimized SFE', q_res_combined)
% title('optimized SFE and SFR', 'FontSize',15);
%% print
sprintf('sigma  = %d ',optimized_sigma_SFE)




