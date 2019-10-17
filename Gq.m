function [res, tau_g] = Gq(q)
global mass_each_link;
global l1;
global l2;
global l3;
global l4;
G_link = mass_each_link*9.8;

tau_g = [G_link*l1/2*cos(q(1)) + G_link*[l1*cos(q(1)) + l2/2*cos(sum(q([1:2])))] + G_link*[l1*cos(q(1)) + l2*cos(sum(q([1:2]))) + l3/2*cos(sum(q([1:3])))]  + G_link*[l1*cos(q(1)) + l2*cos(sum(q([1:2]))) + l3*cos(sum(q([1:3]))) + l4/2 *cos(sum(q))];
                                         G_link*l2/2*cos(sum(q(1:2))) + G_link*(l3/2*cos(sum(q(1:3))) +l2*cos(sum(q(1:2)))) +  G_link*(l4/2 *cos(sum(q))+l3*cos(sum(q(1:3))) +l2*cos(sum(q(1:2))))   ;                                                                                                                                                   ;  
                                         G_link*l3/2*cos(sum(q([1:3]))) + G_link*[l3*cos(sum(q([1:3]))) + l4/2*cos(sum(q))] ;                                                                                                                         
                                         G_link*l4/2 *cos(sum(q)) ];
                                     
res = get_jacob(q)' \ tau_g;

end