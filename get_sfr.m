function [sfr] = get_sfr(q)
global mass_each_link;
global G_link; 
G_link = mass_each_link*9.8;
global m;
global n;

m = 25;   % forward rate
n = 360;  % number of directions

F_des1 = F_des(n);

global B; 
B = F_des1;

X_reality = F_real(q);

delt_x = squeeze(X_reality(1,m,:));
delt_y = squeeze(X_reality(2,m,:));
sfr =[delt_x  delt_y];
% 
% disp(size(delt_x))
% disp(size(delt_y))
end

function [F_des1] = F_des(n)
global kc;
global m;

for j = 1:n
      theta(j) =  pi/(n/2)*(j-1);
    for i = 1:m
        x(j,i) =   0.01*i*cos(theta(j));
        y(j,i) =   0.002*i*sin(theta(j));
    end
end

for i = 1:n
    F_des1(:,:,i) = kc*[x(i,:);y(i,:)]; 
end

end

function [X_reality] = F_real(q_in)
global l1;
global l2;
global l3;
global l4
global kc
global B
global G_link;  
global m;
global n;
for i =1:n
    theta(i) =  pi/(n/2)*(i-1);
    dd(1,i)  =  0.01*cos(theta(i));
    dd(2,i)  =  0.002*sin(theta(i)); 
end 

for j  = 1:n
    q(1) = q_in(1);
    q(2) = q_in(2);
    q(3) = q_in(3);  
    q(4) = q_in(4); 
    F = B(:,:,j);
    
    for i = 1:m
        J1 = get_jacob(q);
        g1(i) = (G_link*l1/2*cos(q(1)) + G_link*[l1*cos(q(1)) + l2/2*cos(sum(q([1:2])))] + G_link*[l1*cos(q(1)) + l2*cos(sum(q([1:2]))) + l3/2*cos(sum(q([1:3])))]  + G_link*[l1*cos(q(1)) + l2*cos(sum(q([1:2]))) + l3*cos(sum(q([1:3]))) + l4/2 *cos(sum(q))]);
        g2(i) = (G_link*l2/2*cos(sum(q(1:2))) + G_link*(l3/2*cos(sum(q(1:3))) +l2*cos(sum(q(1:2)))) +  G_link*(l4/2 *cos(sum(q))+l3*cos(sum(q(1:3))) +l2*cos(sum(q(1:2)))));                                                                                                                                                 
        g3(i) = (G_link*l3/2*cos(q(1)+q(2)+q(3)) + G_link*(l3*cos(q(1)+q(2)+q(3))) + l4/2*cos(q(1)+q(2)+q(3)+q(4)));                                                                                                                     
        g4(i) = (G_link*l4/2 *cos(q(1)+q(2)+q(3)+q(4)));
        
        tauA = J1.'*F +[g1(i);g2(i);g3(i);g4(i)];
        tauA_1 = tauA(1,:);
        tauA_2 = tauA(2,:);
        tauA_3 = tauA(3,:);
        tauA_4 = tauA(4,:);
        
        tauA_1(tauA_1>5) = 5;
        tauA_1(tauA_1<-5) = -5;
        tauA_2(tauA_2>4) = 4;
        tauA_2(tauA_2<-4) = -4;
        tauA_3(tauA_3>3) = 3;
        tauA_3(tauA_3<-3) = -3;  
        tauA_4(tauA_4>2) = 2;
        tauA_4(tauA_4<-2) = -2;
          
        F_reality(:,:,j)  = pinv(J1.')*([tauA_1;tauA_2;tauA_3;tauA_4] -[g1(i);g2(i);g3(i);g4(i)]);   % get reality force
        X_reality(:,:,j)  = inv(kc)*F_reality(:,:,j);
          
        dq(:,i) = pinv(J1)*dd(:,j);
        q(1) = q(1) + dq(1,i);
        q(2) = q(2) + dq(2,i);
        q(3) = q(3) + dq(2,i);  
        q(4) = q(4) + dq(2,i);  
        
        end 
    end
end


