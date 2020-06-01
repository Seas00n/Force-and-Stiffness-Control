function  [] = plot_data2(f, title_name, q_res)
global q_initial;
global points_initial;
figure(f)
title(title_name, 'FontSize',15);
xlabel('x');
ylabel('y');
ellpsoid1 = get_sfe(q_initial);
ellpsoid2 = get_sfe(q_res);
e1 = plot( ellpsoid1(1,:)+points_initial(5,1), ellpsoid1(2,:)+points_initial(5,2), '--', 'Color','#0072BD');
e1.LineWidth = 2;
hold on;
e2=plot(ellpsoid2(1,:)+points_initial(5,1), ellpsoid2(2,:)+points_initial(5,2), 'Color','#D95319');
e2.LineWidth = 2;

 hold on;

sfr1 = get_sfr(q_initial);
s1 = plot(sfr1(:,1)+points_initial(5,1),sfr1(:,2)+points_initial(5,2),'--', 'Color','b');
s1.LineWidth = 2;
sfr2 = get_sfr(q_res);
s2 = plot(sfr2(:,1)+points_initial(5,1),sfr2(:,2)+points_initial(5,2),'r');
s2.LineWidth = 2;

legend([e1 e2 s1 s2],{'origin sfe','optimized sfe','origin sfr','optimized sfr'}, 'Location','northwest', 'FontSize',15);
% hold on;
% legend([s1 s2],{'origin sfr','optimized sfr'}, 'Location','northwest', 'FontSize',15);

end