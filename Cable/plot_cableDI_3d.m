% plots DI as a function of space and time
% use after Dave's plot_cable_3d.m
% written by Trine Krogh-Madsen 2008
%
% I read in a file that has four columns: (1) an integer wave counter (not
% used here); (2) time of occurance of the DI; (3) location of DI
% occurance; (4) DI value

function p=plot_cableDI_3d(filename,dx)

cable = load(filename);

colormap('default')
co=colormap;
z = cable(:,4);
zspread = max(z) - min(z);
% one can also set the spread and the min value manually, e.g., to keep them the
% same for two different files

figure(1)
hold on
for i=1:length(cable),
  plot3(dx*cable(i,3),cable(i,2),-40,'.','Color',co(round((z(i)-min(z))/zspread*63+1),:),'MarkerSize',12,'LineWidth',2)
end


view([90 70])
xlabel('x (cm)','FontSize',14)
ylabel('Time (ms)','FontSize',14)
%zlabel('DI (ms)','FontSize',14)
