function [f3d, TMe, fxyz] = plots(time_i1, factor, time_i2, rg, pt, T_energy, rN)

f3d = figure(1);
for i = time_i1:factor:time_i2      
    plot3(rg{i}(1,:),rg{i}(2,:),rg{i}(3,:),'-','LineWidth',1.5,'Color',[0 0 0]); hold on;   
end
plot3(rg{time_i2}(1,:),rg{time_i2}(2,:),rg{time_i2}(3,:),'-','LineWidth',1.5,'Color',[0 0 0]); hold off;
set(gca,'FontSize',24)
set(gca,'FontSize',24,'TickLabelInterpreter','latex');
ylabel('Y','fontsize',24,'interpreter','latex');
xlabel('X','fontsize',24,'interpreter','latex');
zlabel('Z','fontsize',24,'interpreter','latex');
title('Deformed shapes','fontsize',20,'interpreter','latex');

TMe = figure(2);
plot(pt,T_energy,'-','LineWidth',1.5,'Color',[0 0 0])
set(gca,'FontSize',24)
set(gca,'FontSize',24,'TickLabelInterpreter','latex');
ylabel('Total mechanical energy','fontsize',24,'interpreter','latex');
xlabel('time','fontsize',24,'interpreter','latex');
title('Total mechanical energy vs. time','fontsize',20,'interpreter','latex');

for i = 1:time_i2
    rg_e(i,1) = pt(i);
    rg_e(i,2) = rg{i}(1,rN) - rg{1}(1,rN);
    rg_e(i,3) = rg{i}(2,rN)- rg{1}(2,rN);
    rg_e(i,4) = rg{i}(3,rN)- rg{1}(3,rN);   
end

fxyz = figure(3);
plot((rg_e(:,1)),rg_e(:,2),'-','LineWidth',1.5,'Color',[1 0 0]); hold on;
plot((rg_e(:,1)),rg_e(:,3),'-','LineWidth',1.5,'Color',[0 0 1]); hold on;
plot((rg_e(:,1)),rg_e(:,4),'-','LineWidth',1.5,'Color',[0 0 0]); hold off;

set(gca,'FontSize',24)
set(gca,'FontSize',24,'TickLabelInterpreter','latex');
ylabel('Displacement','fontsize',24,'interpreter','latex');
xlabel('time','fontsize',24,'interpreter','latex');
legend({'$u_{x}$', '$u_{y}$', '$u_{z}$'},'Interpreter','latex','Location','NorthWest','Orientation','horizontal','fontsize',24)
title('Displacement vs. time','fontsize',20,'interpreter','latex');
end

