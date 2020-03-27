clc
clear
close all
%% Input Data %%
omega = 1:1E-2:175;
omega = omega';
Theta = ones(size(omega));
Mt = ones(size(omega));
Jeq = [156146 33834 99530];
Keq = [575043128.23 123869689.09 0];
%% Natural Frequency Calculations Part A (Mteq and Theta 1,2,3 calculation) %%
for i=1:size(omega)
    Mt(i,1) = omega(i)^2*Jeq(1)*Theta(i,1);
end
for i=1:size(omega)
    Theta(i,2) =Theta(i,1)-(Mt(i,1)/Keq(1));
end
for i=1:size(omega)
    Mt(i,2) = Mt(i,1)+omega(i)^2*Jeq(2)*Theta(i,2);
end
for i=1:size(omega)
    Theta(i,3) =Theta(i,2)-(Mt(i,2)/Keq(2));
    Theta(i,4) =omega(i);
end
for i=1:size(omega)
    Mt(i,3) = Mt(i,2)+omega(i)^2*Jeq(3)*Theta(i,3);
    Mt(i,4) = omega(i);
end
%% Natural Frequency Calculations Part B (Natural Frequency tracking) %%
temp = Mt;
for i=size(temp):-1:1
    if temp(i,3)<0
        temp(i,:)=[];
    end
end
minimum1 = min(temp(:,3));
pos1 = find(temp(:,3)==minimum1);
omega1 = temp(pos1,4);
theta_omega1 = ones(1,3);
for i=1:1:3
    theta_omega1(i) = Theta(pos1,i);
end

for i=pos1:-1:1
    temp(i,:)=[];
end
minimum2 = min(temp(:,3));
pos2 = find(temp(:,3)==minimum2);
omega2 = temp(pos2,4);
pos_theta = find(Theta(:,4)==omega2);
theta_omega2 = ones(1,3);
for i=1:1:3
    theta_omega2(i) = Theta(pos_theta,i);
end
%% Natural Frequency Plot %%
figure,hold on,grid on
plot(omega,Mt(:,3),'b','Linewidth',2,'DisplayName','Mt3')
xlabel('omega (rad/s)','Interpreter','latex')
ylabel('Mt3 [kN/m]','Interpreter','latex')
set(gca,'fontsize',16)
plot(omega1,minimum1,'o',...
                'MarkerFaceColor','r',...
                'MarkerEdgeColor','k',...
                'MarkerSize',5,'DisplayName',['\omega1=' num2str(omega1)])
plot(omega2,minimum2,'o',...
                'MarkerFaceColor','g',...
                'MarkerEdgeColor','k',...
                'MarkerSize',5,'DisplayName',['\omega2=' num2str(omega2)])
legend('Location','northwest','Orientation','vertical')
%% Modal Shapes Plot %%
figure,hold on,grid on
plot(theta_omega1(1,:),'b','linewidth',2,'DisplayName',['\omega1=' num2str(omega1)])
plot(theta_omega2(1,:),'r','linewidth',2,'DisplayName',['\omega2=' num2str(omega2)])
plot(1,theta_omega1(1,1),'o',...
                'MarkerFaceColor','g',...
                'MarkerEdgeColor','k',...
                'MarkerSize',5,'DisplayName','M.E.')
plot(2,theta_omega1(1,2),'o',...
                'MarkerFaceColor','y',...
                'MarkerEdgeColor','k',...
                'MarkerSize',5,'DisplayName','S.G.')
plot(3,theta_omega1(1,3),'o',...
                'MarkerFaceColor','k',...
                'MarkerEdgeColor','k',...
                'MarkerSize',5,'DisplayName','P')
plot(1,theta_omega2(1,1),'o',...
                'MarkerFaceColor','g',...
                'MarkerEdgeColor','k',...
                'MarkerSize',5,'DisplayName','M.E.')
plot(2,theta_omega2(1,2),'o',...        
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','y',...
                'MarkerSize',5,'DisplayName','S.G.')
plot(3,theta_omega2(1,3),'o',...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','k',...
                'MarkerSize',5,'DisplayName','P')
legend('Location','southeast','Orientation','Vertical')
set(gca,'fontsize',16)
set(gca,'xlim',[1 3],'xtick',1:1:3)
set(gca,'ylim',[-6 1],'ytick',-6:1:1)





%% Task 3
Nn1 = omega1*30/pi;
Nn2 = omega2*30/pi;
Nmcr = 91;
%Boundaries
N20 = Nmcr * (0.2)^(1/3);
N103 = Nmcr * 1.03;
%Plot
figure
title("Campbell Diagram")
grid on, box on, hold on
ylabel("Natural Frequency [rpm]")
xlabel("Engine Speed [rpm]")
ylim([0 1.1*Nn2])
xlim([0 1.1*N103])
%Plot frequencies
yline(Nn1,'--b','LineWidth',1,'DisplayName',['Nn1=' num2str(Nn1)])
yline(Nn2,'--b','LineWidth',1,'DisplayName',['Nn2=' num2str(Nn2)])
%Plot speeds
xline(N20,'--b','LineWidth',1,'DisplayName',['N20=' num2str(N20)])
xline(N103,'--b','LineWidth',1,'DisplayName',['N103=' num2str(N103)])
%Plot boundaries
plot(N20,Nn1,'r*','DisplayName','Border #1')
plot(N20,Nn2,'r*','DisplayName','Border #2')
plot(N103,Nn1,'r*','DisplayName','Border #3')
plot(N103,Nn2,'r*','DisplayName','Border #4')
%Plot harmonics
a=1.1*N103;
j=1;
for i=1:1:50
    f=i*a;
    if (Nn1/(f/a))>N20 && (Nn1/(f/a)<N103)
        plot([0 a], [0 f],'Color','red','LineWidth',1,'DisplayName',['\lambda =' num2str(i)])
        plot(Nn1/(f/a),Nn1,'-o','MarkerFaceColor','magenta','Color','magenta','DisplayName',['Interception ',num2str(j),'=(', num2str(Nn1/(f/a)),',Nn1)'])
        j=j+1;
    elseif (Nn2/(f/a))>N20 && (Nn2/(f/a)<N103)
        plot([0 a], [0 f],'Color','red','LineWidth',1,'DisplayName',['\lambda =' num2str(i)])
        plot(Nn2/(f/a),Nn2,'-o','MarkerFaceColor','green','Color','green','DisplayName',['Interception ',num2str(j),'=(',num2str(Nn2/(f/a)),',Nn2)']')
        j=j+1;
    else
        plot([0 a], [0 f],'Color','black','DisplayName',['\lambda =' num2str(i)])
    end
end
legend('Location','northwest','Orientation','vertical')
lgd=legend;
lgd.FontSize = 7;
lgd.NumColumns = 3;

clearvars a f i j ans lgd minimum1 minimum2 N103 N20 Nmcr omega Mt pos1 pos2 pos_theta temp Theta
