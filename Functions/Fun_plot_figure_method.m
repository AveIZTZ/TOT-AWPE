function Fun_plot_figure_method(time_vect,PM_data)

CD_origin   = PM_data.CD_origin;
LLR_origin  = PM_data.LLR_origin;
SNR_origin  = PM_data.SNR_origin;
PESQ_origin  = PM_data.PESQ_origin;

CD_matrix   = PM_data.CD_matrix;
LLR_matrix  = PM_data.LLR_matrix;
SNR_matrix  = PM_data.SNR_matrix;
PESQ_matrix = PM_data.PESQ_matrix;

% my_legend = {'original','AWPE','AWPE-RI','AWPE-VR'};
% my_legend = {'original','AWPE','AWPE-RI','AWPE-Mix'};
% my_legend = {'original','AWPE','AWPE-RI','AWPE-MPDR'};
% my_legend = {'original','AWPE','AWPE-RI','AWPE-Switch'};
%my_legend = {'original','AWPE','SKAWPE'};
my_legend = {'original','KAWPE,P=1','KAWPE,P=2','KAWPE,P=4'};

fig=figure('position',[0 0 1100 600]);
linewd          = 1.5;
FontSize        = 12;
subplot(221);
set(gcf,'DefaultTextInterpreter','latex')
set(fig, 'Color', [1, 1, 1]);hold on;
plot(time_vect, CD_origin,'-sb','LineWidth',linewd);grid on;
plot(time_vect, CD_matrix(:,1),'--*m','LineWidth',linewd);grid on;
if length(CD_matrix(1,:)) > 1
plot(time_vect, CD_matrix(:,2),':vk','LineWidth',linewd);grid on;
end
if length(CD_matrix(1,:)) > 2
plot(time_vect, CD_matrix(:,3),'-or','LineWidth',linewd);grid on;
end
legend(my_legend,'Location','SouthEast','Interpreter','latex')
xlabel('Time (s)','fontsize',FontSize,'interpreter', 'latex')
ylabel('CD','fontsize',FontSize,'interpreter', 'latex')
text(4.730,-1,'(a)','fontsize',FontSize,'interpreter', 'latex')
set(gca,'fontsize',FontSize,'fontname','Times new roman');
axis([0 PM_data.x_max PM_data.CD_min PM_data.CD_max]);box on;grid on;

subplot(222);
set(gcf,'DefaultTextInterpreter','latex')
set(fig, 'Color', [1, 1, 1]);hold on;
plot(time_vect, LLR_origin,'-sb','LineWidth',linewd);grid on;
plot(time_vect, LLR_matrix(:,1),'--*m','LineWidth',linewd);grid on;
if length(CD_matrix(1,:)) > 1
plot(time_vect, LLR_matrix(:,2),':vk','LineWidth',linewd);grid on;
end
if length(CD_matrix(1,:)) > 2
plot(time_vect, LLR_matrix(:,3),'-or','LineWidth',linewd);grid on;
end
legend(my_legend,'Location','SouthEast','Interpreter','latex')
xlabel('Time (s)','fontsize',FontSize,'interpreter', 'latex')
ylabel('LLR','fontsize',FontSize,'interpreter', 'latex')
text(4.730,1,'(b)','fontsize',FontSize,'interpreter', 'latex')
set(gca,'fontsize',FontSize,'fontname','Times new roman');
axis([0 PM_data.x_max PM_data.LLR_min PM_data.LLR_max]);box on;grid on;

subplot(223);
set(gcf,'DefaultTextInterpreter','latex')
set(fig, 'Color', [1, 1, 1]);hold on;
plot(time_vect, SNR_origin,'-sb','LineWidth',linewd);grid on;
plot(time_vect, SNR_matrix(:,1),'--*m','LineWidth',linewd);grid on;
if length(CD_matrix(1,:)) > 1
plot(time_vect, SNR_matrix(:,2),':vk','LineWidth',linewd);grid on;
end
if length(CD_matrix(1,:)) > 2
plot(time_vect, SNR_matrix(:,3),'-or','LineWidth',linewd);grid on;
end
legend(my_legend,'Location','NorthEast','Interpreter','latex')
xlabel('Time (s)','fontsize',FontSize,'interpreter', 'latex')
ylabel('fwsegSNR','fontsize',FontSize,'interpreter', 'latex')
text(4.730,-1,'(c)','fontsize',FontSize,'interpreter', 'latex')
set(gca,'fontsize',FontSize,'fontname','Times new roman');
axis([0 PM_data.x_max PM_data.SNR_min PM_data.SNR_max]);box on;grid on;

subplot(224);
set(gcf,'DefaultTextInterpreter','latex')
set(fig, 'Color', [1, 1, 1]);hold on;
plot(time_vect, PESQ_origin,'-sb','LineWidth',linewd);grid on;
plot(time_vect, PESQ_matrix(:,1),'--*m','LineWidth',linewd);grid on;
if length(CD_matrix(1,:)) > 1
plot(time_vect, PESQ_matrix(:,2),':vk','LineWidth',linewd);grid on;
end
if length(CD_matrix(1,:)) > 2
plot(time_vect, PESQ_matrix(:,3),'-or','LineWidth',linewd);grid on;
end
legend(my_legend,'Location','NorthEast','Interpreter','latex')
xlabel('Time (s)','fontsize',FontSize,'interpreter', 'latex')
ylabel('PESQ','fontsize',FontSize,'interpreter', 'latex')
text(4.730,-1,'(d)','fontsize',FontSize,'interpreter', 'latex')
set(gca,'fontsize',FontSize,'fontname','Times new roman');
axis([0 PM_data.x_max PM_data.PESQ_min PM_data.PESQ_max]);box on;grid on;

