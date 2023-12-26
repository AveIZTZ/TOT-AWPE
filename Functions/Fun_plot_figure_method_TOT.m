function Fun_plot_figure_method_TOT(time_vect,PM_data)

P = length(PM_data.CD_matrix(1,:));

CD_origin   = repmat(PM_data.CD_origin, 1, P);
SNR_origin  = repmat(PM_data.SNR_origin, 1, P);
PESQ_origin  = repmat(PM_data.PESQ_origin, 1, P);

CD_matrix   = CD_origin-PM_data.CD_matrix;
SNR_matrix  = PM_data.SNR_matrix-SNR_origin;
PESQ_matrix = PM_data.PESQ_matrix-PESQ_origin;

%my_legend = {'TOT-AWPE,P=1','TOT-AWPE,P=2','TOT-AWPE,P=3'};
my_legend = {'AWPE','TOT-AWPE,P=1','TOT-AWPE,P=2','TOT-AWPE,P=3','TOT-AWPE,P=4'};

x_max = 25;
fig=figure('position',[0 0 500 650]);
linewd          = 1.4;
FontSize        = 11;

subplot('Position', [0.1, 0.67, 0.8, 0.23]);
set(gcf,'DefaultTextInterpreter','latex')
set(fig, 'Color', [1, 1, 1]);hold on;
plot(time_vect, CD_matrix(:,1),'-^b','LineWidth',linewd);grid on;
plot(time_vect, CD_matrix(:,2),'-*c','LineWidth',linewd);grid on;
plot(time_vect, CD_matrix(:,3),'-dg','LineWidth',linewd);grid on;
plot(time_vect, CD_matrix(:,4),'-or','LineWidth',linewd);grid on;
plot(time_vect, CD_matrix(:,5),'-sm','LineWidth',linewd);grid on;
h1 = legend(my_legend, 'Location', 'NorthOutside', 'Interpreter', 'latex', 'NumColumns', 2);
legend_position = get(h1, 'Position');
set(h1, 'Position', [legend_position(1), legend_position(2) + 0.1, legend_position(3:4)]);
xlabel('Time (s)','fontsize',FontSize,'interpreter', 'latex')
ylabel('$\Delta$CD','fontsize',FontSize,'interpreter', 'latex')
text(x_max/2, 2.75,'(a)','fontsize',FontSize,'interpreter', 'latex')
set(gca,'fontsize',FontSize,'fontname','Times new roman');
axis([0 x_max 0 3]);box on;grid on;

subplot('Position', [0.1, 0.37, 0.8, 0.23]);
set(gcf,'DefaultTextInterpreter','latex')
set(fig, 'Color', [1, 1, 1]);hold on;
plot(time_vect, SNR_matrix(:,1),'-^b','LineWidth',linewd);grid on;
plot(time_vect, SNR_matrix(:,2),'-*c','LineWidth',linewd);grid on;
plot(time_vect, SNR_matrix(:,3),'-dg','LineWidth',linewd);grid on;
plot(time_vect, SNR_matrix(:,4),'-or','LineWidth',linewd);grid on;
plot(time_vect, SNR_matrix(:,5),'-sm','LineWidth',linewd);grid on;
%legend(my_legend,'Location','NorthEast','Interpreter','latex')
xlabel('Time (s)','fontsize',FontSize,'interpreter', 'latex')
ylabel('$\Delta$FWSNR (dB)','fontsize',FontSize,'interpreter', 'latex')
text(x_max/2,17.5,'(b)','fontsize',FontSize,'interpreter', 'latex')
set(gca,'fontsize',FontSize,'fontname','Times new roman');
axis([0 x_max 0 20]);box on;grid on;

subplot('Position', [0.1, 0.07, 0.8, 0.23]);
set(gcf,'DefaultTextInterpreter','latex')
set(fig, 'Color', [1, 1, 1]);hold on;
plot(time_vect, PESQ_matrix(:,1),'-^b','LineWidth',linewd);grid on;
plot(time_vect, PESQ_matrix(:,2),'-*c','LineWidth',linewd);grid on;
plot(time_vect, PESQ_matrix(:,3),'-dg','LineWidth',linewd);grid on;
plot(time_vect, PESQ_matrix(:,4),'-or','LineWidth',linewd);grid on;
plot(time_vect, PESQ_matrix(:,5),'-sm','LineWidth',linewd);grid on;
%legend(my_legend,'Location','NorthEast','Interpreter','latex')
xlabel('Time (s)','fontsize',FontSize,'interpreter', 'latex')
ylabel('$\Delta$PESQ','fontsize',FontSize,'interpreter', 'latex')
text(x_max/2,2.25,'(c)','fontsize',FontSize,'interpreter', 'latex')
set(gca,'fontsize',FontSize,'fontname','Times new roman');
axis([0 x_max 0 2.5]);box on;grid on;