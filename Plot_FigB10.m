warning('off', 'MATLAB:MKDIR:DirectoryExists');
mkdir('./Supplementary_figures/');

%defining the parameters
rA=1;
rP=1;
const=0.1;
h=0.1;
x0=3;
y0=2;
kx=400;
ky=300;
al=1;

log_sig_c_an=-1.125;
log_sig_m=-1.25;

log_sig_an=1;
sig_an=log_sig_an;
sig_pl=sig_an;

sig_c_an=exp(log_sig_c_an);
sig_c_pl=sig_c_an;

sig_m=exp(log_sig_m);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%simulate for the case where the selection gradient threshold for
%equilibrium is smaller than the normal one (i.e. equal to 10^(-4))
res=simulation_different_sel_grad_threshold(sig_c_an,sig_m,'mutu',10^(-4));
time_plot=res{1};
trait_plot_an=res{2};
trait_plot_pl=res{3};
pop_plot_an=res{4};
pop_plot_pl=res{5};
sepIndex=res{6};
res2=simulation_invasion_different_threshold_sel_grad(sig_c_an,sig_m,10^(-4),trait_plot_an,trait_plot_pl,pop_plot_an,pop_plot_pl,time_plot,sepIndex);
time_plot=res2{1};
trait_plot_an=res2{2};
trait_plot_pl=res2{3};
pop_plot_an=res2{4};
pop_plot_pl=res2{5};
separating_index=res2{6};
intro_index=res2{7};

trait_plot_an_trans = cellfun(@transpose,trait_plot_an,'un',0);
trait_plot_pl_trans = cellfun(@transpose,trait_plot_pl,'un',0);
pop_plot_an_trans = cellfun(@transpose,pop_plot_an,'un',0);
pop_plot_pl_trans = cellfun(@transpose,pop_plot_pl,'un',0);

% to control the figure size and resolution
width = 7;     % Width in inches
height = 4;    % Height in inches
alw = 0.75;    % AxesLineWidth
fsz = 11;      % Fontsize
lw = 0.5;      % LineWidth
msz = 8;       % MarkerSize

%plotting the trait evolution for each_guild
plot_traits_function_one_guild(intro_index,separating_index,time_plot,trait_plot_an_trans,pop_plot_an_trans,pop_plot_pl_trans,'Animal trait',[2 4])
plot_name=strcat('./Supplementary_figures/FigB10b_animal_traits.png');
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf,'PaperPosition', myfiguresize);
print(plot_name,'-dpng','-r400');
hold off;

plot_traits_function_one_guild(intro_index,separating_index,time_plot,trait_plot_pl_trans,pop_plot_pl_trans,pop_plot_an_trans,'Plant trait',[1 3])
plot_name=strcat('./Supplementary_figures/FigB10d_plant_traits.png');
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf,'PaperPosition', myfiguresize);
print(plot_name,'-dpng','-r400');
hold off;

%ploting the selection gradient through time
all_sel_grad=selection_gradient_function_through_time(pop_plot_an,pop_plot_pl,trait_plot_an,trait_plot_pl,rA,rP,h,const,sig_c_an,sig_c_pl,sig_an,sig_pl,x0,y0,kx,ky,sig_m);
figure('Visible','off');
for i=1:length(time_plot)
    an_side=all_sel_grad{i}(1:length(trait_plot_an{i}));
    plot(time_plot(i),log10(an_side),'bo');
    hold on;
end
set(gca,'ylim',[-21 0]);
set(gca,'Ytick',-21:3:0);
set(gca,'YtickLabel',10.^(-21:3:0))
xlabel('Evolutionary time','FontSize',20);
ylabel ('Selection pressure','FontSize',20);
plot_name=strcat('./Supplementary_figures/FigB10b_selection_pressure.png');
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf,'PaperPosition', myfiguresize);
print(plot_name,'-dpng','-r400');
hold off;

figure('Visible','off');
for i=1:length(time_plot)
    pl_side=all_sel_grad{i}(length(trait_plot_an{i})+1:end);
    plot(time_plot(i),log10(pl_side),'bo');
    hold on;
end
set(gca,'ylim',[-21 0]);
set(gca,'Ytick',-21:3:0);
set(gca,'YtickLabel',10.^(-21:3:0))
xlabel('Evolutionary time','FontSize',20);
ylabel ('Selection pressure','FontSize',20);
plot_name=strcat('./Supplementary_figures/FigB10d_selection_pressure.png');
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf,'PaperPosition', myfiguresize);
print(plot_name,'-dpng','-r400');
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%read the data from file for the case the threshold of the equilibrium is
%the normally used one (10^(-8))
the_folder='./Results_invasion_cross_guild/';
time_plot=read_file_to_array(strcat(the_folder,'time_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
trait_plot_an=read_data_into_cell(strcat(the_folder,'animal_traits_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
trait_plot_pl=read_data_into_cell(strcat(the_folder,'plant_traits_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
pop_plot_an=read_data_into_cell(strcat(the_folder,'animal_abundance_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
pop_plot_pl=read_data_into_cell(strcat(the_folder,'plant_abundance_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
separating_index=read_first_line_to_array(strcat(the_folder,'sepIndex_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
intro_index=read_first_line_to_array(strcat(the_folder,'introIndex_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));

trait_plot_an_trans = cellfun(@transpose,trait_plot_an,'un',0);
trait_plot_pl_trans = cellfun(@transpose,trait_plot_pl,'un',0);
pop_plot_an_trans = cellfun(@transpose,pop_plot_an,'un',0);
pop_plot_pl_trans = cellfun(@transpose,pop_plot_pl,'un',0);

%plotting the trait evolution for each_guild
plot_traits_function_one_guild(intro_index,separating_index,time_plot,trait_plot_an_trans,pop_plot_an_trans,pop_plot_pl_trans,'Animal trait',[2 4])
plot_name=strcat('./Supplementary_figures/FigB10a_animal_traits.png');
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf,'PaperPosition', myfiguresize);
print(plot_name,'-dpng','-r400');
hold off;

plot_traits_function_one_guild(intro_index,separating_index,time_plot,trait_plot_pl_trans,pop_plot_pl_trans,pop_plot_an_trans,'Plant trait',[1 3])
plot_name=strcat('./Supplementary_figures/FigB10c_plant_traits.png');
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf,'PaperPosition', myfiguresize);
print(plot_name,'-dpng','-r400');
hold off;

%ploting the selection gradient through time
all_sel_grad=selection_gradient_function_through_time(pop_plot_an,pop_plot_pl,trait_plot_an,trait_plot_pl,rA,rP,h,const,sig_c_an,sig_c_pl,sig_an,sig_pl,x0,y0,kx,ky,sig_m);
figure('Visible','off');
for i=1:length(time_plot)
    an_side=all_sel_grad{i}(1:length(trait_plot_an{i}));
    plot(time_plot(i),log10(an_side),'bo');
    hold on;
end
set(gca,'ylim',[-21 0]);
set(gca,'Ytick',-21:3:0);
set(gca,'YtickLabel',10.^(-21:3:0))
xlabel('Evolutionary time','FontSize',20);
ylabel ('Selection pressure','FontSize',20);
plot_name=strcat('./Supplementary_figures/FigB10a_selection_pressure.png');
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf,'PaperPosition', myfiguresize);
print(plot_name,'-dpng','-r400');
hold off;

figure('Visible','off');
for i=1:length(time_plot)
    pl_side=all_sel_grad{i}(length(trait_plot_an{i})+1:end);
    plot(time_plot(i),log10(pl_side),'bo');
    hold on;
end
set(gca,'ylim',[-21 0]);
set(gca,'Ytick',-21:3:0);
set(gca,'YtickLabel',10.^(-21:3:0))
xlabel('Evolutionary time','FontSize',20);
ylabel ('Selection pressure','FontSize',20);
plot_name=strcat('./Supplementary_figures/FigB10c_selection_pressure.png');
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf,'PaperPosition', myfiguresize);
print(plot_name,'-dpng','-r400');
hold off;