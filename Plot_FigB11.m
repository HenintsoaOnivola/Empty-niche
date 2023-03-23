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

log_sig_c_an=0.875;
log_sig_m=-1.75;

log_sig_an=1;
sig_an=log_sig_an;
sig_pl=sig_an;

sig_c_an=exp(log_sig_c_an);
sig_c_pl=sig_c_an;

sig_m=exp(log_sig_m);

the_folder='./Results_invasion_intra_guild/';
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

plot_traits_function_no_background(intro_index,separating_index,time_plot,trait_plot_an,trait_plot_pl,pop_plot_an,pop_plot_pl);

% to control the figure size and resolution
width = 7;     % Width in inches
height = 4;    % Height in inches
alw = 0.75;    % AxesLineWidth
fsz = 11;      % Fontsize
lw = 0.5;      % LineWidth
msz = 8;       % MarkerSize

plot_name=strcat('./Supplementary_figures/FigB11.png');
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf,'PaperPosition', myfiguresize);
print(plot_name,'-dpng','-r400');
hold off;

