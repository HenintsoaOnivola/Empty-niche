clear all;
warning('off', 'MATLAB:MKDIR:DirectoryExists');
addpath('Resource_competition_model');
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

parameters={[-1.9375 -2.8125], [-2 -2.625],-1.9375, -2 };
folders={'./Results_intra_guild/','./Results_cross_guild/','./Resource_competition_model/Results_RC/','./Resource_competition_model/Results_RC/'};
invasion_folders={'./Results_invasion_intra_guild/','./Results_invasion_cross_guild/'};
plot_names={'Fig2b','Fig2d','Fig2a','Fig2c'};
limits=[1,4;1.9,3.8;1,4;1.9,3.8];
x_ticks_labels={0:2:6,0:0.5:1.5,0:2:8,0:2:8};

log_sig_an=1;
sig_an=log_sig_an;
sig_pl=sig_an;

for each_pl=1:length(parameters)
    log_sig_c_an=parameters{each_pl}(1);
    sig_c_an=exp(log_sig_c_an);
    sig_c_pl=sig_c_an;
    if (length(parameters{each_pl})==2)
        log_sig_m=parameters{each_pl}(2);
        sig_m=exp(log_sig_m);
        folder_after_invasion=invasion_folders{each_pl};
    end
    
    disp (plot_names{each_pl});
    
    the_folder=folders{each_pl};
    
    if (length(parameters{each_pl})==2)%the case for model with resource competition and mutualism
        tr_an=read_data_into_cell(strcat(the_folder,'/animal_traits_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
        ab_an=read_data_into_cell(strcat(the_folder,'/animal_abundance_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
        time=read_file_to_array(strcat(the_folder,'/time_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
        sep_index=read_first_line_to_array(strcat(the_folder,'/sepIndex_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
        tr_pl=read_data_into_cell(strcat(the_folder,'/plant_traits_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
        ab_pl=read_data_into_cell(strcat(the_folder,'/plant_abundance_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
        
        tr_an = cellfun(@transpose,tr_an,'un',0);
        ab_an = cellfun(@transpose,ab_an,'un',0);
        tr_pl = cellfun(@transpose,tr_pl,'un',0);
        ab_pl = cellfun(@transpose,ab_pl,'un',0);
        
        Pareto_an=sum(read_last_line_to_array(strcat(folder_after_invasion,'/animal_abundance_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat')));
        Pareto_pl=sum(read_last_line_to_array(strcat(folder_after_invasion,'/plant_abundance_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat')));
        
        plot_fitness_background_function(sep_index,time,tr_an,tr_pl,ab_an,ab_pl,Pareto_an,Pareto_pl,rA,h,const,sig_c_an,sig_c_pl,sig_an,sig_pl,x0,y0,kx,ky,sig_m,limits(each_pl,:),x_ticks_labels{each_pl});
    else  % the case for model for a resource-competition model only
        tr_an=read_data_into_cell(strcat(the_folder,'/animal_traits_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'.dat'));
        ab_an=read_data_into_cell(strcat(the_folder,'/animal_abundance_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'.dat'));
        time=read_file_to_array(strcat(the_folder,'/time_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'.dat'));
        sep_index=read_first_line_to_array(strcat(the_folder,'/sepIndex_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'.dat'));
        
        % Cut off the few last trait evolution after convergence in figure 2a and 2c for the
        % sake of aesthetics 
        if each_pl==3
            time=time(time<=7*(10^(5)));
        end
        if each_pl==4
            time=time(time<=6.5*(10^(5)));
        end
        tr_an=tr_an(1:length(time));
        ab_an=ab_an(1:length(time));
        sep_index=[sep_index(sep_index<=length(time)) length(time)];
        
        tr_an = cellfun(@transpose,tr_an,'un',0);
        ab_an = cellfun(@transpose,ab_an,'un',0);
        Pareto_an=sum(ab_an{end});
        plot_fitness_background_function_RC(sep_index,time,tr_an,ab_an,Pareto_an,rA,sig_c_an,sig_an,x0,kx,limits(each_pl,:),x_ticks_labels{each_pl});
    end
    mkdir('./Plots/');
    plot_name=strcat('./Plots/',plot_names{each_pl},'.png');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % to control the figure size and resolution
    width = 7;     % Width in inches
    height = 6;    % Height in inches
    alw = 0.75;    % AxesLineWidth
    fsz = 11;      % Fontsize
    lw = 0.5;      % LineWidth
    msz = 8;       % MarkerSize

    %Here we preserve the size of the image when we save it.
    set(gcf,'InvertHardcopy','on');
    set(gcf,'PaperUnits', 'inches');
    papersize = get(gcf, 'PaperSize');
    left = (papersize(1)- width)/2;
    bottom = (papersize(2)- height)/2;
    myfiguresize = [left, bottom, width, height];
    set(gcf,'PaperPosition', myfiguresize);
    print(plot_name,'-dpng','-r600');
    hold off;

    clf;
   
end


