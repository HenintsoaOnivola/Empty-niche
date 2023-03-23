clear all;
warning('off', 'MATLAB:MKDIR:DirectoryExists');
%defining the parameters
rA=1;
rP=1;
const=0.1;
h=0.1;
x0=3;
y0=2;
kx=400;
ky=300;

log_sig_c_an=-1.5;
log_sig_m=-2.75;
log_sig_an=1;

sig_an=log_sig_an;
sig_pl=sig_an;

input_folders={'./Results_invasion_intra_guild/', '', './Results_invasion_cross_guild/'};
plot_names={'Fig4a','Fig4b','Fig4c'};


sig_c_an=exp(log_sig_c_an);
sig_c_pl=sig_c_an;

sig_m=exp(log_sig_m);

for each_pl=1:length(plot_names)
    disp (plot_names{each_pl});
    
    the_folder=input_folders{each_pl};
    if strcmp(the_folder,'')
        res_scenario2=invasion_different_scenario(rA,rP,h,const,log_sig_c_an,x0,y0,kx,ky,log_sig_m,'intra');
        time_plot=res_scenario2{1};
        trait_plot_an=res_scenario2{2};
        trait_plot_pl=res_scenario2{3};
        pop_plot_an=res_scenario2{4};
        pop_plot_pl=res_scenario2{5};
        separating_index=res_scenario2{6};
        intro_index=res_scenario2{7};
        
    else
        time_file=strcat(the_folder,'time_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat');
        time_plot=read_file_to_array(strcat(the_folder,'time_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
        trait_plot_an=read_data_into_cell(strcat(the_folder,'animal_traits_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
        trait_plot_pl=read_data_into_cell(strcat(the_folder,'plant_traits_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
        pop_plot_an=read_data_into_cell(strcat(the_folder,'animal_abundance_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
        pop_plot_pl=read_data_into_cell(strcat(the_folder,'plant_abundance_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
        separating_index=read_first_line_to_array(strcat(the_folder,'sepIndex_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
        intro_index=read_first_line_to_array(strcat(the_folder,'introIndex_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
    end
    last_tr_an=trait_plot_an{end};
    last_tr_pl=trait_plot_pl{end};
    last_pop_an=pop_plot_an{end};
    last_pop_pl=pop_plot_pl{end};

    trait_plot_an = cellfun(@transpose,trait_plot_an,'un',0);
    trait_plot_pl = cellfun(@transpose,trait_plot_pl,'un',0);
    pop_plot_an = cellfun(@transpose,pop_plot_an,'un',0);
    pop_plot_pl = cellfun(@transpose,pop_plot_pl,'un',0);

    plot_traits_function_with_positive_background_with_fitness(intro_index,separating_index,time_plot,trait_plot_an,trait_plot_pl,pop_plot_an,pop_plot_pl,rA,rP,const,h,x0,y0,kx,ky,sig_an,sig_pl,sig_c_an,sig_c_pl,sig_m,each_pl);
    mkdir('./Plots/');
    plot_name=strcat('./Plots/',plot_names{each_pl},'.png');
    print(plot_name,'-dpng','-r600');
    clf;
end