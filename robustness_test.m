clear all;
warning('off', 'MATLAB:MKDIR:DirectoryExists');
mkdir('./Results_robustness/');

%Defining the parameters at their default values
rA=1;
rP=1;
const=0.1;
h=0.1;
x0=3;
y0=2;
kx=400;
ky=300;
al=1;
log_sig_c_an_range=-2.5:0.125:1;
log_sig_m_range=-3:0.125:-0.5;

log_sig_an=1;
sig_an=log_sig_an;
  
%Vary the width of the intra-specific resource kernel
disp ('Vary the width of the intra-specific resource kernel');
sig_an=0.75;
log_sig_an=sig_an;

all_res_opt={'intra','cross'};

for res_opt_ind=1:length(all_res_opt)
    res_opt=all_res_opt{res_opt_ind};
    the_folder=strcat('./Results_robustness/vary_sigma_a_',res_opt,'/');
    mkdir(the_folder);
    the_inv_folder=strcat('./Results_robustness/vary_sigma_a_',res_opt,'_invasion/');
    mkdir(the_inv_folder);
    for log_sig_c_an=log_sig_c_an_range
        for log_sig_m=log_sig_m_range
            disp (strcat(res_opt,',',num2str(log_sig_c_an),',',num2str(log_sig_m)));
            sig_c_an=exp(log_sig_c_an);
            sig_c_pl=sig_c_an;
            sig_m=exp(log_sig_m);
            res=simulation_for_robustness_test_function(rA,rP,const,x0,y0,kx,ky,sig_an,log_sig_c_an,log_sig_m,res_opt);
            clear_files(strcat(the_folder,'time_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
            clear_files(strcat(the_folder,'animal_traits_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
            clear_files(strcat(the_folder,'plant_traits_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
            clear_files(strcat(the_folder,'animal_abundance_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
            clear_files(strcat(the_folder,'plant_abundance_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
            clear_files(strcat(the_folder,'sepIndex_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
            write_into_files_W(strcat(the_folder,'time_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'),res{1});
            write_into_files_W(strcat(the_folder,'animal_traits_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'),res{2});
            write_into_files_W(strcat(the_folder,'plant_traits_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'),res{3});
            write_into_files_W(strcat(the_folder,'animal_abundance_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'),res{4});
            write_into_files_W(strcat(the_folder,'plant_abundance_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'),res{5});
            write_into_files_W(strcat(the_folder,'sepIndex_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'),res{6});

            %test if there are empty niche
            time_plot=res{1};
            trait_plot_an=res{2};
            trait_plot_pl=res{3};
            pop_plot_an=res{4};
            pop_plot_pl=res{5};
            separating_index=res{6};
            last_tr_an=trait_plot_an{end};
            last_tr_pl=trait_plot_pl{end};
            last_pop_an=pop_plot_an{end};
            last_pop_pl=pop_plot_pl{end};
            EN_value=detect_EN_function_robustness(last_tr_an,last_tr_pl,last_pop_an,last_pop_pl,log_sig_an,log_sig_c_an,log_sig_m,rA,rP,const,h,x0,y0,kx,ky);
            EN_value=EN_value{1};
            %continue with invasion
            if EN_value~=0
                disp ('simulation of invasion after ESS');
                res2=simulation_for_robustness_test_invasion(rA,rP,const,x0,y0,kx,ky,sig_an,log_sig_c_an,log_sig_m,time_plot,trait_plot_an,trait_plot_pl,pop_plot_an,pop_plot_pl,separating_index);
                clear_files(strcat(the_inv_folder,'time_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
                clear_files(strcat(the_inv_folder,'animal_traits_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
                clear_files(strcat(the_inv_folder,'plant_traits_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
                clear_files(strcat(the_inv_folder,'animal_abundance_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
                clear_files(strcat(the_inv_folder,'plant_abundance_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
                clear_files(strcat(the_inv_folder,'sepIndex_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
                clear_files(strcat(the_inv_folder,'introIndex_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
                clear_files(strcat(the_inv_folder,'interval_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));

                write_into_files_W(strcat(the_inv_folder,'time_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'),res2{1});
                write_into_files_W(strcat(the_inv_folder,'animal_traits_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'),res2{2});
                write_into_files_W(strcat(the_inv_folder,'plant_traits_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'),res2{3});
                write_into_files_W(strcat(the_inv_folder,'animal_abundance_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'),res2{4});
                write_into_files_W(strcat(the_inv_folder,'plant_abundance_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'),res2{5});
                write_into_files_W(strcat(the_inv_folder,'sepIndex_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'),res2{6});
                write_into_files_W(strcat(the_inv_folder,'introIndex_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'),res2{7});
                write_into_files_W(strcat(the_inv_folder,'interval_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'),res2{8});
            end
        end
    end

end

%Vary the constant scaling the intensity of mutualistic interactions (c)
disp ('Vary the constant scaling the intensity of mutualistic interactions');
sig_an=1;
log_sig_an=sig_an;
const_range=[0.05 0.2];

all_res_opt={'intra','cross'};

for res_opt_ind=1:length(all_res_opt)
    res_opt=all_res_opt{res_opt_ind};
    for const=const_range
        the_folder=strcat('./Results_robustness/vary_const=',num2str(const),'_',res_opt,'/');
        mkdir(the_folder);
        the_inv_folder=strcat('./Results_robustness/vary_const=',num2str(const),'_',res_opt,'_invasion/');
        mkdir(the_inv_folder);
        for log_sig_c_an=log_sig_c_an_range
            for log_sig_m=log_sig_m_range
                disp (strcat(res_opt,',',num2str(log_sig_c_an),',',num2str(log_sig_m)));
                sig_c_an=exp(log_sig_c_an);
                sig_c_pl=sig_c_an;
                sig_m=exp(log_sig_m);
                res=simulation_for_robustness_test_function(rA,rP,const,x0,y0,kx,ky,sig_an,log_sig_c_an,log_sig_m,res_opt);
                clear_files(strcat(the_folder,'time_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
                clear_files(strcat(the_folder,'animal_traits_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
                clear_files(strcat(the_folder,'plant_traits_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
                clear_files(strcat(the_folder,'animal_abundance_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
                clear_files(strcat(the_folder,'plant_abundance_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
                clear_files(strcat(the_folder,'sepIndex_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
                write_into_files_W(strcat(the_folder,'time_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'),res{1});
                write_into_files_W(strcat(the_folder,'animal_traits_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'),res{2});
                write_into_files_W(strcat(the_folder,'plant_traits_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'),res{3});
                write_into_files_W(strcat(the_folder,'animal_abundance_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'),res{4});
                write_into_files_W(strcat(the_folder,'plant_abundance_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'),res{5});
                write_into_files_W(strcat(the_folder,'sepIndex_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'),res{6});

                %test if there are empty niche
                time_plot=res{1};
                trait_plot_an=res{2};
                trait_plot_pl=res{3};
                pop_plot_an=res{4};
                pop_plot_pl=res{5};
                separating_index=res{6};
                last_tr_an=trait_plot_an{end};
                last_tr_pl=trait_plot_pl{end};
                last_pop_an=pop_plot_an{end};
                last_pop_pl=pop_plot_pl{end};
                EN_value=detect_EN_function_robustness(last_tr_an,last_tr_pl,last_pop_an,last_pop_pl,log_sig_an,log_sig_c_an,log_sig_m,rA,rP,const,h,x0,y0,kx,ky);
                EN_value=EN_value{1};
                %continue with invasion
                if EN_value~=0
                    disp ('simulation of invasion after ESS');
                    res2=simulation_for_robustness_test_invasion(rA,rP,const,x0,y0,kx,ky,sig_an,log_sig_c_an,log_sig_m,time_plot,trait_plot_an,trait_plot_pl,pop_plot_an,pop_plot_pl,separating_index);
                    clear_files(strcat(the_inv_folder,'time_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
                    clear_files(strcat(the_inv_folder,'animal_traits_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
                    clear_files(strcat(the_inv_folder,'plant_traits_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
                    clear_files(strcat(the_inv_folder,'animal_abundance_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
                    clear_files(strcat(the_inv_folder,'plant_abundance_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
                    clear_files(strcat(the_inv_folder,'sepIndex_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
                    clear_files(strcat(the_inv_folder,'introIndex_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
                    clear_files(strcat(the_inv_folder,'interval_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));

                    write_into_files_W(strcat(the_inv_folder,'time_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'),res2{1});
                    write_into_files_W(strcat(the_inv_folder,'animal_traits_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'),res2{2});
                    write_into_files_W(strcat(the_inv_folder,'plant_traits_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'),res2{3});
                    write_into_files_W(strcat(the_inv_folder,'animal_abundance_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'),res2{4});
                    write_into_files_W(strcat(the_inv_folder,'plant_abundance_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'),res2{5});
                    write_into_files_W(strcat(the_inv_folder,'sepIndex_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'),res2{6});
                    write_into_files_W(strcat(the_inv_folder,'introIndex_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'),res2{7});
                    write_into_files_W(strcat(the_inv_folder,'interval_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'),res2{8});
                end
            end
        end
    end
end

%Vary the constant scaling the carrying capacity (K_A)
disp ('Vary the constant scaling the carrying capacity');
const=0.1;
kx=200;
all_res_opt={'intra','cross'};

for res_opt_ind=1:length(all_res_opt)
    res_opt=all_res_opt{res_opt_ind};
    the_folder=strcat('./Results_robustness/vary_kx=',num2str(kx),'_',res_opt,'/');
    mkdir(the_folder);
    the_inv_folder=strcat('./Results_robustness/vary_kx=',num2str(kx),'_',res_opt,'_invasion/');
    mkdir(the_inv_folder);
    for log_sig_c_an=log_sig_c_an_range
        for log_sig_m=log_sig_m_range
            disp (strcat(res_opt,',',num2str(log_sig_c_an),',',num2str(log_sig_m)));
            sig_c_an=exp(log_sig_c_an);
            sig_c_pl=sig_c_an;
            sig_m=exp(log_sig_m);
            res=simulation_for_robustness_test_function(rA,rP,const,x0,y0,kx,ky,sig_an,log_sig_c_an,log_sig_m,res_opt);
            clear_files(strcat(the_folder,'time_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
            clear_files(strcat(the_folder,'animal_traits_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
            clear_files(strcat(the_folder,'plant_traits_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
            clear_files(strcat(the_folder,'animal_abundance_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
            clear_files(strcat(the_folder,'plant_abundance_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
            clear_files(strcat(the_folder,'sepIndex_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
            write_into_files_W(strcat(the_folder,'time_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'),res{1});
            write_into_files_W(strcat(the_folder,'animal_traits_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'),res{2});
            write_into_files_W(strcat(the_folder,'plant_traits_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'),res{3});
            write_into_files_W(strcat(the_folder,'animal_abundance_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'),res{4});
            write_into_files_W(strcat(the_folder,'plant_abundance_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'),res{5});
            write_into_files_W(strcat(the_folder,'sepIndex_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'),res{6});

            %test if there are empty niche
            time_plot=res{1};
            trait_plot_an=res{2};
            trait_plot_pl=res{3};
            pop_plot_an=res{4};
            pop_plot_pl=res{5};
            separating_index=res{6};
            last_tr_an=trait_plot_an{end};
            last_tr_pl=trait_plot_pl{end};
            last_pop_an=pop_plot_an{end};
            last_pop_pl=pop_plot_pl{end};
            EN_value=detect_EN_function_robustness(last_tr_an,last_tr_pl,last_pop_an,last_pop_pl,log_sig_an,log_sig_c_an,log_sig_m,rA,rP,const,h,x0,y0,kx,ky);
            EN_value=EN_value{1};
            %continue with invasion
            if EN_value~=0
                disp ('simulation of invasion after ESS');
                res2=simulation_for_robustness_test_invasion(rA,rP,const,x0,y0,kx,ky,sig_an,log_sig_c_an,log_sig_m,time_plot,trait_plot_an,trait_plot_pl,pop_plot_an,pop_plot_pl,separating_index);
                clear_files(strcat(the_inv_folder,'time_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
                clear_files(strcat(the_inv_folder,'animal_traits_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
                clear_files(strcat(the_inv_folder,'plant_traits_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
                clear_files(strcat(the_inv_folder,'animal_abundance_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
                clear_files(strcat(the_inv_folder,'plant_abundance_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
                clear_files(strcat(the_inv_folder,'sepIndex_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
                clear_files(strcat(the_inv_folder,'introIndex_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
                clear_files(strcat(the_inv_folder,'interval_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));

                write_into_files_W(strcat(the_inv_folder,'time_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'),res2{1});
                write_into_files_W(strcat(the_inv_folder,'animal_traits_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'),res2{2});
                write_into_files_W(strcat(the_inv_folder,'plant_traits_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'),res2{3});
                write_into_files_W(strcat(the_inv_folder,'animal_abundance_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'),res2{4});
                write_into_files_W(strcat(the_inv_folder,'plant_abundance_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'),res2{5});
                write_into_files_W(strcat(the_inv_folder,'sepIndex_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'),res2{6});
                write_into_files_W(strcat(the_inv_folder,'introIndex_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'),res2{7});
                write_into_files_W(strcat(the_inv_folder,'interval_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'),res2{8});
            end
        end
    end
end

%Vary optimum trait value for the resource kernel (X_A_max)
disp ('Vary the constant scaling the carrying capacity');
kx=400;
x0=2.5;
all_res_opt={'intra','cross'};

for res_opt_ind=1:length(all_res_opt)
    res_opt=all_res_opt{res_opt_ind};
    the_folder=strcat('./Results_robustness/vary_x0=',num2str(x0),'_',res_opt,'/');
    mkdir(the_folder);
    the_inv_folder=strcat('./Results_robustness/vary_x0=',num2str(x0),'_',res_opt,'_invasion/');
    mkdir(the_inv_folder);
    for log_sig_c_an=log_sig_c_an_range
        for log_sig_m=log_sig_m_range
            disp (strcat(res_opt,',',num2str(log_sig_c_an),',',num2str(log_sig_m)));
            sig_c_an=exp(log_sig_c_an);
            sig_c_pl=sig_c_an;
            sig_m=exp(log_sig_m);
            res=simulation_for_robustness_test_function(rA,rP,const,x0,y0,kx,ky,sig_an,log_sig_c_an,log_sig_m,res_opt);
            clear_files(strcat(the_folder,'time_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
            clear_files(strcat(the_folder,'animal_traits_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
            clear_files(strcat(the_folder,'plant_traits_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
            clear_files(strcat(the_folder,'animal_abundance_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
            clear_files(strcat(the_folder,'plant_abundance_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
            clear_files(strcat(the_folder,'sepIndex_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
            write_into_files_W(strcat(the_folder,'time_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'),res{1});
            write_into_files_W(strcat(the_folder,'animal_traits_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'),res{2});
            write_into_files_W(strcat(the_folder,'plant_traits_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'),res{3});
            write_into_files_W(strcat(the_folder,'animal_abundance_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'),res{4});
            write_into_files_W(strcat(the_folder,'plant_abundance_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'),res{5});
            write_into_files_W(strcat(the_folder,'sepIndex_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'),res{6});

            %test if there are empty niche
            time_plot=res{1};
            trait_plot_an=res{2};
            trait_plot_pl=res{3};
            pop_plot_an=res{4};
            pop_plot_pl=res{5};
            separating_index=res{6};
            last_tr_an=trait_plot_an{end};
            last_tr_pl=trait_plot_pl{end};
            last_pop_an=pop_plot_an{end};
            last_pop_pl=pop_plot_pl{end};
            EN_value=detect_EN_function_robustness(last_tr_an,last_tr_pl,last_pop_an,last_pop_pl,log_sig_an,log_sig_c_an,log_sig_m,rA,rP,const,h,x0,y0,kx,ky);
            EN_value=EN_value{1};
            %continue with invasion
            if EN_value~=0
                disp ('simulation of invasion after ESS');
                res2=simulation_for_robustness_test_invasion(rA,rP,const,x0,y0,kx,ky,sig_an,log_sig_c_an,log_sig_m,time_plot,trait_plot_an,trait_plot_pl,pop_plot_an,pop_plot_pl,separating_index);
                clear_files(strcat(the_inv_folder,'time_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
                clear_files(strcat(the_inv_folder,'animal_traits_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
                clear_files(strcat(the_inv_folder,'plant_traits_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
                clear_files(strcat(the_inv_folder,'animal_abundance_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
                clear_files(strcat(the_inv_folder,'plant_abundance_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
                clear_files(strcat(the_inv_folder,'sepIndex_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
                clear_files(strcat(the_inv_folder,'introIndex_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
                clear_files(strcat(the_inv_folder,'interval_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));

                write_into_files_W(strcat(the_inv_folder,'time_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'),res2{1});
                write_into_files_W(strcat(the_inv_folder,'animal_traits_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'),res2{2});
                write_into_files_W(strcat(the_inv_folder,'plant_traits_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'),res2{3});
                write_into_files_W(strcat(the_inv_folder,'animal_abundance_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'),res2{4});
                write_into_files_W(strcat(the_inv_folder,'plant_abundance_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'),res2{5});
                write_into_files_W(strcat(the_inv_folder,'sepIndex_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'),res2{6});
                write_into_files_W(strcat(the_inv_folder,'introIndex_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'),res2{7});
                write_into_files_W(strcat(the_inv_folder,'interval_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'),res2{8});
            end
        end
    end
end

