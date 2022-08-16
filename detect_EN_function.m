function EN_value=detect_EN_function(the_folder,log_sig_an,log_sig_c_an,log_sig_m,rA,rP,const,h,x0,y0,kx,ky)

    sig_an=log_sig_an;
    sig_pl=sig_an;
    
    log_sig_c_pl=log_sig_c_an;
    sig_c_an=exp(log_sig_c_an);
    sig_c_pl=sig_c_an;

    sig_m=exp(log_sig_m);

    fil_name_an_tr=strcat(the_folder,'animal_traits_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat');
    fil_name_pl_tr=strcat(the_folder,'plant_traits_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat');

    fil_name_an_pop=strcat(the_folder,'animal_abundance_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat');
    fil_name_pl_pop=strcat(the_folder,'plant_abundance_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat');

    last_tr_an=read_last_line_to_array(fil_name_an_tr);
    last_tr_pl=read_last_line_to_array(fil_name_pl_tr);
    last_pop_an=read_last_line_to_array(fil_name_an_pop);
    last_pop_pl=read_last_line_to_array(fil_name_pl_pop);

    peaks_inter=localize_positive_fitness(last_pop_an,last_pop_pl,last_tr_an,last_tr_pl,rA,rP,h,const,log_sig_c_an,log_sig_c_pl,sig_an,sig_pl,x0,y0,kx,ky,log_sig_m);
    the_interval_an=peaks_inter{2};
    the_interval_pl=peaks_inter{3};

    the_fit_an= fitness_land_function(last_tr_an,last_pop_an,last_pop_pl,last_tr_an,last_tr_pl,rA,h,const,sig_c_an,sig_an,x0,kx,sig_m);
    the_fit_pl= fitness_land_function(last_tr_pl,last_pop_pl,last_pop_an,last_tr_pl,last_tr_an,rP,h,const,sig_c_pl,sig_pl,y0,ky,sig_m);

    if (~isempty(the_interval_an) || ~isempty(the_interval_pl))
        EN_value=1;
        if (sum (abs(the_fit_an)>10^(-3))~=0) ||  (sum (abs(the_fit_pl)>10^(-3))~=0)  %there is oscillation here
            EN_value=2;
        end
    else
        EN_value=0;

    end

end

