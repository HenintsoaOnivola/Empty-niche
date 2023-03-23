function res_array=detect_EN_function_robustness(last_tr_an,last_tr_pl,last_pop_an,last_pop_pl,log_sig_an,log_sig_c_an,log_sig_m,rA,rP,const,h,x0,y0,kx,ky)

    sig_an=log_sig_an;
    sig_pl=sig_an;
    
    log_sig_c_pl=log_sig_c_an;
    sig_c_an=exp(log_sig_c_an);
    sig_c_pl=sig_c_an;

    sig_m=exp(log_sig_m);

    pos_fit_inter=localize_positive_fitness(last_pop_an,last_pop_pl,last_tr_an,last_tr_pl,rA,rP,h,const,log_sig_c_an,log_sig_c_pl,sig_an,sig_pl,x0,y0,kx,ky,log_sig_m);
    the_interval_an=pos_fit_inter{1};
    the_interval_pl=pos_fit_inter{2};

    the_fit_an= fitness_land_function(last_tr_an,last_pop_an,last_pop_pl,last_tr_an,last_tr_pl,rA,h,const,sig_c_an,sig_an,x0,kx,sig_m);
    the_fit_pl= fitness_land_function(last_tr_pl,last_pop_pl,last_pop_an,last_tr_pl,last_tr_an,rP,h,const,sig_c_pl,sig_pl,y0,ky,sig_m);

    if (~isempty(the_interval_an) || ~isempty(the_interval_pl))
        EN_value=1;
        if (sum (abs(the_fit_an)>10^(-3))~=0) ||  (sum (abs(the_fit_pl)>10^(-3))~=0)  %there is oscillation here
            EN_value=2;
        end
    else
        EN_value=0;
        en_type='';

    end
    if any(the_interval_an>max(last_tr_an)) || any (the_interval_pl<min(last_tr_pl)) 
        en_type='intra';
    end
    if ~isempty (the_interval_an)
        if (min(last_tr_an)<min(the_interval_an)) && (max(last_tr_an)>max(the_interval_an)) 
            en_type='cross';
        end
    end
    if ~isempty (the_interval_pl)
        if (min(last_tr_pl)<min(the_interval_pl)) && (max(last_tr_pl)>max(the_interval_pl))
            en_type='cross';
        end
    end
    if any(the_interval_an<min(last_tr_an)) || any (the_interval_pl>max(last_tr_pl))
        en_type='cross';
    end
    if (any(the_interval_an>max(last_tr_an)) || any (the_interval_pl<min(last_tr_pl))) && (any(the_interval_an<min(last_tr_an)) || any (the_interval_pl>max(last_tr_pl)))
        en_type='both';
    end
    
    res_array={EN_value;en_type;the_interval_an;the_interval_pl};

end

