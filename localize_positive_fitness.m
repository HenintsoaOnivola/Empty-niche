function res=localize_positive_fitness(last_pop_an,last_pop_pl,last_tr_an,last_tr_pl,rA,rP,h,const,log_sig_c_an,log_sig_c_pl,sig_an,sig_pl,x0,y0,kx,ky,log_sig_m)

    sig_c_an=exp(log_sig_c_an);
    sig_c_pl=exp(log_sig_c_pl);
    sig_m=exp(log_sig_m);

    min_tr=1.5;
    max_tr=4;
        
    current=min_tr-0.1;
    cont=true;
    while cont==true
        fit_an=fitness_land_function(current,last_pop_an,last_pop_pl,last_tr_an,last_tr_pl,rA,h,const,sig_c_an,sig_an,x0,kx,sig_m);
        fit_pl=fitness_land_function(current,last_pop_pl,last_pop_an,last_tr_pl,last_tr_an,rP,h,const,sig_c_pl,sig_pl,y0,ky,sig_m);
        if (fit_an<=0) && (fit_pl<=0) 
            cont=false;
        else
           cont=true;
           current=current-0.1;
        end
    end
    min_interval=current;            


    current=max_tr+0.1;
    cont=true;
    while cont==true
        fit_an=fitness_land_function(current,last_pop_an,last_pop_pl,last_tr_an,last_tr_pl,rA,h,const,sig_c_an,sig_an,x0,kx,sig_m);
        fit_pl=fitness_land_function(current,last_pop_pl,last_pop_an,last_tr_pl,last_tr_an,rP,h,const,sig_c_pl,sig_pl,y0,ky,sig_m);
        if (fit_an<=0 ) && (fit_pl<=0 ) 
            cont=false;
        else
           cont=true;
           current=current+0.1;
        end
    end
    max_interval=current; 
    trait_interval=min_interval:0.001:max_interval+0.1;
    
   
    all_fit_an=fitness_land_function(trait_interval,last_pop_an,last_pop_pl,last_tr_an,last_tr_pl,rA,h,const,sig_c_an,sig_an,x0,kx,sig_m);
    all_fit_pl=fitness_land_function(trait_interval,last_pop_pl,last_pop_an,last_tr_pl,last_tr_an,rP,h,const,sig_c_pl,sig_pl,y0,ky,sig_m);
    
    interval_an=trait_interval(all_fit_an>0.01);
    interval_pl=trait_interval(all_fit_pl>0.01);
    
    if length(interval_an)==1
        interval_an=[];
    end
    if length(interval_pl)==1
        interval_pl=[];
    end
   
    res={interval_an interval_pl};
        

end