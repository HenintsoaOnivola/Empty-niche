function res=localize_interval_peaks_scaled_param_positive_fitness(last_pop_an,last_pop_pl,last_tr_an,last_tr_pl,rA,rP,h,const,log_sig_c_an,log_sig_c_pl,sig_an,sig_pl,x0,y0,kx,ky,log_sig_m)


    sig_c_an=exp(log_sig_c_an);
    sig_c_pl=exp(log_sig_c_pl);
    sig_m=exp(log_sig_m);

    %[min_tr ind_min]=min([min(last_tr_an) min(last_tr_pl)]);
    %[max_tr ind_max]=max([max(last_tr_an) max(last_tr_pl)]);
    
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
    
%     disp (max(all_fit_an));
%     disp (max(all_fit_pl));
    
  

    an_left_peaks=zeros(size(all_fit_an));
    an_left_peaks(:)=NaN;
    an_left_fit=zeros(size(all_fit_an));
    an_left_fit(:)=NaN;
    count_peak_an_left=1;
    an_right_peaks=zeros(size(all_fit_an));
    an_right_peaks(:)=NaN;
    an_right_fit=zeros(size(all_fit_an));
    an_right_fit(:)=NaN;
    count_peak_an_right=1;

    pl_left_peaks=zeros(size(all_fit_pl));
    pl_left_peaks(:)=NaN;
    pl_left_fit=zeros(size(all_fit_pl));
    pl_left_fit(:)=NaN;
    count_peak_pl_left=1;
    pl_right_peaks=zeros(size(all_fit_pl));
    pl_right_peaks(:)=NaN;
    pl_right_fit=zeros(size(all_fit_pl));
    pl_right_fit(:)=NaN;
    count_peak_pl_right=1;

    for ss=2:length(trait_interval)-1
        the_peak=trait_interval(ss);

        the_fit=all_fit_an(ss);
        if the_fit>10^(-4) &&  all_fit_an(ss-1)<the_fit  && the_fit>all_fit_an(ss+1)
            if  the_peak<min(last_tr_an)-0.01  
                an_left_peaks(count_peak_an_left)=the_peak;
                an_left_fit(count_peak_an_left)=the_fit;
                count_peak_an_left=count_peak_an_left+1;
            end
            if the_peak>max(last_tr_an)+0.01  
                an_right_peaks(count_peak_an_right)=the_peak;
                an_right_fit(count_peak_an_right)=the_fit;
                count_peak_an_right=count_peak_an_right+1;
            end
        end

        the_fit=all_fit_pl(ss);
        if the_fit>=10^(-4)  && all_fit_pl(ss-1)<the_fit  && the_fit>all_fit_pl(ss+1)
            if  the_peak<min(last_tr_pl)-0.01  
                pl_left_peaks(count_peak_pl_left)=the_peak;
                pl_left_fit(count_peak_pl_left)=the_fit;
                count_peak_pl_left=count_peak_pl_left+1;
            end
            if the_peak>max(last_tr_pl)+0.01  
                pl_right_peaks(count_peak_pl_right)=the_peak;
                pl_right_fit(count_peak_pl_right)=the_fit;
                count_peak_pl_right=count_peak_pl_right+1;
            end
        end
    end



    an_left_peaks=an_left_peaks(~isnan(an_left_peaks)); 
    an_right_peaks=an_right_peaks(~isnan(an_right_peaks)); 

    pl_left_peaks=pl_left_peaks(~isnan(pl_left_peaks)); 
    pl_right_peaks=pl_right_peaks(~isnan(pl_right_peaks)); 

    if ~any(abs(an_left_fit-0.1)<10^(-5))
        all_left_fit_an=all_fit_an(trait_interval<min(last_tr_an)-0.01);
        all_left_inter_an=trait_interval(trait_interval<min(last_tr_an)-0.01);
        all_left_fit_an_pos=all_left_fit_an(all_left_fit_an>0);
        all_left_inter_an_pos=all_left_inter_an(all_left_fit_an>0);
        if ~isempty(all_left_fit_an_pos)
            if abs(all_left_fit_an_pos(1)-0.1)<10^(-5)
                an_left_peaks=[an_left_peaks all_left_inter_an_pos(1)];
            end
        end
    end

    if ~any(abs(an_right_fit-0.1)<10^(-5))
        all_right_fit_an=all_fit_an(trait_interval>max(last_tr_an)+0.01);
        all_right_inter_an=trait_interval(trait_interval>max(last_tr_an)+0.01);
        all_right_fit_an_pos=all_right_fit_an(all_right_fit_an>0);
        all_right_inter_an_pos=all_right_inter_an(all_right_fit_an>0);
        if ~isempty(all_right_fit_an_pos)
            if abs(all_right_fit_an_pos(end)-0.1)<10^(-5)
                an_right_peaks=[an_right_peaks all_right_inter_an_pos(end)];
            end
        end
    end

    if ~any(abs(pl_left_fit-0.1)<10^(-5))
        all_left_fit_pl=all_fit_pl(trait_interval<min(last_tr_pl)-0.01);
        all_left_inter_pl=trait_interval(trait_interval<min(last_tr_pl)-0.01);
        all_left_fit_pl_pos=all_left_fit_pl(all_left_fit_pl>0);
        all_left_inter_pl_pos=all_left_inter_pl(all_left_fit_pl>0);
        if ~isempty(all_left_fit_pl_pos)
            if abs(all_left_fit_pl_pos(1)-0.1)<10^(-5)
                pl_left_peaks=[pl_left_peaks all_left_inter_pl_pos(1)];
            end
        end
    end

    if ~any(abs(pl_right_fit-0.1)<10^(-5))
        all_right_fit_pl=all_fit_pl(trait_interval>max(last_tr_pl)+0.01);
        all_right_inter_pl=trait_interval(trait_interval>max(last_tr_pl)+0.01);
        all_right_fit_pl_pos=all_right_fit_pl(all_right_fit_pl>0);
        all_right_inter_pl_pos=all_right_inter_pl(all_right_fit_pl>0);
        if ~isempty(all_right_fit_pl_pos)
            if abs(all_right_fit_pl_pos(end)-0.1)<10^(-5)
                pl_right_peaks=[pl_right_peaks all_right_inter_pl_pos(end)];
            end
        end
    end

    an_peaks=[an_left_peaks an_right_peaks];
    pl_peaks=[pl_left_peaks pl_right_peaks];
    
    interval_an=trait_interval(all_fit_an>0.01);
    %disp (all_fit_an(all_fit_an>10^(-4)));
    interval_pl=trait_interval(all_fit_pl>0.01);
   % disp (all_fit_pl(all_fit_pl>10^(-4)));
    
%     if log_sig_c_an==-3
%         interval_an=trait_interval(all_fit_an>0.3);
%         interval_pl=trait_interval(all_fit_pl>0.3);
%     end
    
    if length(interval_an)==1
        interval_an=[];
    end
    if length(interval_pl)==1
        interval_pl=[];
    end

    peaks={an_peaks; pl_peaks};
    res={peaks interval_an interval_pl};
        

end