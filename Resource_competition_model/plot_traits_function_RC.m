function figg=plot_traits_function(sep_index,time_plot,an_traits,an_pop,r,sig_c,sig_an,x0,kx)
    figg=figure('Visible','on');
    subplot('position',[0.09 0.1 0.4 0.85]);
    %plot the animal
    sol_plot=an_traits;
    count_time_plot=1;
    current_length=length(sol_plot{1});
    i=1;
    while (count_time_plot<=length(time_plot))
        row_each_mat=1;
        temp_array=zeros(length(time_plot),length(sol_plot(end)));
        temp_array(:)=NaN;
        temp_time_plot=zeros(size(time_plot));
        temp_time_plot(:)=NaN;
        the_last_index=sep_index(i);
        while count_time_plot<=length(time_plot) && count_time_plot<=the_last_index
            temp_array(row_each_mat,1:current_length)=sol_plot{count_time_plot};
            temp_time_plot(row_each_mat)=time_plot(count_time_plot);
            row_each_mat=row_each_mat+1;
            count_time_plot=count_time_plot+1;
        end
        temp_array=temp_array(~isnan(temp_array(:,1)),~isnan(temp_array(1,:)));
        temp_time_plot=temp_time_plot(~isnan(temp_time_plot));
        plot(temp_time_plot,temp_array,'b','LineWidth',1.2);
        hold all;
        if count_time_plot<=length(time_plot)
            current_length=length(sol_plot{count_time_plot});
        end
        i=i+1;
        
    end
    
    min_trait=min(cellfun(@(x) min(x(:)),an_traits)) ;
    max_trait=max(cellfun(@(x) max(x(:)),an_traits)) ;
    
    axis([time_plot(1) time_plot(end) min_trait-0.1 max_trait+0.1]);
    
    ylabel('Trait values','FontSize',12);
    xlabel('Evolutionary time','FontSize',12);

    set(gca,'FontSize',12,'LineWidth',1.1);
   
    
    subplot('position',[0.58 0.55 0.4 0.35 ]);
    tr_an_last=(an_traits{end});
    ab_an_last=(an_pop{end});
    
   
    trait_interval=min_trait-6:0.001:max_trait+6;
    plot(tr_an_last,fitness_land_function_RC(tr_an_last,ab_an_last,tr_an_last,r,sig_c,sig_an,x0,kx),'o','MarkerSize',6,'MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7]);
    hold all;
    plot(trait_interval,fitness_land_function_RC(trait_interval,ab_an_last,tr_an_last,r,sig_c,sig_an,x0,kx),'k-','LineWidth',1.2);
    xlabel('Ttrait values','FontSize',11);
    ylabel('fitness','FontSize',11);
    set(gca,'xlim',[1.8 4.2]);
    set(gca,'ylim',[-0.01 0.01]);
    

end