function figg=plot_traits_function_with_fitness(sep_index,time_plot,an_traits,pl_traits,an_pop,pl_pop,rA,rP,h,const,sig_c_an,sig_c_pl,sig_an,sig_pl,x0,y0,kx,ky,sig_m)
    figg=figure('Visible','on');
    subplot('position',[0.1 0.1 0.52 0.85]);
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
        plot((temp_time_plot),temp_array,'b','LineWidth',1.2);
        hold all;
        if count_time_plot<=length(time_plot)
            current_length=length(sol_plot{count_time_plot});
        end
        i=i+1;
        
    end
    %plot the plants
    sol_plot=pl_traits;
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
        plot((temp_time_plot),temp_array,'r','LineWidth',1.2);
        hold all;
        if count_time_plot<=length(time_plot)
            current_length=length(sol_plot{count_time_plot});
        end
        i=i+1;
    end
    
    min_trait=min([min(cellfun(@(x) min(x(:)),an_traits))  min(cellfun(@(x) min(x(:)),pl_traits))]);
    max_trait=max([max(cellfun(@(x) max(x(:)),an_traits))  max(cellfun(@(x) max(x(:)),pl_traits))]);
    
    axis([time_plot(1) time_plot(end) min_trait-0.1 max_trait+0.1]);
    
    ylabel('Trait values','FontSize',15);
    xlabel ('Evolutionary time','FontSize',15);

    set(gca,'FontSize',17,'LineWidth',1.1);
    
    
    %plotting fitness
    subplot('position',[0.725 0.55 0.25 0.35 ]);
    tr_an_last=an_traits{end};
    ab_an_last=an_pop{end};
    tr_pl_last=pl_traits{end};
    ab_pl_last=pl_pop{end};
    trait_interval=1.1:0.001: 4.1;
    plot(tr_an_last,fitness_land_function(tr_an_last,ab_an_last,ab_pl_last,tr_an_last,tr_pl_last,rA,h,const,sig_c_an,sig_an,x0,kx,sig_m),'bo','MarkerSize',8,'MarkerFaceColor','b');
    hold on;
    plot(trait_interval,fitness_land_function(trait_interval,ab_an_last,ab_pl_last,tr_an_last,tr_pl_last,rA,h,const,sig_c_an,sig_an,x0,kx,sig_m),'b-');
    hold on;
    plot(tr_pl_last,fitness_land_function(tr_pl_last,ab_pl_last,ab_an_last,tr_pl_last,tr_an_last,rP,h,const,sig_c_pl,sig_pl,y0,ky,sig_m),'ro','MarkerSize',8,'MarkerFaceColor','r');
    hold on;
    plot(trait_interval,fitness_land_function(trait_interval,ab_pl_last,ab_an_last,tr_pl_last,tr_an_last,rP,h,const,sig_c_pl,sig_pl,y0,ky,sig_m),'r-');
   
    set(gca,'ylim',[-0.2 0.1]);
    set (gca,'xlim',[trait_interval(1) trait_interval(end)]);
    ylabel('fitness','FontSize',12);
    
   
end