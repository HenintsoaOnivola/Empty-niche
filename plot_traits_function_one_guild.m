function figg=plot_traits_function_one_guild(introduction_index,sep_index,time_plot,an_traits,ab_an,ab_pl,the_y_label,trait_limit)
    figg=figure('Visible','off');
    
    the_hot_colormap=hot(110);
    the_hot_colormap=flipud(the_hot_colormap);
    the_hot_colormap=the_hot_colormap(25:105,:);
    colormap(the_hot_colormap);
    
    the_min_ab=min([min(cellfun(@(x) min(x(:)),ab_an))  min(cellfun(@(x) min(x(:)),ab_pl))]);
    the_max_ab=max([max(cellfun(@(x) max(x(:)),ab_an))  max(cellfun(@(x) max(x(:)),ab_pl))]) ;
    
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
        temp_array_abundance=zeros(length(time_plot),length(sol_plot(end)));
        temp_array_abundance(:)=NaN;
        the_last_index=sep_index(i);
        while count_time_plot<=length(time_plot) && count_time_plot<=the_last_index
            temp_array(row_each_mat,1:current_length)=sol_plot{count_time_plot};
            temp_time_plot(row_each_mat)=time_plot(count_time_plot);
            temp_array_abundance(row_each_mat,1:current_length)=ab_an{count_time_plot};
            row_each_mat=row_each_mat+1;
            count_time_plot=count_time_plot+1;
        end
        temp_array=temp_array(~isnan(temp_array(:,1)),~isnan(temp_array(1,:)));
        temp_array_abundance=temp_array_abundance(~isnan(temp_array_abundance(:,1)),~isnan(temp_array_abundance(1,:)));
        temp_time_plot=temp_time_plot(~isnan(temp_time_plot));
        for each_line=1:size(temp_array,2)
            ab_normalized=((temp_array_abundance(:,each_line)-the_min_ab)/(the_max_ab-the_min_ab));
            color_line(temp_time_plot,temp_array(:,each_line),ab_normalized,'LineWidth',1.8);
            hold all;
        end
        if count_time_plot<=length(time_plot)
            current_length=length(sol_plot{count_time_plot});
            if any(count_time_plot==introduction_index)
                if length(sol_plot{count_time_plot})>length(sol_plot{count_time_plot-1})
                    inv_time=time_plot(count_time_plot);
                    inv_tr=sol_plot{count_time_plot}(~ismember(sol_plot{count_time_plot},sol_plot{count_time_plot-1}));
                    plot(inv_time,inv_tr,'bo','MarkerSize',7,'MarkerFaceColor','b');
                    hold all;
                end
            end
        end
        i=i+1;
        
    end
    set(gca,'TickDir','in','TickLength',[0.015 0.015],'YMinorTick','off','XMinorTick','off','layer','top','LineWidth',1.2,'FontSize',12,'Ytick',2.5:0.5:3.5);
    box ('on');
    xlabel('Evolutionary time','FontSize',20);
    ylabel(the_y_label,'FontSize',20);
    set(gca,'ylim',trait_limit);
   
end