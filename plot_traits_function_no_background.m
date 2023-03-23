function figg=plot_traits_function_no_background(introduction_index,sep_index,time_plot,an_traits,pl_traits,ab_an,ab_pl)
    figg=figure('Visible','off');
    
    the_hot_colormap=hot(110);
    the_hot_colormap=flipud(the_hot_colormap);
    the_hot_colormap=the_hot_colormap(25:105,:);
    colormap(the_hot_colormap);
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%plotting the trait evolution %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    h3=subplot(2,2,3);
    the_min_ab=min([min(cellfun(@(x) min(x(:)),ab_an))  min(cellfun(@(x) min(x(:)),ab_pl))]);
    the_max_ab=max([max(cellfun(@(x) max(x(:)),ab_an))  max(cellfun(@(x) max(x(:)),ab_pl))]) ;
    
    min_trait_an=min(cellfun(@(x) min(x(:)),an_traits));
    max_trait_an=max(cellfun(@(x) max(x(:)),an_traits));
    
    min_trait_pl=min(cellfun(@(x) min(x(:)),pl_traits));
    max_trait_pl=max(cellfun(@(x) max(x(:)),pl_traits));
    
    trait_inter_an=[min_trait_an max_trait_an];
    trait_inter_pl=[min_trait_pl max_trait_pl];

    
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
            color_line(temp_array(:,each_line),temp_time_plot,ab_normalized,'LineWidth',1.8);
            hold all;
        end
        if count_time_plot<=length(time_plot)
            current_length=length(sol_plot{count_time_plot});
            if any(count_time_plot==introduction_index)
                if length(sol_plot{count_time_plot})>length(sol_plot{count_time_plot-1})
                    inv_time=time_plot(count_time_plot);
                    inv_tr=sol_plot{count_time_plot}(~ismember(sol_plot{count_time_plot},sol_plot{count_time_plot-1}));
                    plot(inv_tr,inv_time,'bo','MarkerSize',7,'MarkerFaceColor','b');
                    hold all;
                end
            end
        end
        i=i+1;
        
    end

    axis([trait_inter_an(1) trait_inter_an(2) time_plot(1) time_plot(end)]);
    set(gca,'TickDir','in','TickLength',[0.015 0.015],'XMinorTick','off','YMinorTick','off','layer','top','LineWidth',1.2,'FontSize',11,'Xtick',2.5:0.5:3.5);
    box ('on');
    ylabel('Evolutionary time','FontSize',12);
    xlabel('Animal trait','FontSize',12);
    set(gca,'xlim',[2 4]);
    set(h3,'Position',[0.13 0.12 0.3 0.83]);
    
    %plot the plants
    h4=subplot(2,2,4);
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
        temp_array_abundance=zeros(length(time_plot),length(sol_plot(end)));
        temp_array_abundance(:)=NaN;
        the_last_index=sep_index(i);
        while count_time_plot<=length(time_plot) && count_time_plot<=the_last_index
            temp_array(row_each_mat,1:current_length)=sol_plot{count_time_plot};
            temp_time_plot(row_each_mat)=time_plot(count_time_plot);
            temp_array_abundance(row_each_mat,1:current_length)=ab_pl{count_time_plot};
            row_each_mat=row_each_mat+1;
            count_time_plot=count_time_plot+1;
        end
        temp_array=temp_array(~isnan(temp_array(:,1)),~isnan(temp_array(1,:)));
        temp_array_abundance=temp_array_abundance(~isnan(temp_array_abundance(:,1)),~isnan(temp_array_abundance(1,:)));
        temp_time_plot=temp_time_plot(~isnan(temp_time_plot));
        for each_line=1:size(temp_array,2)
            ab_normalized=((temp_array_abundance(:,each_line)-the_min_ab)/(the_max_ab-the_min_ab));
            color_line(temp_array(:,each_line),temp_time_plot,ab_normalized,'LineWidth',1.8);
            hold all;
        end
        if count_time_plot<=length(time_plot)
            current_length=length(sol_plot{count_time_plot});
            if any(count_time_plot==introduction_index)
                if length(sol_plot{count_time_plot})>length(sol_plot{count_time_plot-1})
                    inv_time=time_plot(count_time_plot);
                    inv_tr=sol_plot{count_time_plot}(~ismember(sol_plot{count_time_plot},sol_plot{count_time_plot-1}));
                    plot(inv_tr,inv_time,'bo','MarkerSize',7,'MarkerFaceColor','b');
                    hold all;
                end
            end
        end
        i=i+1;
    end

    axis([trait_inter_pl(1) trait_inter_pl(2) time_plot(1) time_plot(end)]);
    set(gca,'TickDir','in','TickLength',[0.015 0.015],'XMinorTick','off','YMinorTick','off','layer','top','LineWidth',1.2,'FontSize',11,'XTick',1.5:0.5:3,'Yticklabel','');
    xlabel('Plant trait','FontSize',12);
    set(gca,'xlim',[1 3]);
    box ('on');
    
    set(h4,'Position',[0.45 0.12 0.3 0.83]);
    
    cb=colorbar('Position',[0.8 0.10 0.025 0.85]);
    right_ticks=100:100: floor(the_max_ab/100)*100;
    new_ticks=(right_ticks-the_min_ab)/(the_max_ab-the_min_ab);
    set(cb,'XTick',new_ticks);
    ticks_label=cell(1,length(new_ticks));
    for i=1:length(ticks_label)
        ticks_label{i}=num2str(right_ticks(i));
    end
    set(cb,'XTickLabel',ticks_label,'FontSize',11,'TickLength',0.015,'LineWidth',1.2);
    
    clabel=ylabel(cb,'Population density','Rotation',-90,'FontSize',12);
    set(clabel,'Position',get(clabel,'Position')+[1.2 0 0]);
    
    

   
end