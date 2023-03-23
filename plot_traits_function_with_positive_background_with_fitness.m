function figg=plot_traits_function_with_positive_background_with_fitness(introduction_index,sep_index,time_plot,an_traits,pl_traits,ab_an,ab_pl,rA,rP,const,h,x0,y0,kx,ky,sig_an,sig_pl,sig_c_an,sig_c_pl,sig_m,plot_id)
    figg=figure('Visible','off');
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%plotting the fitness %%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ab_an_last=transpose(ab_an{end});
    ab_pl_last=transpose(ab_pl{end});
    tr_an_last=transpose(an_traits{end});
    tr_pl_last=transpose(pl_traits{end});
    trait_inter_an=[1.95 3.85];
    trait_inter_pl=[1.15 3.05];
    
    h1=subplot(2,2,1);
    plot(tr_an_last,fitness_land_function(tr_an_last,ab_an_last,ab_pl_last,tr_an_last,tr_pl_last,rA,h,const,sig_c_an,sig_an,x0,kx,sig_m),'ko','MarkerSize',7,'MarkerFaceColor','k','MarkerEdgeColor','k');
    hold on;
    fit_inter=fitness_land_function(trait_inter_an(1):0.001:trait_inter_an(2),ab_an_last,ab_pl_last,tr_an_last,tr_pl_last,rA,h,const,sig_c_an,sig_an,x0,kx,sig_m);
    plot(trait_inter_an(1):0.001:trait_inter_an(2),fit_inter,'-','Color',[0.5 0.5 0.5],'LineWidth',1.1);
    hold on;
    axis([trait_inter_an(1) trait_inter_an(2) -0.75 0.1]); 
    set(gca,'TickDir','in','TickLength',[0.015 0.015],'XMinorTick','off','YMinorTick','off','layer','top','LineWidth',1.2,'FontSize',14,'Ytick',-0.6:0.4:0.1,'Xtick',2.5:0.5:3.5,'XtickLabel',[]);
    set(h1,'Position',[0.13 0.83 0.3 0.15]);
    ylabel('Fitness','FontSize',14);
    if (plot_id==2 || plot_id==3)
        set(gca,'ylabel',[]);
    end

    h2=subplot (2,2,2);
    plot(tr_pl_last,fitness_land_function(tr_pl_last,ab_pl_last,ab_an_last,tr_pl_last,tr_an_last,rP,h,const,sig_c_pl,sig_pl,y0,ky,sig_m),'ro','MarkerSize',7,'MarkerFaceColor','k','MarkerEdgeColor','k');
    hold on;
    fit_inter=fitness_land_function(trait_inter_pl(1):0.001:trait_inter_pl(2),ab_pl_last,ab_an_last,tr_pl_last,tr_an_last,rP,h,const,sig_c_pl,sig_pl,y0,ky,sig_m);
    plot(trait_inter_pl(1):0.001:trait_inter_pl(2),fit_inter,'-','Color',[0.5 0.5 0.5],'LineWidth',1.1');
    hold on;
    axis([trait_inter_pl(1) trait_inter_pl(2) -0.75 0.1]); 
    set(gca,'TickDir','in','TickLength',[0.015 0.015],'XMinorTick','off','YMinorTick','off','layer','top','LineWidth',1.2,'FontSize',14,'Ytick',-0.6:0.4:0.1,'YtickLabel',[],'XTick',1.5:0.5:3,'XtickLabel',[]);
    set(h2,'Position',[0.45 0.83 0.3 0.15]);
    if (plot_id==2 || plot_id==3)
        set(gca,'ylabel',[]);
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%plotting the trait evolution %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    addpath('freezeColors');
    the_hot_colormap=hot(110);
    the_hot_colormap=flipud(the_hot_colormap);
    the_hot_colormap=the_hot_colormap(25:105,:);
    the_cool_colormap=[repmat([1 1 1],[41,1]); repmat([0.8 0.8 0.8],[40,1])];
    combined_cmap=[the_cool_colormap;the_hot_colormap];
    colormap(combined_cmap);
    
    %plot the fitness landscape as background
    trait_array_an=trait_inter_an(1):0.001:trait_inter_an(2);
    trait_array_pl=trait_inter_pl(1):0.001:trait_inter_pl(2);
    consid_time=time_plot(1:1:end);
    the_matrix_land_an=zeros(length(consid_time)+1,length(trait_array_an)+1);
    the_matrix_land_pl=zeros(length(consid_time)+1,length(trait_array_pl)+1);
    each_col=1;
    for jj=1:1:length(time_plot)
        ff_an=fitness_land_function(trait_array_an,transpose(ab_an{jj}),transpose(ab_pl{jj}),transpose(an_traits{jj}),transpose(pl_traits{jj}),rA,h,const,sig_c_an,sig_an,x0,kx,sig_m);
        the_matrix_land_an(each_col,1:end-1)=transpose(ff_an);
        each_col=each_col+1;
    end

    the_matrix_land_an(the_matrix_land_an>0)=1;
    the_matrix_land_an(the_matrix_land_an<=0)=0;
    the_matrix_land_an(end,:)=NaN ;
    the_matrix_land_an(:,end)=NaN;
    
    
    each_col=1;
    for jj=1:1:length(time_plot)
        ff_pl=fitness_land_function(trait_array_pl,transpose(ab_pl{jj}),transpose(ab_an{jj}),transpose(pl_traits{jj}),transpose(an_traits{jj}),rP,h,const,sig_c_pl,sig_pl,y0,ky,sig_m);
        the_matrix_land_pl(each_col,1:end-1)=transpose(ff_pl);
        each_col=each_col+1;
    end
    
    the_matrix_land_pl(the_matrix_land_pl>0)=1;
    the_matrix_land_pl(the_matrix_land_pl<=0)=0;
    the_matrix_land_pl(end,:)=NaN ;
    the_matrix_land_pl(:,end)=NaN;
    
    h3=subplot(2,2,3);
    
    pp1=pcolor([trait_array_an trait_array_an(end)],transpose([consid_time; consid_time(end)+(consid_time(end)-consid_time(end-1))]),the_matrix_land_an);
    hold on;
    caxis(gca,[0 2.1]);
    set(pp1,'edgecolor','none');
    
    the_min_ab=min([min(cellfun(@(x) min(x(:)),ab_an))  min(cellfun(@(x) min(x(:)),ab_pl))]);
    the_max_ab=max([max(cellfun(@(x) max(x(:)),ab_an))  max(cellfun(@(x) max(x(:)),ab_pl))]) ;
    
    
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
            ab_normalized=1.1+((temp_array_abundance(:,each_line)-the_min_ab)/(the_max_ab-the_min_ab));
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
    set(gca,'TickDir','in','TickLength',[0.015 0.015],'XMinorTick','off','YMinorTick','off','layer','top','LineWidth',1.2,'FontSize',14,'Xtick',2:0.5:3.5);
    box ('on');
    ylabel('Evolutionary time','FontSize',14);
    xlabel('Animal trait','FontSize',14);
    if (plot_id==2 || plot_id==3)
        set(gca,'ylabel',[]);
    end
    if (plot_id==1 || plot_id==2)
        set(gca,'ytickLabel',0:0.5:3);
        
    end
    if (plot_id==3)
        set(gca,'ytickLabel',0:0.5:2);
    end
    l=get(gca,'Ylim');
    text (1.38,(l(2)-0.1*10^(5)),'\times10^{5}','Interpreter','tex','FontSize',14);
    
    set(h3,'Position',[0.13 0.13 0.3 0.67]);
    
    
    h4=subplot(2,2,4);
    
    pp2=pcolor([trait_array_pl trait_array_pl(end)],transpose([consid_time; consid_time(end)+(consid_time(end)-consid_time(end-1))]),the_matrix_land_pl);
    hold on;
    caxis(gca,[0 2.1]);
    set(pp2,'edgecolor','none');
    
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
            ab_normalized=1.1+((temp_array_abundance(:,each_line)-the_min_ab)/(the_max_ab-the_min_ab));
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
    set(gca,'TickDir','in','TickLength',[0.015 0.015],'XMinorTick','off','YMinorTick','off','layer','top','LineWidth',1.2,'FontSize',14,'XTick',1.5:0.5:3,'Yticklabel','');
    xlabel('Plant trait','FontSize',14);
    box ('on');
    if (plot_id==2 || plot_id==3)
        set(gca,'ylabel',[]);
    end
    set(h4,'Position',[0.45 0.13 0.3 0.67]);
    
    cb=colorbar('Position',[0.85 0.13 0.03 0.85]);
    set(cb,'xlim',[1.1 2.1]);
    right_ticks=100:100: floor(the_max_ab/100)*100;
    new_ticks=1.1+(right_ticks-the_min_ab)/(the_max_ab-the_min_ab);
    set(cb,'YTick',new_ticks);
    ticks_label=cell(1,length(new_ticks));
    for i=1:length(ticks_label)
        ticks_label{i}=num2str(right_ticks(i));
    end
    set(cb,'YTickLabel',ticks_label,'FontSize',14,'TickLength',0.015,'LineWidth',1.2);
    
    clabel=ylabel(cb,'Population density','Rotation',-90,'FontSize',14);
    set(clabel,'Position',get(clabel,'Position')+[1.2 0 0]);
  
end