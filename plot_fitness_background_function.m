function plot_fitness_background_function(sep_index,time_plot,an_traits,pl_traits,ab_an,ab_pl,Pareto_an,Pareto_pl,r,h,const,sig_c_an,sig_c_pl,sig_a,sig_p,x0,y0,kx,ky,sig_m,trait_limits,x_tick_lab)
    addpath('freezeColors');
    figure('Visible','off');
    the_hot_colormap=summer(110);
    the_hot_colormap=flipud(the_hot_colormap);
    the_hot_colormap=the_hot_colormap(25:105,:);
    the_cool_colormap=flipud(lbmap(100,'BrownBlue'));
    the_cool_colormap=[the_cool_colormap(11:49,:); the_cool_colormap(55:96,:)];
    combined_cmap=[the_cool_colormap;the_hot_colormap];
    
    colormap(combined_cmap);
    set(gcf,'Renderer','zb');  %this is needed to remove the power 10 of the xlabels
 
    %merge the traits
    for tt=1:length(time_plot)
        an=an_traits{tt};
        pl=pl_traits{tt};
        an_pop=ab_an{tt};
        pl_pop=ab_pl{tt};
        
        an_merged=merge_in_plot(an,an_pop);
        an_tr_merged=an_merged{1};
        an_pop_merged=an_merged{2};

        pl_merged=merge_in_plot(pl,pl_pop);
        pl_tr_merged=pl_merged{1};
        pl_pop_merged=pl_merged{2};
        
        an_traits{tt}=an_tr_merged;
        pl_traits{tt}=pl_tr_merged;
        ab_an{tt}=an_pop_merged;
        ab_pl{tt}=pl_pop_merged;
    end
    add_sep_index=[];
    an_before=length(an_traits{1});
    pl_before=length(pl_traits{1});
    for tt=2:length(time_plot)
        if length(an_traits{tt})~=an_before || length(pl_traits{tt})~= pl_before
            add_sep_index=[add_sep_index tt-1];
            an_before=length(an_traits{tt});
            pl_before=length(pl_traits{tt});
        end
    end
    sep_index=[sep_index add_sep_index];
    sep_index=sort(sep_index);
    sep_index=unique(sep_index);
    
    %%plot the left side (animal side)
    ax=gca;
    pos1=get(ax,'Position');
    pos1(1)=0.08;
    pos1(2)=0.49;
    pos1(3)=0.418;
    pos1(4)=0.35;
    set(ax,'Position',pos1,...
               'XAxisLocation','bottom',...
               'YAxisLocation','right',...
               'Color','none',...
               'XColor','black',...
               'YColor','k',...
               'XLim',[time_plot(1) time_plot(end)]);
       
    %plot the fitness landscape as background
    min_trait=trait_limits(1);
    max_trait=trait_limits(2);
    trait_array=min_trait:0.001:max_trait;
    consid_time=time_plot(1:1:end);
    the_matrix_land_an=zeros(length(trait_array)+1,length(consid_time)+1);
    the_matrix_land_pl=zeros(length(trait_array)+1,length(consid_time)+1);
    each_col=1;
    for jj=1:1:length(time_plot)
        ff_an=fitness_land_function(trait_array,transpose(ab_an{jj}),transpose(ab_pl{jj}),transpose(an_traits{jj}),transpose(pl_traits{jj}),r,h,const,sig_c_an,sig_a,x0,kx,sig_m);
        the_matrix_land_an(1:end-1,each_col)=transpose(ff_an);
        each_col=each_col+1;
    end
    the_matrix_land_an(end,:)=NaN ;
    the_matrix_land_an(:,end)=NaN;
    
    each_col=1;
    for jj=1:1:length(time_plot)
        ff_pl=fitness_land_function(trait_array,transpose(ab_pl{jj}),transpose(ab_an{jj}),transpose(pl_traits{jj}),transpose(an_traits{jj}),r,h,const,sig_c_pl,sig_p,y0,ky,sig_m);
        the_matrix_land_pl(1:end-1,each_col)=transpose(ff_pl);
        each_col=each_col+1;
    end
    the_matrix_land_pl(end,:)=NaN ;
    the_matrix_land_pl(:,end)=NaN;
    
    all_max=[max(max(the_matrix_land_an)) max(max(the_matrix_land_pl))];
    considered_max=max(all_max);
    considered_max=min([considered_max 1.1]);
    all_min=[min(min(the_matrix_land_an)) min(min(the_matrix_land_pl))];
    considered_min=min(all_min);
    considered_min=max([considered_min -1.1]);
    the_limit=min(abs([considered_min considered_max]));
  
    [x_big, y_big]=find(the_matrix_land_an>the_limit);
    for ind=1:length(x_big)
        the_matrix_land_an(x_big(ind),y_big(ind))=the_limit;
    end
    [x_sm, y_sm]=find(the_matrix_land_an<-the_limit);
    for ind=1:length(x_sm)
        the_matrix_land_an(x_sm(ind),y_sm(ind))=-the_limit;
    end
 
    [x_big, y_big]=find(the_matrix_land_pl>the_limit);
    for ind=1:length(x_big)
        the_matrix_land_pl(x_big(ind),y_big(ind))=the_limit;
    end
    [x_sm, y_sm]=find(the_matrix_land_pl<-the_limit);
    for ind=1:length(x_sm)
        the_matrix_land_pl(x_sm(ind),y_sm(ind))=-the_limit;
    end
    the_min=min([min(min(the_matrix_land_an)),min(min(the_matrix_land_pl))]);
    the_max=max([max(max(the_matrix_land_an)),max(max(the_matrix_land_pl))]);
    
    the_matrix_land_an=double((the_matrix_land_an-the_min)/(the_max-the_min));
    the_matrix_land_pl=double((the_matrix_land_pl-the_min)/ (the_max-the_min));
    
    pp1=pcolor(transpose([consid_time; consid_time(end)+(consid_time(end)-consid_time(end-1))]),[trait_array trait_array(end)],the_matrix_land_an,'Parent',ax);
    caxis(ax,[0 2.1]);
    set(pp1,'edgecolor','none');
    
    %plot the trait values and population abundance
    the_min_ab=min([min(cellfun(@(x) min(x(:)),ab_an))  min(cellfun(@(x) min(x(:)),ab_pl))]);
    the_max_ab=max([max(cellfun(@(x) max(x(:)),ab_an))  max(cellfun(@(x) max(x(:)),ab_pl))]) ;
    
    sol_plot=an_traits;
    count_time_plot=1;
    current_length=length(sol_plot{1});
    ind_sep=1;
    
    previous_last_time_plot=time_plot(1);
    
    while (count_time_plot<=length(time_plot))
        row_each_mat=1;
        temp_array=zeros(length(time_plot),length(sol_plot(end)));
        temp_array(:)=NaN;
        temp_array_abundance=zeros(length(time_plot),length(sol_plot(end)));
        temp_array_abundance(:)=NaN;
        temp_time_plot=zeros(size(time_plot));
        temp_time_plot(:)=NaN;
        the_last_index=sep_index(ind_sep);
        while count_time_plot<=length(time_plot) && count_time_plot<=the_last_index
            temp_array(row_each_mat,1:current_length)=sol_plot{count_time_plot};
            temp_array_abundance(row_each_mat,1:current_length)=ab_an{count_time_plot};
            temp_time_plot(row_each_mat)=time_plot(count_time_plot);
            row_each_mat=row_each_mat+1;
            count_time_plot=count_time_plot+1;
        end
        
        temp_array=temp_array(~isnan(temp_array(:,1)),~isnan(temp_array(1,:)));
        temp_array_abundance=temp_array_abundance(~isnan(temp_array_abundance(:,1)),~isnan(temp_array_abundance(1,:)));
        temp_time_plot=temp_time_plot(~isnan(temp_time_plot));
            
        for each_line=1:size(temp_array,2)
            ab_normalized=1.1+((temp_array_abundance(:,each_line)-the_min_ab)/(the_max_ab-the_min_ab));
            color_line([previous_last_time_plot; temp_time_plot],[temp_array(1,each_line); temp_array(:,each_line)],[ab_normalized(1); ab_normalized],'LineWidth',1.6,'Parent', ax);
            hold all;
        end
        previous_last_time_plot=temp_time_plot(end);
        
        if count_time_plot<=length(time_plot)
            current_length=length(sol_plot{count_time_plot});
        end
        ind_sep=ind_sep+1;
    end
    set(ax,'FontSize',12,'XAxisLocation','bottom','YAxisLocation','left','layer','top','Xlim',[time_plot(1) time_plot(end)],'Ylim',[min_trait max_trait],'LineWidth',1.2);
    the_xtick=get(ax,'xtick');
    set(ax,'Xticklabel','');
    ylabel('Animal trait','FontSize',12);
    ax.YLabel.Units = 'normalized';
    ax.YLabel.Position=[-0.11 0.5];
    
  
    
    %Right side (plant)
    %define left axis
    pos2=get(ax,'Position');
    pos2(1)=0.5002;
    ax2 = axes('Position',pos2,...
               'XAxisLocation','bottom',...
               'YAxisLocation','right',...
               'Color','none',...
               'XColor','black',...
               'YColor','k',...
               'XLim',[time_plot(1) time_plot(end)]);
    %plot the fitness landscape as background
    pp2=pcolor(transpose([consid_time; consid_time(end)+(consid_time(end)-consid_time(end-1))]),[trait_array trait_array(end)],the_matrix_land_pl,'Parent',ax2);
    caxis(ax2,[0 2.1]);
    set(pp2,'edgecolor','none');
    
    
    %plot the trait values and population abundance
    sol_plot=pl_traits;
    count_time_plot=1;
    current_length=length(sol_plot{1});
    ind_sep=1;
    previous_last_time_plot=time_plot(1);
    while (count_time_plot<=length(time_plot))
        row_each_mat=1;
        temp_array=zeros(length(time_plot),length(sol_plot(end)));
        temp_array(:)=NaN;
        temp_array_abundance=zeros(length(time_plot),length(sol_plot(end)));
        temp_array_abundance(:)=NaN;
        temp_time_plot=zeros(size(time_plot));
        temp_time_plot(:)=NaN;
        the_last_index=sep_index(ind_sep);
        while count_time_plot<=length(time_plot) && count_time_plot<=the_last_index
            temp_array(row_each_mat,1:current_length)=sol_plot{count_time_plot};
            temp_array_abundance(row_each_mat,1:current_length)=ab_pl{count_time_plot};
            temp_time_plot(row_each_mat)=time_plot(count_time_plot);
            row_each_mat=row_each_mat+1;
            count_time_plot=count_time_plot+1;
        end
        temp_array=temp_array(~isnan(temp_array(:,1)),~isnan(temp_array(1,:)));
        temp_array_abundance=temp_array_abundance(~isnan(temp_array_abundance(:,1)),~isnan(temp_array_abundance(1,:)));
        temp_time_plot=temp_time_plot(~isnan(temp_time_plot));
        for each_line=1:size(temp_array,2)
            ab_normalized=1.1+((temp_array_abundance(:,each_line)-the_min_ab)/(the_max_ab-the_min_ab));
            color_line([previous_last_time_plot; temp_time_plot],[temp_array(1,each_line); temp_array(:,each_line)],[ab_normalized(1); ab_normalized],'LineWidth',1.6,'Parent', ax2);
            hold all;
        end
        
        previous_last_time_plot=temp_time_plot(end);
        
        if count_time_plot<=length(time_plot)
            current_length=length(sol_plot{count_time_plot});
        end
        ind_sep=ind_sep+1;
    end
    set(ax2,'Xdir','reverse','FontSize',12,'XAxisLocation','bottom','YAxisLocation','right','layer','top','Xlim',[time_plot(1) time_plot(end)],'Ylim',[min_trait max_trait],'LineWidth',1.2);
    set(ax2,'xtick',the_xtick);
    set(ax2,'Xticklabel','');
    ylabel('Plant trait','FontSize',12,'Rotation',-90);
    ax2.YLabel.Units = 'normalized';
    ax2.YLabel.Position=[1.18 0.5];
    
    

   %plot the colorbars
    cb1=colorbar('location','northoutside');
    set(cb1,'xlim',[0 1]);
    xxopos=[0.3 0.88 0.398 0.027];
    set(cb1,'Position',xxopos);
    ylabel(cb1,'Invasion fitness');
    right_ticks1=-1:0.5:1;
    new_ticks1=(right_ticks1-the_min)/(the_max-the_min);
    set(cb1,'XTick',new_ticks1);
    ticks_label1=cell(1,length(new_ticks1));
    for i=1:length(ticks_label1)
       ticks_label1{i}=num2str(right_ticks1(i));
    end
    set(cb1,'XTickLabel',ticks_label1);
    set(cb1,'FontSize',12,'LineWidth',1.2);
    
    
    cb2=colorbar('location','southoutside');
    set(cb2,'xlim',[1.1 2.1]);
    xxopos=get(cb2,'Position');
    xxopos=[0.3 0.1 xxopos(3) 0.027];
    set(cb2,'Position',xxopos);
    ylabel(cb2,'Population density');
    right_ticks2=100:100: floor(the_max_ab/100)*100;
    new_ticks2=1.1+(right_ticks2-the_min_ab)/(the_max_ab-the_min_ab);
    set(cb2,'XTick',new_ticks2);
    ticks_label2=cell(1,length(new_ticks2));
    for i=1:length(ticks_label2)
        ticks_label2{i}=num2str(right_ticks2(i));
    end
    set(cb2,'XTickLabel',ticks_label2);
    set(cb2,'FontSize',12,'LineWidth',1.2);

   %Plot the total abundance relative to the Pareto abundance
   %Left side (animal)
    pos1_b=get(ax,'Position');
    pos1_b(2)=0.26;
    pos1_b(4)=0.21;
    
    ax_b = axes('Position',pos1_b,...
               'XAxisLocation','bottom',...
               'YAxisLocation','right',...
               'Color','none',...
               'XColor','black',...
               'YColor','k',...
               'XLim',[time_plot(1) time_plot(end)]);

    plot(time_plot,cellfun(@sum,ab_an)/Pareto_an,'r-','LineWidth',1.8,'Parent',ax_b);
    
    set(ax_b,'FontSize',12,'XAxisLocation','bottom','YAxisLocation','left','layer','top','Xlim',[time_plot(1) time_plot(end)],'Ylim',[0 1.05],'LineWidth',1.2);
    set(ax_b,'xtick',the_xtick);
    set(ax_b,'Xticklabel',x_tick_lab);
    ylabel('Total abundance','FontSize',12);
    ax_b.YLabel.Units = 'normalized';
    ax_b.YLabel.Position=[-0.11 0.5];
    xlabel('Evolutionary time (10^{5})','interpreter','tex','FontSize',12);
  
    %Right side (plant)
    pos2_b=get(ax_b,'Position');
    pos2_b(1)=0.5002;
    
    ax2_b = axes('Position',pos2_b,...
               'XAxisLocation','bottom',...
               'YAxisLocation','right',...
               'Color','none',...
               'XColor','black',...
               'YColor','k',...
               'XLim',[time_plot(1) time_plot(end)]);
           
    plot(time_plot,cellfun(@sum,ab_pl)/Pareto_pl,'r-','LineWidth',1.8,'Parent',ax2_b);
  
    set(ax2_b,'Xdir','reverse','FontSize',12,'XAxisLocation','bottom','YAxisLocation','right','layer','top','Xlim',[time_plot(1) time_plot(end)],'Ylim',[0 1.05],'LineWidth',1.2);
    set(ax2_b,'xtick',the_xtick);
    set(ax2_b,'Xticklabel',x_tick_lab);
    ylabel('Total abundance','FontSize',12,'Rotation',-90);
    ax2_b.YLabel.Units = 'normalized';
    ax2_b.YLabel.Position=[1.18 0.5];
    xlabel('Evolutionary time (10^{5})','Interpreter','tex','FontSize',12);
 
end