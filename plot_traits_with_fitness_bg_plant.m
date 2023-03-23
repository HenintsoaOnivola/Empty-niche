function figg=plot_traits_with_fitness_bg_plant(trait_inter,sep_index,time_plot,an_traits,pl_traits,ab_an,ab_pl,rP,const,h,y0,ky,sig_pl,sig_c_pl,sig_m,the_pl_tr)
    figg=figure('Visible','off');
    
    addpath('freezeColors');
    
    the_cool_colormap=flipud(lbmap(120,'BrownBlue'));
    the_cool_colormap=the_cool_colormap(63:112,:);
    the_negatives=gray(80);
    the_negatives=the_negatives(31:end,:);
    combined_cmap=[the_negatives;the_cool_colormap];
    colormap(combined_cmap);
    
    %plot the fitness landscape as background
    consid_time=time_plot(1:1:end);
    the_matrix_land_pl=zeros(length(consid_time)+1,length(trait_inter)+1);
    each_col=1;
    for jj=1:1:length(time_plot)
        ff_pl=fitness_land_function(trait_inter,transpose(ab_pl{jj}),transpose(ab_an{jj}),transpose(pl_traits{jj}),transpose(an_traits{jj}),rP,h,const,sig_c_pl,sig_pl,y0,ky,sig_m);
        the_matrix_land_pl(each_col,1:end-1)=transpose(ff_pl);
        each_col=each_col+1;
    end
    the_matrix_land_pl(end,:)=NaN ;
    the_matrix_land_pl(:,end)=NaN;
    
    considered_max=max(max(the_matrix_land_pl));
    considered_max=min([considered_max 1.1]);
    considered_min=min(min(the_matrix_land_pl));
    considered_min=max([considered_min -1.1]);
    the_limit=min(abs([considered_min considered_max]));
    
    [x_big, y_big]=find(the_matrix_land_pl>the_limit);
    for ind=1:length(x_big)
        the_matrix_land_pl(x_big(ind),y_big(ind))=the_limit;
    end
    [x_sm, y_sm]=find(the_matrix_land_pl<-the_limit);
    for ind=1:length(x_sm)
        the_matrix_land_pl(x_sm(ind),y_sm(ind))=-the_limit;
    end
    
    the_min=min(min(the_matrix_land_pl));
    the_max=max(max(the_matrix_land_pl));
    
    the_matrix_land_pl=double((the_matrix_land_pl-the_min)/(the_max-the_min));
 
    pp1=pcolor([trait_inter trait_inter(end)],transpose([consid_time; consid_time(end)+(consid_time(end)-consid_time(end-1))]),the_matrix_land_pl);
    hold on;
    %caxis(gca,[0 2.1]);
    set(pp1,'edgecolor','none');
    set(gca,'layer','top');
    
    %plot the animal traits
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
        plot(temp_array,temp_time_plot,'b-','LineWidth',1.5);
        hold all;
        if count_time_plot<=length(time_plot)
            current_length=length(sol_plot{count_time_plot});
        end
        i=i+1;
        
    end
   
    %plot the plant traits
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
        plot(temp_array,temp_time_plot,'r-','LineWidth',1.5);
        hold all;
        
        if count_time_plot<=length(time_plot)
            current_length=length(sol_plot{count_time_plot});
        end
        i=i+1;
    end
    plot(the_pl_tr,0,'kx','Markersize',12);
    
    axis([trait_inter(1) trait_inter(end) time_plot(1) time_plot(end)]);
    ylabel('Evolutionary time','FontSize',15);
    xlabel('Traits','FontSize',15);
    box on;
    
    cb1=colorbar('location','North');
    the_pos=get(cb1,'Position');
    set(cb1,'Position',[the_pos(1) the_pos(2)+0.1 the_pos(3) the_pos(4)-0.02]);
    right_ticks1=-1:0.5:1;
    new_ticks1=(right_ticks1-the_min)/(the_max-the_min);
    set(cb1,'YTick',new_ticks1);
    ticks_label1=cell(1,length(new_ticks1));
    for i=1:length(ticks_label1)
       ticks_label1{i}=num2str(right_ticks1(i));
    end
    set(cb1,'YTickLabel',ticks_label1);
    ylabel(cb1,'Fitness','FontSize',12);
    
    main_pos=get(gca,'Position');
    set(gca,'Position',[main_pos(1) main_pos(2)-0.02 main_pos(3) main_pos(4)-0.05]);
    
    %Here we preserve the size of the image when we save it.
    width = 5;     % Width in inches
    height = 5.5;    % Height in inches
    set(gcf,'InvertHardcopy','on');
    set(gcf,'PaperUnits', 'inches');
    papersize = get(gcf, 'PaperSize');
    left = (papersize(1)- width)/2;
    bottom = ((papersize(2)- height)/2);
    myfiguresize = [left, bottom, width, height];
    set(gcf,'PaperPosition', myfiguresize);
    
end