warning('off', 'MATLAB:MKDIR:DirectoryExists');
%%Define parameters
rA=1;
rP=1;
const=0.1;
h=0.1;
x0=3;
y0=2;
kx=400;
ky=300;

log_sig_c_an_range=-2.5:0.0625:1;
log_sig_m_range=-3:0.0625:-0.5;

log_sig_an=1;
sig_an=1;
sig_pl=sig_an;

all_folders1={'./Results_intra_guild/','./Results_cross_guild/'};
plot_names={'Fig3a','Fig3b'};

%Ploting the first row (Fig3a and Fig3b) (i.e. without invasion)
for ii=1:length(all_folders1)
    the_fold=all_folders1{ii};
    
    % Initializing the matrix containing the number of species
    the_matrix_div=zeros(length(log_sig_c_an_range)+1,length(log_sig_m_range)+1);
    the_matrix_div(end,:)=NaN;
    the_matrix_div(:,end)=NaN;
    % Build the matrix of the number of species
    disp ('Looping through parameters for empty niches')
    for log_sig_c_an=log_sig_c_an_range
        for log_sig_m=log_sig_m_range
            
            disp ([log_sig_an log_sig_c_an log_sig_m]);
    
            log_sig_c_pl=log_sig_c_an;

            sig_c_an=exp(log_sig_c_an);
            sig_c_pl=sig_c_an;

            sig_m=exp(log_sig_m);
            
            fil_name_an_tr=strcat(the_fold,'animal_traits_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat');
            fil_name_pl_tr=strcat(the_fold,'plant_traits_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat');

            last_tr_an=read_last_line_to_array(fil_name_an_tr);
            last_tr_pl=read_last_line_to_array(fil_name_pl_tr);

            the_matrix_div(log_sig_c_an_range==log_sig_c_an,log_sig_m_range==log_sig_m)=length(last_tr_an)+length(last_tr_pl);
            
        end
    end
    
    %plotting the matrix for the number of species (as background color)
    figg=figure('Visible','off');
    box('on');
    all_diversity=2:max(max(the_matrix_div));
    my_colormap=autumn(length(all_diversity)+40);
    my_colormap=flipud(my_colormap);
    first_col=my_colormap(5,:);
    my_colormap=my_colormap(20:length(all_diversity)+39,:);
    my_colormap(1,:)=first_col;
    colormap(my_colormap);

    bg=pcolor([log_sig_m_range -0.4375] ,[log_sig_c_an_range 1.0625],the_matrix_div);
    set(bg, 'EdgeColor', 'none');
    hold on;
    set(gca,'LooseInset',get(gca,'TightInset'))
   
    
    %adding the empty niches to the figure
    disp ('Looping through parameters for empty niches')
    for log_sig_c_an=log_sig_c_an_range
        for log_sig_m=log_sig_m_range
            disp ([log_sig_an log_sig_c_an log_sig_m]);
            log_sig_c_pl=log_sig_c_an;
            sig_c_an=exp(log_sig_c_an);
            sig_c_pl=sig_c_an;

            sig_m=exp(log_sig_m);
            
            re=detect_EN_function(the_fold,log_sig_an,log_sig_c_an,log_sig_m,rA,rP,const,h,x0,y0,kx,ky);
            the_type=re{2};
            if (re{1} ==1)
                if strcmp (the_type,'intra')==1 %case of peripheral empty niches
                    rectangle('Position',[log_sig_m,log_sig_c_an,0.0625 0.0625]);
                    hold on;
                end
                if strcmp (the_type,'cross')==1 % case of central empty niches
                    plot([log_sig_m,log_sig_m+0.0625],[log_sig_c_an+0.0625,log_sig_c_an],'k-');
                    hold on;
                    plot([log_sig_m,log_sig_m+0.0625],[log_sig_c_an,log_sig_c_an+0.0625],'k-');
                end
                if strcmp (the_type,'both')==1  %case where there are both peripheral and central empty niches
                    rectangle('Position',[log_sig_m,log_sig_c_an,0.0625 0.0625]);
                    hold on;
                    plot([log_sig_m,log_sig_m+0.0625],[log_sig_c_an+0.0625,log_sig_c_an],'k-');
                    hold on;
                    plot([log_sig_m,log_sig_m+0.0625],[log_sig_c_an,log_sig_c_an+0.0625],'k-');
                end

            end
            if (re{1} ==2) %case where there are oscillations
                rectangle('Position',[log_sig_m,log_sig_c_an,0.0625 0.0625],'EdgeColor','r');
                hold on;
                plot([log_sig_m,log_sig_m+0.0625],[log_sig_c_an+0.0625,log_sig_c_an],'r-');
                hold on;
                plot([log_sig_m,log_sig_m+0.0625],[log_sig_c_an,log_sig_c_an+0.0625],'r-');
            end
            
        end
            
    end
    %Points needed for the case examples in the figures
    increment=0.0625/2;
    red_m=[-1.9375 -1.7366 -1.5625+0.0625]+increment;
    red_c=-0.375+increment;
    white_m=[-2.5625 -2.25]+increment;
    white_c=-2.3125+increment;
    brown_m=-3+increment+(0.0625);
    brown_c=[-2.0625-(3*0.0625) -1.75 -1.5162 -1.1875 -0.6250]+increment;
    orange_m=[-2.5625 -2.0735 -1.8125 -1.3125]+increment;
    orange_c=-0.1875+increment+(3*0.0625);

    the_m=unique(sort([red_m white_m brown_m orange_m]));
    the_c=unique(sort([red_c white_c brown_c orange_c]));
    
    %Plot the case example points
    for log_sig_c_an=the_c
        for log_sig_m=the_m
            if ii==1 %case of Fig3a
                if (ismember(log_sig_c_an,brown_c) && ismember(log_sig_m,brown_m))
                    disp ('brown');
                    disp ([log_sig_c_an,log_sig_m]);
                    plot(log_sig_m,log_sig_c_an,'o','Color',rgb('Brown'),'MarkerFaceColor',rgb('Brown'),'Markersize',3);
                    hold on;
                end
                if (ismember(log_sig_c_an,orange_c) && ismember(log_sig_m,orange_m))
                    disp ('orange');
                    disp ([log_sig_c_an,log_sig_m]);
                    plot(log_sig_m,log_sig_c_an,'o','Color',rgb('OrangeRed'),'MarkerFaceColor',rgb('OrangeRed'),'Markersize',3);
                    hold on;
                end
            end
            if ii==2 %case of Fig3b
                if (ismember(log_sig_c_an,red_c) && ismember(log_sig_m,red_m))
                    plot(log_sig_m,log_sig_c_an,'o','Color','red','MarkerFaceColor','red','Markersize',3);
                end
                if (ismember(log_sig_c_an,white_c) && ismember(log_sig_m,white_m))
                    plot(log_sig_m,log_sig_c_an,'o','Color','white','MarkerFaceColor','white','Markersize',3);
                end
            end
        end
    end
    %Adjusting the axes
    the_x_ticks=log([0.05,0.1,0.25,0.5:0.5:2]);
    the_x_ticks_label=[0.05,0.1,0.25,0.5:0.5:2];
    set(gca,'TickDir','in','TickLength',[0.015 0.015],'XMinorTick','off','YMinorTick','off','layer','top','LineWidth',1.5);
    set(gca,'Xtick',the_x_ticks,'XtickLabel',the_x_ticks_label,'Ytick',the_x_ticks,'YtickLabel',the_x_ticks_label,'FontSize',15);
    pbaspect([1 1 1]);
    xlab='Tolerance to trait difference, \sigma_{m}';
    ylab='Width of competition kernel, \sigma_{c}';
    xlabel(xlab,'FontSize',15,'Interpreter','tex');
    ylabel(ylab,'FontSize',15,'Interpreter','tex');

    %placing the colorbar
    colormap(my_colormap);
    caxis([min(all_diversity) max(all_diversity)]);
    cc=colorbar();
    cc.Ruler.Scale='log';
    ycb=ylabel(cc,'Number of species','FontSize',15);
    set(ycb,'Rotation',-90);
    set(cc,'YTick',[min(all_diversity) 10:20:max(all_diversity)],'FontSize',15,'TickLength',0.015,'LineWidth',1.5);
    label_pos=cc.Label.Position;
    label_pos(1)=label_pos(1)+1;
    cc.Label.Position=label_pos;
    addpath plotboxpos-pkg/plotboxpos
    fig_pos=plotboxpos(gca);
    cb_pos=get(cc,'Position');
    
    %Adjusting figure position and printing
    set(cc,'Position',[fig_pos(1)+fig_pos(3)+0.05 fig_pos(2)-0.01 cb_pos(3) cb_pos(4)]);
    set(gca,'Position',[fig_pos(1) fig_pos(2)-0.01 fig_pos(3) cb_pos(4)]);
    mkdir('./Plots/');
    name=strcat('./Plots/',plot_names{ii},'.png');
    print(name,'-dpng','-r600');
    clf;
end

%Ploting the second row (Fig3c and Fig3d) (i.e. with invasion)
all_folders2={'./Results_invasion_intra_guild/','./Results_invasion_cross_guild/'};
plot_names={'Fig3c','Fig3d'};
for ii=1:length(all_folders2)
    the_fold=all_folders2{ii}; 
    introduction_matrix=zeros(length(log_sig_c_an_range)+1,length(log_sig_m_range)+1);
    introduction_matrix(end,:)=NaN;
    introduction_matrix(:,end)=NaN;
    for log_sig_c_an=log_sig_c_an_range
        for log_sig_m=log_sig_m_range

            disp ([log_sig_an log_sig_c_an log_sig_m]);

            if log_sig_m<-1.5 
                folder=all_folders2{ii};
            else % case where there is only one evolutionary scenario
                folder=all_folders2{1};
            end
            
            % build the identity (introduced vs native) of each species along the evolutionary
            % time for the specific parameter combination, so as to know
            % the identities at ESS
            an_tr=read_data_into_cell(strcat(folder,'animal_traits_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
            pl_tr=read_data_into_cell(strcat(folder,'plant_traits_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
            sep=read_file_to_array(strcat(folder,'sepIndex_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));

            the_id=define_native_invader_function(an_tr,pl_tr,sep);
            all_an_id=the_id{1};
            all_pl_id=the_id{2};
            
            an_id=all_an_id{end};
            pl_id=all_pl_id{end};
            
            %compute the proportion of invader at the last time step 
            an_invad_number=sum(an_id);
            pl_invad_number=sum(pl_id);

            proportion=(an_invad_number+pl_invad_number)/(length(an_id)+length(pl_id));
            introduction_matrix(log_sig_c_an_range==log_sig_c_an,log_sig_m_range==log_sig_m)=proportion*100;

        end

    end
    
    %Plot the matrix of the proportion of introduced species
    all_proportion=0:0.1:100;
    introduction_matrix(end,end)=100;
    figg=figure('Visible','off');
    box('on');
    my_colormap=flipud(lbmap(80,'RedBlue'));
    colormap(my_colormap);
    h=pcolor([log_sig_m_range -0.4375] ,[log_sig_c_an_range 1.0625],introduction_matrix);
    set(h, 'EdgeColor', 'none');
    pbaspect([1 1 1]);
    set(gca,'LooseInset',get(gca,'TightInset'))
    %Adjust the axes 
    the_x_ticks=log([0.05,0.1,0.25,0.5:0.5:2]);
    the_x_ticks_label=[0.05,0.1,0.25,0.5:0.5:2];
    set(gca,'TickDir','in','TickLength',[0.015 0.015],'XMinorTick','off','YMinorTick','off','layer','top','LineWidth',1.5);
    set(gca,'Xtick',the_x_ticks,'XtickLabel',the_x_ticks_label,'Ytick',the_x_ticks,'YtickLabel',the_x_ticks_label,'FontSize',15);
    xlab='Tolerance to trait difference, \sigma_{m}';
    ylab='Width of competition kernel, \sigma_{c}';
    xlabel(xlab,'FontSize',15,'Interpreter','tex');
    ylabel(ylab,'FontSize',15,'Interpreter','tex');

    %placing the colorbar
    colormap(my_colormap);
    cc=colorbar();
    ycb=ylabel(cc,'Proportion of introduced species (%)','FontSize',15);
    set(ycb,'Rotation',-90);
    set(cc,'YTick',[min(all_proportion) 20:20:max(all_proportion)],'FontSize',15,'TickLength',0.015,'LineWidth',1.5);
    label_pos=cc.Label.Position;
    label_pos(1)=label_pos(1)+1;
    cc.Label.Position=label_pos;
    
    %Adjust plot position and printing
    addpath plotboxpos-pkg/plotboxpos
    fig_pos=plotboxpos(gca);
    cb_pos=get(cc,'Position');
    set(cc,'Position',[fig_pos(1)+fig_pos(3)+0.05 fig_pos(2)-0.01 cb_pos(3) cb_pos(4)]);
    set(gca,'Position',[fig_pos(1) fig_pos(2)-0.01 fig_pos(3) cb_pos(4)]);
    mkdir('./Plots/');
    name=strcat('./Plots/',plot_names{ii},'.png');
    print(name,'-dpng','-r600');
    clf;
end

