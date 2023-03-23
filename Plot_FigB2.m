warning('off', 'MATLAB:MKDIR:DirectoryExists');
%defining the parameters
rA=1;
rP=1;
const=0.1;
h=0.1;
x0=3;
y0=2;
kx=400;
ky=300;


log_sig_an=1;
sig_an=1;
sig_pl=sig_an;

all_folders={'./Results_intra_guild/','Results_intra_guild/','./Results_cross_guild/','./Results_cross_guild/'};
parameters={[-1.3125,-1.3125],[-1.875,-2.6875],[-2.0625,-3],[-1.125,-1.25]};
plot_names={'FigB2a.png','FigB2b.png','FigB2c.png','FigB2d.png'};

for each_pl=1:length(all_folders)

    log_sig_c_an=parameters{each_pl}(1);
    log_sig_m=parameters{each_pl}(2);
    
    sig_c_an=exp(log_sig_c_an);
    sig_c_pl=sig_c_an;
    
    sig_m=exp(log_sig_m);
    disp (plot_names{each_pl});
    the_folder=all_folders{each_pl};

    an_traits=read_data_into_cell(strcat(the_folder,'animal_traits_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
    pl_traits=read_data_into_cell(strcat(the_folder,'plant_traits_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
    an_pop=read_data_into_cell(strcat(the_folder,'animal_abundance_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
    pl_pop=read_data_into_cell(strcat(the_folder,'plant_abundance_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));

    time_plot=read_file_to_array(strcat(the_folder,'time_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
    sep_index=read_file_to_array(strcat(the_folder,'sepIndex_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
    if each_pl==3 %%cut the oscillation time for a better aesthetic of the plot
        time_plot=time_plot(time_plot<=6.5*1e4);
        sep_index=sep_index(sep_index<=length(time_plot));
        sep_index=[sep_index(1:end-1) length(time_plot)];
        an_traits=an_traits(1:length(time_plot));
        pl_traits=pl_traits(1:length(time_plot));
        an_pop=an_pop(1:length(time_plot));
        pl_pop=pl_pop(1:length(time_plot));
    end

    figg=figure('Visible','off');
    subplot('position',[0.1 0.2 0.52 0.7]);
    %plot the animal trait evolution
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
    %plot the plants trait evolution
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
    ylabel('Trait values','FontSize',25);
    xlabel ('Evolutionary time','FontSize',25);
    set(gca,'FontSize',17,'LineWidth',1.1);

    %plotting fitness
    subplot('position',[0.72 0.33 0.25 0.35 ]);
    tr_an_last=an_traits{end};
    ab_an_last=an_pop{end};
    tr_pl_last=pl_traits{end};
    ab_pl_last=pl_pop{end};
    trait_interval=0.9:0.001: 4.1;

    sc=@(x) sign(x).*(log10(1+1000*abs(x)));
    fit_an=fitness_land_function(trait_interval,ab_an_last,ab_pl_last,tr_an_last,tr_pl_last,rA,h,const,sig_c_an,sig_an,x0,kx,sig_m);
    fit_pl=fitness_land_function(trait_interval,ab_pl_last,ab_an_last,tr_pl_last,tr_an_last,rP,h,const,sig_c_pl,sig_pl,y0,ky,sig_m);

    pos_an=sc(fit_an);
    pos_an(pos_an<0)=0;
    fill(trait_interval,pos_an,rgb('LightBlue'));
    hold on; 
    pos_pl=sc(fit_pl);
    pos_pl(pos_pl<0)=0;
    h2=fill(trait_interval,pos_pl,rgb('LightSalmon'));
    set(h2,'facealpha',0.4)
    hold on;

    plot(tr_an_last,sc(fitness_land_function(tr_an_last,ab_an_last,ab_pl_last,tr_an_last,tr_pl_last,rA,h,const,sig_c_an,sig_an,x0,kx,sig_m)),'bo','MarkerSize',2,'MarkerFaceColor','b');
    hold on;
    plot(trait_interval,sc(fit_an),'b-');
    hold on;
    plot(tr_pl_last,sc(fitness_land_function(tr_pl_last,ab_pl_last,ab_an_last,tr_pl_last,tr_an_last,rP,h,const,sig_c_pl,sig_pl,y0,ky,sig_m)),'ro','MarkerSize',2,'MarkerFaceColor','r');
    hold on;
    plot(trait_interval,sc(fit_pl),'r-');
    hold on;
    set(gca,'ylim',sc([-0.6 1.5]));
    set (gca,'ytick',sc(-0.4:0.1:0.8),'yticklabel',{'','','', -0.1,0,0.1, '','', '','','','',''});
    set (gca,'xlim',[trait_interval(1) trait_interval(end)]);
    ylabel('Fitness','FontSize',12);
    box on;
    set(gca,'xlim',[trait_interval(1) trait_interval(end)]);
    xlabel('Trait values','FontSize',12);

    % to control the figure size and resolution
    width = 7;     % Width in inches
    height = 5;    % Height in inches
    alw = 0.75;    % AxesLineWidth
    fsz = 11;      % Fontsize
    lw = 0.5;      % LineWidth
    msz = 8;       % MarkerSize
    %Here we preserve the size of the image when we save it.
    set(gcf,'InvertHardcopy','on');
    set(gcf,'PaperUnits', 'inches');
    papersize = get(gcf, 'PaperSize');
    left = (papersize(1)- width)/2;
    bottom = (papersize(2)- height)/2;
    myfiguresize = [left, bottom, width, height];
    set(gcf,'PaperPosition', myfiguresize);
    %save the figure
    mkdir('./Supplementary_figures/');
    plot_name=strcat('./Supplementary_figures/',plot_names{each_pl},'.png');
    print(plot_name,'-dpng','-r400');
end




