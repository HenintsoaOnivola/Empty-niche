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
al=1;

log_sig_an=1;
sig_an=1;
sig_pl=sig_an;

all_folders={'./Results_intra_guild/','Results_cross_guild/'};
parameters={[-1.9375,-2.8125],[-2,-2.625]};
plot_names={'FigB3a.png','FigB3b.png'};

for each_pl=1:length(all_folders)

    log_sig_c_an=parameters{each_pl}(1);
    log_sig_m=parameters{each_pl}(2);
    
    sig_c_an=exp(log_sig_c_an);
    sig_c_pl=sig_c_an;
    
    sig_m=exp(log_sig_m);
    disp (plot_names{each_pl});
    the_folder=all_folders{each_pl};

    trait_plot_an=read_data_into_cell(strcat(the_folder,'animal_traits_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
    trait_plot_pl=read_data_into_cell(strcat(the_folder,'plant_traits_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
    pop_plot_an=read_data_into_cell(strcat(the_folder,'animal_abundance_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
    pop_plot_pl=read_data_into_cell(strcat(the_folder,'plant_abundance_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
    time_plot=read_file_to_array(strcat(the_folder,'time_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
    separating_index=read_file_to_array(strcat(the_folder,'sepIndex_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
    
    if each_pl==1
        br_points=[separating_index(1:3),separating_index(7),separating_index(end)];
    else
        br_points=[separating_index(1:3),separating_index(end)];
    end

    the_disr_an_min=zeros(1,length(br_points));
    the_disr_comp_an_min=zeros(1,length(br_points));
    the_disr_mutu_an_min=zeros(1,length(br_points));

    the_disr_an_max=zeros(1,length(br_points));
    the_disr_comp_an_max=zeros(1,length(br_points));
    the_disr_mutu_an_max=zeros(1,length(br_points));

    the_disr_pl_min=zeros(1,length(br_points));
    the_disr_comp_pl_min=zeros(1,length(br_points));
    the_disr_mutu_pl_min=zeros(1,length(br_points));

    the_disr_pl_max=zeros(1,length(br_points));
    the_disr_comp_pl_max=zeros(1,length(br_points));
    the_disr_mutu_pl_max=zeros(1,length(br_points));
  
    for t=1:length(br_points)

        tt=br_points(t);

        info=branching_condition_analytical([pop_plot_an{tt} pop_plot_pl{tt} trait_plot_an{tt} trait_plot_pl{tt}],rA,rP,h,const,sig_c_an,sig_c_pl,sig_an,sig_pl,x0,y0,kx,ky,sig_m,length(pop_plot_an{tt}),length(pop_plot_pl{tt}));

        info_an=info{1};
        disr_an=info_an(:,1);
        disr_an_comp=info_an(:,2);
        disr_an_mutu=info_an(:,3);

        the_tr=trait_plot_an{tt};
        [~, i]=max(the_tr);
        [~, j]=min(the_tr);

        the_disr_an_min(t)=disr_an(j);
        the_disr_comp_an_min(t)=disr_an_comp(j);
        the_disr_mutu_an_min(t)=disr_an_mutu(j);

        the_disr_an_max(t)=disr_an(i);
        the_disr_comp_an_max(t)=disr_an_comp(i);
        the_disr_mutu_an_max(t)=disr_an_mutu(i);

        info_pl=info{2};
        disr_pl=info_pl(:,1);
        disr_pl_comp=info_pl(:,2);
        disr_pl_mutu=info_pl(:,3);

        the_tr=trait_plot_pl{tt};
        [the_max_tr, i]=max(the_tr);
        [the_min_tr, j]=min(the_tr);

        the_disr_pl_min(t)=disr_pl(j);
        the_disr_comp_pl_min(t)=disr_pl_comp(j);
        the_disr_mutu_pl_min(t)=disr_pl_mutu(j);

        the_disr_pl_max(t)=disr_pl(i);
        the_disr_comp_pl_max(t)=disr_pl_comp(i);
        the_disr_mutu_pl_max(t)=disr_pl_mutu(i);

    end

    figure('Visible','off');
    font_axis=10;
    font_label=12;

    subplot(2,2,1);
    bar(transpose([the_disr_an_max;the_disr_comp_an_max;the_disr_mutu_an_max]))
    set(gca,'FontSize', font_axis,'xticklabel',[]);
    set(gca,'Position',[0.1 0.52 0.41 0.41]);
    if each_pl==1
        set(gca,'ylim',[-5 50]);
        ll=legend('Total disruptiveness', 'Contribution of competition', 'Contribution of mutualism','Location','Northeast');
        set(ll,'FontSize',7);
    else
        set(gca,'ylim',[-100 100]);
    end

    subplot(2,2,2);
    bar(transpose([the_disr_pl_max;the_disr_comp_pl_max;the_disr_mutu_pl_max]))
    set(gca,'FontSize',font_axis,'xticklabel',[]);
    set(gca,'Position',[0.56 0.52 0.41 0.41]);
    if each_pl==1
        set(gca,'ylim',[-50 100]);
    else
        set(gca,'ylim',[-75 100]);
    end

    subplot(2,2,3);
    bar(transpose([the_disr_an_min;the_disr_comp_an_min;the_disr_mutu_an_min]))
    xlabel('Evolutionary time','FontSize',font_label);
    set(gca,'FontSize', font_axis,'xticklabel',[]);
    set(gca,'ylim',[-75 100]);
    set(gca,'Position',[0.1 0.08 0.41 0.41]);
    yy=ylabel ('Srength of disruptive selection','FontSize',font_label);
    y_pos=get(yy,'Position');
    set(yy,'Position',[y_pos(1) 100 y_pos(3)]);

    subplot(2,2,4);
    bar(transpose([the_disr_pl_min;the_disr_comp_pl_min;the_disr_mutu_pl_min]));
    xlabel('Evolutionary time','FontSize',font_label);
    set(gca,'FontSize', font_axis,'xticklabel',[]);
    set(gca,'Position',[0.56 0.08 0.41 0.41]);
    if each_pl==1
        set(gca,'ylim',[-5 50]);
    else
        set(gca,'ylim',[-100 100]);
    end
    
    width = 18;     % Width in inches
    height = 9;    % Height in inches
    papersize = get(gcf, 'PaperSize');
    left = (papersize(1)- width)/2;
    bottom = (papersize(2)- height)/2;
    myfiguresize = [left, bottom, width, height];
    set(gcf,'PaperPosition', myfiguresize);
    mkdir('./Supplementary_figures/');
    plot_name=strcat('./Supplementary_figures/',plot_names{each_pl},'.png');
    print(plot_name,'-dpng','-r600');
    clf;
end


