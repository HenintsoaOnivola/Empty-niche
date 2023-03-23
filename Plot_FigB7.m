addpath('./addaxis6/');
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

%Plot the last figure (fitness component over tolerance to trait difference
log_sig_m_range_usual=-3:0.0625:-0.5;
small_range=[-1.5625:(0.0625/10):-1.5 -2.0735];
log_sig_m_range=unique(sort([log_sig_m_range_usual small_range]));

the_pl_tr=2.6;

log_sig_c_an=-0.1875;
sig_c_an=exp(log_sig_c_an);
sig_c_pl=sig_c_an;

log_sig_an=1;
sig_an=1;
sig_pl=sig_an;

the_folder_output='./Supplementary_figures/';
mkdir(the_folder_output);

% to control the figure size and resolution
width = 5;     % Width in inches
height = 5;    % Height in inches

intra_res=zeros(1,length(log_sig_m_range));
mutu_res=zeros(1,length(log_sig_m_range));
comp=zeros(1,length(log_sig_m_range));
fit=zeros(1,length(log_sig_m_range));
sp_number=zeros(1,length(log_sig_m_range));
sp_abundance=zeros(1,length(log_sig_m_range));
all_abundance=cell(length(log_sig_m_range),1);
for log_sig_m=log_sig_m_range
    sig_m=exp(log_sig_m);
    if (ismember(log_sig_m,log_sig_m_range_usual)) 
        the_folder='./Results_intra_guild/';
        last_tr_an=read_last_line_to_array(strcat(the_folder,'animal_traits_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an,6),'_logsigM',num2str(log_sig_m,6),'.dat'));
        last_tr_pl=read_last_line_to_array(strcat(the_folder,'plant_traits_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an,6),'_logsigM',num2str(log_sig_m,6),'.dat'));
        last_pop_an=read_last_line_to_array(strcat(the_folder,'animal_abundance_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an,6),'_logsigM',num2str(log_sig_m,6),'.dat'));
        last_pop_pl=read_last_line_to_array(strcat(the_folder,'plant_abundance_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an,6),'_logsigM',num2str(log_sig_m,6),'.dat'));    
    else
        res=more_simulation_for_supp_mat(sig_c_an,sig_m,'res');
        trait_plot_an=res{2};
        trait_plot_pl=res{3};
        pop_plot_an=res{4};
        pop_plot_pl=res{5};
        last_tr_an=trait_plot_an{end};
        last_tr_pl=trait_plot_pl{end};
        last_pop_an=pop_plot_an{end};
        last_pop_pl=pop_plot_pl{end};
    end
    
    intra_res(log_sig_m_range==log_sig_m)=k_x(ky,y0,sig_pl,the_pl_tr);
    mutu_res(log_sig_m_range==log_sig_m)=mutu(the_pl_tr,last_tr_an,last_pop_an,const,sig_m,h);
    comp(log_sig_m_range==log_sig_m)=1-(competition_term_mutant(the_pl_tr,last_tr_pl,last_pop_pl,sig_c_pl)/intra_res(log_sig_m_range==log_sig_m));
    fit(log_sig_m_range==log_sig_m)=fitness_land_function(the_pl_tr,last_pop_pl,last_pop_an,last_tr_pl,last_tr_an,rP,h,const,sig_c_pl,sig_pl,y0,ky,sig_m);
    sp_number(log_sig_m_range==log_sig_m)=length(last_tr_pl);
    sp_abundance(log_sig_m_range==log_sig_m)=mean(last_pop_pl);
    all_abundance{log_sig_m_range==log_sig_m}=[last_pop_an last_pop_pl];
end
figure('visible','off');
plot(log_sig_m_range(log_sig_m_range<=-1.8125),comp(log_sig_m_range<=-1.8125),'-','Color',rgb('magenta'),'LineWidth',1.5);
set(gca,'YColor',rgb('magenta'));
hold on;
plot(log_sig_m_range(log_sig_m_range>=-1.8125 & fit >0),comp(log_sig_m_range>=-1.8125 & fit >0),'-','Color',rgb('magenta'),'LineWidth',1.5);
plot(log_sig_m_range(log_sig_m_range>=-1.8125 & fit<=0),comp(log_sig_m_range>=-1.8125 & fit<=0),'-','Color',rgb('magenta'),'LineWidth',1.5);

the_max=max([fit mutu_res comp])+0.1;
the_min=min([fit mutu_res comp])-0.1;
set(gca,'Ylim',[the_min the_max]);
plot([-2.5625 -2.0735 -1.8125 -1.3125],[the_min the_min the_min the_min],'kx','MarkerSize',12,'LineWidth',1.5);
hold on;
yline(0,'--');
addaxis(log_sig_m_range(log_sig_m_range<=-1.8125),mutu_res(log_sig_m_range<=-1.8125),[the_min the_max],'-','Color',rgb('MediumSeaGreen'),'LineWidth',1.5);
addaxisplot(log_sig_m_range(log_sig_m_range>=-1.8125 & fit >0),mutu_res(log_sig_m_range>=-1.8125 & fit >0),2,'-','Color',rgb('MediumSeaGreen'),'LineWidth',1.5);
addaxisplot(log_sig_m_range(log_sig_m_range>=-1.8125 & fit<=0),mutu_res(log_sig_m_range>=-1.8125 & fit <=0),2,'-','Color',rgb('MediumSeaGreen'),'LineWidth',1.5);
addaxis(log_sig_m_range(log_sig_m_range<=-1.8125),fit(log_sig_m_range<=-1.8125),[the_min the_max],'-','Color',[0 0 0 0.75],'LineWidth',1.5);
addaxisplot(log_sig_m_range(log_sig_m_range>=-1.8125 & fit >0),fit(log_sig_m_range>=-1.8125 & fit >0),3,'-','Color',[0 0 0 0.75],'LineWidth',1.5);
addaxisplot(log_sig_m_range(log_sig_m_range>=-1.8125 & fit<=0),fit(log_sig_m_range>=-1.8125 & fit<=0),3,'-','Color',[0 0 0 0.75],'LineWidth',1.5);
addaxis(log_sig_m_range,sp_number,'o','Color',rgb('orange'),'MarkerFaceColor',rgb('orange'));
addaxis(log_sig_m_range,sp_abundance,[250 600],'o','Color',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6],'LineWidth',1.5);
the_x_ticks=log([0.05,0.1,0.25,0.5:0.5:2]);
the_x_ticks_label=[0.05,0.1,0.25,0.5:0.5:2];
the_pos=get(gca,'Position');
set(gca,'Xtick',the_x_ticks,'XtickLabel',the_x_ticks_label,'FontSize',15,'Position',[the_pos(1) the_pos(2)+0.05 the_pos(3:4)]);
xlab='Tolerance to trait difference, \sigma_{m}';
xlabel(xlab,'FontSize',15,'Interpreter','tex');

ylabels = {'Intra-guild resource benefit','Mutualistic benefit','Fitness','Plant species number','Average population density'};
gx = gca;
hax_all = getaddaxisdata(gx,'axisdata');
hax = [hax_all{1}(1:2) hax_all{2}(1:2) hax_all{3:end}];
hax = [gca hax(1,:)];
first_pos=get(gca,'Position');
for k = 1:numel(hax)
    hax(k).FontSize = 15;
    hax(k).YLabel.String = ylabels{k};
    if rem(k,2)==0
        hax(k).YLabel.Rotation = -90;
        hax(k).YLabel.VerticalAlignment = 'bottom';
    end
    each_pos=get(hax(k),'Position');
    hax(k).Position=[each_pos(1) first_pos(2) each_pos(3:4)];
end

plot_name=strcat(the_folder_output,'/FigB7m.png');

set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width*4, height];
set(gcf,'PaperPosition', myfiguresize);
print(plot_name,'-dpng','-r300');
clf;

%Plot fitness, and traits evo and fitness components for particular cases
all_log_sig_m_range=log_sig_m_range;
log_sig_m_range=[-2.5625 -2.0735 -1.8125 -1.3125];
sc=@(x) sign(x).*(log10(1+1000*abs(x)));
considered_ab=all_abundance(ismember(all_log_sig_m_range,log_sig_m_range));
considered_ab=[considered_ab{:}];
min_ab=min(considered_ab);
max_ab=max(considered_ab);
for log_sig_m=log_sig_m_range
    sig_m=exp(log_sig_m);
    if (ismember(log_sig_m,log_sig_m_range_usual)) 
        the_folder='./Results_intra_guild/';
        trait_plot_an=read_data_into_cell(strcat(the_folder,'animal_traits_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an,6),'_logsigM',num2str(log_sig_m,6),'.dat'));
        trait_plot_pl=read_data_into_cell(strcat(the_folder,'plant_traits_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an,6),'_logsigM',num2str(log_sig_m,6),'.dat'));
        pop_plot_an=read_data_into_cell(strcat(the_folder,'animal_abundance_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an,6),'_logsigM',num2str(log_sig_m,6),'.dat'));
        pop_plot_pl=read_data_into_cell(strcat(the_folder,'plant_abundance_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an,6),'_logsigM',num2str(log_sig_m,6),'.dat'));
        time_plot=read_file_to_array(strcat(the_folder,'time_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an,6),'_logsigM',num2str(log_sig_m,6),'.dat'));
        sepIndex=read_file_to_array(strcat(the_folder,'sepIndex_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an,6),'_logsigM',num2str(log_sig_m,6),'.dat'));
    else
        res=more_simulation_for_supp_mat(sig_c_an,sig_m,'res');
        time_plot=res{1};
        trait_plot_an=res{2};
        trait_plot_pl=res{3};
        pop_plot_an=res{4};
        pop_plot_pl=res{5};
        sepIndex=res{6};
    end

    tr_an_last=trait_plot_an{end};
    ab_an_last=pop_plot_an{end};
    tr_pl_last=trait_plot_pl{end};
    ab_pl_last=pop_plot_pl{end};
   
    trait_interval=1:0.001: 4;

    %%plot fitness functions
    figure('Visible','off');
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
    for an_ind=1:length(tr_an_last)
        if min_ab==max_ab
            ms=7;
        else
            ms=4+((ab_an_last(an_ind)-min_ab)*(10-4)/(max_ab-min_ab));
        end
        plot(tr_an_last(an_ind),sc(fitness_land_function(tr_an_last(an_ind),ab_an_last,ab_pl_last,tr_an_last,tr_pl_last,rA,h,const,sig_c_an,sig_an,x0,kx,sig_m)),'bo','MarkerSize',ms,'MarkerFaceColor','b');
    end
    p1=plot(trait_interval,sc(fit_an),'b-','LineWidth',1.5);
    for pl_ind=1:length(tr_pl_last)
        if min_ab==max_ab
            ms=7;
        else
            ms=4+((ab_pl_last(pl_ind)-min_ab)*(10-4)/(max_ab-min_ab));
        end
        plot(tr_pl_last(pl_ind),sc(fitness_land_function(tr_pl_last(pl_ind),ab_pl_last,ab_an_last,tr_pl_last,tr_an_last,rP,h,const,sig_c_pl,sig_pl,y0,ky,sig_m)),'ro','MarkerSize',ms,'MarkerFaceColor','r');
    end
    p2=plot(trait_interval,sc(fit_pl),'r-','LineWidth',1.5);
    plot(the_pl_tr,sc(-0.6),'kx','Markersize',12);
    xline(the_pl_tr,'--');
    set(gca,'ylim',sc([-0.6 1.5]));
    set (gca,'ytick',sc(-0.4:0.1:0.8),'yticklabel',{'','','', -0.1,0,0.1, '','', '','','','',''});
    set (gca,'xlim',[trait_interval(1) trait_interval(end)]);
    ylabel('Fitness','FontSize',15);
    set(gca,'xlim',[trait_interval(1) trait_interval(end)]);
    axis_pos=get(gca,'Position');
    xlabel('Traits','FontSize',15);
    if log_sig_m==min(log_sig_m_range)
        [~, hobj, ~, ~] =legend([p1,p2],'Animals','Plants','FontSize',12,'Location','NorthWest');
        hl = findobj(hobj,'type','line');
        set(hl,'LineWidth',1.2);
    end
    plot_name=strcat(the_folder_output,'/FigB7_fitness_sigM=',num2str(log_sig_m),'.png');
    %Here we preserve the size of the image when we save it.
    set(gcf,'InvertHardcopy','on');
    set(gcf,'PaperUnits', 'inches');
    papersize = get(gcf, 'PaperSize');
    left = (papersize(1)- width)/2;
    bottom = (papersize(2)- height)/2;
    myfiguresize = [left, bottom, width, height];
    set(gcf,'PaperPosition', myfiguresize);
    print(plot_name,'-dpng','-r300');
    clf;
    
    %%plot components of the fitness
    figure('Visible','off');
    carry_cap=k_x(ky,y0,sig_pl,trait_interval);
    comp_component=1-(competition_term_mutant(trait_interval,tr_pl_last,ab_pl_last,sig_c_pl)./carry_cap);
    comp_component(carry_cap==0)=1-1000;
    mutu_component=mutu(trait_interval,tr_an_last,ab_an_last,const,sig_m,h);
    the_min=min([comp_component mutu_component]);
    the_max=max([comp_component mutu_component]);
    plot(trait_interval,sc(comp_component),'-','Color',rgb('Magenta'),'linewidth',1.5);
    hold on;
    plot(the_pl_tr,sc(-0.6),'kx','Markersize',12);
    xline(the_pl_tr,'--');
    ylabel('Plant intra-guild resource benefit','FontSize',15);
    xlabel('Traits','FontSize',15);
    set (gca,'xlim',[trait_interval(1) trait_interval(end)]);
    set(gca,'ylim',sc([-0.6 1.5]));
    set (gca,'ytick',sc(-0.4:0.1:0.8),'yticklabel',{'','','', -0.1,0,0.1, '','', '','','','',''});
    
    yyaxis right;
    plot(trait_interval,sc(mutu_component),'-','Color',rgb('MediumSeaGreen'),'linewidth',1.5);
    yy=ylabel('Plant mutualistic benefit ','FontSize',15);
    prev_pos=get(yy,'Position');
    set(yy,'Position',[4.4 0 -1],'Rotation',-90);
    
    set(gca,'ylim',sc([-0.6 1.5]));
    set (gca,'ytick',sc(-0.4:0.1:0.8),'yticklabel',{'','','', -0.1,0,0.1, '','', '','','','',''});
    set (gca,'xlim',[trait_interval(1) trait_interval(end)]);
    set(gca,'Position',axis_pos);
    ax = gca;
    ax.YAxis(2).Color = rgb('MediumSeaGreen');
    ax.YAxis(1).Color=rgb('Magenta');
    
    plot_name=strcat(the_folder_output,'/FigB7_fitness_component_sigM=',num2str(log_sig_m),'.png');
    %Here we preserve the size of the image when we save it.
    set(gcf,'InvertHardcopy','on');
    set(gcf,'PaperUnits', 'inches');
    papersize = get(gcf, 'PaperSize');
    left = (papersize(1)- width)/2;
    bottom = (papersize(2)- height)/2;
    myfiguresize = [left, bottom, width, height];
    set(gcf,'PaperPosition', myfiguresize);
    print(plot_name,'-dpng','-r300');
    clf;
    
    %%Plot trait evolution with fitness background
    trait_plot_an = cellfun(@transpose,trait_plot_an,'un',0);
    trait_plot_pl = cellfun(@transpose,trait_plot_pl,'un',0);
    pop_plot_an = cellfun(@transpose,pop_plot_an,'un',0);
    pop_plot_pl = cellfun(@transpose,pop_plot_pl,'un',0);
    the_fig=plot_traits_with_fitness_bg_plant(trait_interval,sepIndex,time_plot,trait_plot_an,trait_plot_pl,pop_plot_an,pop_plot_pl,rP,const,h,y0,ky,sig_pl,sig_c_pl,sig_m,the_pl_tr);
    plot_name=strcat(the_folder_output,'/FigB7_traits_evo_sigM=',num2str(log_sig_m),'.png');
    print(plot_name,'-dpng','-r300');
    clf;
    
end

%Functions
function K=k_x(k0,x0,a,x)
        K=k0*(1-((x-x0)/a).^2);
        K(abs(x-x0)>a)=0;
end

function comp =competition_term_mutant(the_mutant,the_traits, the_pops,sig_c)
    %traits and the_pops should be row vectors, and comp as well
    mat_comp=repmat(the_mutant,[length(the_traits),1]);
    comp=exp(-(1/(2*sig_c^2))*((mat_comp-transpose(the_traits)).^2));
    comp=the_pops*comp;
end

    
function m=mutu(the_traits, the_other_traits, the_other_pops,const,sig_m,h)
    RR=repmat(transpose(the_traits),[1,length(the_other_traits)]);
    SS=repmat(the_other_traits,[length(the_traits),1]);
    benef=exp(-(1/(2*sig_m^2))*((RR-SS).^2));
    int_rateA_num=the_other_pops*transpose(benef);
    int_rateA_denum=1+h*(the_other_pops*transpose(benef));
    m=const*(int_rateA_num./int_rateA_denum);
end

