clear all;
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

parameters={[-1.9375 -2.8125], [-2 -2.625] [-1.875 -2.625]};
folders={'./Results_intra_guild/','./Results_cross_guild/','./Results_cross_guild/'};
positions=[[0.1 0.18 0.28 0.7];[0.4 0.18 0.28 0.7];[0.7 0.18 0.28 0.7]];
titles={'(a)','(b)','(c)'};

log_sig_an=1;
sig_an=1;
sig_pl=sig_an;

figure('visible','off');

for each_pl=1:3
    log_sig_c_an=parameters{each_pl}(1);
    sig_c_an=exp(log_sig_c_an);
    sig_c_pl=sig_c_an;
    log_sig_m=parameters{each_pl}(2);
    sig_m=exp(log_sig_m);
    the_folder=folders{each_pl};
    
    disp (each_pl);

    tr_an_last=read_last_line_to_array(strcat(the_folder,'animal_traits_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
    tr_pl_last=read_last_line_to_array(strcat(the_folder,'plant_traits_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
    ab_an_last=read_last_line_to_array(strcat(the_folder,'animal_abundance_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
    ab_pl_last=read_last_line_to_array(strcat(the_folder,'plant_abundance_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));

    sc=@(x) sign(x).*(log10(1+1000*abs(x)));   
    
    the_sub=subplot(1,3,each_pl);
    trait_interval=0.9:0.00005: 4.1;
    fit_an=fitness_land_function(trait_interval,ab_an_last,ab_pl_last,tr_an_last,tr_pl_last,rA,h,const,sig_c_an,sig_an,x0,kx,sig_m);
    fit_pl=fitness_land_function(trait_interval,ab_pl_last,ab_an_last,tr_pl_last,tr_an_last,rP,h,const,sig_c_pl,sig_pl,y0,ky,sig_m);
    
    pos_an=sc(fit_an);
    pos_an(pos_an<0)=0;
    fill(trait_interval,pos_an,[0.9 0.9 0.9]);
    hold on; 
    pos_pl=sc(fit_pl);
    pos_pl(pos_pl<0)=0;
    fill(trait_interval,pos_pl,[0.9 0.9 0.9]);
    hold on;
    p1=plot(trait_interval,sc(fit_an),'-','Color','b','LineWidth',0.8);
    hold on;
    p2=plot(trait_interval,sc(fit_pl),'-','Color',rgb('SkyBlue'),'LineWidth',0.8);
    hold on;
    plot(tr_an_last,sc(fitness_land_function(tr_an_last,ab_an_last,ab_pl_last,tr_an_last,tr_pl_last,rA,h,const,sig_c_an,sig_an,x0,kx,sig_m)),'o','Color','b','MarkerSize',3,'MarkerFaceColor','b');
    hold on;
    plot(tr_pl_last,sc(fitness_land_function(tr_pl_last,ab_pl_last,ab_an_last,tr_pl_last,tr_an_last,rP,h,const,sig_c_pl,sig_pl,y0,ky,sig_m)),'o','Color',rgb('SkyBlue'),'MarkerSize',3,'MarkerFaceColor',rgb('SkyBlue'));
    hold on;
    
   
    set (gca,'xlim',[trait_interval(1) trait_interval(end)],'FontSize',8);
    set(gca,'ylim',sc([-0.6 1.5]));
    set (gca,'ytick',sc(-0.4:0.1:0.8),'yticklabel',{'','','', -0.1,0,0.1, '','', '','','','',''});
    set(the_sub,'Position',positions(each_pl,:));
    ti=title(titles{each_pl},'FontSize',10);
    set(ti,'Position',[1.05 3.3 0],'FontWeight','Normal');
    if (each_pl==1)
        ylabel('Fitness','FontSize',8);
    else
        set(gca,'yticklabel',[]);
    end
    
    if each_pl==2
        
        text(3.15,sc(0.08),'Animals','FontSize',8,'Color','b');
        text(1.15,sc(0.08),'Plants','FontSize',8,'Color',rgb('SkyBlue'));
        xlabel('Trait','FontSize',8);
    end


end

width = 15;     % Width in inches
height = 5;    % Height in inches
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf,'PaperPosition', myfiguresize);

mkdir('./Plots');
plot_name=strcat('./Plots/Fig1.png');
print(plot_name,'-dpng','-r600');
clf;

