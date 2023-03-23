warning('off', 'MATLAB:MKDIR:DirectoryExists');

%Defining parameters
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

%initilizing dictionary of number, percentage of introduced species and empty niche width
invaders_per=containers.Map;
en_width=containers.Map;
invaders_number=containers.Map;

%looping through all simulated communities
for log_sig_c_an=log_sig_c_an_range
    for log_sig_m=log_sig_m_range
        log_sig_c_pl=log_sig_c_an;

        sig_c_an=exp(log_sig_c_an);
        sig_c_pl=sig_c_an;

        sig_m=exp(log_sig_m);
        
        disp ([log_sig_an log_sig_c_an log_sig_m]);
        
        if log_sig_m<-1.5 % where there are two different evolutionary stable points
            
            re=detect_EN_function('./Results_cross_guild/',log_sig_an,log_sig_c_an,log_sig_m,rA,rP,const,h,x0,y0,kx,ky);
            
            if (re{1}==1)%if there is empty niche when the traits evolve through cross-guild benefits
                an_tr=read_data_into_cell(strcat('./Results_invasion_cross_guild/animal_traits_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
                pl_tr=read_data_into_cell(strcat('./Results_invasion_cross_guild/plant_traits_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
                sep=read_file_to_array(strcat('./Results_invasion_cross_guild/sepIndex_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));

                the_id=define_native_invader_function(an_tr,pl_tr,sep);
                all_an_id=the_id{1};
                all_pl_id=the_id{2};

                an_id=all_an_id{end};
                pl_id=all_pl_id{end};
               
                an_invad_number=sum(an_id);
                pl_invad_number=sum(pl_id);

                proportion=(an_invad_number+pl_invad_number)/(length(an_id)+length(pl_id));
                number=an_invad_number+pl_invad_number;
                invaders_per(strcat(num2str(log_sig_c_an),',',num2str(log_sig_m),',mutu'))=proportion*100;
                invaders_number(strcat(num2str(log_sig_c_an),',',num2str(log_sig_m),',mutu'))=number;
                en_width(strcat(num2str(log_sig_c_an),',',num2str(log_sig_m),',mutu'))=length(re{3})+length(re{4});
            end
            
            re=detect_EN_function('./Results_intra_guild/',log_sig_an,log_sig_c_an,log_sig_m,rA,rP,const,h,x0,y0,kx,ky);
            
            if re{1}==1  %if there is empty niche when the traits evolve through intra-guild benefits
                an_tr=read_data_into_cell(strcat('./Results_invasion_intra_guild/animal_traits_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
                pl_tr=read_data_into_cell(strcat('./Results_invasion_intra_guild/plant_traits_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
                sep=read_file_to_array(strcat('./Results_invasion_intra_guild/sepIndex_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
                the_id=define_native_invader_function(an_tr,pl_tr,sep);
                all_an_id=the_id{1};
                all_pl_id=the_id{2};

                an_id=all_an_id{end};
                pl_id=all_pl_id{end};
                
                an_invad_number=sum(an_id);
                pl_invad_number=sum(pl_id);

                proportion=(an_invad_number+pl_invad_number)/(length(an_id)+length(pl_id));
                number=an_invad_number+pl_invad_number;
                invaders_per(strcat(num2str(log_sig_c_an),',',num2str(log_sig_m),',res'))=proportion*100;
                invaders_number(strcat(num2str(log_sig_c_an),',',num2str(log_sig_m),',res'))=number;
            
                en_width(strcat(num2str(log_sig_c_an),',',num2str(log_sig_m),',res'))=length(re{3})+length(re{4});
            end
                
        else  % when there is only one evolutionary stable point, consider the intra-guild resource case
            re=detect_EN_function('./Results_intra_guild/',log_sig_an,log_sig_c_an,log_sig_m,rA,rP,const,h,x0,y0,kx,ky);
            if re{1}==1  %if there is empty niche when the traits evolve through intra-guild benefits
                an_tr=read_data_into_cell(strcat('./Results_invasion_intra_guild/animal_traits_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
                pl_tr=read_data_into_cell(strcat('./Results_invasion_intra_guild/plant_traits_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
                sep=read_file_to_array(strcat('./Results_invasion_intra_guild/sepIndex_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
                the_id=define_native_invader_function(an_tr,pl_tr,sep);
                all_an_id=the_id{1};
                all_pl_id=the_id{2};

                an_id=all_an_id{end};
                pl_id=all_pl_id{end};
                
                an_invad_number=sum(an_id);
                pl_invad_number=sum(pl_id);

                proportion=(an_invad_number+pl_invad_number)/(length(an_id)+length(pl_id));
                number=an_invad_number+pl_invad_number;
                invaders_per(strcat(num2str(log_sig_c_an),',',num2str(log_sig_m),',all'))=proportion*100;
                invaders_number(strcat(num2str(log_sig_c_an),',',num2str(log_sig_m),',all'))=number;
            
                en_width(strcat(num2str(log_sig_c_an),',',num2str(log_sig_m),',all'))=length(re{3})+length(re{4});
            end
        end
      
    end
end

figg=figure('Visible','off');
box('on');

all_keys=invaders_per.keys();
per_array=zeros(length(all_keys),1);
number_array=zeros(length(all_keys),1);
width_array=zeros(length(all_keys),1);
for each_key_ind=1:length(all_keys)
    each_key=all_keys{each_key_ind};
    per_array(each_key_ind)=invaders_per(each_key);
    number_array(each_key_ind)=invaders_number(each_key);
    width_array(each_key_ind)=(en_width(each_key)-1)*0.001;
end
subplot(1,2,1);
[y_out, low_lim, high_lim, ~]=lowess([width_array per_array],1,1,'lowess1.png');
fill([low_lim(:,1); flipud(high_lim(:,1))],[low_lim(:,2); flipud(high_lim(:,2))],rgb('lightSkyBlue'),'LineStyle','None');
hold on;
plot(width_array,per_array,'o','MarkerSize',2,'MarkerFaceColor','none','MarkerEdgeColor',rgb('SlateGray'),'LineWidth',0.3);
hold on;
plot(y_out(:,1),y_out(:,3),'b-','LineWidth',1);
set(gca,'Xlim',[-0.05 3.65]);
set(gca,'Ylim',[-3 103]);
set(gca,'TickDir','in','TickLength',[0.015 0.015],'XMinorTick','off','YMinorTick','off','layer','top','LineWidth',1.2,'FontSize',9);
pbaspect([1 1 1]);
title('(a)','FontSize',12,'FontWeight','Normal','Units','normalized','Position',[0.04 1.04]);
xlab='Width of empty niches';
ylab='Proportion of introduced species (%)';

xlabel(xlab,'FontSize',9);
ylabel(ylab,'FontSize',9);

subplot(1,2,2);
[y_out, low_lim, high_lim, ~]=lowess([width_array number_array],1,1,'lowess1.png');
fill([low_lim(:,1); flipud(high_lim(:,1))],[low_lim(:,2); flipud(high_lim(:,2))],rgb('lightSkyBlue'),'LineStyle','None');
hold on;
plot(width_array,number_array,'o','MarkerSize',2,'MarkerFaceColor','none','MarkerEdgeColor',rgb('SlateGray'),'LineWidth',0.3);
hold on;
plot(y_out(:,1),y_out(:,3),'b-','LineWidth',1);
set(gca,'Xlim',[-0.05 3.7]);
set(gca,'Ylim',[-0.5 18]);
set(gca,'TickDir','in','TickLength',[0.015 0.015],'XMinorTick','off','YMinorTick','off','layer','top','LineWidth',1.2,'FontSize',9);
set(gca,'YTick',0:4:18);
pbaspect([1 1 1]);
title('(b)','FontSize',12,'FontWeight','Normal','Units','normalized','Position',[0.04 1.04]);
xlab='Width of empty niches';
ylab='Number of introduced species';

xlabel(xlab,'FontSize',9);
ylabel(ylab,'FontSize',9);
mkdir('./Plots/')
name=strcat('./Plots/Fig5.png');

print(name,'-dpng','-r600');
clf;
