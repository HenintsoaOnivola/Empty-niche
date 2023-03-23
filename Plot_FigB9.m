warning('off', 'MATLAB:MKDIR:DirectoryExists');
mkdir('./Supplementary_figures/');
%Define parameters
log_sig_an=1;
log_sig_c_an_range=-2.5:0.0625:1;
log_sig_m_range=-3:0.0625:-0.5;


%%Compute the average over all traits of the selection pressure at each time step, then store the results in the existing selection pressure folders
%to avoid recomputing when plotting
res_folder='./Selection_pressure_intra_guild/';
mutu_folder='./Selection_pressure_cross_guild/';

average_sel_res=cell(100,1);
average_sel_inv=cell(100,1);
time_array=cell(100,1);

count=1;

for log_sig_c_an=log_sig_c_an_range
    for log_sig_m=log_sig_m_range
        disp ([log_sig_an log_sig_c_an log_sig_m]);
        
        file_res=strcat(res_folder,'an_sel_res_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat');
        file_mutu=strcat(mutu_folder,'an_sel_res_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat');
        if log_sig_m<-1.5
            if isfile(file_res)
                an_res=dlmread(file_res);
                pl_res=dlmread(strcat(res_folder,'pl_sel_res_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
                an_res(an_res==0)=NaN;
                pl_res(pl_res==0)=NaN;
                mean_an_res=nanmean(an_res,2);
                mean_pl_res=nanmean(pl_res,2);
                mean_res=(mean_an_res+mean_pl_res)/2;
                average_sel_res{count}=mean_res;
                dlmwrite(strcat(res_folder,'average_sel_res_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'),mean_res,'precision',20);

                an_inv=dlmread(strcat(res_folder,'an_sel_inv_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
                pl_inv=dlmread(strcat(res_folder,'pl_sel_inv_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
                an_inv(an_inv==0)=NaN;
                pl_inv(pl_inv==0)=NaN;
                mean_an_inv=nanmean(an_inv,2);
                mean_pl_inv=nanmean(pl_inv,2);
                mean_inv=(mean_an_inv+mean_pl_inv)/2;
                average_sel_inv{count}=mean_inv;
                dlmwrite(strcat(res_folder,'average_sel_inv_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'),mean_inv,'precision',20);
                time=dlmread(strcat('./Results_invasion_intra_guild/time_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
                time_array{count}=time;
                
                count=count+1;
             end
             if isfile(file_mutu)
                an_res=dlmread(file_mutu);
                pl_res=dlmread(strcat(mutu_folder,'pl_sel_res_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
                an_res(an_res==0)=NaN;
                pl_res(pl_res==0)=NaN;
                mean_an_res=nanmean(an_res,2);
                mean_pl_res=nanmean(pl_res,2);
                mean_res=(mean_an_res+mean_pl_res)/2;
                average_sel_res{count}=mean_res;
                dlmwrite(strcat(mutu_folder,'average_sel_res_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'),mean_res,'precision',20);
                
                an_inv=dlmread(strcat(mutu_folder,'an_sel_inv_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
                pl_inv=dlmread(strcat(mutu_folder,'pl_sel_inv_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
                an_inv(an_inv==0)=NaN;
                pl_inv(pl_inv==0)=NaN;
                mean_an_inv=nanmean(an_inv,2);
                mean_pl_inv=nanmean(pl_inv,2);
                mean_inv=(mean_an_inv+mean_pl_inv)/2;
                average_sel_inv{count}=mean_inv;
                dlmwrite(strcat(mutu_folder,'average_sel_inv_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'),mean_inv,'precision',20);
                time=dlmread(strcat('./Results_invasion_cross_guild/time_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
                time_array{count}=time;
                count=count+1;
            end
        else
          if isfile(file_res)
            an_res=dlmread(file_res);
            pl_res=dlmread(strcat(res_folder,'pl_sel_res_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
            an_res(an_res==0)=NaN;
            pl_res(pl_res==0)=NaN;
            mean_an_res=nanmean(an_res,2);
            mean_pl_res=nanmean(pl_res,2);
            mean_res=(mean_an_res+mean_pl_res)/2;
            average_sel_res{count}=mean_res;
            dlmwrite(strcat(res_folder,'average_sel_res_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'),mean_res,'precision',20);

            an_inv=dlmread(strcat(res_folder,'an_sel_inv_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
            pl_inv=dlmread(strcat(res_folder,'pl_sel_inv_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
            an_inv(an_inv==0)=NaN;
            pl_inv(pl_inv==0)=NaN;
            mean_an_inv=nanmean(an_inv,2);
            mean_pl_inv=nanmean(pl_inv,2);
            mean_inv=(mean_an_inv+mean_pl_inv)/2;
            average_sel_inv{count}=mean_inv;
            dlmwrite(strcat(res_folder,'average_sel_inv_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'),mean_inv,'precision',20);
            time=dlmread (strcat('./Results_invasion_intra_guild/time_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
            time_array{count}=time;
            count=count+1;

          end
        end
              
    end
end

%Read from the selection pressure files, and create a cell array for all the
%simulations
count=1;
average_sel_res=cell(100,1);
average_sel_inv=cell(100,1);
time_array=cell(100,1);
for log_sig_c_an=log_sig_c_an_range
    for log_sig_m=log_sig_m_range
        disp ([log_sig_an log_sig_c_an log_sig_m]);
        mutu_file=strcat(mutu_folder,'average_sel_res_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat');
        res_file=strcat(res_folder,'average_sel_res_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat');
        if log_sig_m<-1.5  %case of two different evolutionary trajectories
            if isfile(res_file)
                mean_res=dlmread(strcat(res_folder,'average_sel_res_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
                mean_inv=dlmread(strcat(res_folder,'average_sel_inv_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
                time=dlmread(strcat('./Results_invasion_intra_guild/time_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
                average_sel_res{count}=mean_res;
                average_sel_inv{count}=mean_inv;
                time_array{count}=time;
                count=count+1;
            end
            if isfile(mutu_file)
                mean_res=dlmread(strcat(mutu_folder,'average_sel_res_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
                mean_inv=dlmread(strcat(mutu_folder,'average_sel_inv_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
                time=dlmread(strcat('./Results_invasion_cross_guild/time_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
                average_sel_res{count}=mean_res;
                average_sel_inv{count}=mean_inv;
                time_array{count}=time;
                count=count+1;
            end
        else %when there is only one evolutionary trajectory
            if isfile(res_file)
                mean_res=dlmread(strcat(res_folder,'average_sel_res_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
                mean_inv=dlmread(strcat(res_folder,'average_sel_inv_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
                time=dlmread(strcat('./Results_invasion_intra_guild/time_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
                average_sel_res{count}=mean_res;
                average_sel_inv{count}=mean_inv;
                time_array{count}=time;
                count=count+1;
            end
                
        end
    end
 
end

for each_sim=1:length(time_array)
    tt=time_array{each_sim};
    new_tt=tt(tt<=(5.5*10^6));
    if length(new_tt)~=length(tt)
        average_sel_res{each_sim}=average_sel_res{each_sim}(1:length(new_tt));
        average_sel_inv{each_sim}=average_sel_inv{each_sim}(1:length(new_tt));
        time_array{each_sim}=new_tt;
    end
end

all_time=[];
for sim=1:length(time_array)
    all_time=[all_time; time_array{sim}];
end
all_time=unique(all_time);
all_time=sort(all_time);


%%time-averaging the selection pressures
time_partition=all_time(end)/400;
considered_time_upper=time_partition:time_partition:all_time(end);
considered_time_lower=0:time_partition:(considered_time_upper(end)-time_partition);
considered_time=(considered_time_upper-considered_time_lower)/2;

res_values=zeros(length(time_array),400);
inv_values=zeros(length(time_array),400);

for each_case=1:length(time_array)
    each_time_array=time_array{each_case};
    each_res=average_sel_res{each_case};
    each_inv=average_sel_inv{each_case};
    previous_time=-1;
    for time_ind=1:length(considered_time_upper)
        disp ([each_case time_ind]);
        the_indexes=(each_time_array>previous_time & each_time_array<=considered_time_upper(time_ind));
        res_values(each_case,time_ind)=nanmean(each_res(the_indexes));
        inv_values(each_case,time_ind)=nanmean(each_inv(the_indexes));
        
        previous_time=considered_time_upper(time_ind);
        
    end

end

log_res_values=log10(res_values);
log_inv_values=log10(inv_values);

the_max=max([max(log_res_values(:)) max(log_inv_values(:))]);
the_min=min([min(log_res_values(:)) min(log_inv_values(:))]);
sel_div=(the_max-the_min)/200;
gr=the_min:sel_div:the_max;

%create the matrices of freqency to be plotted
res_average=nanmean(log_res_values,1);
inv_average=nanmean(log_inv_values,1);
res_std=nanstd(log_res_values,1);
inv_std=nanstd(log_inv_values,1);

res_average=200*(res_average-the_min)/(the_max-the_min);
inv_average=200*(inv_average-the_min)/(the_max-the_min);

res_std=200*(res_std)/(the_max-the_min);
inv_std=200*(inv_std)/(the_max-the_min);

res_matrix=zeros(201,401);
inv_matrix=zeros(201,401);
for t=1:400
    all_res=log_res_values(:,t);
    all_inv=log_inv_values(:,t);
    for j=1:200
        disp ([t j]);
        count_res=length(all_res(all_res>gr(j) & all_res<=gr(j+1)));
        count_inv=length(all_inv(all_inv>gr(j) & all_inv<=gr(j+1)));
        if j==1
            count_res=count_res+length(all_res(all_res==the_min));
            count_inv=count_inv+length(all_inv(all_inv==the_min));
        end
        res_matrix(j,t)=count_res;
        inv_matrix(j,t)=count_inv;
        if count_res==0
            res_matrix(j,t)=NaN;
        end
        if count_inv==0
            inv_matrix(j,t)=NaN;
        end
    end
   
end

res_matrix(:,end)=NaN;
res_matrix(end,:)=NaN;
inv_matrix(:,end)=NaN;
inv_matrix(end,:)=NaN;

all_min=min([min(res_matrix(:)) min(inv_matrix(:))]);
all_max=max([max(res_matrix(:)) max(inv_matrix(:))]);

res_matrix(end,1)=all_min;
res_matrix(end,2)=all_max;
inv_matrix(end,1)=all_min;
inv_matrix(end,2)=all_max;

res_matrix=log10(res_matrix);
inv_matrix=log10(inv_matrix);
log_min=min([min(res_matrix(:)) min(inv_matrix(:))]);
log_max=max([max(res_matrix(:)) max(inv_matrix(:))]);

%Define figure properties (colormap, ticks, labels, positions, ...)
my_colormap=hot(60);
my_colormap=(flipud(my_colormap));
the_first=my_colormap(1,:);
my_colormap=my_colormap (15:end,:);
colormap(my_colormap);

the_x_ticks=0:10^(6):all_time(end);
the_x_ticks=(the_x_ticks*400)/all_time(end);
the_x_ticks_label=0:1:5;
the_x_ticks(1)=1;

relative_sel_pres=[10.^(-9:3:-3) 1];
sel_p=(relative_sel_pres*(10^the_max-10^the_min))+10^the_min;
power_tick=log10(sel_p);
the_y_ticks=((power_tick-the_min)*200)/(the_max-the_min);
the_y_ticks=[1 the_y_ticks];
the_y_ticks(end)=201;
the_y_ticks_label={'0','10^{-9}','10^{-6}','10^{-3}','1'};

%%plot the case for resident species
figure('Visible','off');
colormap(my_colormap);
pp_res=pcolor(res_matrix);
set(pp_res,'edgecolor','none');
hold on;
plot(0.5:1:400,res_average,'k-','LineWidth',1);
hold on;
on_x=0.5:1:400;
errorbar(on_x(1:40:end),res_average(1:40:end),res_std(1:40:end),'ko','LineWidth',1,'MarkerSize',0.2);
set(gca,'TickDir','in','TickLength',[0.015 0.015],'XMinorTick','off','YMinorTick','off','layer','top','LineWidth',1.5);
set(gca,'Xtick',the_x_ticks,'XtickLabel',the_x_ticks_label,'Ytick',the_y_ticks,'YTickLabel',the_y_ticks_label,'FontSize',15);
ylim=get(gca,'ylim');
pbaspect([2 1 1]);

cc=colorbar();
set(cc,'Ylim',[log_min log_max]);
y_t=0:1:3;
y_t_lab=10.^(y_t);
set(cc,'YTick',y_t,'YtickLabel',y_t_lab,'FontSize',15,'TickLength',0.015,'LineWidth',1.5);
ycb=ylabel(cc,'Frequency','FontSize',15);
set(ycb,'Rotation',-90);
label_pos=cc.Label.Position;
label_pos(1)=label_pos(1)+0.7;
cc.Label.Position=label_pos;

xlab='Evolutionary time (10^{6})';
ylab='Selection pressure';
xlabel(xlab,'FontSize',15,'Interpreter','tex');
yy=ylabel(ylab,'FontSize',15,'Interpreter','tex');
set(yy,'Position',get(yy,'Position')+[2 0 0]);
handle=title('Resident');
set(handle,'Position',[345 170 0],'FontSize',15,'FontWeight','Normal'); 
set(gca,'Position',[0.13 0.11 0.67 0.81]);

name=strcat('./Supplementary_figures/FigB9a.png');
print(name,'-dpng','-r600');

%Plot the case for introduced species
figure('Visible','off');
colormap(my_colormap);
pp_inv=pcolor(inv_matrix);
set(pp_inv,'edgecolor','none');
hold on;
plot(0.5:1:400,inv_average,'k-','LineWidth',1);
hold on;
on_x=0.5:1:400;
errorbar(on_x(1:40:end),inv_average(1:40:end),inv_std(1:40:end),'ko','LineWidth',1,'MarkerSize',0.2);
set(gca,'TickDir','in','TickLength',[0.015 0.015],'XMinorTick','off','YMinorTick','off','layer','top','LineWidth',1.5);
set(gca,'Xtick',the_x_ticks,'XtickLabel',the_x_ticks_label,'Ytick',the_y_ticks,'YTickLabel',the_y_ticks_label,'FontSize',15);
pbaspect([2 1 1]);

cc=colorbar();
set(cc,'Ylim',[log_min log_max]);
y_t=0:1:3;
y_t_lab=10.^(y_t);
set(cc,'YTick',y_t,'YtickLabel',y_t_lab,'FontSize',15,'TickLength',0.015,'LineWidth',1.5);
ycb=ylabel(cc,'Frequency','FontSize',15);
set(ycb,'Rotation',-90);
label_pos=cc.Label.Position;
label_pos(1)=label_pos(1)+0.7;
cc.Label.Position=label_pos;

xlab='Evolutionary time (10^{6})';
ylab='Selection pressure';
xlabel(xlab,'FontSize',15,'Interpreter','tex');
yy=ylabel(ylab,'FontSize',15,'Interpreter','tex');
set(yy,'Position',get(yy,'Position')+[2 0 0]);
handle=title('Introduced');
set(handle,'Position',[345 170 0],'FontSize',15,'FontWeight','Normal'); 
set(gca,'Position',[0.13 0.11 0.67 0.81]);

name=strcat('./Supplementary_figures/FigB9b.png');
print(name,'-dpng','-r600');
