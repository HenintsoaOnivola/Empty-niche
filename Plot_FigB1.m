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

% to control the figure size and resolution
width = 3;     % Width in inches
height = 3;    % Height in inches
alw = 0.75;    % AxesLineWidth
fsz = 12;      % Fontsize
lw = 0.5;      % LineWidth
msz = 8;       % MarkerSize

log_sig_an=1;
log_sig_c_an=-1;
log_sig_m=-2.5;

sig_an=log_sig_an;
sig_pl=sig_an;

sig_c_an=exp(log_sig_c_an);
sig_c_pl=sig_c_an;
sig_m=exp(log_sig_m);
disp ([sig_an sig_c_an sig_m]);

generate_traits_phase_portrait(rA,rP,const,h,x0,y0,kx,ky,log_sig_c_an,log_sig_m,log_sig_an);

fil=fopen('data_for_FigB1.dat');
file_content=cell(1000,1);
line=fgetl(fil);
count_index_line=1;
while(ischar(line))
    ll=textscan(line,'%f');
    file_content{count_index_line}=ll{1};
    if isempty(ll{1})==1
        ll=textscan(line,'%s');
        file_content{count_index_line}=ll{1};
    end
    count_index_line=count_index_line+1;
    line=fgetl(fil);
end
fclose(fil);
file_content=file_content(~cellfun('isempty',file_content));
min_an=file_content{1}(1);
max_an=file_content{1}(2);
min_pl=file_content{1}(3);
max_pl=file_content{1}(4);

%set figure properties
figg=figure('Visible','off');
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties

line_ind=2;

all_trait=cell(5,1);
count_pair=1;
count_trait=1;
while line_ind<=length(file_content)
    time_plot=file_content{line_ind+1};
    an_trait=file_content{line_ind+4};
    pl_trait=file_content{line_ind+5};
    anc=file_content{line_ind};
    if strcmp(file_content{line_ind+6}(1),'not')==1
        if rem(count_pair,2)==1
            plot(an_trait,pl_trait,'LineStyle','-','Color',[0.5 0.5 0.5],'LineWidth',0.4);
            hold all;
            all_trait{count_trait}=[an_trait pl_trait];
            count_trait=count_trait+1;
        end
        count_pair=count_pair+1;
    end
    line_ind=line_ind+7;
end


axis([min_pl-0.5 max_pl+0.5 min_an-0.5 max_an+0.5]);

%%%%%%%%%%%%%%%%%%%%%%%%
%localize the population equilibrium
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
equi_pop=cell(4,1);
line_ind=1;
while line_ind<=length(file_content)
    if strcmp(file_content{line_ind+6}(1),'not')==1
        an_pop=file_content{line_ind+2};
        pl_pop=file_content{line_ind+3};
        equi_pop{1}=[an_pop(end) pl_pop(end)];
        break;
    end
    line_ind=line_ind+1;
end
equi_pop_ind=2;
while line_ind<=length(file_content)
    if strcmp(file_content{line_ind+6}(1),'not')==1
        an_pop=file_content{line_ind+2};
        pl_pop=file_content{line_ind+3};
        in_equi=false;
        equi_pop=equi_pop(~cellfun('isempty',equi_pop));
        for each_eq=1:length(equi_pop)
            if abs(equi_pop{each_eq}(1)-an_pop(end))<=0.1 && abs(equi_pop{each_eq}(2)-pl_pop(end))<=0.1
                in_equi=in_equi || true; 
            end
        end
        if in_equi==false
            equi_pop{equi_pop_ind}=[an_pop(end) pl_pop(end)];
            equi_pop_ind=equi_pop_ind+1;
        end
    end
    line_ind=line_ind+7;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%localize the equilibrium
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
equi_points=cell(4,1);
line_ind=1;
while line_ind<=length(file_content)
    if strcmp(file_content{line_ind+6}(1),'not')==1
        an_trait=file_content{line_ind+4};
        pl_trait=file_content{line_ind+5};
        equi_points{1}=[an_trait(end) pl_trait(end)];
        break;
    end
    line_ind=line_ind+1;
end
equi_points_ind=2;
while line_ind<=length(file_content)
    if strcmp(file_content{line_ind+6}(1),'not')==1
        an_trait=file_content{line_ind+4};
        pl_trait=file_content{line_ind+5};
        in_equi=false;
        equi_points=equi_points(~cellfun('isempty',equi_points));
        for each_eq=1:length(equi_points)
            if abs(equi_points{each_eq}(1)-an_trait(end))<=0.1 && abs(equi_points{each_eq}(2)-pl_trait(end))<=0.1
                in_equi=in_equi || true; 
            end
        end
        if in_equi==false
            equi_points{equi_points_ind}=[an_trait(end) pl_trait(end)];
            equi_points_ind=equi_points_ind+1;
        end
    end
    line_ind=line_ind+7;
end


%see the type of equilibrium

%plot diagonal line
x_lim=xlim();
y_lim=ylim();
plot ([x_lim(1) x_lim(2)],[x_lim(1) x_lim(2)],'k--','LineWidth',1);

%plot the arrow head
line_ind=2;
count=1;
for each_line=1:length(all_trait)
    an_trait=all_trait{each_line}(:,1);
    pl_trait=all_trait{each_line}(:,2);
    input_vec=[an_trait pl_trait];
    [~, idxs, ~] = unique(input_vec, 'rows');
    un=input_vec(sort(idxs),:);
    cons_points=zeros(size(un));
    cons_points(:)=NaN;
    cons_points(1,:)=un(1,:);
    for un_ind=2:size(un,1)
        if length(equi_points)==2
            if sqrt((abs(un(un_ind,1)-equi_points{1}(1)))^2 + (abs(un(un_ind,2)-equi_points{1}(2)))^2 )>0.2  && sqrt((abs(un(un_ind,1)-equi_points{2}(1)))^2 + (abs(un(un_ind,2)-equi_points{2}(2)))^2 )>0.2
                cons_points(un_ind,:)=un(un_ind,:);
            else
                break;
            end

        else
            if  sqrt((abs(un(un_ind,1)-equi_points{1}(1)))^2 + (abs(un(un_ind,2)-equi_points{1}(2)))^2 )>0.2
                cons_points(un_ind,:)=un(un_ind,:);
            else
                break;
            end
        end
    end
    cons_points=cons_points(1:(un_ind-1),:);
    if rem(count,3)==1
        if size(cons_points,1)>1
            arrowh(cons_points(:,1),cons_points(:,2),[0.4 0.4 0.4],[100 50],70);
        end
    end
    count=count+1;
end

%plot equi points
for eq=1:length(equi_points)
    plot(equi_points{eq}(1),equi_points{eq}(2),'ko','MarkerFaceColor','k', 'MarkerEdgeColor','k','MarkerSize',5);
    hold on;
end

xlabel('Animal trait value,{\it x}','FontSize',12);
ylabel('Plant trait value,{\it y}','FontSize',12);

axis([2 4 1 3]);

set(gca,'XTick',1:1:4,'FontSize',10);
set(gca,'YTick',1:1:4,'FontSize',10);
set(gca,'LineWidth',1);

%Here we preserve the size of the image when we save it.
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf,'PaperPosition', myfiguresize);

% Save the file as PNG
mkdir('./Supplementary_figures/');
plot_name=strcat('./Supplementary_figures/FigB1.png');
print(plot_name,'-dpng','-r400');

clf;
delete('data_for_FigB1.dat');
