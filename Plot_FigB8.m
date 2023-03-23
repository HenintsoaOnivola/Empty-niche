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

log_sig_c_an_range=-2.5:0.125:1;
log_sig_m_range=-3:0.125:-0.5;

log_sig_an=0.75;
sig_an=log_sig_an;
sig_pl=sig_an;

%This first part of the code only generates all the selection pressure values
%along time in two folders.Saving the results in these folders is to
%avoid re-running the entire code for the purpose of plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp ('Part1: computing the selection pressure values and storing them into files');

res_opt={'intra','cross'};
for each_opt=1:length(res_opt)
    disp (strcat('optimizing for ', res_opt{each_opt},'_guild resources'));
    
    the_folder=strcat('./Results_robustness/vary_sigma_a_',res_opt{each_opt},'_invasion/');
    the_dest=strcat('./Results_robustness/vary_sigma_a_',res_opt{each_opt},'_invasion/Selection_pressure/');
    mkdir(the_dest);
    
    for log_sig_c_an=log_sig_c_an_range
        for log_sig_m=log_sig_m_range

            log_sig_c_pl=log_sig_c_an;
            sig_c_an=exp(log_sig_c_an);
            sig_c_pl=sig_c_an;
            sig_m=exp(log_sig_m);

            disp ([log_sig_an log_sig_c_an log_sig_m]);
            
            intro_file=strcat(the_folder,'introIndex_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat');
            
            if isfile(intro_file)

                time=dlmread(strcat(the_folder,'time_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
                sep_index=dlmread(strcat(the_folder,'sepIndex_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
                an_tr=dlmread(strcat(the_folder,'animal_traits_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
                an_tr=mat2cell(an_tr,ones(1,size(an_tr,1)),size(an_tr,2));
                an_tr=cellfun(@(x) x(x~=0),an_tr,'UniformOutput',0);
                pl_tr=dlmread(strcat(the_folder,'plant_traits_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
                pl_tr=mat2cell(pl_tr,ones(1,size(pl_tr,1)),size(pl_tr,2));
                pl_tr=cellfun(@(x) x(x~=0),pl_tr,'UniformOutput',0);
                an_ab=read_data_into_cell(strcat(the_folder,'animal_abundance_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
                pl_ab=read_data_into_cell(strcat(the_folder,'plant_abundance_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
   
                disp ('done with reading previous files');
                
                the_id=define_native_invader_function(an_tr,pl_tr,sep_index);
                an_id=the_id{1};
                pl_id=the_id{2};
                disp ('indentity of species computed');

                sel_grad_resident_an=cell(length(time),1);
                sel_grad_inv_an=cell(length(time),1);
                sel_grad_resident_pl=cell(length(time),1);
                sel_grad_inv_pl=cell(length(time),1);

                all_sel_grad=selection_gradient_function_different_time(an_ab, pl_ab, an_tr, pl_tr,rA,rP,h,const,sig_c_an,sig_c_pl,sig_an,sig_pl,x0,y0,kx,ky,sig_m);
                disp ('selection gradient computed');

                start=1;
                for ind=sep_index
                    ind_range=start:1:ind;
                    for i=ind_range
                        sel_grad=all_sel_grad{i};

                        sel_grad_an=sel_grad(1:length(an_tr{i}));
                        ind_an_res=find(an_id{i}==0);
                        sel_grad_resident_an{i}=sel_grad_an(ind_an_res);

                        ind_an_inv=find(an_id{i}==1);
                        sel_grad_inv_an{i}=sel_grad_an(ind_an_inv);

                        sel_grad_pl=sel_grad(length(an_tr{i})+1:end);
                        ind_pl_res=find(pl_id{i}==0);
                        sel_grad_resident_pl{i}=sel_grad_pl(ind_pl_res);
                        ind_pl_inv=find(pl_id{i}==1);
                        sel_grad_inv_pl{i}=sel_grad_pl(ind_pl_inv);
                    end
                    start=ind+1;
                end

                the_empty=cellfun(@isempty,sel_grad_resident_an);
                sel_grad_resident_an(the_empty)={NaN};

                the_empty=cellfun(@isempty,sel_grad_inv_an);
                sel_grad_inv_an(the_empty)={NaN};

                the_empty=cellfun(@isempty,sel_grad_resident_pl);
                sel_grad_resident_pl(the_empty)={NaN};

                the_empty=cellfun(@isempty,sel_grad_inv_pl);
                sel_grad_inv_pl(the_empty)={NaN};
                
                mkdir(the_dest);

                clear_files(strcat(the_dest,'an_sel_res_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
                clear_files(strcat(the_dest,'an_sel_inv_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
                clear_files(strcat(the_dest,'pl_sel_res_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
                clear_files(strcat(the_dest,'pl_sel_inv_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
                
                write_into_files_W_sel_pres(strcat(the_dest,'an_sel_res_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'),sel_grad_resident_an);
                write_into_files_W_sel_pres(strcat(the_dest,'an_sel_inv_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'),sel_grad_inv_an);
                write_into_files_W_sel_pres(strcat(the_dest,'pl_sel_res_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'),sel_grad_resident_pl);
                write_into_files_W_sel_pres(strcat(the_dest,'pl_sel_inv_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'),sel_grad_inv_pl);
            end
        end
    end

end

%This second part of the code plots FigB8 by reading data from those files
%generated from part 1. If all the necessary data have already been
%generated, one can only run this second part of the code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sel_folder_res=strcat('./Results_robustness/vary_sigma_a_intra_invasion/Selection_pressure/');
sel_folder_mutu=strcat('./Results_robustness/vary_sigma_a_cross_invasion/Selection_pressure/');
tr_folder_res='./Results_robustness/vary_sigma_a_intra_invasion/';
tr_folder_mutu='./Results_robustness/vary_sigma_a_cross_invasion/';

%initialize all the array containing the selection pressure values
sel_res_early=[];
tr_res_early=[];
sel_inv_early=[];
tr_inv_early=[];
sel_res_int=[];
tr_res_int=[];
sel_inv_int=[];
tr_inv_int=[];
sel_res_late=[];
tr_res_late=[];
sel_inv_late=[];
tr_inv_late=[];

for log_sig_c_an=log_sig_c_an_range
    for log_sig_m=log_sig_m_range
        disp ([log_sig_an log_sig_c_an log_sig_m]);
        
        file_res=strcat(sel_folder_res,'an_sel_res_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat');
        file_mutu=strcat(sel_folder_mutu,'an_sel_res_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat');
        %%intra-guild case
        if isfile(file_res)
            an_res=dlmread_empty(file_res,'\t',0,0,NaN);
            an_res=an_res(:,1:end-1);
            pl_res=dlmread_empty(strcat(sel_folder_res,'pl_sel_res_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'),'\t',0,0,NaN);
            pl_res=pl_res(:,1:end-1);

            an_inv=dlmread_empty(strcat(sel_folder_res,'an_sel_inv_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'),'\t',0,0,NaN);
            an_inv=an_inv(:,1:end-1);
            pl_inv=dlmread_empty(strcat(sel_folder_res,'pl_sel_inv_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'),'\t',0,0,NaN);
            pl_inv=pl_inv(:,1:end-1);

            an_tr=dlmread(strcat(tr_folder_res,'animal_traits_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
            an_tr_for_id=mat2cell(an_tr,ones(1,size(an_tr,1)),size(an_tr,2));
            an_tr_for_id=cellfun(@(x) x(x~=0),an_tr_for_id,'UniformOutput',0);
            pl_tr=dlmread(strcat(tr_folder_res,'plant_traits_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
            pl_tr_for_id=mat2cell(pl_tr,ones(1,size(pl_tr,1)),size(pl_tr,2));
            pl_tr_for_id=cellfun(@(x) x(x~=0),pl_tr_for_id,'UniformOutput',0);
            sep_index=dlmread(strcat(tr_folder_res,'sepIndex_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
            the_id=define_native_invader_function(an_tr_for_id,pl_tr_for_id,sep_index);
            an_id=the_id{1};
            pl_id=the_id{2};
            an_id=cell_to_mat_fill_nan(an_id);
            pl_id=cell_to_mat_fill_nan(pl_id);
            an_tr_res=nan(size(an_res));
            an_tr_inv=nan(size(an_inv));
            pl_tr_res=nan(size(pl_res));
            pl_tr_inv=nan(size(pl_inv));
            for i=1:size(an_tr,1)
                an_trait=an_tr(i,:);
                pl_trait=pl_tr(i,:);
                if sum(an_id(i,:)==0)>0
                    an_tr_res(i,:)=[an_trait(an_id(i,:)==0) nan(1,sum(isnan(an_res(i,:))))];
                end
                if sum(an_id(i,:)==1)>0
                    an_tr_inv(i,:)=[an_trait(an_id(i,:)==1) nan(1,sum(isnan(an_inv(i,:))))] ;
                end
                if sum(pl_id(i,:)==0)>0
                    pl_tr_res(i,:)=[pl_trait(pl_id(i,:)==0) nan(1,sum(isnan(pl_res(i,:))))];
                end
                if sum(pl_id(i,:)==1)>0
                    pl_tr_inv(i,:)=[pl_trait(pl_id(i,:)==1) nan(1,sum(isnan(pl_inv(i,:))))];
                end
            end

            time=dlmread(strcat(tr_folder_res,'time_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));

            %strore the selection pressure values for early
            %times 
            early=find(time<=0.3*10^6);
            early=early(end);

            res_all=[an_res pl_res];
            res_early=res_all(1:early,:);
            res_early=res_early(:);
            sel_res_early=[sel_res_early; res_early];

            res_all=[an_tr_res pl_tr_res];
            res_early=res_all(1:early,:);
            res_early=res_early(:);
            tr_res_early=[tr_res_early; res_early];

            inv_all=[an_inv pl_inv];
            inv_early=inv_all(1:early,:);
            inv_early=inv_early(:);
            sel_inv_early=[sel_inv_early; inv_early];

            inv_all=[an_tr_inv pl_tr_inv];
            inv_early=inv_all(1:early,:);
            inv_early=inv_early(:);
            tr_inv_early=[tr_inv_early; inv_early];

            %strore the selection pressure values for intermediate
            %times 
            int=find(time>0.3*10^6 & time<=2*10^6);
            if ~isempty(int)
                int=int(end);
                res_all=[an_res pl_res];
                res_int=res_all(1:int,:);
                res_int=res_int(:);
                sel_res_int=[sel_res_int; res_int];

                res_all=[an_tr_res pl_tr_res];
                res_int=res_all(1:int,:);
                res_int=res_int(:);
                tr_res_int=[tr_res_int; res_int];

                inv_all=[an_inv pl_inv];
                inv_int=inv_all(1:int,:);
                inv_int=inv_int(:);
                sel_inv_int=[sel_inv_int; inv_int];

                inv_all=[an_tr_inv pl_tr_inv];
                inv_int=inv_all(1:int,:);
                inv_int=inv_int(:);
                tr_inv_int=[tr_inv_int; inv_int];
            end

            %strore the selection pressure values for late
            %times 
            late=find(time>2*10^6);
            if ~isempty(late)
                late=late(end);
                res_all=[an_res pl_res];
                res_late=res_all(1:late,:);
                res_late=res_late(:);
                sel_res_late=[sel_res_late; res_late];

                res_all=[an_tr_res pl_tr_res];
                res_late=res_all(1:late,:);
                res_late=res_late(:);
                tr_res_late=[tr_res_late; res_late];

                inv_all=[an_inv pl_inv];
                inv_late=inv_all(1:late,:);
                inv_late=inv_late(:);
                sel_inv_late=[sel_inv_late; inv_late];

                inv_all=[an_tr_inv pl_tr_inv];
                inv_late=inv_all(1:late,:);
                inv_late=inv_late(:);
                tr_inv_late=[tr_inv_late; inv_late];
            end
        end
        if isfile(file_mutu) %take only those that have introduced species after empty niches, cross-guild case
            an_res=dlmread_empty(file_mutu,'\t',0,0,NaN);
            an_res=an_res(:,1:end-1);
            pl_res=dlmread_empty(strcat(sel_folder_mutu,'pl_sel_res_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'),'\t',0,0,NaN);
            pl_res=pl_res(:,1:end-1);

            an_inv=dlmread_empty(strcat(sel_folder_mutu,'an_sel_inv_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'),'\t',0,0,NaN);
            an_inv=an_inv(:,1:end-1);
            pl_inv=dlmread_empty(strcat(sel_folder_mutu,'pl_sel_inv_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'),'\t',0,0,NaN);
            pl_inv=pl_inv(:,1:end-1);

            an_tr=dlmread(strcat(tr_folder_mutu,'animal_traits_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
            an_tr_for_id=mat2cell(an_tr,ones(1,size(an_tr,1)),size(an_tr,2));
            an_tr_for_id=cellfun(@(x) x(x~=0),an_tr_for_id,'UniformOutput',0);
            pl_tr=dlmread(strcat(tr_folder_mutu,'plant_traits_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
            pl_tr_for_id=mat2cell(pl_tr,ones(1,size(pl_tr,1)),size(pl_tr,2));
            pl_tr_for_id=cellfun(@(x) x(x~=0),pl_tr_for_id,'UniformOutput',0);
            sep_index=dlmread(strcat(tr_folder_mutu,'sepIndex_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
            the_id=define_native_invader_function(an_tr_for_id,pl_tr_for_id,sep_index);
            an_id=the_id{1};
            pl_id=the_id{2};
            an_id=cell_to_mat_fill_nan(an_id);
            pl_id=cell_to_mat_fill_nan(pl_id);

            an_tr_res=nan(size(an_res));
            an_tr_inv=nan(size(an_inv));
            pl_tr_res=nan(size(pl_res));
            pl_tr_inv=nan(size(pl_inv));
            for i=1:size(an_tr,1)
                an_trait=an_tr(i,:);
                pl_trait=pl_tr(i,:);
                if sum(an_id(i,:)==0)>0  %if there is a resident species
                    an_tr_res(i,:)=[an_trait(an_id(i,:)==0) nan(1,sum(isnan(an_res(i,:))))];
                end
                if sum(an_id(i,:)==1)>0  %if there is an introduced species
                    an_tr_inv(i,:)=[an_trait(an_id(i,:)==1) nan(1,sum(isnan(an_inv(i,:))))] ;
                end
                if sum(pl_id(i,:)==0)>0
                    pl_tr_res(i,:)=[pl_trait(pl_id(i,:)==0) nan(1,sum(isnan(pl_res(i,:))))];
                end
                if sum(pl_id(i,:)==1)>0
                    pl_tr_inv(i,:)=[pl_trait(pl_id(i,:)==1) nan(1,sum(isnan(pl_inv(i,:))))];
                end
            end

            time=dlmread(strcat(tr_folder_mutu,'time_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));

            %strore the selection pressure values for intermediate
            %times
            early=find(time<=0.3*10^6);
            early=early(end);
            res_all=[an_res pl_res];
            res_early=res_all(1:early,:);
            res_early=res_early(:);
            sel_res_early=[sel_res_early; res_early];

            res_all=[an_tr_res pl_tr_res];
            res_early=res_all(1:early,:);
            res_early=res_early(:);
            tr_res_early=[tr_res_early; res_early];

            inv_all=[an_inv pl_inv];
            inv_early=inv_all(1:early,:);
            inv_early=inv_early(:);
            sel_inv_early=[sel_inv_early; inv_early];

            inv_all=[an_tr_inv pl_tr_inv];
            inv_early=inv_all(1:early,:);
            inv_early=inv_early(:);
            tr_inv_early=[tr_inv_early; inv_early];

            %strore the selection pressure values for intermediate
            %times
            int=find(time>0.3*10^6 & time<=2*10^6);
            if ~isempty(int)
                int=int(end);
                res_all=[an_res pl_res];
                res_int=res_all(1:int,:);
                res_int=res_int(:);
                sel_res_int=[sel_res_int; res_int];

                res_all=[an_tr_res pl_tr_res];
                res_int=res_all(1:int,:);
                res_int=res_int(:);
                tr_res_int=[tr_res_int; res_int];

                inv_all=[an_inv pl_inv];
                inv_int=inv_all(1:int,:);
                inv_int=inv_int(:);
                sel_inv_int=[sel_inv_int; inv_int];

                inv_all=[an_tr_inv pl_tr_inv];
                inv_int=inv_all(1:int,:);
                inv_int=inv_int(:);
                tr_inv_int=[tr_inv_int; inv_int];
            end
            %strore the selection pressure values for intermediate
            %times
            late=find(time>2*10^6);
            if ~isempty(late)
                late=late(end);
                res_all=[an_res pl_res];
                res_late=res_all(1:late,:);
                res_late=res_late(:);
                sel_res_late=[sel_res_late; res_late];

                res_all=[an_tr_res pl_tr_res];
                res_late=res_all(1:late,:);
                res_late=res_late(:);
                tr_res_late=[tr_res_late; res_late];

                inv_all=[an_inv pl_inv];
                inv_late=inv_all(1:late,:);
                inv_late=inv_late(:);
                sel_inv_late=[sel_inv_late; inv_late];

                inv_all=[an_tr_inv pl_tr_inv];
                inv_late=inv_all(1:late,:);
                inv_late=inv_late(:);
                tr_inv_late=[tr_inv_late; inv_late];
            end
       end
    
    end
end

all_sel=[sel_res_early; sel_res_int; sel_res_late; sel_inv_early; sel_inv_int; sel_inv_late];
all_sel(all_sel==0)=NaN;
all_sel=all_sel(~isnan(all_sel));
sel_min=log10(min(all_sel));
sel_max=log10(max(all_sel));

sel_res_early=log10(sel_res_early);
sel_res_int=log10(sel_res_int);
sel_res_late=log10(sel_res_late);
sel_inv_early=log10(sel_inv_early);
sel_inv_int=log10(sel_inv_int);
sel_inv_late=log10(sel_inv_late);

sel_div=(sel_max-sel_min)/500;
gr_sel=sel_min:sel_div:sel_max;

tr_div=(4-1)/500;
gr_tr=1:tr_div:4;

sel_res_early_ind=discretize(sel_res_early,gr_sel);
tr_res_early_ind=discretize(tr_res_early,gr_tr);
sel_res_int_ind=discretize(sel_res_int,gr_sel);
tr_res_int_ind=discretize(tr_res_int,gr_tr);
sel_res_late_ind=discretize(sel_res_late,gr_sel);
tr_res_late_ind=discretize(tr_res_late,gr_tr);

sel_inv_early_ind=discretize(sel_inv_early,gr_sel);
tr_inv_early_ind=discretize(tr_inv_early,gr_tr);
sel_inv_int_ind=discretize(sel_inv_int,gr_sel);
tr_inv_int_ind=discretize(tr_inv_int,gr_tr);
sel_inv_late_ind=discretize(sel_inv_late,gr_sel);
tr_inv_late_ind=discretize(tr_inv_late,gr_tr);

res_early_matrix=zeros(501,501);
inv_early_matrix=zeros(501,501);

res_int_matrix=zeros(501,501);
inv_int_matrix=zeros(501,501);

res_late_matrix=zeros(501,501);
inv_late_matrix=zeros(501,501);

disp ('Building all the matrices of frequency to be plotted');
for tr=1:500
    tr_res_early_intersect_index=find(tr_res_early_ind==tr);
    tr_res_int_intersect_index=find(tr_res_int_ind==tr);
    tr_res_late_intersect_index=find(tr_res_late_ind==tr);
    tr_inv_early_intersect_index=find(tr_inv_early_ind==tr);
    tr_inv_int_intersect_index=find(tr_inv_int_ind==tr);
    tr_inv_late_intersect_index=find(tr_inv_late_ind==tr);
    for j=1:500
        sel_res_early_ind_considered=sel_res_early_ind(tr_res_early_intersect_index);
        res_early_matrix(tr,j)=sum(sel_res_early_ind_considered==j);
        sel_res_int_ind_considered=sel_res_int_ind(tr_res_int_intersect_index);
        res_int_matrix(tr,j)=sum(sel_res_int_ind_considered==j);
        sel_res_late_ind_considered=sel_res_late_ind(tr_res_late_intersect_index);
        res_late_matrix(tr,j)=sum(sel_res_late_ind_considered==j);
        
        sel_inv_early_ind_considered=sel_inv_early_ind(tr_inv_early_intersect_index);
        inv_early_matrix(tr,j)=sum(sel_inv_early_ind_considered==j);
        sel_inv_int_ind_considered=sel_inv_int_ind(tr_inv_int_intersect_index);
        inv_int_matrix(tr,j)=sum(sel_inv_int_ind_considered==j);
        sel_inv_late_ind_considered=sel_inv_late_ind(tr_inv_late_intersect_index);
        inv_late_matrix(tr,j)=sum(sel_inv_late_ind_considered==j);

    end
end

%Build the frequency matrix of selection pressure
res_early_matrix(res_early_matrix==0)=NaN;
res_int_matrix(res_int_matrix==0)=NaN;
res_late_matrix(res_late_matrix==0)=NaN;
inv_early_matrix(inv_early_matrix==0)=NaN;
inv_int_matrix(inv_int_matrix==0)=NaN;
inv_late_matrix(inv_late_matrix==0)=NaN;

res_early_matrix(:,end)=NaN;
res_early_matrix(end,:)=NaN;
inv_early_matrix(:,end)=NaN;
inv_early_matrix(end,:)=NaN;

res_int_matrix(:,end)=NaN;
res_int_matrix(end,:)=NaN;
inv_int_matrix(:,end)=NaN;
inv_int_matrix(end,:)=NaN;

res_late_matrix(:,end)=NaN;
res_late_matrix(end,:)=NaN;
inv_late_matrix(:,end)=NaN;
inv_late_matrix(end,:)=NaN;

all_min=min([min(res_early_matrix(:)) min(inv_early_matrix(:)) min(res_int_matrix(:)) min(inv_int_matrix(:)) min(res_late_matrix(:)) min(inv_late_matrix(:))]);
all_max=max([max(res_early_matrix(:)) max(inv_early_matrix(:)) max(res_int_matrix(:)) max(inv_int_matrix(:)) max(res_late_matrix(:)) max(inv_late_matrix(:))]);

res_early_matrix(end,1)=all_min;
res_early_matrix(end,2)=all_max;
inv_early_matrix(end,1)=all_min;
inv_early_matrix(end,2)=all_max;

res_int_matrix(end,1)=all_min;
res_int_matrix(end,2)=all_max;
inv_int_matrix(end,1)=all_min;
inv_int_matrix(end,2)=all_max;

res_late_matrix(end,1)=all_min;
res_late_matrix(end,2)=all_max;
inv_late_matrix(end,1)=all_min;
inv_late_matrix(end,2)=all_max;

res_early_matrix=log10(res_early_matrix);
inv_early_matrix=log10(inv_early_matrix);
res_int_matrix=log10(res_int_matrix);
inv_int_matrix=log10(inv_int_matrix);
res_late_matrix=log10(res_late_matrix);
inv_late_matrix=log10(inv_late_matrix);

log_min=min([min(res_early_matrix(:)) min(inv_early_matrix(:)) min(res_int_matrix(:)) min(inv_int_matrix(:)) min(res_late_matrix(:)) min(inv_late_matrix(:))]);
log_max=max([max(res_early_matrix(:)) max(inv_early_matrix(:)) max(res_int_matrix(:)) max(inv_int_matrix(:)) max(res_late_matrix(:)) max(inv_late_matrix(:))]);

%Define all plots and axes properties
my_colormap=hot(60);
my_colormap=(flipud(my_colormap));
the_first=my_colormap(1,:);
my_colormap=my_colormap (15:end,:);
colormap(my_colormap);

relative_sel_pres=[10.^(-15:3:-3) 1];
sel_p=(relative_sel_pres*(10^sel_max-10^sel_min))+10^sel_min;
power_tick=log10(sel_p);
the_x_ticks=((power_tick-sel_min)*500)/(sel_max-sel_min);
the_x_ticks=[1 the_x_ticks];
the_x_ticks(end)=501;
the_x_ticks_label={'0','10^{-15}', '10^{-12}', '10^{-9}','10^{-6}','10^{-3}','1'};

the_y_ticks_label=1:0.5:4;
the_y_ticks=((the_y_ticks_label-1)*500)/(4-1);
the_y_ticks(1)=1;

figure('Visible','off');
colormap(my_colormap);
pp_inv=pcolor(inv_early_matrix);
set(pp_inv,'edgecolor','none');
set(gca,'TickDir','in','TickLength',[0.015 0.015],'XMinorTick','off','YMinorTick','off','layer','top','LineWidth',1.5);
set(gca,'Xtick',the_x_ticks,'XtickLabel',the_x_ticks_label,'Ytick',the_y_ticks,'YTickLabel',the_y_ticks_label,'FontSize',15);
pbaspect([1 1 1]);
box off;
xlab='Selection pressure';
ylab='Trait';
xlabel(xlab,'FontSize',15);
ylabel(ylab,'FontSize',15);

cc=colorbar();
set(cc,'Ylim',[log_min log_max]);
y_t=0:1:3;
y_t_lab=10.^(y_t);
set(cc,'YTick',y_t,'YtickLabel',y_t_lab,'FontSize',15,'TickLength',0.015,'LineWidth',1.5);
ycb=ylabel(cc,'Frequency','FontSize',15);
set(ycb,'Rotation',-90);
label_pos=cc.Label.Position;
label_pos(1)=label_pos(1)+0.5;
cc.Label.Position=label_pos;

addpath plotboxpos-pkg/plotboxpos
fig_pos=plotboxpos(gca);
cb_pos=get(cc,'Position');

set(cc,'Position',[fig_pos(1)+fig_pos(3)+0.05 fig_pos(2) cb_pos(3) cb_pos(4)]);
set(gca,'Position',[fig_pos(1) fig_pos(2) fig_pos(3) cb_pos(4)]);

pos1=[fig_pos(1) fig_pos(2)+0.001 fig_pos(3) cb_pos(4)];
pos2=pos1;
pos2(1)=0.5;
pos2(3)=pos1(3)-(pos2(1)-pos1(1));
ax2 = axes('Position',pos2);
               
dis=nansum(inv_early_matrix,2);
dis=dis(1:end-1);
new_dis=zeros(100,1);
count=1;
for i=1:5:length(dis)
    new_dis(count)=sum(dis(i:1:i+4))/5;
    count=count+1;
end    
new_dis=new_dis/max(new_dis);
plot(new_dis,1:5:500,'-','Color',rgb('ForestGreen'),'LineWidth',1,'Parent',ax2);
set(ax2,'Ytick',the_y_ticks,'YTickLabel',[],'Color','none','XAxisLocation','top','YAxisLocation','right','XColor',rgb('ForestGreen'),'YColor','k');
set(ax2,'TickDir','in','TickLength',[0.015 0.015],'XMinorTick','off','YMinorTick','off','layer','top','LineWidth',1.5,'Xdir','reverse','FontSize',15);
box off;
mkdir('./Supplementary_figures/');
name=strcat('./Supplementary_figures/FigB8a.png');
print(name,'-dpng','-r600');


figure('Visible','off');
colormap(my_colormap);
pp_inv=pcolor(inv_int_matrix);
set(pp_inv,'edgecolor','none');
set(gca,'TickDir','in','TickLength',[0.015 0.015],'XMinorTick','off','YMinorTick','off','layer','top','LineWidth',1.5);
set(gca,'Xtick',the_x_ticks,'XtickLabel',the_x_ticks_label,'Ytick',the_y_ticks,'YTickLabel',the_y_ticks_label,'FontSize',15);
pbaspect([1 1 1]);
box off;
xlab='Selection pressure';
ylab='Trait';
xlabel(xlab,'FontSize',15);
ylabel(ylab,'FontSize',15);

cc=colorbar();
set(cc,'Ylim',[log_min log_max]);
y_t=0:1:3;
y_t_lab=10.^(y_t);
set(cc,'YTick',y_t,'YtickLabel',y_t_lab,'FontSize',15,'TickLength',0.015,'LineWidth',1.5);
ycb=ylabel(cc,'Frequency','FontSize',15);
set(ycb,'Rotation',-90);
label_pos=cc.Label.Position;
label_pos(1)=label_pos(1)+0.5;
cc.Label.Position=label_pos;

addpath plotboxpos-pkg/plotboxpos
fig_pos=plotboxpos(gca);
cb_pos=get(cc,'Position');

set(cc,'Position',[fig_pos(1)+fig_pos(3)+0.05 fig_pos(2) cb_pos(3) cb_pos(4)]);
set(gca,'Position',[fig_pos(1) fig_pos(2) fig_pos(3) cb_pos(4)]);

pos1=[fig_pos(1) fig_pos(2)+0.001 fig_pos(3) cb_pos(4)];
pos2=pos1;
pos2(1)=0.5;
pos2(3)=pos1(3)-(pos2(1)-pos1(1));
ax2 = axes('Position',pos2);
               
dis=nansum(inv_int_matrix,2);
dis=dis(1:end-1);
new_dis=zeros(100,1);
count=1;
for i=1:5:length(dis)
    new_dis(count)=sum(dis(i:1:i+4))/5;
    count=count+1;
end    
new_dis=new_dis/max(new_dis);
plot(new_dis,1:5:500,'-','Color',rgb('ForestGreen'),'LineWidth',1,'Parent',ax2);
set(ax2,'Ytick',the_y_ticks,'YTickLabel',[],'Color','none','XAxisLocation','top','YAxisLocation','right','XColor',rgb('ForestGreen'),'YColor','k');
set(ax2,'TickDir','in','TickLength',[0.015 0.015],'XMinorTick','off','YMinorTick','off','layer','top','LineWidth',1.5,'Xdir','reverse','FontSize',15);
box off;
mkdir('./Supplementary_figures/');
name=strcat('./Supplementary_figures/FigB8b.png');
print(name,'-dpng','-r600');

figure('Visible','off');
colormap(my_colormap);
pp_inv=pcolor(inv_late_matrix);
set(pp_inv,'edgecolor','none');
set(gca,'TickDir','in','TickLength',[0.015 0.015],'XMinorTick','off','YMinorTick','off','layer','top','LineWidth',1.5);
set(gca,'Xtick',the_x_ticks,'XtickLabel',the_x_ticks_label,'Ytick',the_y_ticks,'YTickLabel',the_y_ticks_label,'FontSize',15);
pbaspect([1 1 1]);
box off;
xlab='Selection pressure';
ylab='Trait';
xlabel(xlab,'FontSize',15);
ylabel(ylab,'FontSize',15);


cc=colorbar();
set(cc,'Ylim',[log_min log_max]);
y_t=0:1:3;
y_t_lab=10.^(y_t);
set(cc,'YTick',y_t,'YtickLabel',y_t_lab,'FontSize',15,'TickLength',0.015,'LineWidth',1.5);
ycb=ylabel(cc,'Frequency','FontSize',15);
set(ycb,'Rotation',-90);
label_pos=cc.Label.Position;
label_pos(1)=label_pos(1)+0.5;
cc.Label.Position=label_pos;

addpath plotboxpos-pkg/plotboxpos
fig_pos=plotboxpos(gca);
cb_pos=get(cc,'Position');

set(cc,'Position',[fig_pos(1)+fig_pos(3)+0.05 fig_pos(2) cb_pos(3) cb_pos(4)]);
set(gca,'Position',[fig_pos(1) fig_pos(2) fig_pos(3) cb_pos(4)]);

pos1=[fig_pos(1) fig_pos(2)+0.001 fig_pos(3) cb_pos(4)];
pos2=pos1;
pos2(1)=0.5;
pos2(3)=pos1(3)-(pos2(1)-pos1(1));
ax2 = axes('Position',pos2);
               
dis=nansum(inv_late_matrix,2);
dis=dis(1:end-1);
new_dis=zeros(100,1);
count=1;
for i=1:5:length(dis)
    new_dis(count)=sum(dis(i:1:i+4))/5;
    count=count+1;
end    
new_dis=new_dis/max(new_dis);
plot(new_dis,1:5:500,'-','Color',rgb('ForestGreen'),'LineWidth',1,'Parent',ax2);
set(ax2,'Ytick',the_y_ticks,'YTickLabel',[],'Color','none','XAxisLocation','top','YAxisLocation','right','XColor',rgb('ForestGreen'),'YColor','k');
set(ax2,'TickDir','in','TickLength',[0.015 0.015],'XMinorTick','off','YMinorTick','off','layer','top','LineWidth',1.5,'Xdir','reverse','FontSize',15);
box off;
mkdir('./Supplementary_figures/');
name=strcat('./Supplementary_figures/FigB8c.png');
print(name,'-dpng','-r600');
