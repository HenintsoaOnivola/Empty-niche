%This code simulates the eco-evolutionary dynamics of a
%community in which species are engaged in inter and intra-specific
%competition for resources, and in cross-guild mutualism

clear all;

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

mutant_fraction=0.1;
tolerance_pop_trait=10^(-5);

log_sig_c_an_range=-2.5:0.0625:1;  %this is the range for the width of competition kernel
%log_sig_m_range=-3:0.0625:-0.5;  %this is the range for the width of mutualistic kernel

log_sig_c_an_range=log_sig_c_an_range(16);  
log_sig_m_range=-3;

log_sig_an=1;
sig_an=1;
sig_pl=sig_an;  %this is the width of resource kernel

%%looping over the ranges of competition and mutualistic kernels
for log_sig_c_an=log_sig_c_an_range
    sig_c_an=exp(log_sig_c_an);
    sig_c_pl=sig_c_an;
    for log_sig_m=log_sig_m_range
        sig_m=exp(log_sig_m);
        disp ([sig_an sig_c_an sig_m]);
        
        %initialization
        %%%%%%%%%%%%%%%%
       
       %initial population and trait values when the system optimizes for intra-guild resource intake 
       ini_val=[1 1 3.2 1.8];
       res_file='Results_intra_guild';  %change to'Results_cross-guild'for the second case initial trait values

         
       %initial population and trait values when the system optimizes for cross-guild mutualistic benefits
       %ini_val=[1 1 2.5 2.7];  %to be uncommented if one wants the initial conditions optimizing cross-guild mutualitic benefits
       %res_file='Results_cross_guild'; %to be uncommented if one wants the initial conditions optimizing cross-guild mutualitic benefits

        animal=1;   %initial number of species
        plant=1;

        if log_sig_c_an<=-2.5 
            infinite=1000;
        else
            infinite=200;
        end
        new_sys=true;
        count_into_while=0;
        previous_t_end=0;

        time_plot=zeros(7e005,1); 
        time_plot(:)=NaN;
        pop_plot_an=cell(7e005,1);
        pop_plot_pl=cell(7e005,1);
        trait_plot_an=cell(7e005,1);
        trait_plot_pl=cell(7e005,1);
        time_index=0;
        not_integrated=false;


        separating_index=[];
        %looping over new system
        an_num=[animal];
        pl_num=[plant];
        end_each_sys=[];
        total_time=0;
        while (new_sys==true && count_into_while<15) 
            count_into_while=count_into_while+1;
            disp (['count system: ' num2str(count_into_while)]);
            fprintf ('\n');

            convergence=false;
            count_integrate=0;

            total_number=1;
            new_sys=false;

            %looping over integration until we have convergence
            temporary_time=0;
            while ( (convergence==false) && count_integrate<=infinite )
                tic;
                [time,sol,TE,YE,IE]=solve_eco_evo_subsequent_branching(ini_val,rA,rP,h,const,sig_c_an,sig_c_pl,sig_an,sig_pl,x0,y0,kx,ky,sig_m,animal,plant,temporary_time,count_integrate,count_into_while);
                fprintf ('%d\t',count_integrate);
                temporary_time=temporary_time+time(end);

                if count_integrate>100
                    oscillation=detect_oscillation_using_trait_variation(sol,animal,plant);
                    if oscillation==true
                        an_pop=sol(end,1:animal);
                        pl_pop=sol(end,animal+1:animal+plant);
                        an_trait=sol(end,animal+plant+1:2*animal+plant);
                        pl_trait=sol(end,2*animal+plant+1:2*animal+2*plant);
                        animal=length(an_pop);
                        plant=length(pl_pop);

                        current_time_range=time+previous_t_end;
                        time_plot(time_index+1:time_index+numel(time))=current_time_range;

                        the_pop_an=sol(:,1:animal);
                        solution_pop_an=mat2cell(the_pop_an,ones(1,size(the_pop_an,1)),size(the_pop_an,2));
                        pop_plot_an(time_index+1:time_index+numel(time))=solution_pop_an;
                        the_pop_pl=sol(:,animal+1:animal+plant);
                        solution_pop_pl=mat2cell(the_pop_pl,ones(1,size(the_pop_pl,1)),size(the_pop_pl,2));
                        pop_plot_pl(time_index+1:time_index+numel(time))=solution_pop_pl;
                        the_trait_an=sol(:,animal+plant+1:animal+plant+animal);
                        solution_trait_an=mat2cell(the_trait_an,ones(1,size(the_trait_an,1)),size(the_trait_an,2));
                        trait_plot_an(time_index+1:time_index+numel(time))=solution_trait_an;
                        the_trait_pl=sol(:,animal+plant+animal+1:2*(animal+plant));
                        solution_trait_pl=mat2cell(the_trait_pl,ones(1,size(the_trait_pl,1)),size(the_trait_pl,2));
                        trait_plot_pl(time_index+1:time_index+numel(time))=solution_trait_pl;

                        separating_index=[separating_index time_index+numel(time)];

                        previous_t_end=time_plot(time_index+numel(time));
                        count_integrate=count_integrate+1;
                        time_index=time_index+numel(time);

                        %reinitialization
                        animal=length(an_pop);
                        plant=length(pl_pop);
                        ini_val=[an_pop pl_pop an_trait pl_trait];
                        new_sys=false;
                        disp ('oscillation');
                        separating_index=[separating_index 0];
                        break;
                    end
                end

                if count_integrate==infinite
                    %consider the case like a convergence
                    IE=2;
                end

                if isempty (IE)

                    %storing the traits and populations
                    an_pop=sol(end,1:animal);
                    pl_pop=sol(end,animal+1:animal+plant);
                    an_trait=sol(end,animal+plant+1:2*animal+plant);
                    pl_trait=sol(end,2*animal+plant+1:2*animal+2*plant);
                    animal=length(an_pop);
                    plant=length(pl_pop);

                    current_time_range=time+previous_t_end;
                    time_plot(time_index+1:time_index+numel(time))=current_time_range;

                    the_pop_an=sol(:,1:animal);
                    solution_pop_an=mat2cell(the_pop_an,ones(1,size(the_pop_an,1)),size(the_pop_an,2));
                    pop_plot_an(time_index+1:time_index+numel(time))=solution_pop_an;
                    the_pop_pl=sol(:,animal+1:animal+plant);
                    solution_pop_pl=mat2cell(the_pop_pl,ones(1,size(the_pop_pl,1)),size(the_pop_pl,2));
                    pop_plot_pl(time_index+1:time_index+numel(time))=solution_pop_pl;
                    the_trait_an=sol(:,animal+plant+1:animal+plant+animal);
                    solution_trait_an=mat2cell(the_trait_an,ones(1,size(the_trait_an,1)),size(the_trait_an,2));
                    trait_plot_an(time_index+1:time_index+numel(time))=solution_trait_an;
                    the_trait_pl=sol(:,animal+plant+animal+1:2*(animal+plant));
                    solution_trait_pl=mat2cell(the_trait_pl,ones(1,size(the_trait_pl,1)),size(the_trait_pl,2));
                    trait_plot_pl(time_index+1:time_index+numel(time))=solution_trait_pl;


                    previous_t_end=time_plot(time_index+numel(time));
                    count_integrate=count_integrate+1;
                    time_index=time_index+numel(time);

                    %reinitialization
                    animal=length(an_pop);
                    plant=length(pl_pop);
                    ini_val=[an_pop pl_pop an_trait pl_trait];

                else

                    if any(IE==1)
                        disp ('there is extinction');

                        an_pop_full=sol(end,1:animal);
                        an_pop=an_pop_full(an_pop_full>1e-8);
                        pl_pop_full=sol(end,animal+1:animal+plant);
                        pl_pop=pl_pop_full(pl_pop_full>1e-8);
                        an_trait=sol(end,animal+plant+1:2*animal+plant);
                        an_trait=an_trait(an_pop_full>1e-8);
                        pl_trait=sol(end,2*animal+plant+1:2*animal+2*plant);
                        pl_trait=pl_trait(pl_pop_full>1e-8);

                        %storing the traits and populations
                        current_time_range=time+previous_t_end;
                        time_plot(time_index+1:time_index+numel(time))=current_time_range;

                        the_pop_an=sol(:,1:animal);
                        solution_pop_an=mat2cell(the_pop_an,ones(1,size(the_pop_an,1)),size(the_pop_an,2));
                        pop_plot_an(time_index+1:time_index+numel(time))=solution_pop_an;
                        the_pop_pl=sol(:,animal+1:animal+plant);
                        solution_pop_pl=mat2cell(the_pop_pl,ones(1,size(the_pop_pl,1)),size(the_pop_pl,2));
                        pop_plot_pl(time_index+1:time_index+numel(time))=solution_pop_pl;
                        the_trait_an=sol(:,animal+plant+1:animal+plant+animal);
                        solution_trait_an=mat2cell(the_trait_an,ones(1,size(the_trait_an,1)),size(the_trait_an,2));
                        trait_plot_an(time_index+1:time_index+numel(time))=solution_trait_an;
                        the_trait_pl=sol(:,animal+plant+animal+1:2*(animal+plant));
                        solution_trait_pl=mat2cell(the_trait_pl,ones(1,size(the_trait_pl,1)),size(the_trait_pl,2));
                        trait_plot_pl(time_index+1:time_index+numel(time))=solution_trait_pl;
                        separating_index=[separating_index time_index+numel(time)];


                        previous_t_end=time_plot(time_index+numel(time));
                        count_integrate=count_integrate+1;
                        time_index=time_index+numel(time);


                        if isempty(an_pop)==1 || isempty(pl_pop)==1
                            disp ('one entire side is extinct');
                            separating_index=[separating_index -1];
                            new_sys=false;
                            end_each_sys=[end_each_sys time_index];
                            break;
                        end

                        %reinitialization
                        animal=length(an_pop);
                        plant=length(pl_pop);
                        ini_val=[an_pop pl_pop an_trait pl_trait];

                    else

                        if any (IE==2) 

                            convergence=true;
                            %storing the traits and populations

                            an_pop=sol(end,1:animal);
                            pl_pop=sol(end,animal+1:animal+plant);
                            an_trait=sol(end,animal+plant+1:2*animal+plant);
                            pl_trait=sol(end,2*animal+plant+1:2*animal+2*plant);
                            animal=length(an_pop);
                            plant=length(pl_pop);

                            current_time_range=time+previous_t_end;
                            time_plot(time_index+1:time_index+numel(time))=current_time_range;

                            the_pop_an=sol(:,1:animal);
                            solution_pop_an=mat2cell(the_pop_an,ones(1,size(the_pop_an,1)),size(the_pop_an,2));
                            pop_plot_an(time_index+1:time_index+numel(time))=solution_pop_an;
                            the_pop_pl=sol(:,animal+1:animal+plant);
                            solution_pop_pl=mat2cell(the_pop_pl,ones(1,size(the_pop_pl,1)),size(the_pop_pl,2));
                            pop_plot_pl(time_index+1:time_index+numel(time))=solution_pop_pl;
                            the_trait_an=sol(:,animal+plant+1:animal+plant+animal);
                            solution_trait_an=mat2cell(the_trait_an,ones(1,size(the_trait_an,1)),size(the_trait_an,2));
                            trait_plot_an(time_index+1:time_index+numel(time))=solution_trait_an;
                            the_trait_pl=sol(:,animal+plant+animal+1:2*(animal+plant));
                            solution_trait_pl=mat2cell(the_trait_pl,ones(1,size(the_trait_pl,1)),size(the_trait_pl,2));
                            trait_plot_pl(time_index+1:time_index+numel(time))=solution_trait_pl;

                            separating_index=[separating_index time_index+numel(time)];


                            previous_t_end=time_plot(time_index+numel(time));
                            count_integrate=count_integrate+1;
                            time_index=time_index+numel(time);

                            end_each_sys=[end_each_sys time_index];

                            fprintf ('\n');

                            %merge very similar traits if needed
                            an_before=length(an_trait);
                            pl_before=length(pl_trait);

                            previous_pl_tr=pl_trait;
                            previous_pl_pop=pl_pop;
                            previous_an_tr=an_trait;
                            previous_an_pop=an_pop;

                            an_merged=merge(an_trait,an_pop);
                            an_trait=an_merged{1};
                            an_pop=an_merged{2};

                            pl_merged=merge(pl_trait,pl_pop);
                            pl_trait=pl_merged{1};
                            pl_pop=pl_merged{2};

                            an_after=length(an_trait);
                            pl_after=length(pl_trait);

                            if an_before~=an_after || pl_before~=pl_after %if there were merging, then store the last merged traits

                                disp ('merging');
                                time_plot(time_index+1)=time_plot(time_index);

                                pop_plot_an{time_index+1}=an_pop;
                                pop_plot_pl{time_index+1}=pl_pop;
                                trait_plot_an{time_index+1}=an_trait;
                                trait_plot_pl{time_index+1}=pl_trait;
                                separating_index=[separating_index time_index+1];

                                previous_t_end=time_plot(time_index+1);
                                time_index=time_index+1;
                                end_each_sys(end)=time_index;

                            end


                            %reinitialization
                            animal=length(an_pop);
                            plant=length(pl_pop);
                            ini_val=[an_pop pl_pop an_trait pl_trait];

                            an_num=[an_num animal];
                            pl_num=[pl_num plant];

                            if count_into_while>=2

                                if an_num(end-1)==animal && pl_num(end-1)==plant

                                    an_tr_pr=trait_plot_an{end_each_sys(end-1)};
                                    diff_tr_an=abs(an_tr_pr-an_trait);
                                    
                                    pl_tr_pr=trait_plot_pl{end_each_sys(end-1)};
                                    diff_tr_pl=abs(pl_tr_pr-pl_trait);

                                    if sum(diff_tr_an<10^(-5))==length(diff_tr_an)  && sum(diff_tr_pl<10^(-5))==length(diff_tr_pl) 
                                        new_sys=false;
                                        disp ('no more diversification');
                                        %remove the previous stored
                                        %values
                                        the_prev_end=end_each_sys(end-1);
                                        pop_plot_an(the_prev_end+1:end)={[]};
                                        pop_plot_pl(the_prev_end+1:end)={[]};
                                        trait_plot_an(the_prev_end+1:end)={[]};
                                        trait_plot_pl(the_prev_end+1:end)={[]};
                                        time_plot(the_prev_end+1:end)=NaN;

                                        separating_index(separating_index>the_prev_end)=[];
                                        end_each_sys=end_each_sys(1:end-1);

                                        break;
                                    end
                                end
                            end


                            %check for branching events
                            sel_grad=selection_gradient_function([an_pop pl_pop an_trait pl_trait],rA,rP,h,const,sig_c_an,sig_c_pl,sig_an,sig_pl,x0,y0,kx,ky,sig_m,animal,plant);
                            sel_grad_an=sel_grad{1};
                            sel_grad_pl=sel_grad{2};

                            info=branching_condition_analytical([an_pop pl_pop an_trait pl_trait],rA,rP,h,const,sig_c_an,sig_c_pl,sig_an,sig_pl,x0,y0,kx,ky,sig_m,animal,plant);
                            info_an=info{1};
                            info_pl=info{2};

                            disr_an=info_an(:,1);
                            disr_pl=info_pl(:,1);
                            disr=[disr_an;disr_pl];


                            new_an_pop=zeros(1,2*length(an_pop));
                            new_pl_pop=zeros(1,2*length(pl_pop));
                            new_an_tr=zeros(1,2*length(an_trait));
                            new_pl_tr=zeros(1,2*length(pl_trait));
                            new_an_pop(:)=NaN;
                            new_pl_pop(:)=NaN;
                            new_an_tr(:)=NaN;
                            new_pl_tr(:)=NaN;

                            new_an_pop(1:length(an_pop))=an_pop;
                            new_pl_pop(1:length(pl_pop))=pl_pop;
                            new_an_tr(1:length(an_trait))=an_trait;
                            new_pl_tr(1:length(pl_trait))=pl_trait;
                            count_new_mutant=0;
                            an_br_index=cell(2*animal,1);
                            count_br_index=1;
                            for i=1:length(an_pop)
                                if  disr_an(i)>0
                                    disp ('an is br');
                                    new_sys=true;
                                    count_new_mutant=count_new_mutant+1;
                                    new_an_pop(i)=an_pop(i)*(1-mutant_fraction);
                                    new_an_pop(length(an_pop)+count_new_mutant)=an_pop(i)*(mutant_fraction);
                                    new_an_tr(i)=an_trait(i);
                                    if sel_grad_an(i)<0
                                        new_an_tr(length(an_pop)+count_new_mutant)=an_trait(i)+0.005;
                                    else
                                        new_an_tr(length(an_pop)+count_new_mutant)=an_trait(i)-0.005;
                                    end
                                    an_br_index{count_br_index}=[i length(an_pop)+count_new_mutant];
                                    count_br_index=count_br_index+1;
                                end
         
                            end
                            an_br_index=an_br_index(~cellfun('isempty',an_br_index));
                            count_new_mutant=0;
                            pl_br_index=cell(2*plant,1);
                            count_br_index=1;
                            for i=1:length(pl_pop)
                                if disr_pl(i)>0
                                    disp ('pl is br');
                                    new_sys=true;
                                    count_new_mutant=count_new_mutant+1;
                                    new_pl_pop(i)=pl_pop(i)*(1-mutant_fraction);
                                    new_pl_pop(length(pl_pop)+count_new_mutant)=pl_pop(i)*(mutant_fraction);
                                    new_pl_tr(i)=pl_trait(i);
                                    if sel_grad_pl(i)<0
                                        new_pl_tr(length(pl_pop)+count_new_mutant)=pl_trait(i)+0.005;
                                    else
                                        new_pl_tr(length(pl_pop)+count_new_mutant)=pl_trait(i)-0.005;
                                    end
                                    pl_br_index{count_br_index}=[i length(pl_pop)+count_new_mutant];
                                    count_br_index=count_br_index+1;
                                 end

                            end
                            pl_br_index=pl_br_index(~cellfun('isempty',pl_br_index));

                            new_an_pop=new_an_pop(~isnan(new_an_pop));
                            new_pl_pop=new_pl_pop(~isnan(new_pl_pop));
                            new_an_tr=new_an_tr(~isnan(new_an_tr));
                            new_pl_tr=new_pl_tr(~isnan(new_pl_tr));

                            %re- reinitialization
                            animal=length(new_an_pop);
                            plant=length(new_pl_pop);
                            ini_val=[new_an_pop new_pl_pop new_an_tr new_pl_tr];

                            %check that they are real branching
                            initial_condition=ini_val;

                            continue_check_br=1;
                            temp_total_time=0;
                            count_temporary=0;
                            while(temp_total_time<=10^(6))
                                fprintf ('%d\t',count_temporary);
                                tic;
                                [temp_time,temp_sol,temp_TE, temp_YE, temp_IE]=solve_eco_evo_subsequent_branching(initial_condition,rA,rP,h,const,sig_c_an,sig_c_pl,sig_an,sig_pl,x0,y0,kx,ky,sig_m,animal,plant,temporary_time,count_integrate,count_temporary);
                                temp_total_time=temp_total_time+temp_time(end);
                                count_temporary=count_temporary+1;
                                initial_condition=temp_sol(end,:);
                                if (~isempty(temp_IE) && any (temp_IE==3))  || count_temporary>20
                                    ind_oscillation=detect_oscillation_using_trait_variation(temp_sol,animal,plant);
                                    if ind_oscillation==true
                                        continue_check_br=0;
                                        break;
                                    end
                                end
                                if ~isempty(temp_IE) && any(temp_IE==1)
                                    continue_check_br=0;
                                    break;
                                end


                            end
                            fprintf('\n');

                            if continue_check_br==1
                                temp_an_trait=temp_sol(end,animal+plant+1:2*animal+plant);
                                temp_pl_trait=temp_sol(end,2*animal+plant+1:2*animal+2*plant);
                                if length(an_trait)>1
                                    min_tr_diff=abs(an_trait(2)-an_trait(1));
                                    for tr=3:length(an_trait)
                                       tr_diff=abs(an_trait(tr)-an_trait(tr-1));
                                       if tr_diff<=min_tr_diff
                                           min_tr_diff=tr_diff;
                                       end
                                    end
                                else
                                    min_tr_diff=0.05;
                                end
                                for ind_an=1:length(an_br_index)
                                    the_index_res=an_br_index{ind_an}(1);
                                    the_index_mut=an_br_index{ind_an}(2);

                                    if abs(temp_an_trait(the_index_res)-temp_an_trait(the_index_mut))<=min_tr_diff/2
                                        disp ('an br not confirmed');
                                        new_an_tr(the_index_mut)=NaN;
                                        new_an_pop(the_index_mut)=NaN;
                                        new_an_pop(the_index_res)=an_pop(the_index_res);

                                        ii=the_index_res;
                                        disr_an_comp=info_an(ii,2);
                                        disr_an_mutu=info_an(ii,3);
                                        if sign(disr_an_comp)~=sign(disr_an_mutu)
                                            new_dsr_an_comp=sign(disr_an_comp)*(abs(disr_an_comp-disr_an_mutu)/2);
                                            new_dsr_an_mutu=sign(disr_an_mutu)*(abs(disr_an_comp-disr_an_mutu)/2);
                                        else
                                            new_dsr_an_comp=0;
                                            new_dsr_an_mutu=0;
                                        end
                                        info_an(ii,2)=new_dsr_an_comp;
                                        info_an(ii,3)=new_dsr_an_mutu;
                                        info_an(ii,1)=0;
                                    else
                                        disp ('an branching confirmed');
                                    end

                                end

                                if length(pl_trait)>1
                                    min_tr_diff=abs(pl_trait(2)-pl_trait(1));
                                    for tr=3:length(pl_trait)
                                       tr_diff=abs(pl_trait(tr)-pl_trait(tr-1));
                                       if tr_diff<=min_tr_diff
                                           min_tr_diff=tr_diff;
                                       end
                                    end
                                else
                                    min_tr_diff=0.05;
                                end
                                for ind_pl=1:length(pl_br_index)
                                    the_index_res=pl_br_index{ind_pl}(1);
                                    the_index_mut=pl_br_index{ind_pl}(2);


                                    if abs(temp_pl_trait(the_index_res)-temp_pl_trait(the_index_mut))<=min_tr_diff/2
                                        disp ('pl br not confirmed');
                                        new_pl_tr(the_index_mut)=NaN;
                                        new_pl_pop(the_index_mut)=NaN;
                                        new_pl_pop(the_index_res)=pl_pop(the_index_res);
                                        ii=the_index_res;
                                        disr_pl_comp=info_pl(ii,2);
                                        disr_pl_mutu=info_pl(ii,3);
                                        if sign(disr_pl_comp)~=sign(disr_pl_mutu)
                                            new_dsr_pl_comp=sign(disr_pl_comp)*(abs(disr_pl_comp-disr_pl_mutu)/2);
                                            new_dsr_pl_mutu=sign(disr_pl_mutu)*(abs(disr_pl_comp-disr_pl_mutu)/2);
                                        else
                                            new_dsr_pl_comp=0;
                                            new_dsr_pl_mutu=0;
                                        end
                                        info_pl(ii,2)=new_dsr_pl_comp;
                                        info_pl(ii,3)=new_dsr_pl_mutu;
                                        info_pl(ii,1)=0;
                                    else
                                        disp ('pl branching confirmed');
                                    end
                                end

                                new_an_pop=new_an_pop(~isnan(new_an_pop));
                                new_pl_pop=new_pl_pop(~isnan(new_pl_pop));
                                new_an_tr=new_an_tr(~isnan(new_an_tr));
                                new_pl_tr=new_pl_tr(~isnan(new_pl_tr));

                                if length(new_an_pop)==length(an_pop) && length(new_pl_pop)==length(pl_pop)
                                    
                                    disp ('ESS');
                                    new_sys=false;


                                end
                                ini_val=[new_an_pop new_pl_pop new_an_tr new_pl_tr];
                                animal=length(new_an_pop);
                                plant=length(new_pl_pop);
                            end
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


                            if new_sys==false
                                disp ('ESS');
                            end


                        end

                        if any (IE==3)  
                            %detect oscillation
                            oscillation=detect_oscillation_using_trait_variation(sol,animal,plant);
                            if oscillation==true
                                an_pop=sol(end,1:animal);
                                pl_pop=sol(end,animal+1:animal+plant);
                                an_trait=sol(end,animal+plant+1:2*animal+plant);
                                pl_trait=sol(end,2*animal+plant+1:2*animal+2*plant);
                                animal=length(an_pop);
                                plant=length(pl_pop);

                                current_time_range=time+previous_t_end;
                                time_plot(time_index+1:time_index+numel(time))=current_time_range;

                                the_pop_an=sol(:,1:animal);
                                solution_pop_an=mat2cell(the_pop_an,ones(1,size(the_pop_an,1)),size(the_pop_an,2));
                                pop_plot_an(time_index+1:time_index+numel(time))=solution_pop_an;
                                the_pop_pl=sol(:,animal+1:animal+plant);
                                solution_pop_pl=mat2cell(the_pop_pl,ones(1,size(the_pop_pl,1)),size(the_pop_pl,2));
                                pop_plot_pl(time_index+1:time_index+numel(time))=solution_pop_pl;
                                the_trait_an=sol(:,animal+plant+1:animal+plant+animal);
                                solution_trait_an=mat2cell(the_trait_an,ones(1,size(the_trait_an,1)),size(the_trait_an,2));
                                trait_plot_an(time_index+1:time_index+numel(time))=solution_trait_an;
                                the_trait_pl=sol(:,animal+plant+animal+1:2*(animal+plant));
                                solution_trait_pl=mat2cell(the_trait_pl,ones(1,size(the_trait_pl,1)),size(the_trait_pl,2));
                                trait_plot_pl(time_index+1:time_index+numel(time))=solution_trait_pl;

                                separating_index=[separating_index time_index+numel(time)];

                                previous_t_end=time_plot(time_index+numel(time));
                                count_integrate=count_integrate+1;
                                time_index=time_index+numel(time);

                                %reinitialization
                                animal=length(an_pop);
                                plant=length(pl_pop);
                                ini_val=[an_pop pl_pop an_trait pl_trait];
                                new_sys=false;
                                disp ('oscillation');
                                separating_index=[separating_index 0];
                                break;
                            else
                                %just update the system
                                an_pop=sol(end,1:animal);
                                pl_pop=sol(end,animal+1:animal+plant);
                                an_trait=sol(end,animal+plant+1:2*animal+plant);
                                pl_trait=sol(end,2*animal+plant+1:2*animal+2*plant);
                                animal=length(an_pop);
                                plant=length(pl_pop);

                                current_time_range=time+previous_t_end;
                                time_plot(time_index+1:time_index+numel(time))=current_time_range;

                                the_pop_an=sol(:,1:animal);
                                solution_pop_an=mat2cell(the_pop_an,ones(1,size(the_pop_an,1)),size(the_pop_an,2));
                                pop_plot_an(time_index+1:time_index+numel(time))=solution_pop_an;
                                the_pop_pl=sol(:,animal+1:animal+plant);
                                solution_pop_pl=mat2cell(the_pop_pl,ones(1,size(the_pop_pl,1)),size(the_pop_pl,2));
                                pop_plot_pl(time_index+1:time_index+numel(time))=solution_pop_pl;
                                the_trait_an=sol(:,animal+plant+1:animal+plant+animal);
                                solution_trait_an=mat2cell(the_trait_an,ones(1,size(the_trait_an,1)),size(the_trait_an,2));
                                trait_plot_an(time_index+1:time_index+numel(time))=solution_trait_an;
                                the_trait_pl=sol(:,animal+plant+animal+1:2*(animal+plant));
                                solution_trait_pl=mat2cell(the_trait_pl,ones(1,size(the_trait_pl,1)),size(the_trait_pl,2));
                                trait_plot_pl(time_index+1:time_index+numel(time))=solution_trait_pl;


                                previous_t_end=time_plot(time_index+numel(time));
                                count_integrate=count_integrate+1;
                                time_index=time_index+numel(time);

                                %reinitialization
                                animal=length(an_pop);
                                plant=length(pl_pop);
                                ini_val=[an_pop pl_pop an_trait pl_trait];
                            end
                        end
                    end




                end 


            end
            very_last_time=time_plot(~isnan(time_plot));
            total_time=very_last_time(end);

        end


        time_plot=time_plot(~isnan(time_plot));
        pop_plot_an=pop_plot_an(1:numel(time_plot),:);
        trait_plot_an=trait_plot_an(1:numel(time_plot),:);
        pop_plot_pl=pop_plot_pl(1:numel(time_plot),:);
        trait_plot_pl=trait_plot_pl(1:numel(time_plot),:);
        if length(separating_index)>1
            if separating_index(end)==(separating_index(end-1)+1) %if there were merging, cut the merged trait at the end
                time_plot=time_plot(1:end-1);
                trait_plot_an=trait_plot_an(1:end-1);
                trait_plot_pl=trait_plot_pl(1:end-1);
                pop_plot_an=pop_plot_an(1:end-1);
                pop_plot_pl=pop_plot_pl(1:end-1);
                separating_index=separating_index(1:end-1);
            end
        end

        plot_traits_function_with_fitness(separating_index,time_plot,trait_plot_an,trait_plot_pl,pop_plot_an,pop_plot_pl,rA,rP,h,const,sig_c_an,sig_c_pl,sig_an,sig_pl,x0,y0,kx,ky,sig_m)

        %%saving the simulation results into a folder named 'Results'
        %%%%%%%%%%%%%%%%%%%%%%
        
        

        clear_files(strcat('./',res_file,'/time_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
        clear_files(strcat('./',res_file,'/animal_traits_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
        clear_files(strcat('./',res_file,'/plant_traits_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
        clear_files(strcat('./',res_file,'/animal_abundance_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
        clear_files(strcat('./',res_file,'/plant_abundance_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
        clear_files(strcat('./',res_file,'/sepIndex_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));

        write_into_files_W(strcat('./',res_file,'/time_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'),time_plot);
        write_into_files_W(strcat('./',res_file,'/animal_traits_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'),trait_plot_an);
        write_into_files_W(strcat('./',res_file,'/plant_traits_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'),trait_plot_pl);
        write_into_files_W(strcat('./',res_file,'/animal_abundance_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'),pop_plot_an);
        write_into_files_W(strcat('./',res_file,'/plant_abundance_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'),pop_plot_pl);
        write_into_files_W(strcat('./',res_file,'/sepIndex_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'),separating_index);
    end
end


