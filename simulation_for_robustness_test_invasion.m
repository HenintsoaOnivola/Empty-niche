function res=simulation_for_robustness_test_invasion(rA,rP,const,x0,y0,kx,ky,sig_an,log_sig_c_an,log_sig_m,pre_time,pre_an_tr,pre_pl_tr,pre_an_ab,pre_pl_ab,pre_separating_index,the_fold)
    
    %defining the parameters
    mutant_fraction=0.1;
    invader_fraction=0.01;
    h=0.1;
    sig_pl=sig_an;
    sig_c_an=exp(log_sig_c_an);
    sig_c_pl=sig_c_an;
    sig_m=exp(log_sig_m);
    
    %%get data from the pre-invasion simulation
    last_tr_an=pre_an_tr{end};
    last_tr_pl=pre_pl_tr{end};
    last_pop_an=pre_an_ab{end};
    last_pop_pl=pre_pl_ab{end};

    peaks_inter=localize_positive_fitness(last_pop_an,last_pop_pl,last_tr_an,last_tr_pl,rA,rP,h,const,log(sig_c_an),log(sig_c_pl),sig_an,sig_pl,x0,y0,kx,ky,log(sig_m));
    the_interval_an=peaks_inter{1};
    the_interval_pl=peaks_inter{2};
    
   if ~isempty(the_interval_an)
        inv_an=randsample(length(the_interval_an),1);
        inv_an=the_interval_an(inv_an);
        
    else
        inv_an=[];
   end

    if ~isempty(the_interval_pl)
        inv_pl=randsample(length(the_interval_pl),1);
        inv_pl=the_interval_pl(inv_pl);
    else
        inv_pl=[];
    end
    an_introduced=inv_an;
    pl_introduced=inv_pl;

    oscillation=false;
    if any(pre_separating_index==0)
        oscillation=true;
    end

    separating_index=pre_separating_index;
    intro_index=[];

    count_introduction=0;
    max_introduction=40;
    invade=true;
    repetition=false;
    count_repetition=0;
    while invade==true && count_introduction<max_introduction && repetition==false %&& oscillation==false
        count_introduction=count_introduction+1;
        disp ('number of introduction:')
        disp (count_introduction);

        disp ('animal introduced at:');
        disp (inv_an);
        disp ('plant introduced at:');
        disp (inv_pl);

        an_trait=[last_tr_an inv_an];
        an_pop=[last_pop_an repmat(invader_fraction*min(last_pop_an),[1,length(inv_an)])];

        pl_trait=[last_tr_pl inv_pl];
        pl_pop=[last_pop_pl repmat(invader_fraction*min(last_pop_pl),[1, length(inv_pl)])];

        ini_val=[an_pop pl_pop an_trait pl_trait];
        animal=length(an_trait);
        plant=length(pl_trait);


        %initialization

        infinite=200;
        new_sys=true;
        count_into_while=0;
        previous_t_end=pre_time(end);

        time_plot=zeros(7e005,1);
        time_plot(1:length(pre_time))=pre_time;
        time_plot(length(pre_time)+1:end)=NaN;
        pop_plot_an=cell(7e005,1);
        pop_plot_an(1:length(pre_an_ab))=pre_an_ab;
        pop_plot_pl=cell(7e005,1);
        pop_plot_pl(1:length(pre_pl_ab))=pre_pl_ab;
        trait_plot_an=cell(7e005,1);
        trait_plot_an(1:length(pre_an_tr))=pre_an_tr;
        trait_plot_pl=cell(7e005,1);
        trait_plot_pl(1:length(pre_pl_tr))=pre_pl_tr;
        time_index=length(pre_time);
        not_integrated=false;

        extinction_cell_an=cell(10,1);
        extinction_time_an=zeros(10,1);
        count_extinction_an=0;
        extinction_cell_pl=cell(10,1);
        extinction_time_pl=zeros(10,1);
        count_extinction_pl=0;

        end_each_sys=[];

        intro_index=[intro_index length(pre_time)+1];
        %looping over the system until we get ESS
        an_num=[animal];
        pl_num=[plant];
        maximum_branching=15;
        while (new_sys==true && count_into_while<maximum_branching)
            count_into_while=count_into_while+1;
            disp (['count system: ' num2str(count_into_while)]);
            fprintf ('\n');

            convergence=false;
            count_integrate=0;

            total_number=1;


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
                if count_integrate==infinite  %consider the case like a convergence
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

                            disp (pl_trait);


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

                            disp (pl_trait);

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
                                    diff_tr=abs(an_tr_pr-an_trait);

                                    if sum(diff_tr<10^(-5))==length(diff_tr)
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
                            new_sys=false;
                            if count_into_while<maximum_branching+1
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
                                        cl=localize_closest_empty_niche(an_trait(i),an_pop,pl_pop,an_trait,pl_trait,rA,h,const,sig_c_an,sig_an,x0,kx,sig_m);
                                        if strcmp(cl,'left')
                                            new_an_tr(length(an_pop)+count_new_mutant)=an_trait(i)-0.005;
                                        else
                                            new_an_tr(length(an_pop)+count_new_mutant)=an_trait(i)+0.005;
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
                                        cl=localize_closest_empty_niche(pl_trait(i),pl_pop,an_pop,pl_trait,an_trait,rP,h,const,sig_c_pl,sig_pl,y0,ky,sig_m);
                                        if strcmp(cl,'left')
                                            new_pl_tr(length(pl_pop)+count_new_mutant)=pl_trait(i)-0.005;
                                        else
                                            new_pl_tr(length(pl_pop)+count_new_mutant)=pl_trait(i)+0.005;
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
                                while(temp_total_time<=10^(5))
                                    fprintf ('%d\t',count_temporary);
                                    tic;
                                    [temp_time,temp_sol,temp_TE, temp_YE, temp_IE]=solve_eco_evo_subsequent_branching(initial_condition,rA,rP,h,const,sig_c_an,sig_c_pl,sig_an,sig_pl,x0,y0,kx,ky,sig_m,animal,plant,temporary_time,count_temporary,count_into_while);
                                    temp_total_time=temp_total_time+temp_time(end);
                                    count_temporary=count_temporary+1;
                                    initial_condition=temp_sol(end,:);
                                    if (~isempty(temp_IE) && any (temp_IE==3))  || count_temporary>20
                                        ind_oscillation=detect_oscillation_using_trait_variation(temp_sol,animal,plant);
                                        if ind_oscillation==true
                                            continue_check_br=0;
                                            disp ('there is an oscillation');
                                            break;
                                        end
                                    end
                                    if ~isempty(temp_IE) && any(temp_IE==1)
                                        continue_check_br=0;
                                        disp ('there is an extinction event');
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
                                        disp ('ESS confirmed');
                                        new_sys=false;


                                    end
                                    ini_val=[new_an_pop new_pl_pop new_an_tr new_pl_tr];
                                    animal=length(new_an_pop);
                                    plant=length(new_pl_pop);
                                end
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                %%%%%%%%%%%%%%%%%%%%%%%

                                if new_sys==false
                                    disp ('ESS');
                                end

                            end
                        end

                        if any (IE==3)
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

        end


        time_plot=time_plot(~isnan(time_plot));
        pop_plot_an=pop_plot_an(1:numel(time_plot),:);
        trait_plot_an=trait_plot_an(1:numel(time_plot),:);
        pop_plot_pl=pop_plot_pl(1:numel(time_plot),:);
        trait_plot_pl=trait_plot_pl(1:numel(time_plot),:);

        %check whether there are remaining empty niches
        last_tr_an=trait_plot_an{end};
        last_tr_pl=trait_plot_pl{end};
        last_pop_an=pop_plot_an{end};
        last_pop_pl=pop_plot_pl{end};
        
        peaks_inter=localize_positive_fitness(last_pop_an,last_pop_pl,last_tr_an,last_tr_pl,rA,rP,h,const,log(sig_c_an),log(sig_c_pl),sig_an,sig_pl,x0,y0,kx,ky,log(sig_m));
        the_interval_an=peaks_inter{1};
        the_interval_pl=peaks_inter{2};


        if (~isempty(the_interval_an) || ~isempty(the_interval_pl)) %&& oscillation==false
            invade=true;

            %reinitialize the system
            pre_an_tr=trait_plot_an;
            pre_pl_tr=trait_plot_pl;
            pre_an_ab=pop_plot_an;
            pre_pl_ab=pop_plot_pl;
            pre_time=time_plot;

            if ~isempty(the_interval_an)

                inv_an=randsample(length(the_interval_an),1);
                inv_an=the_interval_an(inv_an);
                disp ([count_introduction count_into_while]);
                if count_introduction==1 && count_into_while==3
                    inv_an=the_interval_an(1);
                end

                diff_an=abs(inv_an-an_introduced);
                count_choice=1;
                while any(diff_an<=10^(-3))  && count_choice<=10

                    inv_an=randsample(length(the_interval_an),1);
                    inv_an=the_interval_an(inv_an);
                    diff_an=abs(inv_an-an_introduced);
                    count_choice=count_choice+1;
                end
                an_introduced=[an_introduced inv_an];
            else
                inv_an=[];
            end

            if ~isempty(the_interval_pl)
                inv_pl=randsample(length(the_interval_pl),1);
                inv_pl=the_interval_pl(inv_pl);
                disp ([count_introduction count_into_while])
                if count_introduction==1 && count_into_while==3
                    inv_pl=the_interval_pl(end);
                end
                diff_pl=abs(inv_pl-pl_introduced);
                count_choice=1;
                while any(diff_pl<=10^(-3)) && count_choice<=10
                    inv_pl=randsample(length(the_interval_pl),1);
                    inv_pl=the_interval_pl(inv_pl);
                    diff_pl=abs(inv_pl-pl_introduced);
                    count_choice=count_choice+1;
                end
                pl_introduced=[pl_introduced inv_pl];
            else
                inv_pl=[];
            end


            %check whether there are system repetitions
            if count_introduction==1
                initial_pop_size_an=length(last_pop_an);
                initial_pop_size_pl=length(last_pop_pl);
                the_initial_pop_an=last_pop_an;
                the_initial_pop_pl=last_pop_pl;
                the_initial_tr_an=last_tr_an;
                the_initial_tr_pl=last_tr_pl;
            end
            if count_introduction>2
                current_pop_size_an=length(last_pop_an);
                current_pop_size_pl=length(last_pop_pl);
                the_current_pop_an=last_pop_an;
                the_current_pop_pl=last_pop_pl;
                the_current_tr_an=last_tr_an;
                the_current_tr_pl=last_tr_pl;

                if current_pop_size_an==initial_pop_size_an && current_pop_size_pl==initial_pop_size_pl
                    the_diff_an=abs(the_initial_pop_an-the_current_pop_an);
                    the_diff_pl=abs(the_initial_pop_pl-the_current_pop_pl);

                    the_diff_an_tr=abs(the_initial_tr_an-the_current_tr_an);
                    the_diff_pl_tr=abs(the_initial_tr_pl-the_current_tr_pl);

                    if ~any(the_diff_an>10^(-1)) && ~any(the_diff_pl>10^(-1))
                        if ~any(the_diff_an_tr>10^(-4)) && ~any(the_diff_pl_tr>10^(-4))
                            count_repetition=count_repetition+1;
                        end
                    end
                    if count_repetition>=2
                        repetition=true;
                    end
                else
                    initial_pop_size_an=current_pop_size_an;
                    initial_pop_size_pl=current_pop_size_pl;
                    the_initial_pop_an=the_current_pop_an;
                    the_initial_pop_pl=the_current_pop_pl;
                    the_initial_tr_an=the_current_tr_an;
                    the_initial_tr_pl=the_current_tr_pl;
                end
            end

        else
            invade=false;
        end

        extinction_cell_an=extinction_cell_an(1:count_extinction_an);
        extinction_time_an=extinction_time_an(1:count_extinction_an);
        extinction_cell_pl=extinction_cell_pl(1:count_extinction_pl);
        extinction_time_pl=extinction_time_pl(1:count_extinction_pl);


    end
    if isempty(the_interval_an)
        the_interval_an=[min(last_tr_an)-0.1 max(last_tr_an)+0.1];
    end
    if isempty(the_interval_pl)
        the_interval_pl=[min(last_tr_pl)-0.1 max(last_tr_pl)+0.1];
    end
    the_inter_min=min([min(the_interval_an) min(the_interval_pl)]);
    the_inter_max=max([max(the_interval_an) max(the_interval_pl)]);
    the_inter=[the_inter_min the_inter_max];

    res={time_plot,trait_plot_an,trait_plot_pl,pop_plot_an,pop_plot_pl,separating_index,intro_index,{the_inter}};
    %plot_traits_function_with_fitness(separating_index,time_plot,trait_plot_an,trait_plot_pl,pop_plot_an,pop_plot_pl,rA,rP,h,const,sig_c_an,sig_c_pl,sig_an,sig_pl,x0,y0,kx,ky,sig_m)

end

