%This code simulates the eco-evolutionary dynamics of a
%community in which species are engaged only in inter and intra-specific competition for resources

clear all;

%defining the parameters
rA=0.1;
x0=3;
kx=400;

mutant_fraction=0.1;

log_sig_an=1;
log_sig_c_an=-0.475;

sig_an=log_sig_an;   %this is the width of resource kernel
sig_c_an=exp(log_sig_c_an);  %this is the width of competition kernel

%initialization
ini_val=[1 2.5]; %initial population abundance and trait value;
animal=1;   %initial number of species

infinite=300;
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
total_time=0;
an_num=[animal];

%looping over new system
while (new_sys==true && count_into_while<15)
    count_into_while=count_into_while+1;
    disp (['count system: ' num2str(count_into_while)]);

    negative_pop=false;
    convergence=false;
    count_integrate=0;
    new_sys=false;

    temporary_time=0;
    while ( (convergence==false) && count_integrate<=infinite )
        tic;

        [time,sol,TE,YE,IE]=solve_eco_evo_RC(ini_val,rA,sig_c_an,sig_an,x0,kx,animal,temporary_time,count_integrate,count_into_while);
        fprintf ('%d\t',count_integrate);
        temporary_time=temporary_time+time(end);
        if count_integrate==infinite  %consider the case like a convergence
            IE=2;
        end

        if isempty (IE)
            %storing the traits and populations
            an_pop=sol(end,1:animal);
            an_trait=sol(end,animal+1:end);
            animal=length(an_pop);

            current_time_range=time+previous_t_end;
            time_plot(time_index+1:time_index+numel(time))=current_time_range;

            the_pop_an=sol(:,1:animal);
            solution_pop_an=mat2cell(the_pop_an,ones(1,size(the_pop_an,1)),size(the_pop_an,2));
            pop_plot_an(time_index+1:time_index+numel(time))=solution_pop_an;
            
            the_trait_an=sol(:,animal+1:end);
            solution_trait_an=mat2cell(the_trait_an,ones(1,size(the_trait_an,1)),size(the_trait_an,2));
            trait_plot_an(time_index+1:time_index+numel(time))=solution_trait_an;


            previous_t_end=time_plot(time_index+numel(time));
            count_integrate=count_integrate+1;
            time_index=time_index+numel(time);

            %reinitialization
            animal=length(an_pop);
            ini_val=[an_pop an_trait];

        else

            if any(IE==1)
                disp ('there is extinction');

                an_pop_full=sol(end,1:animal);
                an_pop=an_pop_full(an_pop_full>1e-8);
                an_trait=sol(end,animal+1:end);
                an_trait=an_trait(an_pop_full>1e-8);

                %storing the traits and populations
                current_time_range=time+previous_t_end;
                time_plot(time_index+1:time_index+numel(time))=current_time_range;

                the_pop_an=sol(:,1:animal);
                solution_pop_an=mat2cell(the_pop_an,ones(1,size(the_pop_an,1)),size(the_pop_an,2));
                pop_plot_an(time_index+1:time_index+numel(time))=solution_pop_an;
                
                the_trait_an=sol(:,animal+1:end);
                solution_trait_an=mat2cell(the_trait_an,ones(1,size(the_trait_an,1)),size(the_trait_an,2));
                trait_plot_an(time_index+1:time_index+numel(time))=solution_trait_an;
    
                separating_index=[separating_index time_index+numel(time)];


                previous_t_end=time_plot(time_index+numel(time));
                count_integrate=count_integrate+1;
                time_index=time_index+numel(time);


                if isempty(an_pop)==1
                    disp ('the entire population is extinct');
                    separating_index=[separating_index -1];
                    new_sys=false;
                    break;
                end

                %reinitialization
                animal=length(an_pop);
                ini_val=[an_pop an_trait];

            else

                if any (IE==2) 

                    convergence=true;
                    %storing the traits and populations
                    an_pop=sol(end,1:animal);
                    an_trait=sol(end,animal+1:end);
                    animal=length(an_pop);

                    current_time_range=time+previous_t_end;
                    time_plot(time_index+1:time_index+numel(time))=current_time_range;


                    the_pop_an=sol(:,1:animal);
                    solution_pop_an=mat2cell(the_pop_an,ones(1,size(the_pop_an,1)),size(the_pop_an,2));
                    pop_plot_an(time_index+1:time_index+numel(time))=solution_pop_an;
                
                    the_trait_an=sol(:,animal+1:end);
                    solution_trait_an=mat2cell(the_trait_an,ones(1,size(the_trait_an,1)),size(the_trait_an,2));
                    trait_plot_an(time_index+1:time_index+numel(time))=solution_trait_an;
                    
                    separating_index=[separating_index time_index+numel(time)];

                    previous_t_end=time_plot(time_index+numel(time));
                    count_integrate=count_integrate+1;
                    time_index=time_index+numel(time);

                    fprintf ('\n');

                    %merge very similar traits if needed
                    an_before=length(an_trait);

                    previous_an_tr=an_trait;
                    previous_an_pop=an_pop;

                    an_merged=merge(an_trait,an_pop);
                    an_trait=an_merged{1};
                    an_pop=an_merged{2};

                    an_after=length(an_trait);

                    %reinitialization
                    animal=length(an_pop);
                    ini_val=[an_pop an_trait];
                    
                    an_num=[an_num animal];
                    if count_into_while>=2

                        if an_num(end-1)==animal
                            new_sys=false;
                            break;
                        end
                    end

                    %check for branching events
                    sel_grad=selection_gradient_function_RC([an_pop an_trait ],rA,sig_c_an,sig_an,x0,kx,animal);
                    disr=branching_condition_RC([an_pop an_trait ],rA,sig_c_an,sig_an,x0,kx,animal);
                    
                    new_an_pop=zeros(1,2*length(an_pop));
                    new_an_tr=zeros(1,2*length(an_trait));
                    new_an_pop(:)=NaN;
                    new_an_tr(:)=NaN;
                    new_an_pop(1:length(an_pop))=an_pop;
                    new_an_tr(1:length(an_trait))=an_trait;

                    count_new_mutant=0;
                    an_br_index=cell(2*animal,1);
                    count_br_index=1;
                    for i=1:length(an_pop)
                        if  disr(i)>0
                            disp ('an is br');
                            new_sys=true;
                            count_new_mutant=count_new_mutant+1;
                            new_an_pop(i)=an_pop(i)*(1-mutant_fraction);
                            new_an_pop(length(an_pop)+count_new_mutant)=an_pop(i)*(mutant_fraction);
                            new_an_tr(i)=an_trait(i);
                            if an_trait(i)>x0
                                new_an_tr(length(an_pop)+count_new_mutant)=an_trait(i)+0.005;
                            else
                                new_an_tr(length(an_pop)+count_new_mutant)=an_trait(i)-0.005;
                            end
                       
                            an_br_index{count_br_index}=[i length(an_pop)+count_new_mutant];
                            count_br_index=count_br_index+1;
                        end
                    end
                    an_br_index=an_br_index(~cellfun('isempty',an_br_index));

                    new_an_pop=new_an_pop(~isnan(new_an_pop));
                    new_an_tr=new_an_tr(~isnan(new_an_tr));

                    %reinitialization
                    animal=length(new_an_pop);
                    ini_val=[new_an_pop new_an_tr];
                    
                    %%%%%%%%%%%%%%%%%%%%%%
                    %check that they are real branching
                    initial_condition=ini_val;
                    continue_check_br=1;
                    temp_total_time=0;
                    count_temporary=0;
                    while(temp_total_time<=10^(6))
                        fprintf ('%d\t',count_temporary);
                        tic;
                        [temp_time,temp_sol,temp_TE, temp_YE, temp_IE]=solve_eco_evo_RC(initial_condition,rA,sig_c_an,sig_an,x0,kx,animal,temporary_time,count_integrate,count_temporary);
                        temp_total_time=temp_total_time+temp_time(end);
                        count_temporary=count_temporary+1;
                        initial_condition=temp_sol(end,:);
                        if (~isempty(temp_IE) && any (temp_IE==3))  || count_temporary>20
                            ind_oscillation=detect_oscillation_using_trait_variation(temp_sol,animal);
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
                        temp_an_trait=temp_sol(end,animal+1:2*animal);
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
                               
                            else
                                disp ('an branching confirmed');
                            end

                        end


                        new_an_pop=new_an_pop(~isnan(new_an_pop));
                        new_an_tr=new_an_tr(~isnan(new_an_tr));

                        if length(new_an_pop)==length(an_pop) 
                            new_sys=false;

                        end
                        ini_val=[new_an_pop  new_an_tr ];
                        animal=length(new_an_pop);
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                 

                    if new_sys==false
                        disp('ESS');
                    end
                    if an_before~=an_after  %if there we merging, then store the last merged traits

                        current_time_range=time+previous_t_end;
                        time_plot(time_index+1)=time_plot(time_index);

                        pop_plot_an{time_index+1}=an_pop;
                        trait_plot_an{time_index+1}=an_trait;
                        separating_index=[separating_index time_index+1];
                       
                    end
                end

                if any (IE==3)
                    %detect oscillation
                    ind_oscillation=detect_oscillation_using_trait_variation(sol,animal);
                    ind_oscillation_an=ind_oscillation{1};
                    if ~isempty(ind_oscillation_an)
                        oscillation=true;
                    else
                        oscillation=false;
                    end
                    if oscillation==true
                        %storing the traits and populations and
                        an_pop=sol(end,1:animal);
                        an_trait=sol(end,animal+1:end);
                        animal=length(an_pop);

                        current_time_range=time+previous_t_end;
                        time_plot(time_index+1:time_index+numel(time))=current_time_range;

                        the_pop_an=sol(:,1:animal);
                        solution_pop_an=mat2cell(the_pop_an,ones(1,size(the_pop_an,1)),size(the_pop_an,2));
                        pop_plot_an(time_index+1:time_index+numel(time))=solution_pop_an;
                        
                        the_trait_an=sol(:,animal+plant+1:animal+plant+animal);
                        solution_trait_an=mat2cell(the_trait_an,ones(1,size(the_trait_an,1)),size(the_trait_an,2));
                        trait_plot_an(time_index+1:time_index+numel(time))=solution_trait_an;
                    
                        separating_index=[separating_index time_index+numel(time)];

                        previous_t_end=time_plot(time_index+numel(time));
                        count_integrate=count_integrate+1;
                        time_index=time_index+numel(time);

                        %reinitialization
                        animal=length(an_pop);
                        ini_val=[an_pop an_trait];
                        new_sys=false;
                        disp ('oscillation');
                        separating_index=[separating_index 0];
                        break;
                    else
                        %just update the system
                        an_pop=sol(end,1:animal);
                        an_trait=sol(end,animal+1:end);
                        
                        animal=length(an_pop);

                        current_time_range=time+previous_t_end;
                        time_plot(time_index+1:time_index+numel(time))=current_time_range;

                        the_pop_an=sol(:,1:animal);
                        solution_pop_an=mat2cell(the_pop_an,ones(1,size(the_pop_an,1)),size(the_pop_an,2));
                        pop_plot_an(time_index+1:time_index+numel(time))=solution_pop_an;
                       
                        the_trait_an=sol(:,animal+plant+1:animal+plant+animal);
                        solution_trait_an=mat2cell(the_trait_an,ones(1,size(the_trait_an,1)),size(the_trait_an,2));
                        trait_plot_an(time_index+1:time_index+numel(time))=solution_trait_an;

                        previous_t_end=time_plot(time_index+numel(time));
                        count_integrate=count_integrate+1;
                        time_index=time_index+numel(time);

                        %reinitialization
                        animal=length(an_pop);
                        ini_val=[an_pop an_trait];
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

plot_traits_function_RC(separating_index,time_plot,trait_plot_an,pop_plot_an,rA,sig_c_an,sig_an,x0,kx);

%%Saving the results of the simulation into a folder named 'Results_RC'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clear_files(strcat('./Results_RC/time_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
% clear_files(strcat('./Results_RC/animal_traits_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
% clear_files(strcat('./Results_RC/animal_abundance_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
% clear_files(strcat('./Results_RC/sepIndex_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'_logsigM',num2str(log_sig_m),'.dat'));
% write_into_files_W(strcat('./Results_RC/time_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'.dat'),time_plot);
% write_into_files_W(strcat('./Results_RC/animal_traits_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'.dat'),trait_plot_an);
% write_into_files_W(strcat('./Results_RC/animal_abundance_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'.dat'),pop_plot_an);
% write_into_files_W(strcat('./Results_RC/sepIndex_logsigA_',num2str(log_sig_an),'_logsigC',num2str(log_sig_c_an),'.dat'),separating_index);

