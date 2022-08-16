function fitness=fitness_land_function(the_mutant,res_pop,other_res_pop,res_trait,other_res_trait,r,h,const,sig_c,sig,x0,kx,sig_m)

    ben=@(sig_m,x1,y1) exp(-(1/(2*sig_m^2))*((x1-y1).^2));   
  
  
    % competition kernal 
    alpha=@(sig_c,x1,x2) exp(-(1/(2*sig_c^2))*((x1-x2).^2));
    
 
    % carrying capacity function for animal
    function K=k_cap(k0,x0,a,x)
        K=k0*(1-((x-x0)/a).^2);
        K(abs(x-x0)>a)=0;
    end
   
    function comp =competition_term_mutant(the_mutant,the_traits, the_pops)
        %traits and the_pops should be row vectors, and comp as well
        mat_comp=repmat(the_mutant,[length(the_traits),1]);
        comp=alpha(sig_c,mat_comp,transpose(the_traits));
        comp=the_pops*comp;
    end

    function benef=benefit(the_traits, the_other_traits)
        %return the matrix Bij
        RR=repmat(transpose(the_traits),[1,length(the_other_traits)]);
        SS=repmat(the_other_traits,[length(the_traits),1]);
        benef=ben(sig_m,RR,SS);
    end

   

    function gr_r=per_cap_growth_rate_mutant(the_mutant,the_traits,the_pops,the_other_pops,benef_mat_mutant)
        int_rate_num=the_other_pops*transpose(benef_mat_mutant);
        int_rate_denum=1+h*(the_other_pops*transpose(benef_mat_mutant));
        int_rate=int_rate_num./int_rate_denum;
        
        if int_rate_num==Inf && int_rate_denum==Inf
            int_rate=1;
        end
        
        comp_t=competition_term_mutant(the_mutant,the_traits, the_pops);
        car_t=k_cap(kx,x0,sig,the_mutant);
      
        the_division=comp_t./car_t;
        
        if car_t==0 
            the_division=1000;
        end
       
%         disp (benef_mat_mutant);
%         disp ([int_rate_num int_rate_denum]);
%         disp ([r  r*the_division  const*int_rate ]);
   
        gr_r=r  - (r*( the_division))   + const* int_rate  ;
        
        
     
    end


    fitness=zeros(1,length(the_mutant));
    for len_mutant=1:length(the_mutant)
        the_mut=the_mutant(len_mutant);
        benef_mat_mutant=benefit(the_mut,other_res_trait);
        the_fit=per_cap_growth_rate_mutant(the_mut,res_trait,res_pop,other_res_pop,benef_mat_mutant);
        fitness(len_mutant)=the_fit;
    end
end