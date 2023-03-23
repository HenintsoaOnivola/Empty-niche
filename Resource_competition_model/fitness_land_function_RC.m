function fitness=fitness_land_function_RC(the_mutant,res_pop,res_trait,r,sig_c,sig,x0,kx)

   
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

    function gr_r=per_cap_growth_rate_mutant(the_mutant,the_traits,the_pops)
        
        num=competition_term_mutant(the_mutant,the_traits, the_pops);
        denum=k_cap(kx,x0,sig,the_mutant);
     
        
        the_division=num./denum;
        
        if num==0 && denum==0
            the_division=10;
        end
       
        

        gr_r=r  - ( r*the_division)    ;
       
       
    end

    fitness=zeros(1,length(the_mutant));
    for len_mutant=1:length(the_mutant)
        the_mut=the_mutant(len_mutant);
        fitness(len_mutant)=per_cap_growth_rate_mutant(the_mut,res_trait,res_pop);
    end
end