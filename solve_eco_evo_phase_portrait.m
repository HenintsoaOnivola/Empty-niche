function [time, solution,TE,YE,IE]=solve_eco_evo_phase_portrait(INI,rA,rP,h,const,sig_c_an,sig_c_pl,sig_an,sig_pl,x0,y0,kx,ky,sig_m)

    
    ben=@(sig_m,x1,y1) exp(-(1/(2*sig_m^2))*((x1-y1).^2));   
  
    ben_prime=@(sig_m,x1,y1) -(1/sig_m^2)*(x1-y1).*ben(sig_m,x1,y1);

    % competition kernal 
    alphaA=@(sig_c_an,x1,x2) exp(-(1/(2*sig_c_an^2))*((x1-x2).^2));
    alphaP=@(sig_c_pl,x1,x2) exp(-(1/(2*sig_c_pl^2))*((x1-x2).^2));
    
    alphaA_prime=@(sig_c_an,x1,x2) -(1/sig_c_an^2)*(x1-x2).*alphaA(sig_c_an,x1,x2);
    alphaP_prime=@(sig_c_pl,x1,x2) -(1/sig_c_pl^2)*(x1-x2).*alphaP(sig_c_pl,x1,x2);
    

     % carrying capacity function for animal
    function K=k_x(k0,x0,a,x)
        K=k0*(1-((x-x0)/a).^2);
        K(abs(x-x0)>a)=0;
    end

    function K=k_x_prime(k0,x0,a,x)
        K=-2*k0*(1/a^2)*(x-x0);
        K(abs(x-x0)>a)=0;
    end

     
    time_scale_and_mut_rate=10^(-6);
    
    function comp =competition_term_animal(the_traits, the_pops)
        %traits and the_pops should be row vectors, and comp as well
        mat_comp=repmat(the_traits,[length(the_traits),1]);
        comp=alphaA(sig_c_an,mat_comp,transpose(mat_comp));
        comp=the_pops*comp;
    end

    function comp =competition_term_plant(the_traits, the_pops)
        %traits and the_pops should be row vectors, and comp as well
        mat_comp=repmat(the_traits,[length(the_traits),1]);
        comp=alphaP(sig_c_pl,mat_comp,transpose(mat_comp));
        comp=the_pops*comp;
    end

    function comp =competition_term_animal_prime(the_traits, the_pops)
        %traits and the_pops should be row vectors, and comp as well
        mat_comp=repmat(the_traits,[length(the_traits),1]);
        comp=alphaA_prime(sig_c_an,mat_comp,transpose(mat_comp));
        comp=the_pops*comp;
    end

    function comp =competition_term_plant_prime(the_traits, the_pops)
        %traits and the_pops should be row vectors, and comp as well
        mat_comp=repmat(the_traits,[length(the_traits),1]);
        comp=alphaP_prime(sig_c_pl,mat_comp,transpose(mat_comp));
        comp=the_pops*comp;
    end


    function benef=benefit(the_traits, the_other_traits)
        %return the matrix Bij
        RR=repmat(transpose(the_traits),[1,length(the_other_traits)]);
        SS=repmat(the_other_traits,[length(the_traits),1]);
        benef=ben(sig_m,RR,SS);
    end

   
    function benef_pr=benefit_prime(the_traits, the_other_traits)
        %return the matrix Bij
        RR=repmat(transpose(the_traits),[1,length(the_other_traits)]);
        PP=repmat(the_other_traits,[length(the_traits),1]);
        benef_pr=ben_prime(sig_m,RR,PP);
    end


    function gr_r=growth_rate_an(the_traits,the_pops,the_other_pops,benef_mat)
        int_rateA_num=the_other_pops*transpose(benef_mat);
        int_rateA_denum=1+h*(the_other_pops*transpose(benef_mat));
        int_rateA=int_rateA_num./int_rateA_denum;
        if (k_x(kx,x0,sig_an,the_traits)==0)
            gr_r=the_pops.*(rA  + const* int_rateA ) ;
        else
            gr_r=the_pops.*(rA  - ( rA*competition_term_animal(the_traits, the_pops)./k_x(kx,x0,sig_an,the_traits))   + const* int_rateA ) ;
        end
    end

    function gr_r=growth_rate_pl(the_traits,the_pops,the_other_pops,benef_mat)
        int_rateP_num=the_other_pops*transpose(benef_mat);
        int_rateP_denum=1+h*(the_other_pops*transpose(benef_mat));
        int_rateP=int_rateP_num./int_rateP_denum;
        if k_x(ky,y0,sig_pl,the_traits)==0
            gr_r=the_pops.*(rP   + const* int_rateP);
        else
            gr_r=the_pops.*(rP  - ( rA*competition_term_plant(the_traits, the_pops)./k_x(ky,y0,sig_pl,the_traits))   + const* int_rateP)  ;
        end
        
    end

    function sel_gr=selection_gradient_an(the_traits,the_pops,the_other_pops,benefit_mat,benefit_mat_pr)
        
        first_term_num=-rA*((competition_term_animal_prime(the_traits, the_pops).*k_x(kx,x0,sig_an,the_traits)) - (k_x_prime(kx,x0,sig_an,the_traits).*competition_term_animal(the_traits, the_pops)));
        first_term_denum=k_x(kx,x0,sig_an,the_traits).^2;
        if first_term_denum==0
            first_term=0;
        else
            first_term=first_term_num./first_term_denum;
        end
        
        second_term_num_1=(the_other_pops*transpose(benefit_mat_pr)).*(1+h*(the_other_pops*transpose(benefit_mat)));
        second_term_num_2=(h*the_other_pops*transpose(benefit_mat_pr)).*(the_other_pops*transpose(benefit_mat));
        second_term_num=second_term_num_1-second_term_num_2;
        second_term_denum=(1+h*(the_other_pops*(transpose(benefit_mat)))).^2;
        second_term=second_term_num./second_term_denum;
        sel_gr=the_pops.*(first_term+const*second_term);
        
    end
    function sel_gr=selection_gradient_pl(the_traits,the_pops,the_other_pops,benefit_mat,benefit_mat_pr)
        
        first_term_num=-rP*((competition_term_plant_prime(the_traits, the_pops).*k_x(ky,y0,sig_pl,the_traits)) - (k_x_prime(ky,y0,sig_pl,the_traits).*competition_term_plant(the_traits, the_pops)));
        first_term_denum=k_x(ky,y0,sig_pl,the_traits).^2;
        if first_term_denum==0
            first_term=0;
        else
            first_term=first_term_num./first_term_denum;
        end
        
        second_term_num_1=(the_other_pops*transpose(benefit_mat_pr)).*(1+h*(the_other_pops*transpose(benefit_mat)));
        second_term_num_2=(h*the_other_pops*transpose(benefit_mat_pr)).*(the_other_pops*transpose(benefit_mat));
        second_term_num=second_term_num_1-second_term_num_2;
        second_term_denum=(1+h*(the_other_pops*(transpose(benefit_mat)))).^2;
        second_term=second_term_num./second_term_denum;
        sel_gr=the_pops.*(first_term+const*second_term);
    end
   options=odeset('Events', @events);%'RelTol',1e-7,'AbsTol',1e-7);
    
    function [value,isterminal,direction]=events(t,y)
        
        value=ones(1,1);
        y=transpose(y);
        %checking for non extinction
        density=y(1:2);
        traits=y(3:4);
        traits_an=traits(1);
        traits_pl=traits(2);
        pop_an=density(1);
        pop_pl=density(2);
      
        %equilibrium event
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        benef_mat_an=benefit(traits_an, traits_pl);
        benef_mat_pl=benefit(traits_pl,traits_an);
        benef_mat_an_pr=benefit_prime(traits_an, traits_pl);
        benef_mat_pl_pr=benefit_prime(traits_pl,traits_an);

        sel_grad_an=selection_gradient_an(traits_an,pop_an,pop_pl,benef_mat_an,benef_mat_an_pr);
        sel_grad_pl=selection_gradient_pl(traits_pl,pop_pl,pop_an,benef_mat_pl,benef_mat_pl_pr);
        sel_gr=[sel_grad_an sel_grad_pl];
        count_conv=0;

        the_threshold=10^(-8);
        
        for i=1:2
            if abs(sel_gr(i))<the_threshold
                count_conv=count_conv+1;
            end
        end
     
        
        if count_conv==2
            value(1)=0;
        else
            value(1)=1;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        isterminal=ones(size(value));
        direction=zeros(size(value));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
   end
 
    function dxdt=ode_function_model(t,IN)
        
        A=IN(1);
        P=IN(2);
        x=IN(3);
        y=IN(4);
        
        
        A=transpose(A);
        P=transpose(P);
        x=transpose(x);
        y=transpose(y);
         
     
        benef_mat_an=benefit(x,y);
        benef_mat_pl=benefit(y,x);
        benef_mat_an_pr=benefit_prime(x, y);
        benef_mat_pl_pr=benefit_prime(y, x);
     
     
        dA=growth_rate_an(x,A,P,benef_mat_an);
        dP=growth_rate_pl(y,P,A,benef_mat_pl);
        
        dx=time_scale_and_mut_rate*selection_gradient_an(x,A,P,benef_mat_an,benef_mat_an_pr);
        dy=time_scale_and_mut_rate*selection_gradient_pl(y,P,A,benef_mat_pl,benef_mat_pl_pr);
        
        dxdt=[dA dP dx dy];
        dxdt=transpose(dxdt);
        
    end
    [time,solution,TE,YE,IE]=ode15s(@ode_function_model,[0,5000],INI,options);      
end