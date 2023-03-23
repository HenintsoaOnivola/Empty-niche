function all_selection=selection_gradient_function_different_time(an_ab,pl_ab,an_tr,pl_tr,rA,rP,h,const,sig_c_an,sig_c_pl,sig_an,sig_pl,x0,y0,kx,ky,sig_m)
    ben=@(sig_m,x1,y1) exp(-(1/(2*sig_m^2))*((x1-y1).^2));   
  
    ben_prime=@(sig_m,x1,y1) -(1/sig_m^2)*(x1-y1).*ben(sig_m,x1,y1);

    % competition kernal 
    alphaA=@(sig_c_an,x1,x2) exp(-(1/(2*sig_c_an^2))*((x1-x2).^2));
    alphaP=@(sig_c_pl,x1,x2) exp(-(1/(2*sig_c_pl^2))*((x1-x2).^2));
    
    alphaA_prime=@(sig_c_an,x1,x2) -(1/sig_c_an^2)*(x1-x2).*alphaA(sig_c_an,x1,x2);
    alphaP_prime=@(sig_c_pl,x1,x2) -(1/sig_c_pl^2)*(x1-x2).*alphaP(sig_c_pl,x1,x2);
    
    time_scale_and_mut_rate=10^(-6);
    

    % carrying capacity function for animal
    function K=k_x(k0,x0,a,x)
        K=k0*(1-((x-x0)/a).^2);
        K(abs(x-x0)>a)=0;
    end

    function K=k_x_prime(k0,x0,a,x)
        K=-2*k0*(1/a^2)*(x-x0);
        K(abs(x-x0)>a)=0;
    end

    
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
    function sel_gr=selection_gradient_an(the_traits,the_pops,the_other_pops,benefit_mat,benefit_mat_pr)
        
        first_term_num=-rA*((competition_term_animal_prime(the_traits, the_pops).*k_x(kx,x0,sig_an,the_traits)) - (k_x_prime(kx,x0,sig_an,the_traits).*competition_term_animal(the_traits, the_pops)));
        first_term_denum=k_x(kx,x0,sig_an,the_traits).^2;
        first_term=first_term_num./first_term_denum;
        first_term(first_term_num==0)=0;
        first_term(first_term_denum==0)=0;
     
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
        first_term=first_term_num./first_term_denum;
        first_term(first_term_num==0)=0;
        second_term_num_1=(the_other_pops*transpose(benefit_mat_pr)).*(1+h*(the_other_pops*transpose(benefit_mat)));
        second_term_num_2=(h*the_other_pops*transpose(benefit_mat_pr)).*(the_other_pops*transpose(benefit_mat));
        second_term_num=second_term_num_1-second_term_num_2;
        second_term_denum=(1+h*(the_other_pops*(transpose(benefit_mat)))).^2;
        second_term=second_term_num./second_term_denum;
        
        sel_gr=the_pops.*(first_term+const*second_term);
    end

    all_selection=cell(length(an_ab),1);

    for t=1:length(an_ab)
        A=an_ab{t};
        P=pl_ab{t};
        x=an_tr{t};
        y=pl_tr{t};

        benef_mat_an=benefit(x,y);
        benef_mat_pl=benefit(y,x);
        benef_mat_an_pr=benefit_prime(x, y);
        benef_mat_pl_pr=benefit_prime(y, x);
  

        sel_an=time_scale_and_mut_rate*selection_gradient_an(x,A,P,benef_mat_an,benef_mat_an_pr);
        sel_pl=time_scale_and_mut_rate*selection_gradient_pl(y,P,A,benef_mat_pl,benef_mat_pl_pr);
        
    


        sg=abs([sel_an sel_pl]);
        all_selection{t}=sg;
        
    end
  

    
end