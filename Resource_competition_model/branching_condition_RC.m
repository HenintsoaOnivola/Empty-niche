function info=branching_condition(INI,rA,sig_c_an,sig_an,x0,kx,animal)
  
    % competition kernal 
    alphaA=@(sig_c_an,x1,x2) exp(-(1/(2*sig_c_an^2))*((x1-x2).^2));
    
    alphaA_prime=@(sig_c_an,x1,x2) -(1/sig_c_an^2)*(x1-x2).*alphaA(sig_c_an,x1,x2);
    
    alphaA_second=@(sig_c_an,x1,x2) -(1/sig_c_an^2)*(alphaA(sig_c_an,x1,x2)+(alphaA_prime(sig_c_an,x1,x2).*(x1-x2)));

    % carrying capacity function for animal
    function K=k_x(k0,x0,a,x)
        K=k0*(1-((x-x0)/a).^2);
        K(abs(x-x0)>a)=0;
    end

    function K=k_x_prime(k0,x0,a,x)
        K=-2*k0*(1/a^2)*(x-x0);
        K(abs(x-x0)>a)=0;
    end

    function K=k_x_second(k0,x0,a,x)
        K=-2*k0*(1/a^2);
        K(abs(x-x0)>a)=0;
    end
    
    
    function comp =competition_term(the_traits, the_pops)
        %traits and the_pops should be row vectors, and comp as well
        mat_comp=repmat(the_traits,[length(the_traits),1]);
        comp=alphaA(sig_c_an,mat_comp,transpose(mat_comp));
        comp=the_pops*comp;
    end

    function comp =competition_term_prime(the_traits, the_pops)
        %traits and the_pops should be row vectors, and comp as well
        mat_comp=repmat(the_traits,[length(the_traits),1]);
        comp=alphaA_prime(sig_c_an,mat_comp,transpose(mat_comp));
        comp=the_pops*comp;
    end
   

    function comp =competition_term_second(the_traits, the_pops)
        %traits and the_pops should be row vectors, and comp as well
        mat_comp=repmat(the_traits,[length(the_traits),1]);
        comp=alphaA_second(sig_c_an,mat_comp,transpose(mat_comp));
        comp=the_pops*comp;
    end
   

    function disr=disruptive_selection(r,carry_cap,carry_cap_pr,carry_cap_sc,competition_term,competition_term_pr,competition_term_sc)
      
        disr_comp_1=(competition_term_sc.*carry_cap - carry_cap_sc.*competition_term);
        disr_comp_2=(carry_cap).^(2);
        disr_comp_3=competition_term_pr.*carry_cap-carry_cap_pr.*competition_term;
        disr_comp_4=2*carry_cap_pr.*carry_cap;
        disr_comp=-r*( ((disr_comp_1.*disr_comp_2)-(disr_comp_3.*disr_comp_4))./((carry_cap).^(4)));
       
        disr=disr_comp;
    end
    
    A=INI(1:animal);
    x=INI(animal+1:end);
  
    carry_cap_an=k_x(kx,x0,sig_an,x);
    carry_cap_an_pr=k_x_prime(kx,x0,sig_an,x);
    carry_cap_an_sc=k_x_second(kx,x0,sig_an,x);
    competition_term_an=competition_term(x,A);
    competition_term_an_pr=competition_term_prime(x,A);
    competition_term_an_sc=competition_term_second(x,A);
  
    disr_an=disruptive_selection(rA,carry_cap_an,carry_cap_an_pr,carry_cap_an_sc,competition_term_an,competition_term_an_pr,competition_term_an_sc);
    
    info=disr_an;

end