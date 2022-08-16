function sg=selection_gradient_function(IN,rA,sig_c_an,sig_an,x0,kx,animal)

    alphaA=@(sig_c_an,x1,x2) exp(-(1/(2*sig_c_an^2))*((x1-x2).^2));
    alphaA_prime=@(sig_c_an,x1,x2) -(1/sig_c_an^2)*(x1-x2).*alphaA(sig_c_an,x1,x2);

    function K=k_x(k0,x0,a,x)
        K=k0*(1-((x-x0)/a).^2);
        K(abs(x-x0)>a)=0;
    end

    function K=k_x_prime(k0,x0,a,x)
        K=-2*k0*(1/a^2)*(x-x0);
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

    function sel_gr=selection_gradient(the_traits,the_pops)
        
        
        first_term_num=-rA*((competition_term_prime(the_traits, the_pops).*k_x(kx,x0,sig_an,the_traits)) - (k_x_prime(kx,x0,sig_an,the_traits).*competition_term(the_traits, the_pops)));
        first_term_denum=k_x(kx,x0,sig_an,the_traits).^2;
        first_term=first_term_num./first_term_denum;
        sel_gr=the_pops.*(first_term);
    end
  

    A=IN(1:animal);
    x=IN(animal+1:end);
    sg=selection_gradient(x,A);

end