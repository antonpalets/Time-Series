function y = sim_arma(A, C, sigma_e2, N)
    if nargin<3
        sigma_e2 = 1;
    end
    if nargin<4
        N = 250;
    else
        N = N+150;
    end
    
    rng(0)
    if nargin==5
        rng(rand_seed)
    end
    e = sqrt(sigma_e2)*randn(N,1);
    y = filter(C, A, e);
    y = y(151:end);
end