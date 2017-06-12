using Distributions

# convert x to a lower triangular matrix
function f1{T}(x::Vector{T})
    n = Int(sqrt(2 * length(x) + 0.25) - 0.5)
    y = zeros(T, n, n)
    k = 1
    for i = 1:n
        for j = i:n
            y[j,i] = x[k]
            k += 1
        end
    end
    y
end


# link function for GLMMs
function h1(link::String, x::Vector)
    if link=="log"
         y = exp(x)
    end
    if link=="logit"
         y = log(1+exp(x))
    end
    return(y)
end


# differentiation of link function for GLMMs
function h2(link::String, x::Vector)
    if link=="log"
        y = exp(x)
    end
    if link=="logit"
        y = 1./(1+exp(-x))
    end
    return(y)
end


# compute log h(theta) for GLMM #
function logh_GLMM(link::String, theta::Vector, y::Array{Array{Float64,1},1}, 
    X::Array{Array{Float64,2},1}, Z::Array{Array{Float64,2},1},  
    Idiag::Array{Bool,1}, vbeta0::Float64, vzeta0::Float64)
    
    n = length(y)                           # no. of subjects
    k = size(X[1],2)                        # length of beta
    p = size(Z[1],2)                        # length of b_i
    beta = theta[(p*n)+1:(p*n)+k]
    zeta = theta[p*n+k+1:end]
    W = f1(zeta)
    W[diagind(W)] = exp(diag(W))
    L = 0.0
    for i in 1:n
        b_i = theta[(i-1)*p+1:i*p]
        Wb = W\b_i
        eps = X[i]*beta + Z[i]*b_i
        L += vecdot(y[i],eps) - sum(h1(link,eps)) - 0.5sum(Wb.^2)
    end
    L += -n*sum(zeta[Idiag]) - 0.5sum(beta.^2)/vbeta0 - 0.5sum(zeta.^2)/vzeta0
    return(L)  
end

# logh_GLMM(link, theta, y, X, Z, Idiag, vbeta0, vzeta0)
    
  

# compute gradient of log h(theta) for GLMM #
function gradlogh_GLMM(link::String, theta::Vector, y::Array{Array{Float64,1},1}, 
    X::Array{Array{Float64,2},1}, Z::Array{Array{Float64,2},1}, id::Array{Int64,1}, 
    Idiag::Array{Bool,1}, vbeta0::Float64, vzeta0::Float64)

    n = length(y)                  # no. of subjects
    k = size(X[1],2)               # length of beta
    p = size(Z[1],2)               # length of b_i
    d = Int(n*p + k + 0.5p*(p+1))  # length of theta = [b_1, ...,b_n, beta, zeta]
    g = Array(Float64,d)
    beta = theta[(p*n)+1:(p*n)+k]
    zeta = theta[p*n+k+1:end]
    W = f1(zeta)
    W[diagind(W)] = exp(diag(W))
    S1 = zeros(k)
    S2 = zeros(p,p)
    
    # gradient for b_i
    for i in 1:n
        b_i = theta[(i-1)*p+1:i*p]
        Wb = W\b_i
        s = y[i] - h2(link,X[i]*beta + Z[i]*b_i)
        S1 += X[i]'*s
        S2 += Wb * Wb'
        g[(i-1)*p+1:i*p] = Z[i]'*s - W'\Wb    
    end
    
    # gradient for beta
    g[(p*n)+1:(p*n)+k] = S1 - beta/vbeta0

    # gradient for zeta
    A = W'\S2
    A = A[id]
    A[Idiag] = A[Idiag].*diag(W)
    g[p*n+k+1:end] = - n*Idiag + A - zeta/vzeta0

    return g
end


# compute log h(theta) for SMM #
function logh_SSM(theta::Array{Float64,1}, y::Array{Float64,1}, valpha::Float64,
    vlambda::Float64, vpsi::Float64)
    n = length(y)
    x = theta[1:n]
    alpha = theta[n+1]
    lambda = theta[n+2]
    psi = theta[n+3]
    phi = (exp(psi)-1)/(exp(psi)+1)
    L = (0.5log(1-phi^2) - 0.5(1-phi^2)*x[1]^2 - 0.5sum((x[2:n] - phi*x[1:(n-1)]).^2) 
      - 0.5n*lambda - 0.5exp(alpha)*sum(x) + 0.5sum(y) - 0.5sum(exp(y - lambda - exp(alpha)*x)) 
      - 0.5lambda^2/vlambda - 0.5alpha^2/valpha - 0.5psi^2/vpsi)
return L
end

# logh_SSM(theta, y, valpha, vlambda, vpsi)

# compute gradient of log h(theta) for SMM (restrict phi as positive)
function gradlogh_SSM(theta::Array{Float64,1}, y::Array{Float64,1}, valpha::Float64, 
    vlambda::Float64, vpsi::Float64)
    n = length(y)
    x = theta[1:n]
    alpha = theta[n+1]
    lambda = theta[n+2]
    psi = theta[n+3]
    phi = 1/(exp(-psi)+1)
    g = zeros(n+3)
    s1 = exp(alpha)
    s2 = x[2:n] - phi*x[1:(n-1)] 
    s3 = exp(y - lambda - s1*x)
    g[1:n] = 0.5s1*(s3 - 1)
    g[2:n] -= s2
    g[1:(n-1)] += phi*s2
    g[1] -= x[1]*(1-phi^2)
    g[n+1] = -0.5s1*sum(x) + 0.5s1*sum(x.*s3) - alpha/valpha
    g[n+2] = -0.5n +0.5sum(s3) - lambda/vlambda
    g[n+3] = (sum(s2.*x[1:(n-1)]) - phi/(1-phi^2) + phi*x[1]^2)*exp(psi)/(exp(psi)+1)^2 - psi/vpsi
    return g 
end


# GLMM: I, J give indices of non-zero elements in T
function createIJ(n::Int64, p::Int64, k::Int64)
    f = Int(0.5p*(p+1))                                # no. of non-zero elements in each T_{ii} for i =1,...,n
    e = Int(0.5(k+f)*(k+f+1))                          # no. of non-zero elements in T_{n+1,n+1}
    d = Int(n*p + k + f)                               # length of theta = [b_1, ...,b_n, beta, zeta]
    nz = Int(f*n + e + (k+f)*p*n)                      # total no. of non-zero elements
    I = Array(Int64,nz)                                # row index
    J = Array(Int64,nz)                                # column index
    (rb, cb) = findn(tril(ones(p,p)).==1)              # row indices and column indices  
    (reta, ceta) = findn(tril(ones(k+f,k+f)).==1)      # row indices and column indices  
    for i in 1:n
        I[(i-1)*f+1:i*f] = rb + p*(i-1)
        J[(i-1)*f+1:i*f] = cb + p*(i-1)
    end
    I[f*n+1:f*n+e] = reta + p*n
    J[f*n+1:f*n+e] = ceta + p*n
    I[f*n+e+1:end] = repeat( n*p+1:d, outer=n*p)
    J[f*n+e+1:end] = repeat( 1:n*p, inner=k+f)
    return I, J
end

# SSM: I, J give indices of non-zero elements in T (n: length of time series)
function createIJ_SSM(n::Int64)
    nz = 5(n+1)
    I = Array(Int64,nz)
    J = Array(Int64,nz)
    I[1:(2n-1)] = [1:n; 2:n]
    J[1:(2n-1)] = [1:n; 1:(n-1)]
    I[(2n):end] = [repeat([n+1], outer=[n+1]); repeat([n+2], outer=[n+2]); repeat([n+3], outer=[n+3])]
    J[(2n):end] = [1:(n+1); 1:(n+2); 1:(n+3)]
    return I, J
end

# I, J give indices of all lower triangle elements
function createIJ_ltri(d::Int64)
    nz = Int(d*(d+1)/2)                                # total no. of non-zero elements
    I = Array(Int64,nz)
    J = Array(Int64,nz)
    l = 0
    for i in 1:d
        for j in 1:i
            l += 1
            I[l] = i
            J[l] = j
        end
    end
    return I, J
end
  

# Algorithm 2a (diagonals of T positive) stepsize: adadelta
function GLMMAlg2a(link::String, y::Array{Array{Float64,1},1}, X::Array{Array{Float64,2},1}, 
    Z::Array{Array{Float64,2},1}, rho::Float64, F::Int64, tol::Int64)

    N = 1e5                                  # max no. of iterations
    n = length(y)                            # no. of subjects
    k = size(X[1],2)                         # length of beta
    p = size(Z[1],2)                         # length of b_i
    d = Int(n*p + k + 0.5p*(p+1))            # length of theta = [b_1, ...,b_n, beta, zeta]
    id = find(tril(ones(p,p)).==1)           # indices of lower-triangle elements
    Idiag = diagm(ones(p))[id]               # p*(p+1)/2 x 1 vector, 1 if diagonal element
    Idiag = round(Int, Idiag)
    Idiag = round(Bool, Idiag)
    vbeta0 = 100.0
    vzeta0 = 100.0
    I, J = createIJ(n, p, k)
    nz = length(I)
    Tdiag = find(I.==J)
    mu = zeros(d)
    T = speye(d)                                     # print(full(T)) to see the dense matrix
    Tp = zeros(Float64,nz)                           # take log of diagonals of T
    par = Array(Float64, Int(N/F), 2*k+p*(p+1)+1)
    it_save = 0
    srand(1)
    epsilon = 1.0e-6
    Egmu = zeros(d)
    Ecmu = zeros(d)
    EgTp = zeros(nz)
    EcTp = zeros(nz)
    LBavg = 0.0
    it = 0
    cnt = 0
    maxLB = -Inf

    while (it<N) & (cnt<tol)
        it += 1     
        s = randn(d)
        g = T*s
        Ts =  A_ldiv_B!(factorize(T'), s)
        theta = mu + Ts
        LBavg += (logh_GLMM(link, theta, y, X, Z, Idiag, vbeta0, vzeta0) + 0.5d*log(2pi) - log(det(T)) + 0.5vecdot(s,s))/F

        g += gradlogh_GLMM(link, theta, y, X, Z, id, Idiag, vbeta0, vzeta0)
        Egmu = rho*Egmu + (1-rho)*g.^2
        cmu = sqrt(Ecmu + epsilon) ./ sqrt(Egmu + epsilon) .* g
        Ecmu = rho*Ecmu + (1-rho)*cmu.^2
        mu += cmu
        
        Tg = A_ldiv_B!(factorize(T), g)
        Tpgrad = -Ts[I].*Tg[J]
        Tpgrad[Tdiag] = diag(T).*Tpgrad[Tdiag]     
        EgTp = rho*EgTp + (1-rho)*Tpgrad.^2
        cTp = sqrt(EcTp + epsilon) ./ sqrt(EgTp + epsilon) .* Tpgrad
        EcTp = rho*EcTp + (1-rho)*cTp.^2
        Tp += cTp
        T = sparse(I,J,Tp,d,d)
        T[diagind(T)] = exp(Tp[Tdiag])
        
        if mod(it,F)==0
            if (LBavg < maxLB) 
                cnt +=1
            elseif (LBavg > maxLB)
                cnt = 0
            end
            maxLB = max(maxLB, LBavg)
            println("Iteration:",it," LB=", round(LBavg,3), " maxLB=", round(maxLB,3), " cnt=", cnt) 
            it_save += 1
            par[it_save,:] = [mu[n*p+1:end]; diag(T)[n*p+1:end]; LBavg] 
            LBavg = 0.0 
        end 
        
        if mod(it,10000)==0
            println(round([mu[n*p+1:end]; diag(T)[n*p+1:end]],2))  
        end             
    end 
    par = par[1:it_save,:] 
    VA = Dict("T" => T, "mu" => mu, "par" => par)
    return VA  
end





# Algorithm 1 (L diagonal positive, closed form) stepsize: adadelta
function GLMMAlg1diagb(link::String, y::Array{Array{Float64,1},1}, X::Array{Array{Float64,2},1}, 
    Z::Array{Array{Float64,2},1}, rho::Float64, F::Int64, tol::Int64)

    N = 1e5                                  # max no. of iterations
    n = length(y)                            # no. of subjects
    k = size(X[1],2)                         # length of beta
    p = size(Z[1],2)                         # length of b_i
    d = Int(n*p + k + 0.5p*(p+1))            # length of theta = [b_1, ...,b_n, beta, zeta]
    id = find(tril(ones(p,p)).==1)           # indices of lower-triangle elements
    Idiag = diagm(ones(p))[id]               # p*(p+1)/2 x 1 vector, 1 if diagonal element
    Idiag = round(Int, Idiag)
    Idiag = round(Bool, Idiag)
    vbeta0 = 100.0
    vzeta0 = 100.0
    mu = zeros(d)
    L = ones(d)
    Lp = zeros(d)                             # L is a diagonal matrix, only store diagonals
    par = Array(Float64, Int(N/F), 2*k+p*(p+1)+1)
    it_save = 0
    srand(1)
    epsilon = 1.0e-6
    Egmu = zeros(d)
    Ecmu = zeros(d)
    EgLp = zeros(d)
    EcLp = zeros(d)
    LBavg = 0.0
    it = 0
    cnt = 0
    maxLB = -Inf

    while (it<N) & (cnt<tol)
        it += 1          
        s = randn(d)
        theta = mu + L.*s
        LBavg += (logh_GLMM(link, theta, y, X, Z, Idiag, vbeta0, vzeta0) + 0.5d*log(2pi) + log(sum(L)) + 0.5vecdot(s,s))/F

        g = gradlogh_GLMM(link, theta, y, X, Z, id, Idiag, vbeta0, vzeta0)
        Egmu = rho*Egmu + (1-rho)*g.^2
        cmu = sqrt(Ecmu + epsilon) ./ sqrt(Egmu + epsilon) .* g
        Ecmu = rho*Ecmu + (1-rho)*cmu.^2
        mu += cmu
        
        Lpgrad = L.*g.*s + 1
        EgLp = rho*EgLp + (1-rho)*Lpgrad.^2
        cLp = sqrt(EcLp + epsilon) ./ sqrt(EgLp + epsilon) .* Lpgrad
        EcLp = rho*EcLp + (1-rho)*cLp.^2
        Lp += cLp
        L = exp(Lp)
        
        if mod(it,F)==0
            if (LBavg < maxLB) 
                cnt +=1
            elseif (LBavg > maxLB)
                cnt = 0
            end
            maxLB = max(maxLB, LBavg)
            it_save += 1
            par[it_save,:] = [mu[n*p+1:end]; L[n*p+1:end]; LBavg] 
            println("Iteration:",it," LB=", round(LBavg,3), " maxLB=", round(maxLB,3), " cnt=", cnt)  
            LBavg = 0.0
        end 

        if mod(it,1e4)==0
            println(round([mu[n*p+1:end]; L[n*p+1:end]],2))  
        end   
    end  
    par = par[1:it_save,:]
    VA = Dict("L" => L, "mu" => mu, "par" => par)
    return VA  
end


# Algorithm 1b (L full matrix, diagonal positive) closed form, stepsize: adadelta
function GLMMAlg1b(link::String, y::Array{Array{Float64,1},1}, X::Array{Array{Float64,2},1}, 
    Z::Array{Array{Float64,2},1}, rho::Float64, F::Int64, tol::Int64)

    N = 1e5                                  # max no. of iterations
    n = length(y)                            # no. of subjects
    k = size(X[1],2)                         # length of beta
    p = size(Z[1],2)                         # length of b_i
    d = Int(n*p + k + 0.5p*(p+1))            # length of theta = [b_1, ...,b_n, beta, zeta]
    id = find(tril(ones(p,p)).==1)           # indices of lower-triangle elements
    Idiag = diagm(ones(p))[id]               # p*(p+1)/2 x 1 vector, 1 if diagonal element
    Idiag = round(Int, Idiag)
    Idiag = round(Bool, Idiag)
    vbeta0 = 100.0
    vzeta0 = 100.0
    I, J = createIJ_ltri(d)
    nz = length(I)
    Ldiag = find(I.==J)
    mu = zeros(d)
    L = speye(d)                             # L is a lower triangular matrix
    Lp = zeros(nz)
    par = Array(Float64, Int(N/F), 2*k+p*(p+1)+1)
    it_save = 0
    srand(1)
    epsilon = 1.0e-6
    Egmu = zeros(d)
    Ecmu = zeros(d)
    EgLp = zeros(nz)
    EcLp = zeros(nz)
    LBavg = 0.0
    it = 0
    cnt = 0
    maxLB = -Inf

    while (it< N) & (cnt<tol)
        it += 1          
        s = randn(d)
        theta = mu + L*s
        LBavg += (logh_GLMM(link, theta, y, X, Z, Idiag, vbeta0, vzeta0) + 0.5d*log(2pi) + log(det(L)) + 0.5vecdot(s,s))/F

        g = gradlogh_GLMM(link, theta, y, X, Z, id, Idiag, vbeta0, vzeta0)
        Egmu = rho*Egmu + (1-rho)*g.^2
        cmu = sqrt(Ecmu + epsilon) ./ sqrt(Egmu + epsilon) .* g
        Ecmu = rho*Ecmu + (1-rho)*cmu.^2
        mu += cmu
        
        Lpgrad = g[I].*s[J] 
        Lpgrad[Ldiag] = diag(L).*Lpgrad[Ldiag] + 1
        EgLp = rho*EgLp + (1-rho)*Lpgrad.^2
        cLp = sqrt(EcLp + epsilon) ./ sqrt(EgLp + epsilon) .* Lpgrad
        EcLp = rho*EcLp + (1-rho)*cLp.^2
        Lp += cLp
        L = sparse(I,J,Lp,d,d)
        L[diagind(L)] = exp(Lp[Ldiag])
       
        if mod(it,F)==0
           if (LBavg < maxLB) 
                cnt +=1
            elseif (LBavg > maxLB)
                cnt = 0
            end
            maxLB = max(maxLB, LBavg)
            println("Iteration:",it," LB=", round(LBavg,3), " maxLB=", round(maxLB,3), " cnt=", cnt) 
            it_save += 1
            par[it_save,:] = [mu[n*p+1:end]; diag(L)[n*p+1:end]; LBavg] 
            LBavg = 0.0 
        end 

        if mod(it,1e4)==0
            println(round([mu[n*p+1:end]; diag(L)[n*p+1:end]],2))  
        end        
    end  
    par = par[1:it_save,:]
    VA = Dict("L" => L, "mu" => mu, "par" => par)
    return VA  
end


# Algorithm 2 (diagonals of T positive) stepsize: adadelta
function SSMAlg2a(y::Array{Float64,1}, rho::Float64, F::Int64, tol::Int64)

    N = 1e5                             # max no. of iterations
    n = length(y)                       # no. of subjects
    d = Int(n+3)                        # length of theta = [x_1, ...,x_n, alpha,lambda,psi]
    valpha = 10.0
    vlambda = 10.0
    vpsi = 10.0
    I, J = createIJ_SSM(n)
    nz = length(I)
    Tdiag = find(I.==J)
    mu = zeros(d)
    # T = speye(d)                        # print(full(T)) to see the dense matrix
    T = eye(d)                        # print(full(T)) to see the dense matrix
    Tp = zeros(nz)                      # take log of diagonals of T
    par = Array(Float64, Int(N/F), 7)
    it_save = 0
    srand(1)
    epsilon = 1.0e-6
    Egmu = zeros(Float64,d)
    Ecmu = zeros(d)
    EgTp = zeros(nz)
    EcTp = zeros(nz)
    LBavg = 0.0
    it = 0
    cnt = 0
    maxLB = -Inf

    while (it < N) & (cnt<tol)
        it += 1             
        s = randn(d)
        g = T*s
        Ts =  A_ldiv_B!(factorize(T'), s)
        theta = mu+Ts
        LBavg += (logh_SSM(theta, y, valpha, vlambda, vpsi) + 0.5d*log(2pi) - log(det(T)) + 0.5vecdot(s,s))/F

        g += gradlogh_SSM(theta, y, valpha, vlambda, vpsi)
        Egmu = rho*Egmu + (1-rho)*g.^2
        cmu = sqrt(Ecmu + epsilon) ./ sqrt(Egmu + epsilon) .* g
        Ecmu = rho*Ecmu + (1-rho)*cmu.^2
        mu += cmu
        
        Tg = A_ldiv_B!(factorize(T), g)
        Tpgrad = -Ts[I].*Tg[J]
        Tpgrad[Tdiag] = diag(T).*Tpgrad[Tdiag]    
        EgTp = rho*EgTp + (1-rho)*Tpgrad.^2
        cTp = sqrt(EcTp + epsilon) ./ sqrt(EgTp + epsilon) .* Tpgrad
        EcTp = rho*EcTp + (1-rho)*cTp.^2
        Tp += cTp
        T = sparse(I,J,Tp,d,d)
        T[diagind(T)] = exp(Tp[Tdiag])
 
        if mod(it,F)==0
           if (LBavg < maxLB) 
                cnt +=1
            elseif (LBavg > maxLB)
                cnt = 0
            end
            maxLB = max(maxLB, LBavg)
            println("Iteration:",it," LB=", round(LBavg,3), " maxLB=", round(maxLB,3), " cnt=", cnt) 
            it_save += 1
            par[it_save,:] = [mu[n+1:end]; diag(T)[n+1:end]; LBavg] 
            LBavg = 0.0 
        end 
        if mod(it,1e4)==0
            println(round([mu[n+1:end]; diag(T)[n+1:end]],2))  
        end      
    end  
    par = par[1:it_save,:]
    VA = Dict("T" => T, "mu" => mu, "par" => par)
    return VA  
end



# Algorithm 1diagb (L diagonal positive) stepsize: adadelta
function SSMAlg1diagb(y::Array{Float64,1}, rho::Float64, F::Int64, tol::Int64)

    N = 1e5                             # max no. of iterations
    n = length(y)                       # no. of subjects
    d = Int(n+3)                        # length of theta = [x_1, ...,x_n, alpha,phi,psi]
    valpha = 10.0
    vlambda = 10.0
    vpsi = 10.0
    mu = zeros(d)
    L = ones(d)
    Lp = zeros(d)                       # L is a diagonal matrix, only store diagonals
    par = Array(Float64, Int(N/F), 7)
    it_save = 0
    srand(1)
    epsilon = 1.0e-6
    Egmu = zeros(d)
    Ecmu = zeros(d)
    EgLp = zeros(d)
    EcLp = zeros(d)
    LBavg = 0.0
    it = 0
    cnt = 0
    maxLB = -Inf

    while (it < N) & (cnt<tol)
        it += 1      
        s = randn(d)
        theta = mu + L.*s
        LBavg += (logh_SSM(theta, y, valpha, vlambda, vpsi) + 0.5d*log(2pi) + log(sum(L)) + 0.5vecdot(s,s))/F

        g = gradlogh_SSM(theta, y, valpha, vlambda, vpsi)
        Egmu = rho*Egmu + (1-rho)*g.^2
        cmu = sqrt(Ecmu + epsilon) ./ sqrt(Egmu + epsilon) .* g
        Ecmu = rho*Ecmu + (1-rho)*cmu.^2
        mu += cmu
        
        Lpgrad = L.*g.*s + 1
        EgLp = rho*EgLp + (1-rho)*Lpgrad.^2
        cLp = sqrt(EcLp + epsilon) ./ sqrt(EgLp + epsilon) .* Lpgrad
        EcLp = rho*EcLp + (1-rho)*cLp.^2
        Lp += cLp
        L = exp(Lp)
 
        if mod(it,F)==0
           if (LBavg < maxLB) 
                cnt +=1
            elseif (LBavg > maxLB)
                cnt = 0
            end
            maxLB = max(maxLB, LBavg)
            println("Iteration:",it," LB=", round(LBavg,3), " maxLB=", round(maxLB,3), " cnt=", cnt) 
            it_save += 1
            par[it_save,:] = [mu[n+1:end]; L[n+1:end]; LBavg] 
            LBavg = 0.0 
        end     
        if mod(it,1e4)==0
            println(round([mu[n+1:end]; L[n+1:end]],2))  
        end
    end  
    par = par[1:it_save,:]
    VA = Dict("L" => L, "mu" => mu, "par" => par)
    return VA  
end



# Algorithm 1b (L full matrix, diagonal positive) stepsize: adadelta
function SSMAlg1b(y::Array{Float64,1}, rho::Float64, F::Int64, tol::Int64)

    N = 1e5                           # max no. of iterations
    n = length(y)                     # no. of subjects
    d = Int(n+3)                      # length of theta = [b_1, ...,b_n, beta, zeta]
    valpha = 10.0
    vlambda = 10.0
    vpsi = 10.0
    I, J = createIJ_ltri(d)
    nz = length(I)
    Ldiag = find(I.==J)
    mu = zeros(d)
    L = speye(d)                      # L is a lower triangular matrix
    Lp = zeros(nz)
    par = Array(Float64, Int(N/F), 7)
    it_save = 0
    srand(1)
    epsilon = 1.0e-6
    Egmu = zeros(d)
    Ecmu = zeros(d)
    EgLp = zeros(nz)
    EcLp = zeros(nz)
    LBavg = 0.0
    it = 0
    cnt = 0
    maxLB = -Inf

    while (it < N) & (cnt < tol)
        it += 1      
        s = randn(d)
        theta = mu + L*s
        LBavg += (logh_SSM(theta, y, valpha, vlambda, vpsi) + 0.5d*log(2pi) + log(det(L)) + 0.5vecdot(s,s))/F

        g = gradlogh_SSM(theta, y, valpha, vlambda, vpsi)
        Egmu = rho*Egmu + (1-rho)*g.^2
        cmu = sqrt(Ecmu + epsilon) ./ sqrt(Egmu + epsilon) .* g
        Ecmu = rho*Ecmu + (1-rho)*cmu.^2
        mu += cmu
        
        Lpgrad = g[I].*s[J] 
        Lpgrad[Ldiag] = diag(L).*Lpgrad[Ldiag] + 1
        EgLp = rho*EgLp + (1-rho)*Lpgrad.^2
        cLp = sqrt(EcLp + epsilon) ./ sqrt(EgLp + epsilon) .* Lpgrad
        EcLp = rho*EcLp + (1-rho)*cLp.^2
        Lp += cLp
        L = sparse(I,J,Lp,d,d)
        L[diagind(L)] = exp(Lp[Ldiag])
 
        if mod(it,F)==0
           if (LBavg < maxLB) 
                cnt +=1
            elseif (LBavg > maxLB)
                cnt = 0
            end
            maxLB = max(maxLB, LBavg)
            println("Iteration:",it," LB=", round(LBavg,3), " maxLB=", round(maxLB,3), " cnt=", cnt) 
            it_save += 1
            par[it_save,:] = [mu[n+1:end]; diag(L)[n+1:end]; LBavg] 
            LBavg = 0.0 
        end     
        if mod(it,1e4)==0
            println(round([mu[n+1:end]; diag(L)[n+1:end]],2))  
        end
    end 
    par = par[1:it_save,:] 
    VA = Dict("L" => L, "mu" => mu, "par" => par)
    return VA  
end

