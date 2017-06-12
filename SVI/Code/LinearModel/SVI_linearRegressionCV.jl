using Distributions
function vech(X::Array{Float64,2})
    # @assert issymmetric(X) "X is not symmetric"
    n = size(X,1)
    v = zeros(div(n*(n+1),2))
    for j = 1:n
        for i = j:n
            v[n*(j-1)-div((j-1)*(j-2),2)+i-j+1] = X[i,j]
        end
    end
    return v
end

# revert a vector to a lower-triangular matrix
function xpnd(v::Array{Float64,1},p::Int64)
    m = zeros(p,p)
    for j = 1:p
        for i = j:p
            m[i,j] = v[p*(j-1)-div((j-1)*(j-2),2)+i-j+1]
        end
    end
    return m
end



function duplicationMatrix(n::Int64)
    @assert n>=2 "argument n is less than 2"
    k = zeros(Int64,n,n)
    for j = 1:n
        for i = j:n
            k[i,j] = (j-1)*n+i-div(j*(j-1),2)
        end
    end
    index = 1:(n*(n+1)/2)
    k = Symmetric(k,:L)
    tmpVec = vec(k)
    t = length(tmpVec)
    s = length(index)
    res = zeros(t,s)
    for k = 1:t
        for u = 1:s
            if tmpVec[k] == index[u]
                res[k,u] = 1
            else
                res[k,u] = 0
            end
        end
    end
    return res
end


function logh(n::Int64,p::Int64,y::Array{Float64,1},X::Array{Float64,2},β::Array{Float64,1},
              mu0::Array{Float64,1},sig0::Array{Float64,2},sigma2::Float64,A::Float64,B::Float64)
    # likelihood   = MvNormal(X*β,sigma2*eye(n))
    # βprior       = MvNormal(mu0,sig0)
    # σ2prior      = InverseGamma(A,B)
    # v = logpdf(likelihood,y)+logpdf(βprior,β)+logpdf(σ2prior,sigma2)
    v = -0.5*n*log(2*pi*sigma2)-0.5/sigma2*sum((y-X*β).^2)-0.5*logdet(sig0)-0.5*p*log(2*pi)-0.5*sum((β-mu0).*(sig0\(β-mu0)))-(A+1)*log(sigma2)-B/sigma2+A*log(B)-lgamma(A)
    return v
end

function logq(p::Int64,β::Array{Float64,1},muq::Array{Float64,1},invsigq::Array{Float64,2},sigma2::Float64,Aq::Float64,Bq::Float64)
    # βvariational  = MvNormal(muq,sigq)
    # σ2variational = InverseGamma(Aq,Bq)
    # v = logpdf(βvariational,β)+logpdf(σ2variational,sigma2)
    v = -p*0.5*log(2*pi)+0.5*logdet(invsigq)-0.5*sum((β-muq).*(invsigq*(β-muq)))-(Aq+1)*log(sigma2)-Bq/sigma2+Aq*log(Bq)-lgamma(Aq)    
    return v
end

function grad(βSample::Array{Float64,2},muq::Array{Float64,1},sigq::Array{Float64,2},invsigq::Array{Float64,2},
              igSample::Array{Float64},Aq::Float64,Bq::Float64,gaussM::Array{Float64,1},igM::Array{Float64,1},
              mu0::Array{Float64,1},sig0::Array{Float64,2},A::Float64,B::Float64,n::Int64,p::Int64,y::Array{Float64,1},X::Array{Float64,2})
    S = length(igSample)
    gaussGrad = zeros(Int(p+p*(p+1)/2))
    igGrad    = zeros(2)
    lb        = 0
    for j = 1:S
        β   = βSample[:,j]
        sigma2 = igSample[j]
        logh_v = logh(n,p,y,X,β,mu0,sig0,sigma2,A,B)
        logq_v = logq(p,β,muq,invsigq,sigma2,Aq,Bq)
        f_v    = logh_v-logq_v
        lb    += f_v
        gaussGrad += ([β-muq;vech(A_mul_Bt(β,β)-sigq-A_mul_Bt(muq,muq))]-gaussM).*f_v
        igGrad    += ([-log(sigma2)+log(Bq)-digamma(Aq);-1/sigma2+Aq/Bq]-igM).*f_v
    end
    gaussGrad /= (S-1)
    igGrad    /= (S-1)

    return Dict("gaussGrad"=>gaussGrad,"igGrad"=>igGrad)
end

# function rMNorm(n::Int64,mu::Array{Float64,1},Sigma::Array{Float64,2})
#     p = length(mu)
#     X = randn(n,p)

# end

function BVI(y::Array{Float64,1},X::Array{Float64,2})
    # initialize hyperparameters for prior distributions
    n,p  = size(X)
    sig0 = eye(p)
    mu0  = zeros(p)
    A    = 10.
    B    = 10.

    # initialize variational parameters
    sigq    = eye(p)
    origsigq = copy(sigq)
    invsigq = inv(sigq)
    muq     = zeros(p)
    Aq      = 10.
    Bq      = 10.

    # initialize auxiliary components
    Dmat  = duplicationMatrix(p)
    Dplus = pinv(Dmat)
    S     = 200
    nIter = 50
    aux   = zeros(Int(p*(p+1)/2),p)
    lb    = 0.
    lbc   = zeros(0)
    muqc  = zeros(p,nIter)

    # natural parameters
    gaussLambda1 = sigq\muq
    gaussLambda2 = -0.5*At_mul_B(Dmat,vec(invsigq))
    gaussLambda  = [gaussLambda1;gaussLambda2]
    invgamLambda = [Aq;Bq]

    for i = 1:nIter
        # initialize samplers
        mvnrnd          = MvNormal(muq,sigq)
        igrnd           = InverseGamma(Aq,Bq)
        # store old values
        gaussLambdaOld  = copy(gaussLambda)
        muqOld          = copy(muq)
        sigqOld         = copy(sigq)
        invgamLambdaOld = copy(invgamLambda)
        AqOld           = Aq
        BqOld           = Bq
        rho             = 1/(5+i)
        Mmat            = 2*Dplus*(kron(muq,eye(p)))
        Smat            = 2*Dplus*A_mul_Bt(kron(sigq,sigq),Dplus)
        # aux             = Smat\Mmat
        A_ldiv_B!(aux,factorize(Smat),Mmat)
        gaussInvFisher  = vcat(hcat(invsigq+At_mul_B(Mmat,aux),-aux.'),hcat(-aux,inv(Smat)))
        # gaussInvFisher  = [invsigq+At_mul_B(Mmat,aux) -aux';-aux inv(Smat)]
        invgamFisher    = [trigamma(Aq) -1/Bq;-1/Bq Aq/(Bq*Bq)]

        # get samples
        βSample  = rand(mvnrnd,S)
        igSample = rand(igrnd,S)

        # compute control variate
        gaussM = zeros(length(gaussLambda))
        igM    = zeros(2)
        for k = 1:S
            β = βSample[:,k]
            gaussM += [β-muq;vech(A_mul_Bt(β,β)-sigq-A_mul_Bt(muq,muq))]
            sigma2 = igSample[k]
            igM    += [-log(sigma2)+log(Bq)-digamma(Aq);-1/sigma2+Aq/Bq]
        end
        gaussM /= S
        igM    /= S

        # compute gradients
        G = grad(βSample,muq,sigq,invsigq,igSample,Aq,Bq,gaussM,igM,mu0,sig0,A,B,n,p,y,X)
        global G
        gaussGrad = G["gaussGrad"]
        igGrad    = G["igGrad"]
        # gaussGrad = zeros(length(gaussLambda))
        # igGrad    = zeros(2)
        # lb        = 0
        # for j = 1:S
        #     β   = βSample[:,j]
        #     sigma2 = igSample[j]
        #     logh_v = logh(n,p,y,X,β,mu0,sig0,sigma2,A,B)
        #     logq_v = logq(p,β,muq,invsigq,sigma2,Aq,Bq)
        #     f_v    = logh_v-logq_v
        #     lb    += f_v
        #     gaussGrad += ([β-muq;vech(A_mul_Bt(β,β)-sigq-A_mul_Bt(muq,muq))]-gaussM).*f_v
        #     igGrad    += ([-log(sigma2)+log(Bq)-digamma(Aq);-1/sigma2+Aq/Bq]-igM).*f_v
        # end
        # gaussGrad /= (S-1)
        # igGrad    /= (S-1)

        # update
        gaussLambda  = (1-rho)*gaussLambda+rho*(gaussInvFisher*gaussGrad)
        invgamLambda = (1-rho)*invgamLambda+rho*(invgamFisher\igGrad)
        gaussLambda1 = gaussLambda[1:p]
        gaussLambda2 = gaussLambda[(p+1):end]
        sigq         = -0.5*inv(reshape(At_mul_B(Dplus,gaussLambda2),p,p))
        muq          = sigq*gaussLambda1
        println(sigq)
        println("sigqsym:",issymmetric(sigq))
        if !isposdef(sigq)
            gaussLambda = copy(gaussLambdaOld)
            sigq        = copy(sigqOld)
            muq         = copy(muqOld)
        else 
            invsigq     = inv(sigq)
        end
        Aq = invgamLambda[1]
        Bq = invgamLambda[2]
        if Aq < 0
            Aq = AqOld
        end
        if Bq < 0
            Bq = BqOld
        end
        invgamLambda = [Aq;Bq]
        lb /= S
        push!(lbc,lb)
        muqc[:,i] = muq
        println("Iteration:",i,", LB=",round(lb,3))
    end
    return Dict("muq"=>muq,"sigq"=>sigq,"lbc"=>lbc,"Aq"=>Aq,"Bq"=>Bq,"muqc"=>muqc,"origsigq"=>origsigq,"igSample"=>igSample,"βSample"=>βSample,"G"=>G)
end


