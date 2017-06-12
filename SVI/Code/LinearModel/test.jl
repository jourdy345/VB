include("/Users/daeyounglim/Desktop/Github/VB/SVI/Code/julia/SVI_linearRegressionCV.jl")
cd("/Users/daeyounglim/Desktop/Github/VB/SVI/Code/julia/")
# writedlm("ig.txt",x)

n        = 100;
p        = 4;
betaTrue = [0.3,10.,2.,6.];
X        = randn(100,4)
y        = X*betaTrue+randn(100)*1.3
# VA       = BVI(y,X)


n,p  = size(X)
sig0 = eye(p)
mu0  = zeros(p)
A    = 10.
B    = 10.

# initialize variational parameters
sigq    = randn(p,p)
sigq    = At_mul_B(sigq,sigq)
origsigq = copy(sigq)
invsigq = inv(sigq)
muq     = zeros(p)
Aq      = 10.
Bq      = 10.

# initialize auxiliary components
Dmat   = duplicationMatrix(p)
Dplus = pinv(Dmat)
S      = 40
nIter  = 50
it     = 0
cnt    = 0
lb     = 0.
lbc    = zeros(0)
muqc   = zeros(p,nIter)

# natural parameters
gaussLambda1 = zeros(p)
A_ldiv_B!(gaussLambda1,factorize(sigq),muq)
# gaussLambda1 = sigq\muq
gaussLambda2 = -0.5*At_mul_B(Dmat,invsigq[:])
gaussLambda  = [gaussLambda1;gaussLambda2]
invgamLambda = [Aq;Bq]

# intialize auxiliary components
gaussM       = zeros(length(gaussLambda))
igM          = zeros(2)
gaussGrad    = zeros(length(gaussLambda))
igGrad       = zeros(2)
rho          = 0.
gaussLambdaOld = zeros(length(gaussLambda))
muqOld         = zeros(p)
sigqOld        = zeros(p,p)
βSample        = zeros(p,S)
igSample       = zeros(S)
β              = zeros(p)
sigma2         = 1.
f_v            = 0.
logh_v         = 0.
logq_v         = 0.

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
    aux             = Smat\Mmat
    gaussInvFisher  = vcat(hcat(invsigq+At_mul_B(Mmat,aux),-aux.'),hcat(-aux,inv(Smat)))
    # gaussInvFisher  = [invsigq+At_mul_B(Mmat,aux) -aux';-aux inv(Smat)]
    invgamFisher    = [trigamma(Aq) -1/Bq;-1/Bq Aq/(Bq*Bq)]

    # get samples
    βSample  = rand(mvnrnd,S)
    igSample = rand(igrnd,S)
    # igSample = quantile(igrnd,rand(S))
    # igSample = rig(Aq,Bq,S)
    # for k = 1:S
    #     igSample[k] = rand_inverse_gamma(Aq,Bq)
    # end
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
    gaussGrad = zeros(length(gaussLambda))
    igGrad    = zeros(2)
    lb        = 0.
    for j = 1:S
        β   = βSample[:,j]
        sigma2 = igSample[j]
        logh_v = logh(n,p,y,X,β,mu0,sig0,sigma2,A,B)
        logq_v = logq(p,β,muq,invsigq,sigma2,Aq,Bq)
        f_v    = logh_v-logq_v
        lb    += f_v
        gaussGrad += ([β-muq;vech(A_mul_Bt(β,β)-sigq-A_mul_Bt(muq,muq))]-gaussM)*f_v
        igGrad    += ([-log(sigma2)+log(Bq)-digamma(Aq);-1/sigma2+Aq/Bq]-igM)*f_v
    end
    gaussGrad /= (S-1)
    igGrad    /= (S-1)

    # update
    gaussLambda  = (1-rho)*gaussLambda+rho*(gaussInvFisher*gaussGrad)
    invgamLambda = (1-rho)*invgamLambda+rho*(invgamFisher\igGrad)
    gaussLambda1 = gaussLambda[1:p]
    gaussLambda2 = gaussLambda[(p+1):end]
    sigq         = -0.5*inv(reshape(At_mul_B(Dplus,gaussLambda2),p,p))
    muq          = sigq*gaussLambda1
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


# x = map(igSample) do ig
#     return [-log(ig)+log(Bq)-digamma(Aq);-1/ig+Aq/Bq]
# end