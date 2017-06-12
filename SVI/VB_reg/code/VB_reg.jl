using Distributions
using AutoDiffSource,ReverseDiffSource
using ReverseDiff
function logpdfn(y::Float64,mu::Float64,sigma2::Float64)
    return -0.5*log(2.*pi*sigma2)-0.5/sigma2*(y-mu)*(y-mu)
end

function logmvnpdf(y::Array{Float64,1},mu::Array{Float64,1},Sig::Array{Float64,2})
    p = size(mu,1)
    return -p/2.*log(2.*pi)-0.5*logdet(Sig)-0.5*sum((y-mu).*(Sig\(y-mu)))
end


function logpdfIG(x::Float64,A::Float64,B::Float64)
    return A*log(B)-lgamma(A)-(A+1.)*log(x)-B/x
end

function logh(parameters::Array{Float64,1},hyperparameters::Array{Float64,1},y::Array{Float64,1},vphi::Array{Float64,2},n::Int64,J::Int64)
    # extract parameters
    theta = parameters[1:J]
    psi   = parameters[J+1]
    sig2  = parameters[J+2]
    tau2  = parameters[J+3]
    gam   = parameters[J+4]

    # extract hyperparameters
    r_sig = hyperparameters[1]
    s_sig = hyperparameters[2]
    r_tau = hyperparameters[3]
    s_tau = hyperparameters[4]
    a     = hyperparameters[5]
    b     = hyperparameters[6]
    w0    = hyperparameters[7]

    res = logmvnpdf(y,vphi*theta,sig2*eye(n))+logpdfIG(sig2,r_sig/2.,s_sig/2.)+logpdfIG(tau2,r_tau/2.,s_tau/2.)+logpdfIG(psi,a,b)+log(w0)-w0*gam
    for j = 1:J
        res = res+logpdfn(theta[j],0.,sig2*tau2*exp(-j*gam))
    end
    return res
end
input = (Ω,hyperparameters,y,vphi,n,J)
∇logh = ReverseDiff.gradient(logh,input)

[:hyperparameters,:y,:vphi,:n,:J]
∇logh = rdiff(logh,init=(Array{Float64,1},Array{Float64,1},Array{Float64,1},Array{Float64,2},Int64,Int64),order=1,ignore=[:hyperparameters;:y;:vphi;:n;:J;])


fun1(x) = 2.*pi*x+sin(2.*pi*x)
n = 200;
x = rand(n);
fx = fun1.(x);
sig2_true = 3.4;
xmin = minimum(x);
xmax = maximum(x);
x = (x-xmin)./(xmax-xmin);
y = fx+randn(n)*sqrt(sig2_true);

J = 20;
ψ = 2.4;

vphi = sqrt(2.)*cos(A_mul_Bt(x,pi*linspace(1,J,J)));
r0_sig,s0_sig,r0_tau,s0_tau,a,b,w0 = 10.,10.,10.,10.,10.,10.,2.;
hyperparameters = (r0_sig,s0_sig,r0_tau,s0_tau,a,b,w0);
μ_post = zeros(J+4);
L_post  = eye(J+4);
LB = 0.;
LBC = [];
S = 20;
ϵ = 1.0e-16;
s_k = ones(J+4);
ρ = zeros(J+4);
α = 0.3;
t = 1.;
iterMax = 4000;
thetaContainer = zeros(J+4,iterMax);
for nIter = 1:iterMax
    z_η = randn(S,J+4);
    grad_μ = zeros(J+4);
    grad_L = zeros(J+4,J+4);
    for s = 1:S
        z = z_η[s,:];
        ζ = L_post*z+μ_post;

        θ = ζ[1:J];
        ψ = log(exp(ζ[J+1])+1.);
        σ2 = log(exp(ζ[J+2])+1.);
        τ2 = log(exp(ζ[J+3])+1.);
        γ  = log(exp(ζ[J+4])+1.);
        Ω  = [θ;ψ;σ2;τ2;γ];
        fx,∂Ω = δlogh(Ω,hyperparameters,y,vphi,n,J)
        ∇Ω = ∂Ω()[1];
        ∇_Tinv = [ones(J);exp(ζ[J+1])/(1.+exp(ζ[J+1]));exp(ζ[J+2])/(1.+exp(ζ[J+2]));exp(ζ[J+3])/(1.+exp(ζ[J+3]));exp(ζ[J+4])/(1.+exp(ζ[J+4]))];
        ∇_logdet = [zeros(J);1./(1.+exp(ζ[J+1]));1./(1.+exp(ζ[J+2]));1./(1.+exp(ζ[J+3]));1./(1.+exp(ζ[J+4]))];
        grad_μ += ∇Ω/S.*∇_Tinv+∇_logdet/S;
        grad_L += A_mul_Bt(grad_μ,z)+transpose(inv(L_post))/S;
    end
    if nIter == 0
        s_k = grad_μ .* grad_μ;
    else
        s_k = α*grad_μ+(1.-α)*s_k;
    end
    ρ = 0.1*nIter^(-0.5+ϵ)./(t+sqrt.(s_k));
    μ_post += ρ.*grad_μ;
    L_post += diagm(ρ)*grad_L;
    LB = 0.;
    Sig = A_mul_Bt(L_post,L_post);
    for s = 1:S
        ζ = L_post*randn(J+4)+μ_post;
        θ = ζ[1:J];
        ψ = log(exp(ζ[J+1])+1.);
        σ2 = log(exp(ζ[J+2])+1.);
        τ2 = log(exp(ζ[J+3])+1.);
        γ  = log(exp(ζ[J+4])+1.);
        Ω  = [θ;ψ;σ2;τ2;γ];
        LB += logh(Ω,hyperparameters,y,vphi,n,J)+sum(ζ[(J+1):end])-sum(log(1.+exp(ζ[(J+1):end])))-logmvnpdf(ζ,μ_post,Sig);
    end
    LB /= S;
    println("Iteration:",nIter,", LB=",round(LB,3))
    LBC = [LBC;LB];
    thetaContainer[:,nIter] = μ_post;
end

parameters::Array{Float64,1},hyperparameters::Array{Float64,1},y::Array{Float64,1},vphi::Array{Float64,2},n::Int64,J::Int64
# ∇logh = rdiff(logh,order=1,parameters=(Vector{Float64},),ignore=[:hyperparameters;:y;:vphi;:n;:J])
∇logh = rdiff(logh,order=1,(Vector{Float64},Vector{Float64},Vector{Float64},Matrix{Float64},Int64,Int64))


