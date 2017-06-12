function [Z] = Zmat(X,S,lambda,n,m,d)
    Z = zeros(n,m);
    for i = 1:n
       for j = 1:m
          Z(i,j) = sum(S(j,:).*X(i,:).*lambda');
       end
    end
    Z = [cos(Z) sin(Z)];
end