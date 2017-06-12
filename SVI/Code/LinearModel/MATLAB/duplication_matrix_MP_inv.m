function D_plus = duplication_matrix_MP_inv(n)
% The Moore-Penrose inverse of the duplication matrix of order n.
D = full(DuplicationM(n));
D_plus = (D'*D)\D';
end
