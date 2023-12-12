function W_null = compute_ranksum_null(x,z,num_perm)

nx = length(x);
nz = length(z);

for i = 1: num_perm
    W_null(i) = sum(randperm(nx+nz,nx));
end

W_null = W_null';

end
