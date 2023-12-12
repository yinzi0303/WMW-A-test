function [p1,p2] = wmwa_all5(x,y,z,num_perm)
nx = length(x);
ny = length(y);
nz = length(z);
x = x(:);   % ensure columns
y = y(:);
z = z(:);
alpha1 = nx/ny;

Wx_obs = compute_ranksum_obs(x,[y;z]);
Wy_obs = compute_ranksum_obs(y,[x;z]);
W_obs1 = Wx_obs - alpha1* Wy_obs;

% Wx_obs2 = compute_ranksum_obs(x,z);
% Wy_obs2 = compute_ranksum_obs(y,z);
% W_obs2 = Wx_obs2 - alpha1* Wy_obs2;


all = [x;y;z];
label(1:nx) = 1;
label(nx+1:nx+ny) = 2;
label(nx+ny+1:nx+ny+nz) = 3;


for k = 1:num_perm
    plabel(k,:) = label(randperm(nx+ny+nz));
    px = all(plabel(k,:)==1);
    py = all(plabel(k,:)==2);
    pz = all(plabel(k,:)==3);
    
    Wx_null(k) = compute_ranksum_obs(px,[py;pz]);
    Wy_null(k) = compute_ranksum_obs(py,[px;pz]);
    W_null_1(k) = Wx_null(k) - alpha1* Wy_null(k);
    
%     Wx_null2(k) = compute_ranksum_obs(px,pz);
%     Wy_null2(k) = compute_ranksum_obs(py,pz);
%     W_null_2(k) = Wx_null2(k) - alpha1*Wy_null2(k);
end

plo = sum( W_null_1 <= W_obs1) / num_perm ;
phi = sum( W_null_1 >= W_obs1) / num_perm ;
p_tail = min(plo,phi);
p1 = min(2*p_tail, 1);

p2 = sum( abs(W_null_1) >= abs(W_obs1) ) / num_perm ;

% plo = sum( W_null_2 <= W_obs2) / num_perm ;
% phi = sum( W_null_2 >= W_obs2) / num_perm ;
% p_tail = min(plo,phi);
% p3 = min(2*p_tail, 1);
% plo = sum( W_null_2 <= W_obs2) / num_perm ;
% phi = sum( W_null_2 >= W_obs2) / num_perm ;
% p_tail = min(plo,phi);
% p3 = min(2*p_tail, 1);
% 
% p4 = sum( abs(W_null_2) >= abs(W_obs2) ) / num_perm ;

%
% p5 = min(p1,p3);
% p6 = min(p2,p4);

% figure(1)
% histogram(W_null_1)
% figure(2)
% histogram(W_null_2)

end

