function [power,mu_1,sigma_1] = cal_power_2(X,Y,Z)
alpha = 0.05;
X = X(:);
Y = Y(:);
Z = Z(:);
XZ = [X;Z];
YZ = [Y;Z];

m = length(X);
n = length(Y);
l = length(Z);
N = m + n + l;
s = n + l;
t = m + l;

sum_1 = 0;
for i =1:m
    for j = 1:s
        sum_1 = sum_1 + (X(i) > YZ(j));
    end
end
theta_1 = sum_1/m/s;

sum_2 = 0;
for i =1:n
    for j = 1:t
        sum_2 = sum_2 + (Y(i) > XZ(j));
    end
end
theta_2 = sum_2/n/t;

[F_X,F_YZ]=my_empirical_dist(X,YZ);
[G_Y,G_XZ]=my_empirical_dist(Y,XZ);
[U_YZ,U_X]=my_empirical_dist(YZ,X);
[V_XZ,V_Y]=my_empirical_dist(XZ,Y);

Var_F_YZ = var(F_YZ);
Var_G_XZ = var(G_XZ);
Var_U_X = var(U_X);
Var_V_Y = var(V_Y);

Cov_UX_GX = cov(U_X,G_XZ(1:m));
Cov_FY_VY = cov(V_Y,F_YZ(1:n));
Cov_FZ_VZ = cov(F_YZ(n+1:end),G_XZ(m+1:end));
Cov_1_2 = - ((n+l) * Cov_UX_GX + (m+l) * Cov_FY_VY - l * Cov_FZ_VZ);

Var_1 = (n+l)^2 / m * Var_U_X + (n+l) * Var_F_YZ;
Var_2 = (m+l)^2 / n * Var_V_Y + (m+l) * Var_G_XZ;

mu_0 = 0;
all = [X; Y; Z];
if length(unique(all)) == (m + n + l)
    sigma_0 = (m + n + l + 1) / 12 * (( n + l) / m + ( m + l) / n + 2);
else
    tall = tabulate(all);
    ind = find(tall(:,2) > 1);
    sm = 0;
    for h = 1:length(ind)
        sm = sm + tall(ind(h),2)^3 - tall(ind(h),2);
    end
    sigma_0 = ((N + 1) / 12 - sm / 12 / N / (N - 1)) * ((n + l) / m + (m + l) / n + 2);
end

mu_1 = (n + l) * theta_1 + (m + 1) / 2 - (m + l) * theta_2 - (n + 1) / 2;
sigma_1 = Var_1 + Var_2 - 2 * Cov_1_2(1,2);

delta = (mu_1 - mu_0) / sqrt(sigma_0);
sigma = sqrt(sigma_1 / sigma_0);    

u_1 = norminv(1 - alpha /2);
Phi_1 = normcdf((u_1 - delta) / sigma);

u_2 = norminv(alpha /2);
Phi_2 = normcdf((u_2 - delta) / sigma);

% power_1 = 1 - 2 * min(1 - Phi_1, Phi_2);
power = 1 - (Phi_1 - Phi_2);



