clear;
N = [2 4 5 8 10 12 15 18 20 25];
T = 1000;
D = [0.48,0.45,0.42,0.4];
alpha = 0.5;
beta = 100;
for d = 1:length(D)
    for i = 1:length(N)
        for j = 1:T
            n = N(i);
            X{d,i,j} = nbinrnd(beta,alpha,[1,n]);
            Y{d,i,j} = nbinrnd(beta,D(d),[1,n]);
        end
    end
end

save data_k_eq X Y N D;



