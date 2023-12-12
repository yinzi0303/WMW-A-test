

clear;
N = [10 15 20 25 30 40 50 60 80 100];
T = 1000;
D = [0.5,1,1.5,2,2.5];
alpha = 0.5;
beta = 1;
for d = 1:length(D)
    for i = 1:length(N)
        for j = 1:T
            n = N(i);
            X{d,i,j} = gamrnd(alpha,beta,[1,5]);
            Y{d,i,j} = gamrnd(D(d),beta,[1,n]);
        end
    end
end

save data_k_ueq X Y N D;

