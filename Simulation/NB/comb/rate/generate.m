
clear;
N = 20;
rate1 = 0.1:0.1:0.9;
rate2 = 0.9:-0.1:0.1;
K = [1,2];
n_case = fix(N*rate1);
n_control = fix(N*rate2);
T = 1000;
D = [0.48,0.45];
alpha = 0.5;
beta = 100;
for d = 1:length(D)
    for i = 1:length(n_control)
        for j = 1:T
            X{d,i,j} = nbinrnd(beta,alpha,[1,n_case(i)]);
            Y{d,i,j} = nbinrnd(beta,D(d),[1,n_control(i)]);
        end
    end
end

save('data_rate_1.mat', 'X', 'Y', 'N', 'n_case', 'n_control', 'D');


clear;
N = 20;
rate1 = 0.1:0.1:0.9;
rate2 = 0.9:-0.1:0.1;
K = [1,2];
n_case = fix(N*rate1);
n_control = fix(N*rate2);
T = 1000;
D = [95,90];
alpha = 100;
beta = 0.5;
for d = 1:length(D)
    for i = 1:length(n_control)
        for j = 1:T
            X{d,i,j} = nbinrnd(alpha,beta,[1,n_case(i)]);
            Y{d,i,j} = nbinrnd(D(d),beta,[1,n_control(i)]);
        end
    end
end

save('data_rate_2.mat', 'X', 'Y', 'N', 'n_case', 'n_control', 'D');