



clear;
N = 30;
rate1 = 0.1:0.1:0.9;
rate2 = 0.9:-0.1:0.1;
K = [1,2];
n_case = fix(N*rate1);
n_control = fix(N*rate2);
T = 1000;
D = [0.5,1];
alpha = 0;
beta = 1;
for d = 1:length(D)
    for i = 1:length(n_control)
        for j = 1:T
            X{d,i,j} = normrnd(alpha,beta,[1,n_case(i)]);
            Y{d,i,j} = normrnd(D(d),beta,[1,n_control(i)]);
        end
    end
end

save('data_rate.mat', 'X', 'Y', 'N', 'n_case', 'n_control', 'D');



