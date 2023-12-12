

clear;
N = 30;
rate1 = 0.1:0.1:0.9;
rate2 = 0.9:-0.1:0.1;
n_case = fix(N*rate1);
n_control = fix(N*rate2);
T = 1000;
D = [1,1.5,2,2.5];
alpha = 0.5;
beta = 1;
for d = 1:length(D)
    for i = 1:length(n_control)
        for j = 1:T
            X{d,i,j} = gamrnd(alpha,beta,[1,n_case(i)]);
            Y{d,i,j} = gamrnd(D(d),beta,[1,n_control(i)]);
        end
    end
end

save data_rate X Y N n_case n_control D;

