
clear;

load data_asy_m30;
nperm = 1000;
T = size(X,3); %trial times
% K = [2 5 10]; %auxiliary sample size
K = 2;
for d = 2:2%1:length(D) %difference

    for i = 1:length(N) %case/control sample
        [d i]
        n = N(i);
        for j = 1:T %trial times
            x = X{d,i,j};y = Y{d,i,j};
%             pwmw(i,j) = ranksum(x,y);
%             [~,pt(i,j)] = ttest2(x,y,'tail','both');
            
            allmin = min([x y]);
            allmax = max([x y]);
            
            for k = 1:length(K) %auxiliary sample size
%                 z = [gamrnd(0.5,1,[1,K(k)*n]) gamrnd(D(d),1,[1,K(k)*n])];
%                 z = [gamrnd(0.5,1,[1,K(k)*30]) gamrnd(D(d),1,[1,K(k)*n])];
                z = [nbinrnd(100,0.5,[1,K(k)*30]) nbinrnd(100,D(d),[1,K(k)*n])];
                idx = intersect(find(z>=allmin),find(z<=allmax));
                z = z(idx);
                x = x'; y = y'; z= z';
%                 pi_u = length(intersect(1:K(k)*n,idx))/(length(z)+n);
%                 pi_v = (length(intersect(1:K(k)*n,idx))+n)/(length(z)+n);
                [power(d,i,k,j),~,~]=cal_power_2(x,y,z);
%                 [power_2(d,i,k,j),~,~]=cal_power_3(x,y,z,pi_u,pi_v,d);
                [pwmwa1(i,k,j),pwmwa2(i,k,j)]=wmwa(x,y,z,nperm);
            end
        end
        
        for k = 1:length(K)
            power_wmwa_asy(d,i,k) = sum(power(d,i,k,:))/T;
%             power_wmwa_asy_2(d,i,k) = sum(power_2(d,i,k,:))/T;
            power_wmwa1(d,i,k) = length(find(pwmwa1(i,k,:)<=0.05))/T;
            power_wmwa2(d,i,k) = length(find(pwmwa2(i,k,:)<=0.05))/T;
        end

    end

end

save("result_asy_m30.mat", 'power_wmwa_asy', 'power_wmwa1', 'power_wmwa2', 'N', 'T', 'D', 'K');