
clear;

load data_eq;

nperm = 1000;
T = size(X,3); %trial times
K = 1:0.5:10; %auxiliary sample size

for d = 3:3%1:length(D) %difference
    for i = 3:3%1:length(N) %case/control sample
        [d i]
        n = N(i);
        for j = 1:T %trial times
            x = X{d,i,j};y = Y{d,i,j};
            x = x(:);   % ensure columns
            y = y(:);
%             pwmw(i,j) = ranksum(x,y);
%             [~,pt(i,j)] = ttest2(x,y,'tail','both');
            
            allmin = min([x;y]);
            allmax = max([x;y]);

            for k = 1:length(K) %auxiliary sample size
                %                 [d i k]
                z = [normrnd(0,1,[1,floor(K(k)*n)]) normrnd(D(d),1,[1,floor(K(k)*n)])];
                idx = intersect(find(z>=allmin),find(z<=allmax));
                z = z(idx);
                x = x'; y = y'; z= z';
                [pwmwa1(i,k,j),pwmwa2(i,k,j)]=wmwa(x,y,z,nperm);
            end
        end
        
        for k = 1:length(K)
            power_wmwa1(d,i,k) = length(find(pwmwa1(i,k,:)<=0.05))/T;
            power_wmwa2(d,i,k) = length(find(pwmwa2(i,k,:)<=0.05))/T;

        end
    end
end

save result_k_eq power_wmwa1 power_wmwa2;