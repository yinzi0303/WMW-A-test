
clear;

load data_k_eq;

nperm = 1000;
T = size(X,3); %trial times
K = 1:0.5:10; %auxiliary sample size
sum = [];
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
            sum(i,j) = 0;
            
            for k = 1:length(K)%1:length(K) %auxiliary sample size
                
                z = [gamrnd(0.5,1,[1,floor(K(k)*n)]) gamrnd(D(d),1,[1,floor(K(k)*n)])];
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
        %         power_wmw(d,i) = length(find(pwmw(i,:)<=0.05))/T;
        %         power_t(d,i) = length(find(pt(i,:)<=0.05))/T;
    end
    
end
save result_k_equal power_wmwa1 power_wmwa2 N T D K;