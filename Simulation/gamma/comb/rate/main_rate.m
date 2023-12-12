
clear;

load data_rate;

nperm = 1000;
T = size(X,3); %trial times
K = [1,2]; %auxiliary sample size

for d = 1:length(D) %difference
    tic
    for i = 1:length(n_control) %case/control sample
        [d i]
        m = n_case(i);
        n = n_control(i);
        for j = 1:T %trial times
            x = X{d,i,j};y = Y{d,i,j};
            x = x(:);   % ensure columns
            y = y(:);
            nx = length(x); ny = length(y);
            V=tiedrank([x;y]); % rank transformation of concatenated vectors
            pwmw(i,j) = ranksum(x,y);
            [~,pt(i,j)] = ttest2(x,y);
            [~,ptR(i,j)]=ttest2(V(1:nx),V(nx+1:end)); % two sample t-test after rank transformation(t-testR)
            [~,pwc(i,j)]=ttest2(x,y,[],[],'unequal'); % two-sample t-test using unequal variances option (Welch test)
            
            allmin = min([x;y]);
            allmax = max([x;y]);
            
            for k = 1:length(K) %auxiliary sample size
                for g = 1:11%1:11
                    z = [gamrnd(0.5,1,[1,fix(K(k)*N*(g*0.1-0.1))]) gamrnd(D(d),1,[1,K(k)*N-fix(K(k)*N*(g*0.1-0.1))])];
                    idx = intersect(find(z>=allmin),find(z<=allmax));
                    z = z(idx);
                    z = z(:);
                    [pwmwa1(d,i,k,g,j),pwmwa2(d,i,k,g,j)]=wmwa_all5(x,y,z,nperm);
                end
            end
        end
        
        for k = 1:length(K)
            for g = 1:11
                power_wmwa1(d,i,k,g) = length(find(pwmwa1(d,i,k,g,:)<=0.05))/T;
                power_wmwa2(d,i,k,g) = length(find(pwmwa2(d,i,k,g,:)<=0.05))/T;
            end
        end
        power_wmw(d,i) = length(find(pwmw(i,:)<=0.05))/T;
        power_t(d,i) = length(find(pt(i,:)<=0.05))/T;
        power_tR(d,i) = length(find(ptR(i,:)<=0.05))/T;
        power_wc(d,i) = length(find(pwc(i,:)<=0.05))/T;
    end
    toc
end


save('rate_comb.mat', 'power_wmwa1', 'power_wmwa2', 'power_wmw', 'power_t', 'power_wc', 'power_tR', 'pwmw', 'pwmwa1', 'pwmwa2', 'pt', 'N', 'T', 'D', 'K');

