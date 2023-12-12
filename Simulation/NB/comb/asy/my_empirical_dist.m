function [cumpr,cumpr_2]=my_empirical_dist(data, data_2)
if min(size(data))~=1
    error('data must be a vector')
end
n=length(data);
data_sort=sort(data);
[x,a]=unique(data_sort);
frequency=[a(1);diff(a)];
cumpr=cumsum(frequency)/n;

for i = 1:length(data_2)
    temp_2 = [];
    temp_2 = find(data_2(i) >= x);
    if ~isempty(temp_2)
        ind_2(i) = temp_2(end);
        cumpr_2(i) = cumpr(ind_2(i));
    else
        cumpr_2(i) = 0; 
    end
end

for j = 1:length(data)
    ind(j) = find(x == data(j));
    temp = cumpr(ind(j));
    cumpr_temp(j) = temp;
end
cumpr = cumpr_temp;