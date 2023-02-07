function res = f_function(object,Q,V)


N = object.N;

tmp1 = 0;
for i=2:N
    tmp1 = tmp1+G(object,V,1,i);
end

tmp2 = zeros(N-1,1);

for j=1:N-1
    for i=1:N
        tmp2(j) = tmp2(j)+G(object,V,j+1,i);
    end
end
res = -Q(1)*tmp1+Q(2:end).*tmp2;

end


