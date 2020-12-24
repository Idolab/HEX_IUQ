function y_interp = func_interp(u,x_samp,y_samp)

for i=1:size(u,1)
if u(i)>=x_samp(1) & u(i)<x_samp(end)
ind1=max(find(x_samp<=u(i)));
ind2=min(find(x_samp>u(i)));
elseif u(i)<x_samp(1)
    ind1=1;ind2=2;
elseif u(i)>=x_samp(end)
    ind1=size(x_samp,1)-1;ind2=size(x_samp,1);
end

y1=y_samp(ind1);
y2=y_samp(ind2);
x1=x_samp(ind1);
x2=x_samp(ind2);

y_interp(i,:)=(y2-y1)/(x2-x1)*(u(i)-x1)+y1;
end