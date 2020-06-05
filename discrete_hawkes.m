function [Nt]=discrete_hawkes(t,T)

Nt=zeros(T,1);

for i=1:T
    Nt(i)=sum(t>=(i-1)&t<i);
end

end