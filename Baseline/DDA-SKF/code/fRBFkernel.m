function [B]=fRBFkernel(A,sigma)

[m,n]=size(A);

for i=1:m
    for j=1:i
        B(i,j)= exp(-norm(A(i,:)-A(j,:))^2/(2*sigma^2));
        B(j,i)=B(i,j);
    end
end

end
