function B=sharematrix(A)
[l1,l2]=size(A);
% B=zeros(l1,l1);
for i=1:l1
    for j=i:l1
        
  ai= find(A(i,:)==1); 
  aj=find(A(j,:)==1);     
  B(i,j)=length(intersect(ai,aj));
  B(j,i)=B(i,j);
  
    end
end

