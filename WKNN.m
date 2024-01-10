function [DR_mat_new] = WKNN( DR_mat, D_mat, R_mat, K, r )

[rows,cols]=size(DR_mat);
y_d=zeros(rows,cols);  
y_r=zeros(rows,cols);  

knn_network_d = KNN( D_mat, K );  
for i = 1 : rows   
         w=zeros(1,K);
        [sort_d,idx_d]=sort(knn_network_d(i,:),2,'descend'); 
        sum_d=sum(sort_d(1,1:K));   
        for j = 1 : K
            w(1,j)=r^(j-1)*sort_d(1,j); 
            y_d(i,:) =  y_d(i,:)+ w(1,j)* DR_mat(idx_d(1,j),:); 
        end                      
            y_d(i,:)=y_d(i,:)/sum_d;              
end

knn_network_r = KNN( R_mat , K );  
for i = 1 : cols   
        w=zeros(1,K);
        [sort_r,idx_r]=sort(knn_network_r(i,:),2,'descend');
        sum_r=sum(sort_r(1,1:K));
        for j = 1 : K
            w(1,j)=r^(j-1)*sort_r(1,j);
            y_r(:,i) =  y_r(:,i)+ w(1,j)* DR_mat(:,idx_r(1,j)); 
        end                      
            y_r(:,i)=y_r(:,i)/sum_r;               
end

y_dr=(y_d+y_r)/2;  

 for i = 1 : rows
     for j = 1 : cols
         DR_mat_new(i,j)=max(DR_mat(i,j),y_dr(i,j));
     end    
 end

end

function [ knn_network ] = KNN( network , k )
    [rows, cols] = size( network );
    network= network-diag(diag(network)); 
    knn_network = zeros(rows, cols);
    [sort_network,idx]=sort(network,2,'descend');
    for i = 1 : rows
        knn_network(i,idx(i,1:k))=sort_network(i,1:k);
    end
end


