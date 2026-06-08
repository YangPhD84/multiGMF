function [ result ] =myFunNorm( M )
    [num1,num2] = size(M);
    nM = zeros(num1,num2);
    result = zeros(num1,num2);
    
    for ii = 1:num1
        rnM(ii) = sum(M(ii,:));
    end
    
    for ij = 1:num2
        cnM(ij) = sum(M(:,ij));
    end

    for i = 1:num1
        rsum = rnM(i);
        for j = 1:num2
            csum = cnM(j);
            if((rsum==0)||(csum==0))
                result(i,j) = 0;
            else
                result(i,j) = M(i,j)/sqrt(rsum*csum);
            end
        end
    end
end