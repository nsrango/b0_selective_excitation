function fields = evalFields(srcPts, srcW, evalPts)
% srcPts{evalPts} = 3 x n_srcs{n_evals} matrix whose columns are [x,y,z]
% srcW = n-length vector of source charges, or 3xn matrix of src currents
% evalPts = 3 x n_evals matrix of field evaluation points

nSrc = size(srcPts,2);
nEval = size(evalPts,2);

srcXYZ = ones(2,nSrc);
evalXYZ = -ones(2,nEval);
mats2e = cell(3);
Rsq = zeros(nSrc,nEval);
% Form x-x', y-y' and z-z' annd (dx^2+dy^2+dz^2)^1.5
for i = 1:3
    srcXYZ(1,:) = srcPts(i,:);
    evalXYZ(2,:) = evalPts(i,:);
    mats2e{i} = srcXYZ' * evalXYZ;
    Rsq = Rsq + mats2e{i}.^2;
end
Rcubed = Rsq.*sqrt(Rsq);
fields = zeros(3,nEval);
if size(srcW,1) == 1  % Efields from src charges
    for i = 1:3
       fields{i} = srcW'*(mats2e{i}./Rcubed);
    end
else % B fields from src currents
    for i = 1:3
        ip1 = 1+mod(i,3);
        ip2 = 1+mod(i+1,3);
        fields(i,:) = srcW(ip1,:)*(mats2e{ip2}./Rcubed) - srcW(ip2,:)*(mats2e{ip1}./Rcubed);
    end
end
end