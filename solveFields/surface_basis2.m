function [bfull] = surface_basis2(rk,X,Y,Z)
    nx = size(X,1);
    ny = size(X,2);
    nz = size(X,3);
    n = 8;
    n = 8;
    rk_t = padarray(rk,[n-mod(size(rk,1),n),0],'post');
    bfull = zeros(nx*ny*nz,size(rk_t,1));
    for i = 1:size(rk_t,1)/n
        mc = rk_t((1:n)+(i-1)*n,1:3);
        mv = -rk_t((1:n)+(i-1)*n,4:6);
        evalPts = [X(:), Y(:), Z(:)];
        srcXYZ = gpuArray(ones(2,size(mc,1)));
        evalXYZ = gpuArray(-ones(2,numel(X)));
        Rsq = gpuArray(zeros(size(mc,1),numel(X)));
        r = gpuArray(zeros(size(mc,1),numel(X),3));
        for j = 1:3
            srcXYZ(1,:) = mc(:,j);
            evalXYZ(2,:) = evalPts(:,j);
            r(:,:,j) = (srcXYZ' * evalXYZ);
            Rsq = Rsq + r(:,:,j).^2;
        end
        Rcubed = Rsq.*sqrt(Rsq);
        Rfifth = Rcubed.*Rsq;

        evalXYZ(1,:) = 1;
        mdr = gpuArray(zeros(size(mc,1),numel(X)));
        for j = 1:3
            mdr = mdr + mv(:,j).*r(:,:,j);
        end
%         Bx = 3*rk_t(:,1).*mdr./Rfifth-X(:)'./Rcubed;
%         By = 3*rk_t(:,2).*mdr./Rfifth-Y(:)'./Rcubed;
        Bz = 3*r(:,:,3).*mdr./Rfifth-mv(:,3)./Rcubed;
        bfull(:,(1:n)+(i-1)*n) = gather(Bz');

    end
    bfull = bfull(:,1:size(rk,1));
end

