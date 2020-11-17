function lhs = coil_vector_product(x,rk,X,Y,Z,mask,transpose)
nx = size(X,1);
ny = size(X,2);
nz = size(X,3);
n = 14;
rku = padarray(permute(rk(:,1:6),[3,2,1]),[0,0,n-mod(size(rk,1),n)],'post');
if ~strcmp(transpose,'transp')
    x = padarray(x,n-mod(size(rk,1),n),'post');
end

Xv = repmat(gpuArray(X(:)),1,1,n);
Yv = repmat(gpuArray(Y(:)),1,1,n);
Zv = repmat(gpuArray(Z(:)),1,1,n);
if ~strcmp(transpose,'transp')
    lhs = gpuArray(zeros((nx-1)*(ny-1)*nz,1));
    for i = 1:size(rku,3)/n
        mc = rku(:,1:3,(1:n)+(i-1)*n);
        mv = rku(:,4:6,(1:n)+(i-1)*n);

        r = [Xv-(mc(:,1,:)) Yv-mc(:,2,:) Zv-mc(:,3,:)];
        rm = repmat(sqrt(r(:,1,:).^2+r(:,2,:).^2+r(:,3,:).^2),1,3,1);
%             mcross = gpuArray(zeros(3,3,n));
%             mcross(1,2,:) = -mv(1,3,:);
%             mcross(1,3,:) = mv(1,2,:);
%             mcross(2,1,:) = mv(1,3,:);
%             mcross(2,3,:) = -mv(1,1,:);
%             mcross(3,1,:) = -mv(1,2,:);
%             mcross(3,2,:) = mv(1,1,:);
%                 mtx = blkdiag(mcross(:,:,1),mcross(:,:,2),mcross(:,:,3),mcross(:,:,4),...
%                     mcross(:,:,5),mcross(:,:,6),mcross(:,:,7),mcross(:,:,8));
%                 
%         rtmp = reshape(r,[],n*3).';
%         rmtmp = reshape(rm,[],n*3).^3;
        
        A = cat(2,mv(:,2,:).*r(:,3,:) - mv(:,3,:).*r(:,2,:),...
                  mv(:,1,:).*r(:,3,:) - mv(:,3,:).*r(:,1,:),...
                  mv(:,1,:).*r(:,2,:) - mv(:,2,:).*r(:,1,:))./rm;
        t1 = diff(reshape(A(:,2,:),nx,ny,nz,n),1,1);
        t2 = diff(reshape(A(:,1,:),nx,ny,nz,n),1,2);
        b = reshape(t1(:,1:end-1,:,:)-t2(1:end-1,:,:,:),(nx-1)*(ny-1)*nz,n);
        lhs = lhs + b*x((1:n)+(i-1)*n);
    end
    lhs = gather(lhs(mask));
else
    lhs = gpuArray(zeros(1,size(rku,3)));
    for i = 1:size(rku,3)/n%             mcross = gpuArray(zeros(3,3,n));
        mc = rku(:,1:3,(1:n)+(i-1)*n);
        mv = rku(:,4:6,(1:n)+(i-1)*n);
        r = [Xv-(mc(:,1,:)) Yv-mc(:,2,:) Zv-mc(:,3,:)];
        rm = (sqrt(r(:,1,:).^2+r(:,2,:).^2+r(:,3,:).^2));
%         mcross = gpuArray(zeros(3,3,n));
%         mcross(1,2,:) = -mv(1,3,:);
%         mcross(1,3,:) = mv(1,2,:);
%         mcross(2,1,:) = mv(1,3,:);
%         mcross(2,3,:) = -mv(1,1,:);
%         mcross(3,1,:) = -mv(1,2,:);
%         mcross(3,2,:) = mv(1,1,:);
%         mtx = blkdiag(mcross(:,:,1),mcross(:,:,2),mcross(:,:,3),mcross(:,:,4),...
%             mcross(:,:,5),mcross(:,:,6),mcross(:,:,7),mcross(:,:,8));
%         rtmp = reshape(r,[],n*3).';
%         rmtmp = reshape(rm,[],n*3).^3;
%         A = (mtx*rtmp).'./rmtmp;

        
        a1 = (mv(:,2,:).*r(:,3,:) - mv(:,3,:).*r(:,2,:))./rm;
        a2 = (mv(:,1,:).*r(:,3,:) - mv(:,3,:).*r(:,1,:))./rm;
        a3 = (mv(:,1,:).*r(:,2,:) - mv(:,2,:).*r(:,1,:))./rm;

        
        
        
        A = cat(2,a1,...
                  a2,...
                  a3);
              
        A = reshape(A,[],3,n);
        t1 = diff(reshape(A(:,2,:),nx,ny,nz,n),1,1);
        t2 = diff(reshape(A(:,1,:),nx,ny,nz,n),1,2);
        b = reshape(t1(:,1:end-1,:,:)-t2(1:end-1,:,:,:),(nx-1)*(ny-1)*nz,n);
        lhs((1:n)+(i-1)*n) = x'*b(mask,:);
    end
    lhs = gather(lhs(1:size(rk,1))');
end

