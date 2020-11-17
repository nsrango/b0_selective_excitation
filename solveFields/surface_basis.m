function [bfull_out] = surface_basis(rk,X,Y,Z)
nx = size(X,1);
ny = size(X,2);
nz = size(X,3);
n = 8;
rku = padarray(permute(rk(:,:),[3,2,1]),[0,0,n-mod(size(rk,1),n)],'post');

Xv = repmat(gpuArray(X(:)),1,1,n);
Yv = repmat(gpuArray(Y(:)),1,1,n);
Zv = repmat(gpuArray(Z(:)),1,1,n);

bfull = zeros((nx-1)*(ny-1)*nz,size(rku,3));
for i = 1:size(rku,3)/n
    mc = rku(:,1:3,(1:n)+(i-1)*n);
    mv = -rku(:,4:6,(1:n)+(i-1)*n);

    r = [Xv-(mc(:,1,:)) Yv-mc(:,2,:) Zv-mc(:,3,:)];
    rm = repmat(sqrt(r(:,1,:).^2+r(:,2,:).^2+r(:,3,:).^2),1,3,1);
    mcross = gpuArray(zeros(3,3,n));
    mcross(1,2,:) = -mv(1,3,:);
    mcross(1,3,:) = mv(1,2,:);
    mcross(2,1,:) = mv(1,3,:);
    mcross(2,3,:) = -mv(1,1,:);
    mcross(3,1,:) = -mv(1,2,:);
    mcross(3,2,:) = mv(1,1,:);
    mtx = blkdiag(mcross(:,:,1),mcross(:,:,2),mcross(:,:,3),mcross(:,:,4),...
        mcross(:,:,5),mcross(:,:,6),mcross(:,:,7),mcross(:,:,8));
    A = reshape((mtx*reshape(r,[],n*3).').'./reshape(rm,[],n*3).^3,[],3,n);
%     A = ((mtx*reshape(r,[],n*3).').'./reshape(rm,[],n*3).^3);
    t1 = diff(reshape(A(:,2,:),nx,ny,nz,n),1,1);
    t2 = diff(reshape(A(:,1,:),nx,ny,nz,n),1,2);
    b = reshape(t1(:,1:end-1,:,:)-t2(1:end-1,:,:,:),(nx-1)*(ny-1)*nz,n);
    bfull(:,(1:n)+(i-1)*n) = gather(b);
end

bfull_out = bfull(:,1:size(rk,1));
end

