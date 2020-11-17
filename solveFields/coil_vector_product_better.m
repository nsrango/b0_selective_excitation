function lhs = coil_vector_product_better(x,rk,X,Y,Z,mask,transpose)
nx = size(X,1);
ny = size(X,2);
nz = size(X,3);
n = 8;
rku = padarray(permute(rk(:,1:6),[3,2,1]),[0,0,n-mod(size(rk,1),n)],'post');
if ~strcmp(transpose,'transp')
    x = padarray(x,n-mod(size(rk,1),n),'post');
end

Xv = repmat(gpuArray(X(:)),1,1,n);
Yv = repmat(gpuArray(Y(:)),1,1,n);
Zv = repmat(gpuArray(Z(:)),1,1,n);

coords = [Xv Yv Zv];

if ~strcmp(transpose,'transp')
    mag = @(x,y,z) sqrt(x^2+y^2+z^2);
    lhs = gpuArray(zeros((nx-1)*(ny-1)*nz,1));
    for i = 1:size(rku,3)/n
        r = arrayfun(@minus,coords,rku(:,1:3,(1:n)+(i-1)*n));
        rm = arrayfun(mag,r(:,1,:),r(:,2,:),r(:,3,:));
        A = arrayfun(@(x1,x2,x3,x4,rm) (x1*x2-x3*x4)/rm,...
                rku(:,[5,4,4],(1:n)+(i-1)*n),...
                r(:,[3,3,2],:),...
                rku(:,[6,6,5],(1:n)+(i-1)*n),...
                r(:,[2,1,1],:),rm);

        t1 = diff(reshape(A(:,2,:),nx,ny,nz,n),1,1);
        t2 = diff(reshape(A(:,1,:),nx,ny,nz,n),1,2);
        b = reshape(t1(:,1:end-1,:,:)-t2(1:end-1,:,:,:),(nx-1)*(ny-1)*nz,n);
        lhs = lhs + b*x((1:n)+(i-1)*n);
    end
    lhs = gather(lhs(mask));
else
    mag = @(x,y,z) sqrt(x^2+y^2+z^2);
    lhs = gpuArray(zeros(1,size(rku,3)));
    for i = 1:size(rku,3)/n%             mcross = gpuArray(zeros(3,3,n));
        r = arrayfun(@minus,coords,rku(:,1:3,(1:n)+(i-1)*n));
        rm = arrayfun(mag,r(:,1,:),r(:,2,:),r(:,3,:));

        A = arrayfun(@(x1,x2,x3,x4,rm) (x1*x2-x3*x4)/rm,...
                rku(:,[5,4,4],(1:n)+(i-1)*n),...
                r(:,[3,3,2],:),...
                rku(:,[6,6,5],(1:n)+(i-1)*n),...
                r(:,[2,1,1],:),rm);
        t1 = diff(reshape(A(:,2,:),nx,ny,nz,n),1,1);
        t2 = diff(reshape(A(:,1,:),nx,ny,nz,n),1,2);
        b = reshape(t1(:,1:end-1,:,:)-t2(1:end-1,:,:,:),(nx-1)*(ny-1)*nz,n);
        lhs((1:n)+(i-1)*n) = x'*b(mask,:);
    end
    lhs = gather(lhs(1:size(rk,1))');
end

