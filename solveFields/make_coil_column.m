function b = make_coil_column(rk,X,Y,Z,i)
    nx = size(X,1);
    ny = size(X,2);
    nz = size(X,3);
    n = 1;
    rku = permute(rk(i,1:6),[3,2,1]);
    
    Xv = repmat(gpuArray(X(:)),1,1,n);
    Yv = repmat(gpuArray(Y(:)),1,1,n);
    Zv = repmat(gpuArray(Z(:)),1,1,n);

    coords = [Xv Yv Zv];

    mag = @(x,y,z) sqrt(x^2+y^2+z^2);

    r = arrayfun(@minus,coords,rku(:,1:3,(1:n)));
    rm = arrayfun(mag,r(:,1,:),r(:,2,:),r(:,3,:));
    A = arrayfun(@(x1,x2,x3,x4,rm) (x1*x2-x3*x4)/rm,...
            rku(:,[5,4,4],(1:n)),...
            r(:,[3,3,2],:),...
            rku(:,[6,6,5],(1:n)),...
            r(:,[2,1,1],:),rm);

    t1 = diff(reshape(A(:,2,:),nx,ny,nz,n),1,1);
    t2 = diff(reshape(A(:,1,:),nx,ny,nz,n),1,2);
    b = gather(reshape(t1(:,1:end-1,:,:)-t2(1:end-1,:,:,:),(nx-1)*(ny-1)*nz,n));
    
end


