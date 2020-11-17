function b = make_coil_row(rk,X,Y,Z,i)
    nx = size(X,1);
    ny = size(X,2);
    nz = size(X,3);
    
    [xind, yind, zind] = ind2sub([nx-1,ny-1,nz], i);
    
    
    n = 1;
    rku = permute(rk(:,1:6),[3,2,1]);
    
    Xv = repmat(gpuArray(X(xind:xind+1,yind:yind+1,zind)),1,1,n);
    Yv = repmat(gpuArray(Y(xind:xind+1,yind:yind+1,zind)),1,1,n);
    Zv = repmat(gpuArray(Z(xind:xind+1,yind:yind+1,zind)),1,1,n);

    coords = [Xv(:) Yv(:) Zv(:)];

    mag = @(x,y,z) sqrt(x^2+y^2+z^2);

    r = arrayfun(@minus,coords,rku(:,1:3,:));
    rm = arrayfun(mag,r(:,1,:),r(:,2,:),r(:,3,:));
    A = arrayfun(@(x1,x2,x3,x4,rm) (x1*x2-x3*x4)/rm,...
            rku(:,[5,4,4],:),...
            r(:,[3,3,2],:),...
            rku(:,[6,6,5],:),...
            r(:,[2,1,1],:),rm);

    t1 = diff(reshape(A(:,2,:),2,2,1,size(rku,3)),1,1);
    t2 = diff(reshape(A(:,1,:),2,2,1,size(rku,3)),1,2);
    b = gather(reshape(t1(:,1:end-1,:,:)-t2(1:end-1,:,:,:),(2-1)*(2-1)*1,size(rku,3)));
    
end


