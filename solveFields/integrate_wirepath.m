function [Bx,By,Bz] = integrate_wirepath(Gamma,X,Y,Z,I)
mu0 = 4*pi*10^-7;
pts = cat(3,Gamma(1:end-1,:),Gamma(2:end,:));
npts = size(pts,1);
ntgts = numel(X);
l = pts(:,:,2)-pts(:,:,1);
normL = vecnorm(l,2,2);

v1x = [pts(:,1,1) -ones( npts,1)]*[ones(1,ntgts); X(:)'];
v1y = [pts(:,2,1) -ones( npts,1)]*[ones(1,ntgts); Y(:)'];
v1z = [pts(:,3,1) -ones( npts,1)]*[ones(1,ntgts); Z(:)'];

v1norm2 = v1x.^2+v1y.^2+v1z.^2;

v2x = [pts(:,1,2) -ones( npts,1)]*[ones(1,ntgts); X(:)'];
v2y = [pts(:,2,2) -ones( npts,1)]*[ones(1,ntgts); Y(:)'];
v2z = [pts(:,3,2) -ones( npts,1)]*[ones(1,ntgts); Z(:)'];

v2norm2 = (v2x.^2+v2y.^2+v2z.^2);


spv1 = bsxfun(@times,v1x,l(:,1)./normL)+...
       bsxfun(@times,v1y,l(:,2)./normL)+...
       bsxfun(@times,v1z,l(:,3)./normL);

ax = v1x-bsxfun(@times,spv1,l(:,1)./normL);
ay = v1y-bsxfun(@times,spv1,l(:,2)./normL);
az = v1z-bsxfun(@times,spv1,l(:,3)./normL);

spv2 = bsxfun(@times,v2x,l(:,1)./normL)+...
       bsxfun(@times,v2y,l(:,2)./normL)+...
       bsxfun(@times,v2z,l(:,3)./normL);

ax2 = v2x-bsxfun(@times,spv2,-l(:,1)./normL);
ay2 = v2y-bsxfun(@times,spv2,-l(:,2)./normL);
az2 = v2z-bsxfun(@times,spv2,-l(:,3)./normL);


normA2 = ax.^2+ay.^2+az.^2;

cos_theta1 = sqrt(v1norm2 -normA2)./sqrt(v1norm2);
cos_theta2 = sqrt(v2norm2 -normA2)./sqrt(v2norm2);

Bm =  mu0*I./(4*pi*sqrt(normA2)).*(cos_theta1+cos_theta2);

ctermx = bsxfun(@times,v1z,l(:,2))-bsxfun(@times,v1y,l(:,3));
ctermy = bsxfun(@times,v1x,l(:,3))-bsxfun(@times,v1z,l(:,1));
ctermz = bsxfun(@times,v1y,l(:,1))-bsxfun(@times,v1x,l(:,2));
cmag = sqrt(ctermx.^2+ctermy.^2+ctermz.^2);

Bx = sum(Bm.*(ctermx./cmag),1);
By = sum(Bm.*(ctermy./cmag),1);
Bz = sum(Bm.*(ctermz./cmag),1);

% t1 = acos(  sum   )

end

