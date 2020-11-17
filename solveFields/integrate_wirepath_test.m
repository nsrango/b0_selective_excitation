L = 1;
shift = 0;
Gamma1 = [0 0 -L;0 0 L];
Gamma2 = [0 0 -L; 0 0 0; 0 0 L];
[X, Y, Z] = meshgrid(linspace(-20,20,100),linspace(-20,20,100),[-L,0,L]);

[Bx1,By1,Bz1] = integrate_wirepath(Gamma1,X,Y,Z,1);
B1 = reshape(sqrt(Bx1.^2+By1.^2+Bz1.^2),size(X));
Bx1i = reshape(Bx1,size(X)); By1i = reshape(By1,size(X)); Bz1i = reshape(Bz1,size(X));
[Bx2,By2,Bz2] = integrate_wirepath(Gamma2,X,Y,Z,1);
B2 = reshape(sqrt(Bx2.^2+By2.^2+Bz2.^2),size(X));
Bx2i = reshape(Bx2,size(X)); By2i = reshape(By2,size(X)); Bz2i = reshape(Bz2,size(X));



a = sqrt(X.^2+Y.^2);
mu0 = 4*pi*10^-7;
B_gt = mu0*1./(2*pi*a)*L./sqrt(L^2+a.^2);

s =  50;
imagesc([B2(:,:,2);B2(:,:,3)],1e-8*[-1,1]); set(gca,'colormap',colorcet('D1')); axis image;
%%
[Xt, Yt, Zt] = meshgrid(linspace(-20,20,100),20,0);

[Bxt,Byt,Bzt] = integrate_wirepath(Gamma1,Xt,Yt,Zt,1);

%%

fields = evalFields([0;0;0],[0;0;1],[X(:)'; Y(:)'; Z(:)']);
bx = reshape(fields(1,:),size(X));

imagesc(bx(:,:,3)-bx(:,:,2),1e-2*[-1,1]);

%%


L = 2;
Gammaz = linspace(-L,L,20);
Gamma1 = [0*Gammaz; 0*Gammaz; Gammaz];
[X, Y, Z] = meshgrid(linspace(-20,20,100),linspace(-20,20,100),20);
fields = integrate_wirepath2(Gamma1,X,Y,Z,1);
bx = reshape(fields(1,:),size(X));
imagesc(bx(:,:,1),1e-9*[-1,1]);





