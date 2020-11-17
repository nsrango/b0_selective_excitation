function B = integrate_wirepath3(LS,LD,X,Y,Z,I)

mpt = [];
l = [];
mu0 = 1.256e-6;
    l = cat(2,l,LS(:,:)-LD(:,:));
    mpt = cat(2,mpt,(LS(:,:)+LD(:,:,1))/2);


B = mu0/(4*pi)*evalFields(mpt, l*I, [X(:)'; Y(:)'; Z(:)']);

end

