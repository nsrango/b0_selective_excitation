function B = integrate_wirepath2(Gamma,X,Y,Z,I)

mpt = [];
l = [];
mu0 = 1.256e-6;
for i = 1:length(Gamma)
    pts = cat(3,Gamma{i}(:,1:end-1),Gamma{i}(:,2:end));
    l = cat(2,l,pts(:,:,2)-pts(:,:,1));
    mpt = cat(2,mpt,(pts(:,:,2)+pts(:,:,1))/2);
end


B = mu0/(4*pi)*evalFields(mpt, l*I, [X(:)'; Y(:)'; Z(:)']);

end

