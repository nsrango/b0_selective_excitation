function A = construct_unitTransforms(coil)
A = zeros(4,4,size(coil.ConnectivityList,1));
    for i = 1:size(coil.ConnectivityList,1)
        p = coil.ConnectivityList(i,:);
        A(:,4,i) = [coil.Points(p(3),:)';1];
        A(1:3,1,i) = coil.Points(p(1),:)'-A(1:3,4,i);
        A(1:3,2,i) = coil.Points(p(2),:)'-A(1:3,4,i);
        A(1:3,3,i) = (coil.Points(p(2),:)+coil.Points(p(3),:))'/2 - A(1:3,4,i)- sqrt(2)*A(1:3,2,i);
    end
end