function [c1,c2] = centroid(q1,q2)
    q1 = [q1;q1(1)];
    q2 = [q2;q2(1)];
    
    A = 1/2*sum(q1(1:end-1).*q2(2:end)-q1(2:end).*q2(1:end-1));
    
    c1 = 1/(6*A)*sum(...
        ((q1(1:end-1)+q1(2:end))...
            .*(q1(1:end-1).*q2(2:end)-q1(2:end).*q2(1:end-1))));
    c2 = 1/(6*A)*sum(...
        ((q2(1:end-1)+q2(2:end))...
            .*(q1(1:end-1).*q2(2:end)-q1(2:end).*q2(1:end-1))));
end