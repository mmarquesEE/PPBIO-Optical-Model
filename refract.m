function [k2,P2] = refract(k1,P1,R,O,n1,n2)

n  = R(:,3);

if transpose(n)*k1 < 0
    n = -n;
end

P2 = P1 + transpose(n)*(O - P1)./(transpose(n)*k1).*k1;

k2p = (n1/n2)*(k1 - (transpose(n)*k1).*n); 
k2s = sqrt(1 - sum(k2p.^2,1)).*n;

k2 = k2p + k2s;

end

