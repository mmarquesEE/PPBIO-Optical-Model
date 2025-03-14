function [k2,P2] = reflect(k1,P1,R,O)

n  = R(:,3);

if transpose(n)*k1 > 0
    n = -n;
end

P2 = P1 + transpose(n)*(O - P1)./(transpose(n)*k1).*k1;
k2 = k1 - 2*(transpose(n)*k1).*n;

end

