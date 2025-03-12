function [k2,P2] = lens(k1,P1,R,O,fx,fy)

n  = R(:,3);

if transpose(n)*k1 < 0
    n = -n;
end

P2 = P1 + transpose(n)*(O - P1)./(transpose(n)*k1).*k1;
P2_ = transpose(R)*(P2 - O);
k1_ = transpose(R)*k1;

tx1 = asin(k1_(1,:)./sqrt(1 - k1_(2,:).^2));
ty1 = asin(k1_(2,:)./sqrt(1 - k1_(1,:).^2));

tx2 = -P2_(1,:)/fx + tx1;
ty2 = -P2_(2,:)/fy + ty1;

k2_ = zeros(3,size(k1,2));

k2_(1,:) = (sin(tx2).*cos(ty2))./sqrt(1 - (sin(tx2).*sin(ty2)).^2);
k2_(2,:) = (sin(ty2).*cos(tx2))./sqrt(1 - (sin(tx2).*sin(ty2)).^2);
k2_(3,:) = sqrt(1 - k2_(1,:).^2 - k2_(2,:).^2);

k2 = R*k2_;

end

