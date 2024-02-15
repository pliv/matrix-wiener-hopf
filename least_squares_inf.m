function c = least_squares_inf(f, t, s)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here



    
if s==1
    h_1 = 1./t;
    h_3 = 1./t.^3;
    hh = [h_1; h_3];
    n = 2;
elseif s==0
    h_2 = 1./t.^2;
    h_4 = 1./ t.^4;
    hh = [h_2; h_4];
    n = 2;
elseif s==2
    h_2 = 1./t.^2;
    hh = h_2;
    n=1;
elseif s==3
    h_3 = 1./t.^3;
    hh = h_3;
    n=1;
elseif s==4
    h_3 = 1./t.^4;
    hh = h_3;
    n=1;    
elseif s==5
    h_1 = 1./t;
    hh = h_1;
    n=1;
elseif s==100 %testing
    h_3 = 1./t;
    hh = h_3;
    n=1;
end
    H = zeros(n,n);
    B = zeros(n,1);
    for j=1:n
        for q=1:n
            H(j,q) = sum(hh(j,:).*hh(q,:));
        end
        B(j) = sum(hh(j,:).*f);
    end
    c = H\B;
    
c = c.';
c = double(c);
end


