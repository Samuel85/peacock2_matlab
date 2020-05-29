function [pValue, KS] = peackock2(x, y)
    [m_x,n_x] = size(x);
    [m_y,n_y] = size(y);
    xx = x;
    yy = y;
    %%
    n1 = m_x;
    n2 = m_y;
    n = n1 + n2;
    %%
    dd = gcd(n1, n2);
    L = n1/dd*n2;
    d1 = L/n1;
    d2 = L/n2;    
    %%
    dim1 = n_x;
    dim2 = n_y;
    dmin = min([dim1 dim2]);
    
    if dmin < 2
        disp('The dimensions of both samples should be at least 2'); 
        return;
    end
    
    xy1 = [xx(:,1); yy(:,1)];
    xy2 = [xx(:,2); yy(:,2)];
    
    [B1,I1] = sort( xy1 );
    [B2,I2] = sort( xy2 );      
    
    max_hnn = 0;
    max_hpn = 0;
    max_hnp = 0;
    max_hpp = 0;
    
    for zu=xy1(I1)
        
        hnn = 0;
        hpn = 0;
        for t=1:n
            v = I2(t);
            if xy1(v) <= zu
                if v <= n1
                    hnn = hnn + d1;
                else
                    hnn = hnn - d2;
                end
                max_hnn = max([max_hnn abs(hnn)]);
            else
                if v <= n1
                    hpn = hpn + d1;
                else
                    hpn = hpn - d2;
                end
                max_hpn = max([max_hpn abs(hpn)]);
            end
        end
        
        hnp = 0;
        hpp = 0;        
        for t=n:-1:1
            v = I2(t);
            if xy1(v) <= zu
                if v <= n1                    
                    hnp = hnp + d1;
                else
                    hnp = hnp - d2;
                end
                max_hnp = max([max_hnp, abs(hnp)]);
            else
                if v <= n1
                    hpp = hpp + d1;
                else
                    hpp = hpp - d2;
                end
                max_hpp = max([max_hpp, abs(hpp)]);
            end
        end                                        
    end
    KS = max([max_hnn, max_hpn, max_hnp, max_hpp])/L;      
    
    % Average correlation coefficients
    r1 = corrcoef(x); r1 = r1(1,2);
    r2 = corrcoef(y); r2 = r2(1,2);
    rr = 0.5*(r1*r1 + r2*r2);

    pValue = probks(n1,n2,KS,rr);   
end

function p = probks(n1,n2,D,rr)
%Numerical Recipes in C, section 14.7
N = (n1*n2)/(n1+n2);
lambda = (sqrt(N)*D) / (1 + sqrt(1 - rr)*(.25 - .75/sqrt(N)));

j = (1:101)';
p = 2 * sum((-1).^(j-1).*exp(-2*lambda*lambda*j.^2));
p = min(max(p,0),1);

end

