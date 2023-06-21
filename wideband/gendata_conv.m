function X = gendata_conv(s,P,N,sigma)

% Construct h(t-k) and multiply sk to get x(t)
t_new = linspace(0,N-1/P,N*P);
k=0:N-1;
for n=1:N
    for m=1:length(t_new)
        t_k=t_new(m)-k(n);
        if t_k<0.25 && t_k>=0
            h=1;
            x(m)=h*s(n);
        end
        if t_k>=0.25 && t_k<0.5
            h=-1;
            x(m)=h*s(n);
        end
        if t_k>=0.5 && t_k<0.75
            h=1;
            x(m)=h*s(n);
        end
        if t_k>=0.75 && t_k<1
            h=-1;
            x(m)=h*s(n);
        end
    end
end
%Create complex noise
real_part = sqrt(2)/2*sigma * randn(1, N*P);
imaginary_part = sqrt(2)/2*sigma * randn(1, N*P);
n = complex(real_part, imaginary_part);

% 
% Construct the sequence X(0:N-1/P)
% for i=1:N
%     X(i)=x((i-1)*200/P+1);
% end
X=x+n;

