function output = k2h_h2k(input,interval,flag)
% interval=4:(h=1,k=1,2,3,4)(h=2,k=5,6,7,8)
if flag == 'k2h'
    k = input;
    h = ceil(k/interval);
    output = h;
elseif flag == 'h2k'
    h = input;
    k = (h-1)*4+3;% 取第三个
    output = k;   
end