%% 2. Fibonacci number generator.
function f = fibnum(n)
    if n < 1
        error('n must be at least 1');
    elseif n == 1
        f = 1; % f = [1];
        return;
    elseif n == 2
        f = [1; 2]; % ?? f = [1, 2];
        return;
    end
    
    f = zeros(n,1);
    f(1) = 1; f(2) = 2;
    for k= 3:n
        f(k) = f(k-1) + f(k-2);
    end
end