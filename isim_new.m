function [isk,nc]=isim_new(x,y)
% ISIM_NEW Intersection similarity with ranking vectors


if ~isvector(x)||~isvector(y), error('isim:input','x and y must be vectors'); end
if numel(x)~=numel(y), error('isim:size','arguments must have equal length'); end

n = numel(x); X=x(:); Y=y(:); % if x,y are row vectors, this line transforms them in columns  
nc = zeros(n,1); nc(1) = NaN;

    % the case without ties is simple
    isk = zeros(n,1);
    
    for k=1:n
        sumsymdiff = 0; 
        
        for i = 1:k %summation index
            x = sort(X(1:i),'descend');
            y = sort(Y(1:i),'descend');
            for j = 1:i
                if isempty(find(y(1:i) == x(j),1)), 
                    sumsymdiff = sumsymdiff + 1/(i); 
                end
            end
        end
        isk(k) = sumsymdiff/k; 
        if k > 1
        nc(k) = k*isk(k) - (k-1)*isk(k-1);
        end
    end
    
end

