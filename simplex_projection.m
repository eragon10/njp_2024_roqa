function [x, y] = simplex_projection(xin, varargin)

    niter = 100;
    rho   = 0.5;
    scale = 1;
    
    while ~isempty(varargin)
        switch lower(varargin{1})
            case 'iterations'
                niter = varargin{2};
            case 'rho'
                rho = varargin{2};
            case 'scale'
                scale = varargin{2};
            otherwise
                error(['Unexpected option: ' varargin{1}])
         end
         varargin(1:2) = [];
    end

    
    u = xin / scale;
    z = u;
    y = ones(size(xin,1),1);
    
    M = eye(size(xin,1)) - 1 / size(xin,1);
    
    for i = 1:niter
        x = 1 / (1+rho) * M * (u + rho * z - y) + 1 / size(xin,1);
        z = max(0.0, x + y / rho);
        y = y + rho * (x - z);
    end
    
    x = scale * z;
end