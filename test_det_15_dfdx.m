syms t h V;
rank([1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0;
    0, 0, 0, t, t, t, 2*t, 2*t, 2*t, -1, -1, -1, -1, -1, -1;
    0, 0, 0, t^2, t^2, t^2, 4*t^2, 4*t^2, 4*t^2, 0, 0, 0, -2*t, -2*t, -2*t;
    0, 0, 0, t, t, t, 8*t, 8*t, 8*t, 0, 0, 0, -3, -3, -3;
    0, 0, 0, t, t, t, 16*t, 16*t, 16*t, 0, 0, 0, -4 ,-4, -4;
    -h, 0, h, -h, 0, h, -h, 0, h,-V, -V, -V, -V, -V, -V;
    h, 0, h, h, 0, h, h,0,  h, 2*V, 0, -2*V, 2*V, 0, -2*V;
    -h, 0, h, -h, 0, h, -h, 0, h,-3*V, 0, -3*V, -3*V, 0, -3*V;
    h, 0, h, h, 0, h, h,0,  h, 4*V, 0, -4*V, 4*V, 0, -4*V; 
    0, 0, 0, -t*h, 0, t*h, -2*t*h, 0, 2*t*h,h, 0, -h, h-V*t, -V*t, -h-V*t;
    0, 0, 0, t*h^2, 0, t*h^2, 2*t*h^2, 0, 2*t*h^2, -h^2, 0, -h^2, 2*V*t*h-h^2, 0, -h^2-2*V*t*h;
    0, 0, 0, -t*h^3, 0, t*h^3, -2*t*h^3, 0, 2*t*h^3, h^3, 0, -h^3, -3*V*t*h^2+h^3, 0, -3*V*t*h^2-h^3;
    0, 0, 0 ,-t^2*h, 0, t^2*h, -4*t^2*h, 0, 4*t^2*h, 0, 0, 0, -V*t^2+2*t*h, -V*t^2, -V*t^2-2*t*h;
    0, 0, 0, t^2*h^2, 0, t^2*h^2, 4*t^2*h^2, 0, 4*t^2*h^2, 0, 0, 0, -2*t*h^2+2*V*t^2*h, 0, -2*t*h^2-2*V*t^2*h;
    0, 0, 0, -t^3*h, 0, t^3*h, -8*t^3*h, 0, 8*t^3*h, 0, 0, 0, -V*t^3+3*t^2*h, -V*t^3, -V*t^3-3*t^2*h])