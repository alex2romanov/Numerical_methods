syms v
syms w
syms z

z = (v*exp(-i*w)-v*exp(i*w)+4)/(v*exp(i*w)-v*exp(-i*w)+4);

S30 = simplify(abs(z),'Steps',30)
%S50 = simplify(im,'Steps',50)
%S100 = simplify(im,'Steps',100)
%S500 = simplify(im,'Steps', 500)
%S1000 = simplify(im,'Steps', 1000)
