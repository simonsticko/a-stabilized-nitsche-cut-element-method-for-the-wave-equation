%Returns the constant of the hat function for the element with nodes at x
%and y. phi=a+b*x+c*y
function[a,b,c,area]=getPhi(x,y)
area=polyarea(x,y);
a=.5*(x([2 3 1]).*y([3 1 2])-x([3 1 2]).*y([2 3 1]))/area;
b=.5*(y([2 3 1])-y([3 1 2]))/area;
c=.5*(x([3 1 2])-x([2 3 1]))/area;
end