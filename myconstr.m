function c = myconstr(x)
nC=length(x);
c(1:nC-2)=0;
c(nC-1) =  x(nC-1)/(x(nC-1)+x(nC))-1;
c(4) =  x(nC)/(x(nC-1)+x(nC))-1;
end
