function c = myconstr2(x)
nC=length(x);
%c(1:nC-2)=0;
c(1) =  x(1)/(x(1)+x(2))-1;
c(2) =  x(2)/(x(1)+x(2))-1;
end
