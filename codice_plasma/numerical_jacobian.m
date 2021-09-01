function numJac=numerical_jacobian(fun,TT,dT)
numJac=zeros(length(TT),length(TT));
for jj=1:length(TT)
    TT_up=TT;
    TT_up(jj)=TT_up(jj)+dT/2;
    TT_dw=TT;
    TT_dw(jj)=TT_dw(jj)-dT/2;
    f_up=fun(TT_up);
    f_dw=fun(TT_dw);
    numJac(:,jj)=(f_up-f_dw)/dT;   
end
end