function[xx]=myNewton(funct,dfunct,x0,toll)
    err=10*toll;
    residual=10*toll;
    fun=funct(x0);
    dfun=dfunct(x0);
    xx=x0;
    n_iter=0;

    while err>toll | residual>toll
        deltax=-fun/dfun;
        xx=xx+deltax;
        n_iter=n_iter+1;
        fun=funct(xx);
        dfun=dfunct(xx);
        err(n_iter)=norm(deltax);
        residual(n_iter)=norm(fun);
    end
end