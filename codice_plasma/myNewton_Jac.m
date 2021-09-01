function [TT,err,residual,niter]=myNewton_Jac(funct,jfunct,TTguess,toll,nmax)

err=10*toll;
niter=0;
fun=funct(TTguess);
jac=jfunct(TTguess);
residual=norm(fun);
TT=TTguess;

while (err(end)>toll || residual(end)>toll) && niter<nmax   
    if condest(jac)>1e15            %se il numero di condizionamento è molto elevato (spannometricamente > 1e15) si innescano errori numerici molto rilevanti, dunque conviene procedere con il "preconditioning"
         [P,R,C]=equilibrate(jac);  %funzione di Matlab per il "preconditioning" di una matrice
         B=R*P*jac*C;
         d=R*P*fun;
         y=B\d;         
         deltaT=-C*y;
    else
         deltaT=-jac\fun;
    end
    
    TT=TT+deltaT; 
    niter=niter+1;
   
    %calcolo parametri di arresto
    
    err(niter)=norm(deltaT);
      
    %aggiornamento
    
    fun=funct(TT);
    jac=jfunct(TT);
    
    residual(niter)=norm(fun);
  
end

end
    