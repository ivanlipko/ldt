function f = rkStepPeterbFullSystem(t, h, X, right, A,B,G, cntrl, ptb)

h2 = 0.5*h;
h6 = 0.166666666*h;
Fs = right(t,X, A,B,G, cntrl,ptb);
t = t + h2;

Xr = X + h2 * Fs;

F = right(t,Xr, A,B,G, cntrl,ptb);

s=F;
Fs=Fs+s+s;
Xr=X+h2*s;

F = right(t,Xr, A,B,G, cntrl,ptb);
t=t+h2;

s=F;
Fs=Fs+s+s;
Xr=X+h*s;

F = right(t,Xr, A,B,G, cntrl,ptb);

X=X+h6*(Fs+F);

f=X;
