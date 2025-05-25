function dydt = wkadydt(t,y,sys)
% implementation of ode
% sys holds system parameters
% get coordinates

% equations of motion */
%
%
%       (M[O] + cos(Ω t)*M[C] + sin(Ω t)*M[S])*Q'' + 
%  Ω*   (G[O] + cos(Ω t)*G[C] + sin(Ω t)*G[S])*Q'  + 
%        K[K]                                 *Q   + 
%  Ω^2* (K[O] + cos(Ω t)*K[C] + sin(Ω t)*K[S])*Q    = 0 

P = y( 1:14,1);
Q = y(15:28,1);

Oga = sys.Oga(sys.i);
M =                 sys.MO + cos(Oga*t)*sys.MC + sin(Oga*t)*sys.MS;
G =          Oga  *(sys.GO + cos(Oga*t)*sys.GC + sin(Oga*t)*sys.GS);
K = sys.KK + Oga^2*(sys.KO + cos(Oga*t)*sys.KC + sin(Oga*t)*sys.KS);

%disp(size(M))

dydt(15:28,1) = -inv(M)*(G*P+K*Q);

%% velocities
dydt(1:14,1) = y(15:28,1);

end