% method of Casumano et al (1987)(rooting interpol. conditions)
clear
pom=[4, 13, -35, 22]; 

G=tf([1, -2, 2],[1, 1])*tf([1, -4, 5],pom);
n=4; % order of the plant
m=1; % number of inputs/outputs

% coprime factorization of G
eps=1;  % try also 0.1;
bQ=tf(pom,[1, 2*eps, eps^2])*tf([1,1], [1, 2*eps, eps^2]);
bP=tf([1, -2, 2],[1, 2*eps, eps^2])*tf([1, -4, 5],[1, 2*eps, eps^2]);

[A,B,C,D]=ssdata(G);
[n,n1]=size(A);
[m1,m]=size(B);
[p,p1]=size(C);

 s=tzero(bP)
 
 % interpolation points
 s1=s(1);
 s2=s(2);
 s3=s(3);
 s4=s(4);

 ell=4; % number of interpol. points

 [AY,BY,CY,DY]=ssdata(bQ);
 [nY,nY1]=size(AY);
 % interpolation data/conditions
 y1=(DY+CY*(s1*eye(nY)-AY)^(-1)*BY); 
 y2=(DY+CY*(s2*eye(nY)-AY)^(-1)*BY); 
 y3=(DY+CY*(s3*eye(nY)-AY)^(-1)*BY); 
 y4=(DY+CY*(s4*eye(nY)-AY)^(-1)*BY); 
 
Api=diag([s1,s2,s3,s4]);
Cmin=[1 1 1 1];
Cpl=[y1,y2,y3,y4];

mm=4; % minimal exponent

CCpl=[y1^(1/mm), y2^(1/mm), y3^(1/mm), y4^(1/mm) ], % rooting the interpol. cond.

% Checking the compex Pick matrix X
X=lyap(Api',-2*(Cmin'*CCpl+CCpl'*Cmin))
eig(X)
% another check
for i=1:4
    for k=1:4
Lam(i,k)=(CCpl(i)'+CCpl(k))/(s(i)'+s(k));
end
end
2*eig(Lam)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% transformation to a problem with real coefficients
T0=[1 -complex(0,1); 1 complex(0,1)]/sqrt(2);
T=[T0 zeros(2,2); zeros(2,2), T0];
hApi=real(T'*Api*T);
hCmin=real(Cmin*T);
hCCpl=real(CCpl*T);

% finding SPR RF Hm that satisfies the rooted interpol. conditions in CCpl
tCmin=hCmin+hCCpl;
tCpl=hCmin-hCCpl;

XX=lyap(hApi',-(tCmin'*tCmin-tCpl'*tCpl)) % real Pick matrix
% checking the eig. of the real Pick matrix
e=eig(XX)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Hm=dss(XX*hApi-tCmin'*(tCmin+tCpl), tCmin', -2*tCpl, eye(m), XX); % SPR interpolant
[AHm,BHm,CHm,DHm]=ssdata(Hm);
[nHm,hHm1]=size(AHm);

% checking Hm
pole(Hm)
tzero(Hm)
% checking the int. cond. on Hm
(DHm+CHm*(s1*eye(nHm)-AHm)^(-1)*BHm)-CCpl(1)
(DHm+CHm*(s2*eye(nHm)-AHm)^(-1)*BHm)-CCpl(2)
(DHm+CHm*(s3*eye(nHm)-AHm)^(-1)*BHm)-CCpl(3)
(DHm+CHm*(s4*eye(nHm)-AHm)^(-1)*BHm)-CCpl(4)
% checking the SPR property of Hm
P=ispr(AHm,BHm,CHm,DHm);
eig(P)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

H=Hm^mm; % interpolant
[AH,BH,CH,DH,EH]=dssdata(H);
[nH,nH1]=size(AH);

% checking of the interpolant H
p=pole(H)
z=tzero(H) 
% checking the int. cond. on H
(DH+CH*(s1*EH-AH)^(-1)*BH)-y1
(DH+CH*(s2*EH-AH)^(-1)*BH)-y2
(DH+CH*(s3*EH-AH)^(-1)*BH)-y3
(DH+CH*(s4*EH-AH)^(-1)*BH)-y4
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Cons=(bQ-H)*bP^(-1); % controller
[Cons1,unst]=stabsep(Cons); % canceling its unstable modes
Conr=reduce(Cons1); % reduced controller

% Checking
pp=pole(Cons1)  %  stability of the cont.
zz=tzero([ 1 G; Conr 1]) % stability of the CLS
pole(Conr) % stability of the reduced cont.
tzero([ 1 G; Conr 1]) % stability of the CLS with reduced contr.

function Popt=ispr(A,B,C,D)
% checking the SPR property
[n,n1]=size(A);

setlmis([])
P=lmivar(1,[n,1]);

lmiterm([1 1 1 P],1,A,'s');
lmiterm([1 1 2 P],1,B);
lmiterm([1 1 2 0],-C');
lmiterm([1 2 2 0],-D-D');
lmiterm([-2 1 1 P],1,1);

LMIs=getlmis;

[val,xfeas]=feasp(LMIs);
Popt=dec2mat(LMIs,xfeas,P);
end