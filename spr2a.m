% Finding strong controller for the plant G
clear
addpath('C:\Users\jovans\Desktop\hanso2_0');
addpath('C:\Users\jovans\Desktop\HIFOO3.501');

A=-[2 -3 -2 1 -1 0 1 2 -1 3 ; 
   0 0 -2 1 -1 -2 3 1 1 0  ;
   -2 1 3 2 -3 0 -1 1 3 3  ;
   -3 1 0 -2 0 0 -1 2 3 1  ;
    1 0 -1 0 1 -3 -1 0 0 -1  ;
    -2 3 2 0 -1 0 1 2 -1 -3  ; 
    2 -1 -3 -2 -3 1 0 -1 3 -3  ;
    3 -1 0 -2 0 1 -1 2 3 1  ; 
    2 0 1 3 -1 3 -1 1 0 -1 ;
    1 0 -1 0 1 -3 1 -1 0 1 ];
B=-[ -2 -1 ; 1 3 ; -1 0 ; -2 3 ; 2 0 ; 
    -1 2 ; -2 -3 ; 0 1 ; -1 0 ; -1 3 ];
    
C=-[ -3 0 -1 3 2 1 -2 0 1 -3  ; 
    -2 3 0 -1 -1 -2 -1 -3 1 0 ];
D=[0 -1 ; 0 1];

tic
% application of HIFOO for finding strong controller
PP.A=A;
PP.B=B;
PP.C=C;
PP.D=D;

Psys={PP,'K'};

ORDER=18; % 25; % order of the required controller
INIT=[];
FUN='ss';
UPPERBND=[-0.00,-0.00];
OPTIONS.nrand=100; % how many initial controllers are considered
[K,val,viol]=hifoo(Psys,ORDER,INIT,FUN,UPPERBND,OPTIONS);

AK=K.a;
BK=K.b;
CK=K.c;
DK=K.d;

pom=[A+B*DK*C, B*CK; BK*C, AK];

eig(AK) % poles of the controller
eig(pom) % modes of the closed-loop system
toc