% strong stabilization: SISO plant
clear

pom=[4, 13, -35, 22]; 
G=tf([1, -2, 2],[1, 1])*tf([1, -4, 5],pom);

% %%%%%%%%%%%%%%%%%%%%
n=4;
m=1;

[A,B,C,D]=ssdata(G);
[m1,m]=size(B);
%[p,p1]=size(C);

% transformation of the complex argument
bet=-0.03;
B0=B;
A0=(A-bet*eye(n))*(eye(n)-bet*A)^(-1);
C0=(1-bet^2)*C*(eye(n)-bet*A)^(-2);
D0=D+bet*C*(eye(n)-bet*A)^(-1)*B;
G0=ss(A0,B0,C0,D0); % transformed plant

 [X,L,F] = care(A0,B0,1*eye(n),10*eye(m)); % needed for param.stab.contr.
 % Checking
% eig(A0-B0*F) % stable
 % %%%%%%%%%%%%%%%%%%%
 
 bQ=ss(A0-B0*F,B0,-F,eye(m));
 bP=ss(A0-B0*F,B0,C0-D0*F,D0);
 
 s=tzero(bP)
 s1=s(1);
 s2=s(2);
 s3=s(3);
 s4=s(4);
 
 ell=4;
  
 [AM,BM,CM,DM]=ssdata(bP);
 x1=null1(DM+CM*(s1*eye(n)-AM)^(-1)*BM);
 x2=null1(DM+CM*(s2*eye(n)-AM)^(-1)*BM);
 x3=null1(DM+CM*(s3*eye(n)-AM)^(-1)*BM);
 x4=null1(DM+CM*(s4*eye(n)-AM)^(-1)*BM);
 

 [AY,BY,CY,DY]=ssdata(bQ);
 y1=(DY+CY*(s1*eye(n)-AY)^(-1)*BY)*x1;
 y2=(DY+CY*(s2*eye(n)-AY)^(-1)*BY)*x2;
 y3=(DY+CY*(s3*eye(n)-AY)^(-1)*BY)*x3;
 y4=(DY+CY*(s4*eye(n)-AY)^(-1)*BY)*x4; 

Api=diag([s1,s2,s3,s4]);
Cmin=[x1,x2,x3,x4];
Cpl=[y1,y2,y3,y4];
T0=[1 -complex(0,1); 1 complex(0,1)]/sqrt(2);
T=[T0 zeros(2,2); zeros(2,2) T0];
hApi=real(T'*Api*T);
hCmin=real(Cmin*T);
hCpl=real(Cpl*T);

[Qaux]=nult(hApi,hCmin);

Q(:,:,1)=hCpl;
for i=1:15
    i
[al(i),ind,Q(:,:,i+1)]=prv(hApi,Q(:,:,i), Qaux, hCmin); % Qaux (third argument) or hCmin
    if ind==1
    break
    end
end

H=1;
for ii=1:i+1
    Cplii=Q(:,:,ii);
    
    if ii==i+1
    Cminii=hCmin;
    else
    Cminii=Q(:,:,ii+1);
    end

XX=lyap(hApi',-(Cminii'*Cplii+Cplii'*Cminii)); % Pick matrix

% checking
e=eig(XX)
% %%%%%%%%

% SPR interpolant
% finding vector xi needed for reduced-order interpolant
      xi1=null((-s1'*eye(ell)+hApi')*XX);
      vr1=((Cminii-Cplii)*xi1) / ((Cminii+Cplii)*xi1);
      norm(vr1)
      
      xi2=null((-s2'*eye(ell)+hApi')*XX);   
      vr2=((Cminii-Cplii)*xi2) / ((Cminii+Cplii)*xi2);
      norm(vr2)
      
      xi3=null((-s3'*eye(ell)+hApi')*XX);
      vr3=((Cminii-Cplii)*xi3) / ((Cminii+Cplii)*xi3);
      norm(vr3)
      
      xi4=null((-s4'*eye(ell)+hApi')*XX);
      vr4=((Cminii-Cplii)*xi4) / ((Cminii+Cplii)*xi4);
      norm(vr4)
      
      if ii==7 || ii==8
      vr=vr3;
      sn=s3;
      xi=real(vr);
      eta=imag(vr);
      sigma=real(-sn');
      omega=imag(-sn');
      zeta=xi/eta;
      betaa=sqrt(sigma^2+omega^2+omega^2*zeta^2)-omega*zeta;
      alfa=vr*(-sn'+betaa)/(-sn'-betaa);
      alfa=real(alfa)
      u=alfa*tf([1,-betaa],[1,betaa]);
    
% Checking the int. condition on u      
      su1=alfa*(-sn'-betaa)/(-sn'+betaa);
      su1-vr
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
      Th=ss(hApi,XX^(-1)*[Cplii'-Cminii',Cplii'+Cminii'], [Cplii; Cminii], [-eye(m), eye(m); eye(m), eye(m)]);
        Th11=Th(1:m,1:m);
         Th12=Th(1:m,m+1:2*m);
          Th21=Th(m+1:2*m,1:m);
           Th22=Th(m+1:2*m,m+1:2*m);
      
      Hm=(Th11*u+Th12)*(Th21*u+Th22)^(-1);
      Hm=balreal(minreal(Hm));     
      else
Hm=dss(XX*hApi-(Cminii+Cplii)'*Cminii,(Cminii+Cplii)' , Cplii-Cminii, eye(m), XX);
      end

     % checking Hm
[AH,BH,CH,DH]=ssdata(Hm);
[AH,BH,CH,DH]=ssdata(Hm);
[eH,eH1]=size(AH);
pole(Hm)
tzero(Hm)
% checking interpol. condit. on Hm
pom1=Cplii*T';
pom2=Cminii*T';
(DH+CH*(s1*eye(eH)-AH)^(-1)*BH)*pom2(:,1)-pom1(:,1)
(DH+CH*(s2*eye(eH)-AH)^(-1)*BH)*pom2(:,2)-pom1(:,2)
(DH+CH*(s3*eye(eH)-AH)^(-1)*BH)*pom2(:,3)-pom1(:,3)
(DH+CH*(s4*eye(eH)-AH)^(-1)*BH)*pom2(:,4)-pom1(:,4)

% checking the SPR property of Hm
PP=ispr(AH,BH,CH,DH);
ePP=eig(PP)
% %%%%%%%%%%%%%%%%%%%%%%%%

H=balreal(H*Hm); % H is product of i+1 RMs Hm
end

% checking H
[AH,BH,CH,DH,EH]=dssdata(H);
[eH,eH1]=size(AH);
pole(H)
tzero(H)
% checking interpol. condit. on Hm
(DH+CH*(s1*EH-AH)^(-1)*BH)*Cmin(:,1)-Cpl(:,1)
(DH+CH*(s2*EH-AH)^(-1)*BH)*Cmin(:,2)-Cpl(:,2)
(DH+CH*(s3*EH-AH)^(-1)*BH)*Cmin(:,3)-Cpl(:,3)
(DH+CH*(s4*EH-AH)^(-1)*BH)*Cmin(:,4)-Cpl(:,4)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Cons=(bQ-H)*bP^(-1);
%Contr=balreal(minreal(minreal(Cons,1e-8)));
Contr=balreal((minreal(stabsep(minreal(Cons))))); %same result with minreal

[AA, BB, CC, DD, EE]=dssdata(Contr);

pp=pole(Contr)  %  kontr. is stable
zz=tzero([ eye(m) G0; Contr eye(m)]) % i CLS is stable

% checking
[pomn,xyz]=size(zz);
pom=zeros(pomn,1);
for iii=1:pomn
pom(iii)=(zz(iii)+bet)/(1+bet*zz(iii));
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% inverse argument transformation of the contr.
[AC,BC,CC,DC]=ssdata(Contr);
[nC,nC1]=size(AC);
betm=-bet;
B00=BC;
A00=(AC-betm*eye(nC))*(eye(nC)-betm*AC)^(-1);
C00=(1-betm^2)*CC*(eye(nC)-betm*AC)^(-2);
D00=DC+betm*CC*(eye(nC)-betm*AC)^(-1)*BC;

% Controller
Con=ss(A00,B00,C00,D00); 

% checking the stab. of contr. and of CLS with Con
pp1=pole(Con)
zz1=tzero([ eye(m) G; Con eye(m)])

% reductiom of controller
Conr=reduce(Con);

% checking the stab. of contr. and of CLS with Conr
ppr=pole(Conr)  %  kontr. is stable
zzr=tzero([ eye(m) G; Conr eye(m)]) % CLS is stable


function [TH,ind,Qopt]=prv(A,Q,Qaux,hCmin)

[m,el]=size(hCmin);

setlmis([])
X=lmivar(1,[el,1]);
X1=lmivar(1,[el,1]);
Y=lmivar(2,[m,el]);
al=lmivar(1,[1,1]);

lmiterm([1 1 1 X],1,A,'s');
lmiterm([1 1 1 Y],Q',-1,'s');

lmiterm([-2 1 1 X],1,1);

lmiterm([3 1 1 X1],1,A,'s');
lmiterm([3 1 1 Y],hCmin',-1,'s');

lmiterm([-4 1 1 X1],1,1);

LMIs=getlmis;

Nn = decnbr(LMIs); 
c = zeros(Nn,1);

for jj=1:Nn, 
	[tauj] = defcx(LMIs,jj,al); 
	c(jj) = tauj; 
end

options=[0,0,0,0,0];
output1 = evalc('[TH,xfeas] = mincx(LMIs,c,options);');
%Xopt=dec2mat(LMIs,xfeas,X); 

% Check feasibility
if isempty(TH)
    ind=0; % disp('No feasible solution found.');
elseif isinf(TH) || isnan(TH)
    ind=0; % disp('Optimization failed or returned an invalid cost.');
else
    ind=1; % disp('Feasible solution found.');
    Qopt=dec2mat(LMIs,xfeas,Y);
    Qopt=Qopt/norm(Qopt);
    return
end

% if the previous code was not succesfull (ind=0):

setlmis([])
X=lmivar(1,[el,1]);
Y=lmivar(2,[m,el]);
al=lmivar(1,[1,1]);

lmiterm([1 1 1 X],1,A,'s');
lmiterm([1 1 1 Y],Q',-1,'s');

lmiterm([-2 1 1 X],1,1);

lmiterm([-3,1,1,al],1,1);
lmiterm([-3,1,2,Y],1,1);
lmiterm([-3,1,2,0],-Qaux);
lmiterm([-3,2,2,al],1,1);

LMIs=getlmis;

Nn = decnbr(LMIs); 
c = zeros(Nn,1);

for jj=1:Nn, 
	[tauj] = defcx(LMIs,jj,al); 
	c(jj) = tauj; 
end

options=[0,0,0,0,0];
output1 = evalc('[TH,xfeas] = mincx(LMIs,c,options);');
Xopt=dec2mat(LMIs,xfeas,X);

Qopt=dec2mat(LMIs,xfeas,Y);

end


function [Qaux]=nult(A,hCmin)

[m,el]=size(hCmin);

setlmis([])
X=lmivar(1,[el,1]);
Y=lmivar(2,[m,el]);
al=lmivar(1,[1,1]);

lmiterm([1 1 1 X],1,A,'s');
lmiterm([1 1 1 Y],hCmin',-1,'s');

lmiterm([-2 1 1 X],1,1);
lmiterm([-2,1,1,al],1,1);

lmiterm([-3,1,1,0],1);
lmiterm([-3,1,2,Y],1,1);
lmiterm([-3,2,2,0],1);

LMIs=getlmis;

Nn = decnbr(LMIs); 
c = zeros(Nn,1);

for jj=1:Nn, 
	[tauj] = defcx(LMIs,jj,al); 
	c(jj) = tauj; 
end

options=[0,0,0,0,0];
output1 = evalc('[TH,xfeas] = mincx(LMIs,c,options);');
Qaux=dec2mat(LMIs,xfeas,Y);

% checking the eig(Xopt)

%Xopt=dec2mat(LMIs,xfeas,X);
%e=eig(Xopt)
%comparison with Qaux=hCmin
%XX=lyap(A',-2*hCmin'*hCmin); % Pick matrix
%eXX=eig(XX)

end


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

output1 = evalc('[val,xfeas]=feasp(LMIs);');
Popt=dec2mat(LMIs,xfeas,P);
end