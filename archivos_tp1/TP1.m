%definimos dH/dt = f
% x1=H  x1'=H'=f


syms la lb L Qi S u g x1

optionss = bodeoptions;
optionss.PhaseMatching = 'on';
optionss.PhaseMatchingFreq = 1;
optionss.PhaseMatchingValue = -180;
optionss.Grid = 'on';

f= (Qi - S*u*sqrt(2*g*x1))/(lb + (la - lb)*x1/L);
x=x1;
y=x1;

% definimos valores conocidos

Qi=8/(1000*60);  %metros^3/segundo
d_tuberia= 10.65e-3; % metros
lb = 10e-2;  % metros
la= 40e-2; % metros
S= pi*((d_tuberia/2)^2); % metros^2
L=0.9; % metros
g=9.81; % metros/seg^2

f= (Qi - S*u*sqrt(2*g*x1))/(lb + (la - lb)*x1/L);

% valores equilibrio x1e=0.45;
x1e=0.45;
ue=solve(Qi - S*u*sqrt(2*g*x1e)==0,u);
ue=double(ue);


%%
% valores de equilibrio
x1e=0.45;
ue=0.5037;
ye=0.45;

A=jacobian(f,x);
B=jacobian(f,u);
C=jacobian(y,x);
D=jacobian(y,u);

A=subs(A,{'x1','u','y'},{x1e,ue,ye});
B=subs(B,{'x1','u','y'},{x1e,ue,ye});
C=subs(C,{'x1','u','y'},{x1e,ue,ye});
D=subs(D,{'x1','u','y'},{x1e  ,ue,ye});

Ass=double(A);
Bss=double(B);
Css=double(C);
Dss=double(D);

G= zpk(ss(Ass,Bss,Css,Dss));


%% Linealización a distintas alturas
clc;
Qi=8/(1000*60);  %metros^3/segundo

% he=0.1;
G1= linealizar_Fun(0.1,valores_Equilibrio(0.1,Qi),0.1,Qi);

%he = 0.2
G2= linealizar_Fun(0.2,valores_Equilibrio(0.2,Qi),0.2,Qi);

% he= 0.3
G3= linealizar_Fun(0.3,valores_Equilibrio(0.3,Qi),0.3,Qi);

% he= 0.4
G4= linealizar_Fun(0.4,valores_Equilibrio(0.4,Qi),0.4,Qi);

% he= 0.5
G5= linealizar_Fun(0.5,valores_Equilibrio(0.5,Qi),0.5,Qi);

% he= 0.6
G6= linealizar_Fun(0.6,valores_Equilibrio(0.6,Qi),0.6,Qi);

% he= 0.7
G7= linealizar_Fun(0.7,valores_Equilibrio(0.7,Qi),0.7,Qi);

% he= 0.8
G8= linealizar_Fun(0.8,valores_Equilibrio(0.8,Qi),0.8,Qi);

% BODE de las linealizaciones a lazo abierto
figure ()
hold on 
bode(G1);
bode(G2);
bode(G3);
bode(G4);
bode(G5);
bode(G6);
bode(G7);
bode(G8);
legend();
title('Bode de las G')
hold off


%% Diseño del controlador 
k=db2mag(-140);
%Probé con un polo en -0.01, y fui desplazandolo a la derecha ligeramente. 
%Se utilizó un PI. + red de atraso de fase
C=(1/G)* zpk([-1000],[0, -0.016],k);
L=C*G; %Lazo abierto;
H=L/(1+L); %Lazo cerrado;
S=1/(L+1);
T=1-S;


% GRAFICOS
subplot(2,2,1)
bode(L,optionss);
title('Lazo abierto')

subplot(2,2,2)
step(H)
stepinfo(H)
title('Step del lazo cerrado')

subplot(2,2,3)
bode(T,optionss);
title('Sensibilidad complementaria')


%% BODE de lazo abierto, con controlador.

H1=C*G1;
H2=C*G2;
H3=C*G3;
H4=C*G4;
H5=C*G5;
H6=C*G6;
H7=C*G7;
H8=C*G8;

% % GRAFICOS
% figure ()
% hold on 
% bode(H1);
% bode(H2);
% bode(H3);
% bode(H4);
% bode(H5);
% bode(H6);
% bode(H7);
% bode(H8);
% legend();
% title("Bode de las H=G*C");
% hold off

%% BODE de la sensibilidad complementaria (siempre a lazo cerrado)
%Se aprecia que aumenta el ancho de banda a medida que el punto de
%equilibrio es una altura más pequeña, lo que implica una respuesta más
%rápida del sistema a esas alturas.
T1=1 - 1/(1+H1);
T2=1 - 1/(1+H2);
T3=1 - 1/(1+H3);
T4=1 - 1/(1+H4);
T5=1 - 1/(1+H5);
T6=1 - 1/(1+H6);
T7=1 - 1/(1+H7);
T8=1 - 1/(1+H8);

% % GRAFICOS
% figure ()
% title("Bode de las T");
% hold on 
% bode(T1);
% bode(T2);
% bode(T3);
% bode(T4);
% bode(T5);
% bode(T6);
% bode(T7);
% bode(T8);
% legend();
% hold off



function ue = valores_Equilibrio(h1e,Qi)
    syms u
    d_tuberia= 10.65e-3; % metros
    lb = 10e-2;  % metros
    la= 40e-2; % metros
    S= pi*((d_tuberia/2)^2); % metros^2
    L=0.9; % metros
    g=9.81; % metros/seg^2

    ue=solve(Qi - S*u*sqrt(2*g*h1e)==0,u);
    ue=double(ue);
end

function G=linealizar_Fun(x1e,ue,ye,Qi)
    syms x1 u
    d_tuberia= 10.65e-3; % metros
    lb = 10e-2;  % metros
    la= 40e-2; % metros
    S= pi*((d_tuberia/2)^2); % metros^2
    L=0.9; % metros
    g=9.81; % metros/seg^2
    f= (Qi - S*u*sqrt(2*g*x1))/(lb + (la - lb)*x1/L);

    x=x1;
    y=x1;
    A=jacobian(f,x);
    B=jacobian(f,u);
    C=jacobian(y,x);
    D=jacobian(y,u);

    A=subs(A,{'x1','u','y'},{x1e,ue,ye});
    B=subs(B,{'x1','u','y'},{x1e,ue,ye});
    C=subs(C,{'x1','u','y'},{x1e,ue,ye});
    D=subs(D,{'x1','u','y'},{x1e,ue,ye});

    Ass=double(A);
    Bss=double(B);
    Css=double(C);
    Dss=double(D);

    G= zpk(ss(Ass,Bss,Css,Dss));
end



