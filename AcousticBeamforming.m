clc; clear; close all;

c = 343; % Velocidad del Sonido
f = 1000; % Frecuencia 
lambda = c/f; % Longitud de Onda
res = lambda/32; % Resolución del Dominio
w = 2*pi*f; % Frecuencia Angular
k = w/c; % Numero de Onda
p0 = 1.21; % Densidad del Aire
q = 1; % Amplitud Compleja

yx = 6;
room = [yx,yx,0]; % Dimensiones del Dominio
Fuente = [room(1)/2,0,0];  % Posición de la Fuente

x_room = 0:res:room(1);
y_room = 0:res:room(2);
z_room = 0:res:room(3);

[X,Y,Z] = meshgrid(x_room-Fuente(1),y_room-Fuente(2),z_room-Fuente(3));

%% Receptores

N_r = 37; % Numero de Receptores
Radio = room(1)/2;  % Radio de los Receptores
Recs = (linspace(pi,0,N_r)).'; % Receptores en cordenadas polares

theta_beam = deg2rad(115); % Angulo de Dirección del Beam
theta_dif = abs(Recs-theta_beam);
[theta_min,~] = find(theta_dif == min(theta_dif));
p_m = zeros(N_r,1);
recs_Active = 1;
p_m(theta_min-recs_Active:1:theta_min+recs_Active) = 1;

idx0_pm = find(p_m == 0);
idx1_pm = find(p_m == 1);

[x_rec,y_rec,z_rec] = pol2cart(Recs,Radio,zeros(N_r,1));
Recs_final = [x_rec, y_rec, z_rec, p_m];

%% Fuentes

N_f = 16; % Numero de fuentes
F_dx = 0.0215; % Distancia entre Fuentes
M = zeros(N_f,3);
lentgh_array = N_f*F_dx;
M(1:N_f,1) = (((-lentgh_array+F_dx)/2) : F_dx : (lentgh_array/2));

r = zeros(N_r,N_f);
H = zeros(N_r,N_f);
for i = 1:N_r  
    x_src = Recs_final(i,1)-M(:,1);
    y_src = Recs_final(i,2)-M(:,2);
    z_src = Recs_final(i,3)-M(:,3);
    r(i,:) = sqrt(x_src.^2 + y_src.^2 + z_src.^2);
    H(i,:) = exp(-1j*k*r(i,:))./(4*pi*r(i,:));
end

%% Problema Inverso 

H_inv = pinv(H);
Q = H_inv * Recs_final(:,4);

Pt = 0;
rr = zeros(length(x_room),length(y_room),length(z_room));
for i = 1:N_f
    for iii = 1:length(Z(1,1,:))
        for i2 = 1:length(X(1,:,1))
            for ii = 1:length(Y(:,1,1))
                X_src = X(ii,i2,iii)-M(i,1);
                Y_src = Y(ii,i2,iii)-M(i,2);
                Z_src = Z(ii,i2,iii)-M(i,3);
                rr(ii,i2,iii) = sqrt(X_src.^2 + Y_src.^2 + Z_src.^2);
            end       
        end
    end
    Pt = Pt + ((1j*w*p0*q*exp(-1j*k*rr))./(4*pi*rr)) * Q(i);
end

% p_ref = 2e-5;
p_ref = p0;
Spl_Pt = 20*log10(abs(Pt)/p_ref);

%% Plot 

X_size=35;
Y_size=60;
width=1300;
height=600;

z_point = 1;

figure(1)
pcolor(X(:,:,z_point),Y(:,:,z_point),Spl_Pt(:,:,z_point))
shading interp;
c = colorbar;
c.Label.String = 'SPL (dB)';
colormap jet;
axis equal tight;
grid on;
grid minor;
set(gca,'FontSize',16,'LineWidth',2)
xlabel('x[m]')
ylabel('y[m]')
title(strcat('Beamforming',' (',num2str(f),'Hz',',',num2str(rad2deg(theta_beam)),'°',')'))
caxis([40 120])
hold on
h1 = scatter(M(:,1),M(:,2),60,'kd','filled');
h2 = scatter(Recs_final(idx1_pm,1),Recs_final(idx1_pm,2),60,'k*');
h3 = scatter(Recs_final(idx0_pm,1),Recs_final(idx0_pm,2),60,'ko','filled');
legend([h1 h2 h3],strcat('Sources',' (',num2str(N_f),')'),strcat('Bright point',' (',num2str(length(idx1_pm)),')'),strcat('Dark point',' (',num2str(length(idx0_pm)),')'));
yticks(0:0.99:room(2))
yticklabels(0:1:room(2))
set(gcf,'position',[X_size,Y_size,width,height])

%%
% figure(2)
% subplot(211),plot(abs(Q),'r','linewidth',2.5),axis tight,grid on,subplot(212),plot(unwrap(atan(imag(Q)/real(Q))),'b','linewidth',2.5),axis tight,grid on