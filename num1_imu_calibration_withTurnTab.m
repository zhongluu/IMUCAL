close all;
clear;
load("1#inner.mat");

refRot = [0, 5, 10, 15, 20, 25, -5, -10, -15, -20, -25]; % deg/s
refRot1 = [0, 5, 10, 15, 20, 25, 25/2 - 2.5, -5, -10, -15, -20, -25]; % deg/s

refRot2 = [5,10,15,20,25];


%% 内框转动 gyroY 标定Ky、Exy、Ezy
Yaxis = gyro_y(1:157431);
Xaxis = gyro_x(1:157431);
Zaxis = gyro_z(1:157431);

innerAccXaxis = acc_x(157431:end);
innerAccYaxis = acc_y(157431:end);
innerAccZaxis = acc_z(157431:end);

Index = zeros(1,length(refRot1) + 1);

deltaInt = 500;
j=1;
for i =1:length(Yaxis)
    if j <= 6
    if Yaxis(i) > (refRot1(j) - 2.5)
        Index(j) = i;
        j = j + 1;
        if j > length(refRot1)
            break;
        end
    end
    else
    if Yaxis(i) < (refRot1(j) + 2.5)
        Index(j) = i;
        j = j + 1;
        if j > length(refRot1)
            break;
        end
    end
    end
end
Index = Index(1:end-1);
Lindex = zeros(1,12);
Rindex = zeros(1,12);
j = 1;
for i = 2 : length(Index)
   Rindex(i) = Index(i) - 200;
   Lindex(i) = Index(i) + 200;
end
Lindex(1) = 1;
Rindex(1:end -1) =  Rindex(2:end);
Rindex(end) = length(Yaxis);
actRotY = zeros(1, length(Lindex));
actRotX = zeros(1, length(Lindex));
actRotZ = zeros(1, length(Lindex));
for i =1:length(Lindex)
    actRotY(i) = mean(Yaxis(Lindex(i):Rindex(i)));
    actRotX(i) = mean(Xaxis(Lindex(i):Rindex(i)));
    actRotZ(i) = mean(Zaxis(Lindex(i):Rindex(i)));
end
% Dx = actRot';
Dy = [actRotY(1:6)';actRotY(8:end)'];
wy = refRot';
Omega = [ones(size(wy)),wy];
tmpM = (Omega'*Omega)\Omega'*Dy;
Dy0 = tmpM(1); Ky = tmpM(2);
wy = refRot2';
Dxp = actRotX(2:6)';
Dxm = actRotX(8:end)';
Dzp = actRotZ(2:6)';
Dzm = actRotZ(8:end)';

Exy = (Dxp - Dxm)./(2*Ky*wy);
Exy = mean(Exy);

Ezy = (Dzp - Dzm)./(2*Ky*wy);
Ezy = mean(Ezy);



clearvars -except Ky Exy Ezy innerAccXaxis innerAccYaxis innerAccZaxis Dy0
%% 中框转动 gyroX 标定Kx、Eyx、Ezx
load("1#mid.mat");

% refRot = [0, 5, 10, 15, 20, 25, -5, -10, -15, -20, -25]; % deg/s
% refRot1 = [0, 5, 10, 15, 20, 25, 25/2 - 2.5, -5, -10, -15, -20, -25]; % deg/s
 refRot = [0, 5, 10, 15, 20, 25, -5, -10, -15, -20, -25]; % deg/s
refRot1 = [0, 5, 10, 15, 20, 25, 25/2 - 2.5, -5, -10, -15, -20, -25]; % deg/s
refRot = -refRot;
refRot1 = -refRot1;
refRot2 = [5,10,15,20,25];

Yaxis = gyro_y(1:110795);
Xaxis = gyro_x(1:110795);
Zaxis = gyro_z(1:110795);

midAccXaxis = acc_x(110795:end);
midAccYaxis = acc_y(110795:end);
midAccZaxis = acc_z(110795:end);

Index = zeros(1,length(refRot1) + 1);

deltaInt = 500;
j=1;
for i =1:length(Xaxis)
    if j <= 6
    if Xaxis(i) < (refRot1(j) + 2.5)
        Index(j) = i;
        j = j + 1;
        if j > length(refRot1)
            break;
        end
    end
    else
    if Xaxis(i) > (refRot1(j) - 2.5)
        Index(j) = i;
        j = j + 1;
        if j > length(refRot1)
            break;
        end
    end
    end
end
Index = Index(1:end-1);
Lindex = zeros(1,12);
Rindex = zeros(1,12);
j = 1;
for i = 2 : length(Index)
   Rindex(i) = Index(i) - 200;
   Lindex(i) = Index(i) + 200;
end
Lindex(1) = 1;
Rindex(1:end -1) =  Rindex(2:end);
Rindex(end) = length(Xaxis);
actRotY = zeros(1, length(Lindex));
actRotX = zeros(1, length(Lindex));
actRotZ = zeros(1, length(Lindex));
for i =1:length(Lindex)
    actRotY(i) = mean(Yaxis(Lindex(i):Rindex(i)));
    actRotX(i) = mean(Xaxis(Lindex(i):Rindex(i)));
    actRotZ(i) = mean(Zaxis(Lindex(i):Rindex(i)));
end
% Dx = actRot';
Dx = [actRotX(1:6)';actRotX(8:end)'];
wx = refRot';
Omega = [ones(size(wx)),wx];
tmpM = (Omega'*Omega)\Omega'*Dx;
Dx0 = tmpM(1); Kx = tmpM(2);
wx = refRot2';

Dyp = actRotY(8:end)';
Dym = actRotY(2:6)';
Dzp = actRotZ(8:end)';
Dzm = actRotZ(2:6)';

Eyx = (Dyp - Dym)./(2*Kx*wx);
Eyx = mean(Eyx);

Ezx = (Dzp - Dzm)./(2*Kx*wx);
Ezx = mean(Ezx);

% 标定 Ka
g = 9.7803267715;
Ax1 = mean(midAccXaxis(40506:69078)); Ay1 = mean(midAccYaxis(40506:69078)); Az1 = mean(midAccZaxis(40506:69078));
Ax2 = mean(innerAccXaxis(71326:98388)); Ay2 = mean(innerAccYaxis(71326:98388)); Az2 = mean(innerAccZaxis(71326:98388));
Ax3 = mean(innerAccXaxis(103951:130337)); Ay3 = mean(innerAccYaxis(103951:130337)); Az3 = mean(innerAccZaxis(103951:130337)); 
Ax4 = mean(midAccXaxis(8701:37603)); Ay4 = mean(midAccYaxis(8701:37603)); Az4 = mean(midAccZaxis(8701:37603));
Ax5 = mean(innerAccXaxis(7620:37233)); Ay5 = mean(innerAccYaxis(7620:37233)); Az5 = mean(innerAccZaxis(7620:37233)); 
Ax6 = mean(innerAccXaxis(42420:68219)); Ay6 = mean(innerAccYaxis(42420:68219)); Az6 = mean(innerAccZaxis(42420:68219)); 

Ax0 = (Ax1 + Ax2 + Ax3 + Ax4) / 4;
Ay0 = (Ay2 + Ay3 + Ay5 + Ay6) / 4;
Az0 = (Az1 + Az4 + Az5 + Az6) / 4;

AKxp = (Ax5 - Ax0) / g; AKxm = (-Ax6 + Ax0) / g; Kax = g*(AKxp + AKxm)/2;
AKyp = (Ay4 - Ay0) / g; AKym = (-Ay1 + Ay0) / g; Kay = g*(AKyp + AKym)/2;
AKzp = (Az3 - Az0) / g; AKzm = (-Az2 + Az0) / g; Kaz = g*(AKzp + AKzm)/2;

Fyx = (Ay5 - Ay6)/((AKyp + AKym) * g); Fzx = (Az5 - Az6)/((AKzp + AKzm) * g); 
Fxy = (Ax4 - Ax1)/((AKxp + AKxm) * g); Fzy = (Az4 - Az1)/((AKzp + AKzm) * g); 
Fxz = (Ax3 - Ax2)/((AKxp + AKxm) * g); Fyz = (Ay3 - Ay2)/((AKyp + AKym) * g); 

Ka = [Kax,     Kax*Fxy, Kax*Fxz;
      Kay*Fyx, Kay,     Kay*Fyz;
      Kaz*Fzx, Kaz*Fzy, Kaz];


offsetAcc = [Ax0; Ay0; Az0];

clearvars -except Ky Exy Ezy Kx Eyx Ezx Dy0 Dx0 Ka offsetAcc

%% 外框转动 gyroZ 标定Kz、Exz、Eyz

load("1#outer.mat");

% refRot = [0, 5, 10, 15, 20, 25, -5, -10, -15, -20, -25]; % deg/s
% refRot1 = [0, 5, 10, 15, 20, 25, 25/2 - 2.5, -5, -10, -15, -20, -25]; % deg/s
 refRot = [0, 5, 10, 15, 20, 25, -5, -10, -15, -20, -25]; % deg/s
refRot1 = [0, 5, 10, 15, 20, 25, 25/2 - 2.5, -5, -10, -15, -20, -25]; % deg/s
refRot = -refRot;
refRot1 = -refRot1;
refRot2 = [5,10,15,20,25];

Yaxis = gyro_y(1:105348);
Xaxis = gyro_x(1:105348);
Zaxis = gyro_z(1:105348);

Index = zeros(1,length(refRot1) + 1);

deltaInt = 500;
j=1;
for i =1:length(Zaxis)
    if j <= 6
    if Zaxis(i) < (refRot1(j) + 2.5)
        Index(j) = i;
        j = j + 1;
        if j > length(refRot1)
            break;
        end
    end
    else
    if Zaxis(i) > (refRot1(j) - 2.5)
        Index(j) = i;
        j = j + 1;
        if j > length(refRot1)
            break;
        end
    end
    end
end
Index = Index(1:end-1);
Lindex = zeros(1,12);
Rindex = zeros(1,12);
j = 1;
for i = 2 : length(Index)
   Rindex(i) = Index(i) - 200;
   Lindex(i) = Index(i) + 200;
end
Lindex(1) = 1;
Rindex(1:end -1) =  Rindex(2:end);
Rindex(end) = length(Xaxis);
actRotY = zeros(1, length(Lindex));
actRotX = zeros(1, length(Lindex));
actRotZ = zeros(1, length(Lindex));
for i =1:length(Lindex)
    actRotY(i) = mean(Yaxis(Lindex(i):Rindex(i)));
    actRotX(i) = mean(Xaxis(Lindex(i):Rindex(i)));
    actRotZ(i) = mean(Zaxis(Lindex(i):Rindex(i)));
end
% Dx = actRot';
Dz = [actRotZ(1:6)';actRotZ(8:end)'];
wz = refRot';
Omega = [ones(size(wz)),wz];
tmpM = (Omega'*Omega)\Omega'*Dz;
Dz0 = tmpM(1); Kz = tmpM(2);
wz = refRot2';

Dxp = actRotX(8:end)';
Dxm = actRotX(2:6)';
Dyp = actRotY(8:end)';
Dym = actRotY(2:6)';

Exz = (Dxp - Dxm)./(2*Kz*wz);
Exz = mean(Exz);

Eyz = (Dyp - Dym)./(2*Kz*wz);
Eyz = mean(Eyz);


Kg = [Kx,     Kx*Exy,  Kx*Exz;
      Ky*Eyx, Ky,      Ky*Eyz;
      Kz*Ezx, Kz*Ezy,  Kz];
offsetGyro = [Dx0; Dy0; Dz0];
% Kg = inv(Kg);
time = 1/100*(1:length(gyro_x));
figure; 
subplot(311); plot(time, gyro_x); title('gyro_x');
subplot(312); plot(time, gyro_y); title('gyro_y');
subplot(313); plot(time, gyro_z); title('gyro_z');
sgtitle('before calibration');

tmpGyro = [gyro_x, gyro_y,gyro_z];
tmpGyro = Kg \ tmpGyro' - offsetGyro;
tmpGyro = tmpGyro';
figure; 
subplot(311); plot(time, tmpGyro(:,1)); title('gyro_x');
subplot(312); plot(time, tmpGyro(:,2)); title('gyro_y');
subplot(313); plot(time, tmpGyro(:,3)); title('gyro_z');
sgtitle('after calibration');

figure; 
subplot(311); plot(time, acc_x); title('acc_x');
subplot(312); plot(time, acc_y); title('acc_y');
subplot(313); plot(time, acc_z); title('acc_z');
sgtitle('before calibration');

tmpAcc = [acc_x, acc_y, acc_z];
tmpAcc = Ka \ tmpAcc' - offsetAcc;
tmpAcc = tmpAcc';
figure; 
subplot(311); plot(time, tmpAcc(:,1)); title('acc_x');
subplot(312); plot(time, tmpAcc(:,2)); title('acc_y');
subplot(313); plot(time, tmpAcc(:,3)); title('acc_z');
sgtitle('after calibration');

disp('Gyro. Scale-factor matrix:');
disp(inv(Kg));
disp('Gyro. Offset:');
disp(-offsetGyro);

disp('Acc. Scale-factor matrix:');
disp(inv(Ka));
disp('Acc. Offset:');
disp(-offsetAcc);

fd = fopen('num1.txt','w');
fprintf(fd,'Gyro. Scale-factor matrix:\n');
tmpK = inv(Kg);
[m,n] = size(tmpK);
for i = 1:m
    for j = 1:n
        if j==n
            fprintf(fd,'%ld\n',tmpK(i,j));
        else
            fprintf(fd,'%ld ',tmpK(i,j));
        end
    end
end
fprintf(fd,'Gyro. Offset:\n');
tmpD = -offsetGyro;
[m,n] = size(tmpD);
for i = 1:m
    for j = 1:n
        if j==n
            fprintf(fd,'%ld\n',tmpD(i,j));
        else
            fprintf(fd,'%ld ',tmpD(i,j));
        end
    end
end
fprintf(fd,'\n');

fprintf(fd,'Acc. Scale-factor matrix:\n');
tmpK = inv(Ka);
[m,n] = size(tmpK);
for i = 1:m
    for j = 1:n
        if j==n
            fprintf(fd,'%ld\n',tmpK(i,j));
        else
            fprintf(fd,'%ld ',tmpK(i,j));
        end
    end
end
fprintf(fd,'Acc. Offset:\n');
tmpD = -offsetAcc;
[m,n] = size(tmpD);
for i = 1:m
    for j = 1:n
        if j==n
            fprintf(fd,'%ld\n',tmpD(i,j));
        else
            fprintf(fd,'%ld ',tmpD(i,j));
        end
    end
end
