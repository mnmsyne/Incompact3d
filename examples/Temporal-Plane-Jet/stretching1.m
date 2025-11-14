clear;clc;

%%
yly = 10;
ny = 601;
istret = 1;
beta = 2.25;
ncly1 = 1;
nclyn = 1;

%%
yinf = -yly/2;
den = 2*beta*yinf;
xnum = -yinf-sqrt(pi*pi*beta*beta+yinf*yinf);
alpha = abs(xnum/den);
xcx = 1/beta/alpha;

%%
nym = ny-1;
yp = zeros(ny,1);
yeta = zeros(ny,1);

for j = 2:ny
    yeta(j) = (j-1)*(1/nym);
    den1 = sqrt(alpha*beta+1);
    xnum = den1/sqrt(alpha/pi)/sqrt(beta)/sqrt(pi);
    den = 2*sqrt(alpha/pi)*sqrt(beta)*pi*sqrt(pi);
    den3 = ((sin(pi*yeta(j)))*(sin(pi*yeta(j)))/beta/pi)+alpha/pi;
    den4 = 2*alpha*beta-cos(2*pi*yeta(j))+1;
    xnum1 = (atan(xnum*tan(pi*yeta(j))))*den4/den1/den3/den;
    cst = sqrt(beta)*pi/(2*sqrt(alpha)*sqrt(alpha*beta+1));

    if (istret == 1)
       if (yeta(j) < 0.5) 
           yp(j) = xnum1-cst-yinf;
       end
       if (yeta(j) == 0.5) 
           yp(j) = 0-yinf;
       end
       if (yeta(j) > 0.5) 
           yp(j) = xnum1+cst-yinf;
       end
    end
end

%%
ypi = zeros(ny,1);
yetai = zeros(ny,1);

for j = 1:ny
    yetai(j) = (j-0.5)*(1/nym);
    den1 = sqrt(alpha*beta+1);
    xnum = den1/sqrt(alpha/pi)/sqrt(beta)/sqrt(pi);
    den = 2*sqrt(alpha/pi)*sqrt(beta)*pi*sqrt(pi);
    den3 = ((sin(pi*yetai(j)))*(sin(pi*yetai(j)))/beta/pi)+alpha/pi;
    den4 = 2*alpha*beta-cos(2*pi*yetai(j))+1;
    xnum1 = (atan(xnum*tan(pi*yetai(j))))*den4/den1/den3/den;
    cst = sqrt(beta)*pi/(2*sqrt(alpha)*sqrt(alpha*beta+1));

    if (istret == 1)
       if (yetai(j) < 0.5) 
           ypi(j) = xnum1-cst-yinf;
       end
       if (yetai(j) == 0.5) 
           ypi(j) = 0-yinf;
       end
       if (yetai(j) > 0.5) 
           ypi(j) = xnum1+cst-yinf;
       end
    end
end

% plot(1:j,yp,1:j,ypi)
ypi = ypi(1:nym);

%%
ppy = zeros(ny,1);
pp2y = zeros(ny,1);
pp4y = zeros(ny,1);

for j = 1:ny
    ppy(j) = yly*(alpha/pi+(1/pi/beta)*sin(pi*yeta(j))*sin(pi*yeta(j)));
    pp2y(j) = ppy(j)*ppy(j);
    pp4y(j) = (-2/beta*cos(pi*yeta(j))*sin(pi*yeta(j)));
end

%%
% yp2 = importdata('yp.dat');
% ypi2 = importdata('ypi.dat');

%% Hayashi et al., J. Fluid Mech. (2021)
alphay = 1.5;
j = 1:ny;
y = -yinf+yinf/alphay*atanh(tanh(alphay)*(1-2*(j-1)/(ny-1)));
yi = -yinf+yinf/alphay*atanh(tanh(alphay)*(1-2*(j-0.5)/(ny-1)));

plot(j,yp,j,y)
