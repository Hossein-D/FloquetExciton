function [eigVal] = FloquetExciton(m, v, omega, A0, coulU, alpha, pmax, Nk, Gamma, NB, changeVar)
% In this function we consider infinite momentum summation to reproduce the
% results of the 2D hydrogen atom. To do the summation we apply a
% transformation of variables to the momentum to create a finite range for
% the transformed momenta. 


deltaP = 2*pmax/(Nk-1);

pvec = -pmax:deltaP:pmax;
pxVec = -pmax:deltaP:pmax;
pyVec = -pmax:deltaP:pmax;

%kres = 1/(2*sqrt(v^2))*(omega^2 - 4*m^2)^(1/2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[epsilond, zetak] = SemiCondDisp(m, v, omega, kxVec, kyVec, A0, NB, Gamma);
%display(epsilond);

%omegad = omega - 2*m;

gammav = NB/(1 + 2*NB);
gammac = (1+NB)/(1 + 2*NB);
%parameters = [kres, gammav, gammac];

syms x
roots = vpasolve(legendreP(Nk,x) == 0);
%roots(1:10)
%roots(end-10:end)

if Nk>=10
    tanRoot = eval(tan(pi/2*roots(1:10)));
else
    tanRoot = eval(tan(pi/2*roots(1:end)));
end

weightsLgndr = zeros(Nk, 1);
weightsFun = zeros(Nk, 1);
xRoots = zeros(Nk, 1);
rootFun = zeros(Nk, 1);
fVec = zeros(Nk, 1);
diffLeg = zeros(Nk, 1);
a=-1;b=1;

for n=1:Nk
    xRoots(n) = (a+b)/2 + (b-a)/2*(roots(n));
    diffLeg(n) = vpa(subs(diff(legendreP(Nk, x),x,1),x,xRoots(n)));
    weightsLgndr(n) = 2/((1-xRoots(n)^2)*diffLeg(n)^2);
    %fVec(n) = f(xRoots(n));
end

xRootsTr = xRoots';
weightsTr = weightsLgndr';

excitonMat = zeros(Nk*Nk, Nk*Nk);
%excitonMatSmooth = zeros(Nk*Nk, Nk*Nk);
epsilondOfk = @(kx, ky) (tan(pi/2*kx).^2+tan(pi/2*ky).^2);
epsilond = zeros(Nk, Nk);
a=-1;b=1;

%excitonMat = zeros((Nk+1)*(Nk+1), (Nk+1)*(Nk+1));

Rabi = v*A0;
derivVar = zeros(Nk, 1);
xRoots = xRoots;

if changeVar
    % qx is equivalent to the rootFun=xroots
    momFun = @(qx) tan(pi/2*qx);
    weightsFun = weightsLgndr;
    delMomCoeff = pi/2;
    rootFun = xRoots;     
    derivVar = @(qx) 1+(tan(pi/2*qx)).^2;
else
    momFun = @(kx) kx; 
    weightsFun(n) = 1;
    delMomCoeff = deltaP;
    rootFun = pvec;
    derivVar = @(kx) 1;
end

dOfk = @(qx, qy) 2*m + v^2*(tan(pi/2*qx).^2 + tan(pi/2*qy).^2)/(2*m);

epsilondOfk = @(qxVec, qyVec) (dOfk(qxVec, qyVec) - omega);

%epsilond = epsilondOfk(kxVec, kyVec);

RabiRWASqOfk = @(qxVec, qyVec) (1 - v^2*(momFun(qxVec).^2 + momFun(qyVec).^2)/(4*m^2))*Rabi^2;

RabiRWASqOfk = @(qx, qy) (1 - v^2*(tan(pi/2*qx).^2 + tan(pi/2*qy).^2)/(4*m^2))*Rabi^2;

%RabiRWASq = RabiRWASqOfk(kxVec, kyVec);


zetak = @(qxVec, qyVec) RabiRWASqOfk(qxVec, qyVec)./(epsilondOfk(qxVec, qyVec).^2 ...
    + (1/2 + NB)^2*Gamma^2);
epsilond = zeros((Nk), (Nk));


deltaNcNv = zeros(Nk, Nk);
delMomCoeff = delMomCoeff
for nx = 1:Nk
    for ny = 1:Nk
        epsilond(nx, ny) = epsilondOfk(rootFun(nx), rootFun(ny));
        excitonMat(ny+(Nk)*(nx-1), ny+(Nk)*(nx-1)) = epsilondOfk(rootFun(nx),rootFun(ny)) ...
            +1i*Gamma*(1/2+NB);
        %excitonMatSmooth(ny+(Nk)*(nx-1), ny+(Nk)*(nx-1)) = epsilondOfk(xRoots(nx),xRoots(ny)) ...
        %    +1i*Gamma*(1/2+NB) + C*sqrt(momFun(xRoots(nx))^2 + momFun(xRoots(ny))^2);
        deltaNcNv(nx, ny) = (gammac-gammav) / (2*zetak(momFun(rootFun(nx)), momFun(rootFun(ny))) ...
            + gammac + gammav);
        
        for nxp = 1:Nk
            for nyp = 1:Nk
                if pxVec(nx)~=pxVec(nxp) || pyVec(ny)~=pyVec(nyp)
                    potential = coulU/sqrt((momFun(rootFun(nx))-momFun(rootFun(nxp)))^2 + ...
                        (momFun(rootFun(ny))-momFun(rootFun(nyp)) )^2 + alpha^2);
                    
                    excitonMat(ny+Nk*(nx-1), nyp+Nk*(nxp-1)) = -1/(2*pi) * ... %deltaP^2 * ...
                        delMomCoeff^2 * weightsFun(nxp) * weightsFun(nyp) * ...
                        derivVar(rootFun(nxp)) * derivVar(rootFun(nyp)) * ...
                        (gammac-gammav) / (2*zetak((rootFun(nx)), (rootFun(ny))) + gammac + gammav) ...
                        * potential;
                end
                
                %{
                if nx ~= nxp || ny ~= nyp
                    excitonMat(ny+(Nk)*(nx-1), nyp+(Nk)*(nxp-1)) = -1/(2*pi) * ...
                        (pi/2)^2 * weights(nxp) * weights(nyp) * ...
                        (1+momFun(pi/2*xRoots(nxp))^2) * (1+momFun(pi/2*xRoots(nyp))^2) * ...
                        (gammac-gammav)/(2*zetak(kxVec(nx), kyVec(ny))+gammac+gammav) * ...
                        2/sqrt( (momFun(pi/2*xRoots(nx))-momFun(pi/2*xRoots(nxp)))^2 + ...
                        (momFun(pi/2*xRoots(ny))-momFun(pi/2*xRoots(nyp)))^2 + alpha^2);
                end
                %}
            end
        end
    end
end

Rydberg = coulU^4;
eyeMatrix = 2*Rydberg*eye(Nk*Nk, Nk*Nk);


largestEigs=20;

[eigVec, eigVal] = eig(excitonMat);                  
eigVal = sort(diag(real(eigVal)));

eigVal1to20 = eigVal(1:largestEigs)';
%eigValSmooth1to20 = eigValSmooth(1:largestEigs)'
end





function [epsilond, zetak] = SemiCondDisp(m, v, omega, kxVec, kyVec, A0, NB, Gamma)
Rabi = v*A0;
dOfk = @(kxVec, kyVec) sqrt(v^2*(kxVec.^2 + kyVec.^2) + m^2);
dzOfk = @(kxVec, kyVec) m;
epsilondOfk = @(kxVec, kyVec) abs(2*dOfk(kxVec, kyVec) - omega);

epsilond = epsilondOfk(kxVec, kyVec);

RabiRWASqOfk = @(kxVec, kyVec) (1 - v^2*(kxVec.^2 + kyVec.^2)/(4*m^2))*Rabi^2;
RabiRWASq = RabiRWASqOfk(kxVec, kyVec);

zetak = RabiRWASq./(epsilond.^2 + (1/2 + NB)^2*Gamma^2);

end



%kxVec = ;
%{
%cutoff = 0.01;
vsqr = v^2;


kmax = 1/sqrt(v^2-2*m*kappa)*sqrt(1/4*(omega + cutoff)^2 - m^2);%sqrt(m*(omegad + A0*v))./sqrt(v^2 - m*rabiCurv);
%kmax = 1/v*(m*(A0 + omegad)).^(1/2);
%kzero = (m*omegad).^(1/2)./sqrt(v^2 - m*rabiCurv);
kmin = 1/sqrt(v^2-2*m*kappa)*sqrt(1/4*(omega - cutoff)^2 - m^2); %sqrt(m*(abs(omegad - A0*v)))./sqrt(v^2 - m*rabiCurv);
%{
kmax = sqrt(m/(v^2 - 2*m*kappa)*(omega - 2*m + cutoff));%sqrt(m*(omegad + A0*v))./sqrt(v^2 - m*rabiCurv);
kmin = sqrt(m/(v^2 - 2*m*kappa)*(omega - 2*m - cutoff));
%}
%{
kmin = 0;
kmax = 6;
%}
%v^2 - 2*m*kappa

Beta = kres*(v^2 - 2*m*kappa)/m;

if imag(kmin)~=0
    kmin = 0;
end

kminmax = [kmin, kmax, kres];
%
%OmegaCW = A0*v*(4*kappa*m - v^2)*(kx.^2)/(8*m^2);
%OmegaCCW = A0*v*(1 + k.^2*(4*m*kappa-v^2)/(4*m^2));
%{
A0=.5; v=0.4; m=1; kappa=0.03;omega=3;
OmegaK = A0*v * sqrt( (1 + kappa*kx.^2/m - v^2/(4*m^2)*3*kx.^2).^2 + ...
    (v/m*kx).^2 );
OmegaKPrime = A0*v^2/(4*m^2)*kx.^2;
dk = sqrt(v^2*kx.^2 + (m-kappa*kx.^2).^2);

HtotK = sqrt((dk - omega/2).^2 + OmegaK.^2);
HtotKPrime = sqrt((dk - omega/2).^2 + 0.1*OmegaKPrime.^2);
figure; plot(kx, HtotK); hold on; plot(kx, HtotKPrime, 'r')
%}
gammavBos = NB/(1 + 2*NB);
gammacBos = (1+NB)/(1 + 2*NB);

effectiveCurv = 4*m*kappa-v^2;
[v, effectiveCurv];
effectiveCurv = -v^2;

% Things to Modify in case, kappa
rabiTildeSqBosDimlessA0 = @(k) (1 + k.^2/(4*m^2)*(v^2-4*m*kappa))*DimlessRabi^2./(epsilond(k).^2/GammaBos^2 + (1/2 + NB)^2);
rabiTildeSqMinKBosDimlessA0 = @(k) 1*(DimlessRabi./(4*m^2).*k.^2*(4*kappa*m*0 - v^2)).^2./(epsilond(k).^2/GammaBos^2 + (1/2 + NB)^2);


nvBosPlusKDimlessA0 = @(k) 1/2 + (gammacBos - gammavBos)./(2 * (2*rabiTildeSqBosDimlessA0(k) + gammavBos + gammacBos));
ncBosPlusKDimlessA0 = @(k) 1/2 - (gammacBos - gammavBos)./(2 * (2*rabiTildeSqBosDimlessA0(k) + gammavBos + gammacBos));
%nvBosMinusKDimlessA0 = @(k) 1/2 + (gammacBos - gammavBos)./(2 * (2*rabiTildeSqMinKBosDimlessA0(k) + gammavBos + gammacBos));
%ncBosMinusKDimlessA0 = @(k) 1/2 - (gammacBos - gammavBos)./(2 * (2*rabiTildeSqMinKBosDimlessA0(k) + gammavBos + gammacBos));
NdeltaBosDimlessA0 = @(k)(nvBosPlusKDimlessA0(k) + ncBosMinusKDimlessA0(k) - 1) .* mu ./ ( mu^2 + GammaBos^2*(1/2+NB).^2 );




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
GammaBos=0.0001;NB=0.001;mu=0.1;v=0.4;kappa=0.03;m=1;cutoff=0.5;omega=2.5;DimlessA0=500;
DimlessRabi = v*DimlessA0;


%{
NU=8
U = logspace(-1, 0, NU)
eigvalvsAGamma0p0001v2 = zeros(NU, 1);
parfor n = 1:NU
[eig] = FloquetExcitonNonlinear(1, 1, 2., 0.00100, U(n), 0, 1, 30, .0, 0, 1);
%[eig] = FloquetExciton(1, 1, 2., .1, U(n), .0, 1, 40, .0, 0, 1);
eigvalvsAGamma0p0001v2(n) = eig(1);
end

NA=16
A = logspace(-4, 0, NA)
eigvalvsAU1 = zeros(NA, 1);
parfor n = 1:NA
[eig] = FloquetExcitonNonlinear(1, 1, 2.0, A(n), 1, 0, 1, 30, .001, 0, 1);

eigvalvsAU1(n) = eig(1);
end

figure;plot(log10(A), log10(-eigvalvsA), '-o' ,'LineWidth', 1.5)

figure;hold on; plot(log10(A), [-eig0p01], 'b-o', 'linewidth', 1.5);plot(log10(A), [-eig0p001], 'k-d', 'linewidth', 1.5);plot(log10(A), [-eig0p0001], 'r-*', 'linewidth', 1.5);
xlabel('$log(A_0)$','FontSize',32,'FontName','Times New Roman', 'interpreter', 'latex')
ylabel('$|E_{GS}|$','FontSize',32,'FontName','Times New Roman', 'interpreter', 'latex')
ax=gca; set(gca,'LineWidth', 1.5,'FontSize',24, 'FontName','Times New Roman', 'TickLength',[0.012, 0.01]); box on;
title('$m=1; v=1; e=1; \omega=2$', 'FontName','Times New Roman', 'interpreter', 'latex', 'FontSize', 20)

NA=16
omega = linspace(1, 3, NA);
eigvalvsOmA0p01Gamma0p01 = zeros(NA, 1);
parfor n = 1:NA
[eig] = FloquetExcitonNonlinear(1, 1, omega(n), .1, 1, .0, 1, 30, .01, 0, 1);
%[eig] = FloquetExciton(1, 1, omega(n), .1, 1, .0, 1, 40, .01, 0, 1);
eigvalvsOmA0p01Gamma0p01(n) = eig(1);
end

figure;hold on; plot(omega, [-eigvalvsOmA0p01Gamma0p01], 'b-o', 'linewidth', 1.5);
%plot(log10(A), [-eig0p001], 'k-d', 'linewidth', 1.5);plot(log10(A), [-eig0p0001], 'r-*', 'linewidth', 1.5);
xlabel('$\omega$','FontSize',32,'FontName','Times New Roman', 'interpreter', 'latex')
ylabel('$|E_{GS}|$','FontSize',32,'FontName','Times New Roman', 'interpreter', 'latex')
ax=gca; set(gca,'LineWidth', 1.5,'FontSize',24, 'FontName','Times New Roman', 'TickLength',[0.012, 0.01]); box on;
title('$m=1; v=1; e=1; \omega=2$', 'FontName','Times New Roman', 'interpreter', 'latex', 'FontSize', 20)
%}


nvBosPlusKDimlessA0 = @(k) 1/2 + (gammacBos - gammavBos)./(2 * (2*rabiTildeSqBosDimlessA0(k) + gammavBos + gammacBos));
ncBosPlusKDimlessA0 = @(k) 1/2 - (gammacBos - gammavBos)./(2 * (2*rabiTildeSqBosDimlessA0(k) + gammavBos + gammacBos));
nvBosMinusKDimlessA0 = @(k) 1/2 + (gammacBos - gammavBos)./(2 * (2*rabiTildeSqMinKBosDimlessA0(k) + gammavBos + gammacBos));
ncBosMinusKDimlessA0 = @(k) 1/2 - (gammacBos - gammavBos)./(2 * (2*rabiTildeSqMinKBosDimlessA0(k) + gammavBos + gammacBos));
NdeltaBosDimlessA0 = @(k)(nvBosPlusKDimlessA0(k) + ncBosMinusKDimlessA0(k) - 1);% .* mu ./ ( mu^2 + GammaBos^2*(1/2+NB).^2 );
%figure;
set(gca,'LineWidth', 1,'FontSize',16, 'FontName','Times New Roman', 'TickLength',[0.012, 0.01])%, 'XAxisLocation', 'top'
hold on; k = 0:0.01:5; %plot(k,abs(NdeltaBosDimlessA0(k)), 'LineWidth',2, 'Color', 'k')
plot(k,(ncBosPlusKDimlessA0(k)), 'LineWidth',2, 'Color', 'r')
plot(k,(ncBosMinusKDimlessA0(k)), 'LineWidth',2, 'Color', 'b')
%}
