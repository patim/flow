% Del^2 psi = ksi
% Del^2 ksi = 1/v*(dpsi/dy*dksi/dx - dpsi/dx*dksi/dy)
function FluidFlow2(eps)
    Lc = 1; Lb = 2; La =20; Le = 5; Lf = 9; %scaling
    v0 = 1; nu = 0.5;
    
    n = 4; % number of zones per unit L
    Nx = (Le+Lc+La)*n + 1; % number of points along x
    Ny = Lf*n + 1; % number of points along y
    dx = 1.0/n; 
    dy = dx;
    ksi = zeros(Nx, Ny);
    psi = zeros(Nx, Ny);

    % Boundary conditions
    ksi(:,Ny) = 0;
    psi(:,Ny) = v0*dy + psi(:,Ny-1);
    
    Ne = Le*n+1;
    Nc = Lc*n+1;
    Nb = Lb*n+1;
    
    x = 0:dx:(Le+Lc+La);
    y = 0:dy:Lf;
    
    % SOR Loop
    omega = 1.6;
    for l = 1:1000
        rho1 = ksi;
        rho2 = ksi_rho(nu, psi, ksi,dx,dy);
        psi = SORpsi(v0,psi,rho1,dx,dy,omega,Ne,Nc,Nb);
        
        ksi = SORksi(ksi, rho2, dx, dy, omega, Ne, Nc, Nb);
        
        % updating boundary on D
        for j=1:Nb
            ksi(Ne,j) = 2*psi(Ne-1,j)/(dx*dx);
            ksi(Ne+Nc-1,j) = 2*psi(Ne+Nc,j)/(dx*dx);
        end
        
        for i=Ne+1:Ne+Nc-1
            ksi(i,Nb) = 2*psi(i,Nb+1)/(dx*dx);
        end
        
        res_psi = residual(psi,rho1,dx,dy,Nb,Nc,Ne);
        res_ksi = residual(ksi,rho2,dx,dy,Nb,Nc,Ne);
        resnorm = norm(res_psi)/sqrt(dx*dy);
        %resnorm = norm(res_ksi)/sqrt(dx*dy);
        %fprintf('|res| = %g\n',resnorm);
        figure(1)
        cla
        contour(x,y,psi',100);
        %view(-35,45)
         figure(2)
%         cla
%         contour(x,y,res_psi',100);
        
%         figure(3)
        cla
        contour(x,y,ksi',20);
%         figure(4)
%         cla
%         contour(x,y,res_ksi',20);
        
        %pause(0.01)
        if resnorm < eps
            fprintf('Stopped at l=%d, |res| = %g\n',l,resnorm)
            break
        end
    end
    
end

%left hand side for the second equation
function rho = ksi_rho(nu,psi,ksi,dx,dy)
    [Nx,Ny] = size(psi);
    rho = zeros(Nx,Ny);
    for i=1:Nx-1
        for j=1:Ny-1
            dpsidx = (psi(i+1,j)-psi(i,j))/dx;
            dpsidy = (psi(i,j+1)-psi(i,j))/dy;
            dksidx = (ksi(i+1,j)-ksi(i,j))/dx;
            dksidy = (ksi(i,j+1)-ksi(i,j))/dy;
            rho(i,j) = (dpsidy*dksidx - dpsidx*dksidy)/nu;
        end
    end
    
    rho(Nx,:)=0;
    rho(:,Ny)=rho(:,Ny-1);
    
end

function res = residual(v,rho,dx,dy,Nb,Nc,Ne)
    vol = dx*dy;
    rx = dy/dx;
    ry = dx/dy;
    res = - 2.0*(rx+ry)*v - vol*rho;
    [nx,ny] = size(v);
    res(:,1) = 0.0; res(:,ny)=0.0;res(1,:)=0.0;res(nx,:)=0.0;
    
    for i=2:nx-1
        for j=2:ny-1
            res(i,j) = res(i,j) + rx*(v(i+1,j)+v(i-1,j)) + ry*(v(i,j+1)+v(i,j-1));
        end
    end
    
    for i=1:Nb
        res(Ne,i) = 0.0;
        res(Ne+Nc-1,i) = 0.0;
    end
    
    for i=Ne:Ne+Nc-1
        res(i,Nb) = 0;
    end
end

function v = SORksi(v,rho,dx,dy,omega,Ne,Nc,Nb)
    vol = dx*dy;
    rx = dy/dx;
    ry = dx/dy;
    [Nx,Ny] = size(v);
    
    % region (i)
    for i=2:Ne-1
        for j=2:Ny-1
            delta = rx*(v(i+1,j)-2.0*v(i,j)+v(i-1,j)) + ry*(v(i,j+1)-2.0*v(i,j)+v(i,j-1)) - vol*rho(i,j);
            v(i,j) = v(i,j) + omega*delta/(2.0*(rx+ry));
        end
    end
    
    % region (ii)
    for i=Ne:Ne+Nc-1
        for j=Nb+1:Ny-1
            delta = rx*(v(i+1,j)-2.0*v(i,j)+v(i-1,j)) + ry*(v(i,j+1)-2.0*v(i,j)+v(i,j-1)) - vol*rho(i,j);
            v(i,j) = v(i,j) + omega*delta/(2.0*(rx+ry));
        end
    end    
    
    % region (iii)
    for i=Ne+Nc:Nx-1
        for j=2:Ny-1
            delta = rx*(v(i+1,j)-2.0*v(i,j)+v(i-1,j)) + ry*(v(i,j+1)-2.0*v(i,j)+v(i,j-1)) - vol*rho(i,j);
            v(i,j) = v(i,j) + omega*delta/(2.0*(rx+ry));
        end
    end 
    
    for i=2:Ny-1
        v(Nx,i)=v(Nx-1,i);
    end
end

function v = SORpsi(v0,v,rho,dx,dy,omega,Ne,Nc,Nb)
    vol = dx*dy;
    rx = dy/dx;
    ry = dx/dy;
    [Nx,Ny] = size(v);
    
    % region (i)
    for i=2:Ne-1
        for j=2:Ny-1
            delta = rx*(v(i+1,j)-2.0*v(i,j)+v(i-1,j)) + ry*(v(i,j+1)-2.0*v(i,j)+v(i,j-1)) - vol*rho(i,j);
            v(i,j) = v(i,j) + omega*delta/(2.0*(rx+ry));
        end
    end
    
    % region (ii)
    for i=Ne:Ne+Nc-1
        for j=Nb+1:Ny-1
            delta = rx*(v(i+1,j)-2.0*v(i,j)+v(i-1,j)) + ry*(v(i,j+1)-2.0*v(i,j)+v(i,j-1)) - vol*rho(i,j);
            v(i,j) = v(i,j) + omega*delta/(2.0*(rx+ry));
        end
    end    
    
    % region (iii)
    for i=Ne+Nc:Nx-1
        for j=2:Ny-1
            delta = rx*(v(i+1,j)-2.0*v(i,j)+v(i-1,j)) + ry*(v(i,j+1)-2.0*v(i,j)+v(i,j-1)) - vol*rho(i,j);
            v(i,j) = v(i,j) + omega*delta/(2.0*(rx+ry));
        end
    end 
    
    % updating boundaries
    for i=2:Ny-1
        v(1,i) = v(2,i);
        v(Nx,i) = v(Nx-1,i);
    end
    
    for i=2:Nx-1
        v(i,Ny) = v(i,Ny-1) + v0*dy;
    end
    
    % updating upper boundary
end