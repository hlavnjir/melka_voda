clear all; close all; clc;
%% skript pro pocitani vysky hladiny melke vody
% Radim Dvorak
% 18.11.2020

%____________________________PREPROCESSOR__________________________________
% vstupy
    ndx = 100; % pocet uzlu na x
    Termi = 1; % doba reseni
    Cour = 0.1; % courantovo c. ulohy (1 = mez stability pro linearni ulohu, zde je stabilita cca na 0,1)
    h_prehrady = 2; % pocatecni vyska hladiny prehrady
    h_reky = 1; % pocatecni vyska hladiny reky za prehradou
    l_prehrady = 10; % delka prehrady
    l_reky = 10; % delka reky za prehradou
    u_prehrady = 0; % horizontalni rychlost prehrady
    u_reky = 1; % horizontalni rychlost reky

% konstanty
    g = 9.81; % tihove zrychleni

%______________________________SOLVER______________________________________
% urceni casoveho kroku s ohledem na stabilitu
    dx = (l_prehrady + l_reky)./ndx;
    dt = Cour .* dx;
% pocet casovych kroku
    ndt = round(Termi./dt);
% pole reseni konzerv. promennych
    [fi, m] = deal(zeros(ndt,ndx));
% pocatecni podminky
    ndx_prehrady = round((l_prehrady)./(l_prehrady + l_reky) .* ndx);
    ndx_reky = ndx - ndx_prehrady;
    
    fi(1,1:ndx_prehrady) = g.*h_prehrady;
    fi(1,ndx_prehrady+1:end) = g.*h_reky;
    m(1,1:ndx_prehrady) = g.*h_prehrady.*u_prehrady;
    m(1,ndx_prehrady+1:end) = g.*h_reky.*u_reky;
    
%reseni
    % cykly
        for n = 1:ndt-1 % cyklus pres cas
            % numericke okrajove podminky (lze volit ruzne tvary...)
                fi(n+1,1) = fi(n,1);
                fi(n+1,end) = fi(n,end);
                m(n+1,1) = m(n,1);
                m(n+1,end) = m(n,end);

            for k = 2:ndx-1 % cyklus pres x (krajni hodnoty se pocitaji extra)
                % Lax - Friedrichsovo schema pro fi
                    fi(n+1,k) = (1./2).*(fi(n,k+1)+fi(n,k-1)) - (Cour./2).*(m(n,k+1)-m(n,k-1));
                % Lax - Friedrichsovo schema pro m
                    Fp = ((m(n,k+1).^2)./(fi(n,k+1)) + (fi(n,k+1).^2)./(2));
                    Fl = ((m(n,k-1).^2)./(fi(n,k-1)) + (fi(n,k-1).^2)./(2));
                    m(n+1,k) = (1./2).*(m(n,k+1)+m(n,k-1)) - (Cour./2).*(Fp-Fl);
            end
        end
    % prevod z konzervativnich promennych
        H = fi./g;
        U = m./H./g;
    
%_____________________________POSTPROCESSOR________________________________
% pracovni animace, prepisu podle toho, v jakym formatu bude vystup c++
    % limity os
        ylimH = [min(min(H))-0.1.*abs(min(min(H))-max(max(H))) max(max(H))+0.1.*abs(min(min(H))-max(max(H)))];
        ylimU = [min(min(U))-0.1.*abs(min(min(U))-max(max(U))) max(max(U))+0.1.*abs(min(min(U))-max(max(U)))];
    % plot
        for n = 1:ndt
            subplot(2,1,1)
            plot(1:ndx,H(n,:))
            title(['vyska hladiny h ' num2str(round(n./ndt.*100)) '%'])
            ylim(ylimH)
            subplot(2,1,2)
            plot(1:ndx,U(n,:))
            title(['horizontalni rychlost u ' num2str(round(n./ndt.*100)) '%'])
            ylim(ylimU)
            pause(0.1)
        end

        
    