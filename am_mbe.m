classdef am_mbe

    properties (Constant)
        tiny      = 1E-4;          % precision of atomic coordinates
        r_0       = 2.81794032E-6; % [nm]       classical electron radius
        N_A       = 6.022141E23;   % [mol]      Avogadro's number
    end
    
    % pressure related things

    methods (Static)
        
        function [P]     = flux2pressure(J,T,m,th)
            %
            % P = flux2pressure(J,T,m,th) 
            % 
            % J     [atoms/cm^2-s]     flux 
            % P     [Torr]             pressure
            % T     [K]                temperature
            % m     [amu]              mass of atomic species
            % th    [deg]              molecular beam incident angle
            % 
            % Conversion factor 
            %   = 1 / ( cm^2 * s ) * sqrt( K * amu *  boltzman's constant) / Torr
            %   = 1.135674×10^-23 
            %
            % Eq 1 from C. D. Theis, J. Yeh, D. G. Schlom, M. E. Hawley,
            % and G. W. Brown, Thin Solid Films 325, 107 (1998). 
            %
            
            if nargin < 4
                th = 1; % assume normal flux if theta is not supplied
            end
            
            P = J / cosd(th) * sqrt(pi*m*T/8) * 1.135674E-23;
            
        end
        
        function [J]     = pressure2flux(P,T,m,th)
            %
            % J = pressure2flux(P,T,m,th) 
            % 
            % P     [Torr]             pressure
            % J     [atoms/cm^2-s]     flux 
            % T     [K]                temperature
            % m     [amu]              mass of atomic species
            % th    [deg]              molecular beam incident angle
            % 
            % Conversion factor 
            %   = 1 / ( cm^2 * s ) * sqrt( K * amu *  boltzman's constant) / Torr
            %   = 1.135674×10^-23 
            %
            % Eq 1 from C. D. Theis, J. Yeh, D. G. Schlom, M. E. Hawley,
            % and G. W. Brown, Thin Solid Films 325, 107 (1998). 
            %
            
            if nargin < 4
                th = 1; % assume normal flux if theta is not supplied
            end
            
            J = P * cosd(th) / sqrt(pi*m*T/8) / 1.135674E-23;
            
        end
        
        function [mfp]   = pressure2mfp(P,T,sigma)
            %
            % mfp = pressure2mfp(P,T,sigma)
            % 
            % P     [Torr]             pressure
            % T     [K]                temperature
            % sigma [cm^2]             scattering cross-section
            % mfp   [cm]               mean free path
            % 
            % Conversion factor
            %	= Boltzmann's constant * K / ( cm^2 * Torr ) / cm
            %   = 1.03557E-19 
            %
            % Eq 10.1 in A. Rockett, The Materials Science of
            % Semiconductors (Springer Science & Business Media, 2007). 
            %
            
            if nargin < 3 
                sigma = 1E-16; % assume a typical values if sigma is not supplied
            end
            
            mfp = T / ( sqrt(2) * sigma * P ) * 1.03557E-19;
            
        end
        
    end
    
    % x-ray related things
    
    methods (Static)
        
        function           demo_simulate_xrr()
            %
            % R = simulate_xray_reflectivity(layer,th,hv,thickness,filling,roughness)
            % 
            % layers [struct]           layers
            %    layer(i).mass_density      [g/cm^3]        mass density
            %    layer(i).Z                 [Z]             atomic species
            %    layer(i).stoichiometry     [atoms/u.c.]    number of atoms per formula unit
            % th        [deg]           angles
            % thickness [nm]            thickness
            % filling   [unitless]      multiplies density
            % roughness [nm]            inteface roughness
            % 
            
            import am_mbe.*
            
            % define photon energy and angles
            hv = get_atomic_emission_line_energy(get_atomic_number('Cu'),'kalpha1');
            th = [0:0.005:4].';

            % define material (species, stoichiometry [per f.u.], density)
%             layer(1) = define_material({'Ge'},[1],5.323);
%             layer(1) = define_material({'W'},[1],19.3);
            layer(1) = define_material({'Ta','O'},[2,5],8.2);
            layer(2) = define_material({'Si'},[1],2.53);

            thickness = [50 1E3];
            roughness = [0 0];
            filling   = [1 1];

            % get dielectric contributions
            nlayers = numel(layer);
            for i = 1:nlayers
                layer(i).xray_refractive_index = get_xray_refractive_index( ...
                    layer(i).Z,layer(i).atomic_density.*layer(i).stoichiometry/sum(layer(i).stoichiometry),hv,th);
            end

            % simulate and analyze
            % x = [thickness,filling,roughness,I0,background]
            nlayers = numel(layer);
            simulate_ = @(x,method) simulate_xrr(layer,th,hv,...
                x(0*nlayers+[1:nlayers]),...
                x(1*nlayers+[1:nlayers]),...
                x(2*nlayers+[1:nlayers]),...
                x(3*nlayers+1),method)*x(end-1)+x(end);
            R_transfer = simulate_([thickness,filling,roughness,10,2,1E-3],1);
            R_parratt  = simulate_([thickness,filling,roughness,10,2,1E-3],2);

            % exclude above noise floor
            ex_ = th<1.8; ex_ = th>0;
            analyze_xrr_with_fft(th(ex_),R_parratt(ex_),hv)
            subplot(3,1,1); hold on; plot(th(ex_),R_transfer(ex_),'.r');
        end
 
        function           demo_xrd_fit()
            %
            clear; clc; import am_mbe.*
            hv = get_atomic_emission_line_energy(get_atomic_number('Cu'),'kalpha1');
            d  = xrdmlread('~/Downloads/ZW628_70minPbZr0p3Ti0p7O3_onSrTiO3_tth.xrdml');
            % semilogy(d.Theta2,d.data)
            domain = [18,25];
            [d,th2] = analyze_xrd_with_fit(d.Theta2,d.data,hv,domain);
        end
        
        function           demo_xrr_fit()
            
            clear; clc; import am_mbe.*
           
            % define photon energy and angles
            hv = get_atomic_emission_line_energy(get_atomic_number('Cu'),'kalpha1');
            
            % load data
            fname = './test_data/0.8W_on_Si_refs.txt';
%             fname = './test_data/0.8Ta2O5_on_Si_refs.txt';
            fid = fopen(fname,'r'); x = fscanf(fid,'%f %f\n'); fclose(fid);
            th = x(1:2:end); intensity = x(2:2:end);
           
            % define material (species, stoichiometry, density)
            layer(1) = define_material({'W'},[1],19.3);
%             layer(1) = define_material({'Ta','O'},[2,5],8.2);
            layer(2) = define_material({'Si'},[1],2.53);
            
            isfixed = [  0   1   0   0   0   0  0   0    0];
            lb      = [  0 1E8 0.6 0.6 0.1 0.1  7   1 1E-8];
            ub      = [100 1E8 1.4 1.4   5   5 13 1E8 1E+2];
            x0      = mean([lb;ub]);
            
            x = analyze_xrr_with_fit(th,intensity,hv,layer,x0,ub,lb,isfixed);
        end
        
        function           analyze_xrr_with_fft(th,intensity,hv,thc)
            %
            % beta = analyze_xrr_with_fft(th,xrr,hv,thc)
            % 
            % th         [deg]          incident x-ray angle
            % xintensity [a.u.]         xrr intensity
            % hv         [eV]           photon energy
            % thc        [deg]          (optional) critical angle
            %
            % Example:
            % ex_ = th<1.8;
            % analyze_xrr_with_fft(th(ex_),intensity(ex_),hv,0.36)
            % analyze_xrr_with_fft(th(ex_),intensity(ex_),hv)
            %

            import am_mbe.* am_lib.*

            if nargin < 4 || isempty(thc)
                qc = [];
            else
                qc = get_qz(thc,hv);
            end

            % fft analysis
            F = intensity(:);
            % get wavevector
            qz = get_qz(th,hv);
            % subtract critical wavevector
            if isempty(qc); [~,c]=max(-diff(F)./diff(th)); qc = qz(c); end
            q = sqrt(qz.^2 - qc.^2); Fq=F(q>0); q=q(q>0); nqs = numel(q);
            % resample Fq,q on equidistant intervals
            qe = linspace(min(q),max(q),nqs).'; Fq = interp1(q,Fq,qe); q = qe; 
            % apply background correction
            Fq= (Fq.*q.^4-mean(Fq.*q.^4)) .* tukey_(nqs,0.90);
            % evalute DFT
            x  = [0:0.5:200]; 
            K  = exp(-1i*x(:).*q(:).');
            Fx = K*Fq(:);
            % evaluate FFT
            nxs=2^10;
            x_fft = [0:nxs-1]/(nxs)/(q(2)-q(1))*pi*2; x_fft = x_fft([1:floor(nxs/2)]);
            F_fft = fft(Fq,2^10);                     F_fft = F_fft([1:floor(nxs/2)]);
            % plot th, Fx, and Fq
            subplot(3,1,1); semilogy(th, F); xlabel('\theta [deg]'); ylabel('Intensity [a.u.]'); axis tight;
            subplot(3,1,2); plot(q, Fq); xlabel('(q^2 - q_c^2)^{1/2} [1/nm]'); ylabel('Intensity [a.u.]'); axis tight;
            subplot(3,1,3); semilogy(x, abs(Fx).^2,'-', x_fft, abs(F_fft).^2,'.'); xlabel('x [nm]'); ylabel('|FFT|^2 [a.u.]'); axis tight; xlim([0 200])

        end

        function [x]     = analyze_xrr_with_fit(th,intensity,hv,layer,x0,ub,lb,isfixed,domain)
            import am_mbe.*
            import am_lib.*
            
            % if window is defined, crop and filter
            if nargin > 8
                ex_ = and( th < max(domain) , th > min(domain) );
                th = th(ex_); intensity = intensity(ex_);
            end
            
            % get dielectric contributions
            nlayers = numel(layer);
            for i = 1:nlayers
                layer(i).xray_refractive_index = get_xray_refractive_index( ...
                    layer(i).Z,layer(i).atomic_density.*layer(i).stoichiometry/sum(layer(i).stoichiometry),hv,th);
            end
            % x = [thickness,filling,roughness,I0,background]
            nlayers = numel(layer);
            simulate_ = @(x) simulate_xrr(layer,th,hv,...
                x(0*nlayers+[1:nlayers]),...
                x(1*nlayers+[1:nlayers]),...
                x(2*nlayers+[1:nlayers]),...
                x(3*nlayers+1))*x(end-1)+x(end);
            % define rescaling
            fscale_= @(x)   [x(1:(3*nlayers+1)),log(x((3*nlayers+2):end))];
            rscale_= @(x)   [x(1:(3*nlayers+1)),exp(x((3*nlayers+2):end))];
            % scale lower and upper bounds
            lb = fscale_(lb); ub = fscale_(ub);
            % define objective function: robust + pearson's correlation coefficient
            cost_  = @(x) sum(abs( log(simulate_(rscale_(x))).' - log(intensity(:)) )) ;

            switch 'GA'
                case 'GA'
                    % optimization options
                    opts_hybrid_ = optimoptions(...
                        @fmincon,'Display','off',...
                                 'MaxFunctionEvaluations',1E3,...
                                 'MaxIterations',1E3,...
                                 'StepTolerance',1E-18,...
                                 'FunctionTolerance',1E-18,...
                                 'Algorithm','active-set');
                    opts_ = optimoptions(@ga,'PopulationSize',200, ...
                                             'InitialPopulationMatrix',x0(~isfixed), ...
                                             'EliteCount',1, ...
                                             'FunctionTolerance',1, ...
                                             'Display','off',...
                                             'PlotFcns',{@plot_func_ga_},...
                                             'HybridFcn',{@fmincon,opts_hybrid_});
                    % perform optimization
                    [x,~] = ga_(cost_,x0,isfixed,[],[],[],[],lb(~isfixed),ub(~isfixed),[],opts_); x = rscale_(x);
                case 'SA'
                    % optimization options
                    opts_hybrid_ = optimoptions(@fmincon,...
                                 'Display','off',...
                                 'MaxFunctionEvaluations',1E3,...
                                 'MaxIterations',1E3,...
                                 'StepTolerance',1E-18,...
                                 'FunctionTolerance',1E-18,...
                                 'Algorithm','active-set');
                    opts_ = optimoptions(@simulannealbnd,...
                                 'Display','off',...
                                 'MaxIterations',100,...
                                 'InitialTemperature',1.5*abs(ub(~isfixed)-lb(~isfixed)),...
                                 'ReannealInterval',10,...
                                 'FunctionTolerance',1E-18,...
                                 'PlotFcns',{@plot_func_sa_},...
                                 'HybridFcn',{@fmincon,opts_hybrid_});
                    % perform optimization
                    [x,~] = sa_(cost_,x0,isfixed,lb(~isfixed),ub(~isfixed),opts_);  x = rscale_(x);
            end
            
            function state = plot_func_ga_(~,state,flag,~)
                switch flag
                    % Plot initialization
                    case 'init'
                        clf; figure(1); set(gcf,'color','w');
                    case 'iter'
                        % find best population
                        [~,j]=min(state.Score);
                        % reconstruct full vector
                        f = find( isfixed); xf = x0(f);
                        l = find(~isfixed); xl = state.Population(j,:);
                        x([f,l]) = [xf,xl];
                        % plot it
                        semilogy(th,simulate_(rscale_(x)),th,intensity); 
                        % title
                        title(sprintf('(generation = %i)',state.Generation)); 
                        % state.Population(i,3), 2*asind(1239.842/hv/(4*pi)*state.Population(i,2))
                        axis tight; xlabel('\theta [deg]'); ylabel('Intensity [a.u.]');
                        ylim([min(intensity),max(intensity)]);
                    case 'done'
                        % initialize
                        figure(1); clf; set(gcf,'color','w');
                        % plot
                        semilogy(th,intensity,'.',th,simulate_(rscale_(x)),'--');
                        axis tight; xlabel('\theta [deg]'); ylabel('Intensity [a.u.]');
                        % title
                        title(sprintf('thickness = %.2f nm; film density %.2f g/cm3; substrate density %.2f g/cm3; \n film roughness = %.2f nm, substrate roughness = %.2f nm',...
                            x(1),x(3)*layer(1).mass_density,x(4)*layer(2).mass_density,x(5),x(6)));
                        % define labels manually
                        switch nlayers
                            case 2; label = {'T1','T2','D1','D2','R1','R2','W','I','B'};
                            case 3; label = {'T1','T2','T3','D1','D2','D3','R1','R2','R3','W','I','B'};
                        end
                        % bounds
                        axes('position',[0.4 0.6 0.45 0.25]);
                        a = 1:sum(~isfixed); b = (x(~isfixed)-lb(~isfixed))./(ub(~isfixed)-lb(~isfixed)); plot(a,b,'.-','markersize',20); 
                        h=text(a,b+0.05,strread(sprintf('%0.3g\n',x(~isfixed)),'%s'),'HorizontalAlignment','left'); set(h,'rotation',90);
                        ylim([0 1]); set(gca,'YTick',[0 1],'YTickLabel',{'LB','UB'}); set(gca,'XTickLabel',label(~isfixed)); 
                end
            end
            function state = plot_func_sa_(~,state,flag,~)
                switch flag
                    % Plot initialization
                    case 'init'
                        clf; figure(1); set(gcf,'color','w');
                    case 'iter'
                        % reconstruct full vector
                        f = find( isfixed); xf = x0(f);
                        l = find(~isfixed); xl = state.bestx;
                        x([f,l]) = [xf,xl];
                        % plot it
                        semilogy(th,simulate_(rscale_(x)),th,intensity); 
                        % title
                        title(sprintf('(iteration = %i)',state.iteration)); 
                        axis tight; xlabel('\theta [deg]'); ylabel('Intensity [a.u.]');
                        ylim([min(intensity),max(intensity)]);
                    case 'done'
                        % initialize
                        figure(1); clf; set(gcf,'color','w');
                        % plot
                        semilogy(th,intensity,'.',th,simulate_(rscale_(x)),'--');
                        axis tight; xlabel('\theta [deg]'); ylabel('Intensity [a.u.]');
                        % title
                        title(sprintf('thickness = %.2f nm; film density %.2f g/cm3; substrate density %.2f g/cm3; \n film roughness = %.2f nm, substrate roughness = %.2f nm',...
                            x(1),x(3)*layer(1).mass_density,x(4)*layer(2).mass_density,x(5),x(6)));
                        % define labels manually
                        switch nlayers
                            case 2; label = {'T1','T2','D1','D2','R1','R2','W','I','B'};
                            case 3; label = {'T1','T2','T3','D1','D2','D3','R1','R2','R3','W','I','B'};
                        end
                        % bounds
                        axes('position',[0.4 0.6 0.45 0.25]);
                        a = 1:sum(~isfixed); b = (x(~isfixed)-lb(~isfixed))./(ub(~isfixed)-lb(~isfixed)); plot(a,b,'.-','markersize',20); 
                        h=text(a,b+0.05,strread(sprintf('%0.3g\n',x(~isfixed)),'%s'),'HorizontalAlignment','left'); set(h,'rotation',90);
                        ylim([0 1]); set(gca,'YTick',[0 1],'YTickLabel',{'LB','UB'}); set(gca,'XTickLabel',label(~isfixed)); 
                end
                state = false;
            end
        end
        
        function           analyze_xrd_with_fft(th2,intensity,hv,domain)
            %
            % beta = analyze_xrd_with_fft(th,xrr,hv,domain)
            % 
            % th         [deg]          incident x-ray angle
            % xintensity [a.u.]         xrr intensity
            % hv         [eV]           photon energy
            % domain     [deg]          (optional) [min max]
            %
            % Example:
            % ex_ = th<1.8;
            % analyze_xrr_with_fft(th(ex_),intensity(ex_),hv,0.36)
            % analyze_xrr_with_fft(th(ex_),intensity(ex_),hv)
            %

            import am_mbe.*
            
            % resample intensity and q on equidistant intervals
            q = get_qz(th2(:)/2,hv); qe = linspace(min(q),max(q),numel(q)).'; Fq = interp1(q,intensity,qe); q = qe; th2 = get_th(q,hv)*2;
            % if window is defined, crop and filter
            if nargin > 3; if ~isempty(domain)
                ex_ = and( th2 < max(domain) , th2 > min(domain) );
                th2 = th2(ex_); q = q(ex_); Fq = Fq(ex_); Fq = Fq .* tukey_(numel(Fq),0.90); 
            end; end
            % evalute DFT
            x  = [0:0.5:200]; 
            K  = exp(-1i*x(:).*q(:).');
            Fx = K*Fq(:);
            % evaluate FFT
            nxs=2^10;
            x_fft = [0:nxs-1]/(nxs)/(q(2)-q(1))*pi*2; x_fft = x_fft([1:floor(nxs/2)]);
            F_fft = fft(Fq,2^10);                     F_fft = F_fft([1:floor(nxs/2)]);
            % plot th, Fx, and Fq
            subplot(2,1,1); semilogy(th2,Fq,'-'); xlabel('2\theta [deg]'); ylabel('Intensity [a.u.]'); axis tight;
            subplot(2,1,2); semilogy(x, abs(Fx).^2,'-', x_fft, abs(F_fft).^2,'.'); xlabel('x [nm]'); ylabel('|FFT|^2 [a.u.]'); axis tight; xlim([0 200]);
        end
        
        function [d,th2] = analyze_xrd_with_fit(th2,intensity,hv,domain,profile)
            %
            % beta = analyze_xrd_with_fit(th,xrr,hv)
            % 
            % th         [deg]          incident x-ray angle
            % xintensity [a.u.]         xrr intensity
            % hv         [eV]           photon energy
            % domain     [deg]          (optional) [min max]
            %

            import am_mbe.* am_lib.*
            
            % resample intensity and q on equidistant intervals
            q = get_qz(th2(:)/2,hv); qe = linspace(min(q),max(q),numel(q)).'; 
            intensity = interp1(q,intensity,qe); intensity(intensity==0)=min(intensity(intensity~=0)); %intensity = intensity./max(intensity);
            q = qe; th2 = get_th(q,hv)*2;
            % if window is defined, crop and filter
            if nargin > 3
                ex_ = and( th2 < max(domain) , th2 > min(domain) );
                th2 = th2(ex_); q = q(ex_); intensity = intensity(ex_);
            else 
                ex_ = true;
            end
            if nargin < 4
                profile='f';
            end
            if      contains(profile,'fs')
                % define rescaling
                fscale_= @(x) [log(x(1)),x(2:3),log(x(4)),x(5:6),log(x(7))];
                rscale_= @(x) [exp(x(1)),x(2:3),exp(x(4)),x(5:6),exp(x(7))];
                % define gaussian, sinc, peak, and objective functions
                gaus_ = @(x) x(1).*exp(-((q-x(2))./x(3)).^2);
                sinc_ = @(x) x(1).*sinc_( (q-x(2)).*x(3)).^2 ;
                func_ = @(x) sinc_(x(1:3))+gaus_(x(4:6)) + x(7);
                cost_ = @(x) sum(abs(log(func_(rscale_(x))) - log(intensity(:)) ));
                % guess substrate position
                [~,i]=max(intensity); q_sub = q(i);
                % lower and upper bounds; starting condition
                isfixed = [0 0 0 0 0 0 0];
                lb = [1E-6 min(q)  0 1E-1 q_sub-0.1 0.1 1E-6]; lb = fscale_(lb); 
                ub = [1E+6 max(q) 10 1E+9 q_sub+0.1 1.0 1E-4]; ub = fscale_(ub);
                x0 = mean([lb;ub]);
            elseif contains(profile,'f')
                % define rescaling
                fscale_= @(x) [log10(x(1)),x(2:3),log10(x(4))];
                rscale_= @(x) [ 10.^(x(1)),x(2:3), 10.^(x(4))];
                % define gaussian, sinc, peak, and objective functions
                sinc_ = @(x) x(1).*sinc( (q-x(2)).*x(3)).^2;
                func_ = @(x) sinc_(x(1:3)) + x(4);
                cost_ = @(x) sum(abs(log(func_(rscale_(x))) - log(intensity(:)) ));
                % lower and upper bounds; starting condition
                isfixed = [0 0 0 0];
                lb = [1E-6 min(q)  5/2/pi 1E-6]; lb = fscale_(lb); 
                ub = [1E6  max(q) 50/2/pi 1E-4]; ub = fscale_(ub);
                x0 = mean([lb;ub]);
            elseif contains(profile,'pv')
                % define rescaling
                fscale_= @(x) [log10(x(1)),x(2:4),log10(x(5))];
                rscale_= @(x) [ 10.^(x(1)),x(2:4), 10.^(x(5))];
                % define pseudo-voigt
                func_ = @(x) x(1).*pvoigt_((q-x(2))./x(3),x(4)).*x(3) + x(5);
                cost_ = @(x) sum(abs(log10(func_(rscale_(x))) - log10(intensity(:)) ));
                % lower and upper bounds; starting condition
                isfixed = [0 0 0 0 0];
                lb = [    min(intensity) min(q) 0.1 0 min(intensity)]; lb = fscale_(lb); 
                ub = [110*max(intensity) max(q) 1.0 1 max(intensity)]; ub = fscale_(ub);
                x0 = mean([lb;ub]);
            end
            % optimization options
            opts_hybrid_ = optimoptions(...
                @fmincon,'Display','off',...
                         'MaxFunctionEvaluations',1E3,...
                         'MaxIterations',1E3,...
                         'StepTolerance',1E-18,...
                         'FunctionTolerance',1E-18,...
                         'Algorithm','active-set');
            opts_ = optimoptions(@ga,'PopulationSize',100, ...
                                     'InitialPopulationMatrix',x0(~isfixed), ...
                                     'MutationFcn',{@mutationadaptfeasible}, ...
                                     'EliteCount',1, ...
                                     'FunctionTolerance',1E-5, ...
                                     'Display','off',...
                                     'PlotFcns',{@plot_func_},...
                                     'HybridFcn',{@fmincon,opts_hybrid_});
            % perform optimization
            [x,~] = ga_(cost_,x0,isfixed,[],[],[],[],lb(~isfixed),ub(~isfixed),[],opts_); x = rscale_(x);
            
            figure(1); clf; set(gcf,'color','w');
            semilogy(th2,intensity,'.',th2,func_(x),'--');
            axis tight; xlabel('2\theta [deg]'); ylabel('Intensity [a.u.]');
            
            % output
            d = 2*pi*x(3); % thickness [nm]
            th2 = 2*get_th(x(2),hv); % 2theta peak position [deg]
            
            function state = plot_func_(~,state,flag,~)
                
                switch flag
                    % Plot initialization
                    case 'init'
                        clf; figure(1); set(gcf,'color','w');
                    case 'iter'
                        % find best population
                        [~,j]=min(state.Score);
                        % reconstruct full vector
                        f = find( isfixed); xf = x0(f);
                        l = find(~isfixed); xl = state.Population(j,:);
                        x([f,l]) = [xf,xl];
                        % plot it
                        semilogy(th2,func_(rscale_(x)),th2,intensity); 
                        % title
                        title(sprintf('(generation = %i)',state.Generation)); 
                        % state.Population(i,3), 2*asind(1239.842/hv/(4*pi)*state.Population(i,2))
                        title(sprintf('thickness = %2.2f nm   2th = %2.3f deg   (generation = %i)',...
                            2*pi*state.Population(j,3), 2*asind(1239.842/hv/(4*pi)*state.Population(j,2)),state.Generation )); 
                        axis tight; xlabel('2\theta [deg]'); ylabel('intensity [a.u.]');
                        %ylim([min(intensity),max(intensity)]);
                end
            end
        end
        
        function [intensity] = simulate_xrr(layer,th,hv,thickness,filling,roughness,sample_length,method)
            %
            % R = simulate_xray_reflectivity(layer,th,hv,thickness,filling,roughness)
            % 
            % layers [struct]           layers
            %    layer(i).mass_density      [g/cm^3]        mass density
            %    layer(i).Z                 [Z]             atomic species
            %    layer(i).stoichiometry     [atoms/u.c.]    number of atoms per formula unit
            % th        [deg]           angles
            % thickness [nm]            thickness
            % filling   [unitless]      multiplies density
            % roughness [nm]            inteface roughness
            % method    [1,2]           transfer matrix (explicit,slow), recursive parratt (default,fast)
            
            import am_mbe.*
            
            if nargin < 8; method = 2; end
            if nargin < 7; sample_length = []; end

            switch method
                case 1 % explicit transfer matrix (slower)                    
                    % get kx and kz
                    % Note: k are incident wavevectors and not diffraction vectors
                    lambda = get_photon_wavelength(hv);
                    get_kz = @(th,hv) 2*pi/lambda*sind(th);
                    get_kx = @(th,hv) 2*pi/lambda*cosd(th);

                    % [nths,nlayers] solve wavevector boundary conditions to get
                    % out-of-plane wavevector component in layer kz [nths,nlayers]
                    % 1. k2  = (n2/n1) k1  is  snell's law
                    % 2. k1x = k2x        
                    nlayers = numel(layer); nths = numel(th); kz = zeros(nths,nlayers);
                    get_kz_in_layer = @(n1,n2,k1z,k1x) sqrt( (n2./n1).^2 .* k1z.^2 + ((n2./n1).^2-1) .* k1x.^2 );
                    for i = 1:nlayers
                        kz(:,i) = get_kz_in_layer(1,filling(i)*(layer(i).xray_refractive_index(:)-1)+1,get_kz(th,hv),get_kx(th,hv));
                    end

                    % get transfer matricies R and T which describe propagation
                    % across interfaces and media
                    % syms a b c d k1 k2 z
                    % % boundary conditions for plane waves at an interface
                    % eq(1) = a * exp(1i*k1*z) + b * exp(-1i*k1*z) == c*exp(1i*k2*z) + d*exp(-1i*k2*z);
                    % eq(2) = diff(eq(1),z);
                    % solution = solve(eq,a,b);
                    % % transfer matrix
                    % T(1,1:2) = simplify(equationsToMatrix(solution.a,c,d));
                    % T(2,1:2) = simplify(equationsToMatrix(solution.b,c,d));
                    % % set interface at 0 for simplicity
                    % T = subs(T,z,0)
                    % simplify
                    % T = simplify(expand(T));
                    get_transfer_matrix_at_interface = @(k1,k2,sigma) [ ...
                        [ (exp(-(k1 - k2).^2*sigma.^2/2) * (k1 + k2))/(2*k1), (exp(-(k1 + k2).^2*sigma.^2/2) * (k1 - k2))/(2*k1)]
                        [ (exp(-(k1 + k2).^2*sigma.^2/2) * (k1 - k2))/(2*k1), (exp(-(k1 - k2).^2*sigma.^2/2) * (k1 + k2))/(2*k1)] ];

                    get_transfer_matrix_in_medium = @(k,l) [ ...
                            [ exp(-1i*k*l), 0]
                            [ 0, exp(+1i*k*l)]];

                    % get reflection at each theta value
                    intensity = zeros(1,nths);
                    for j = 1:nths
                        % vacuum/first layer
                        M =     get_transfer_matrix_at_interface(get_kz(th(j),hv)    , kz(j,1), roughness(1) );
                        M = M * get_transfer_matrix_in_medium(                         kz(j,1), thickness(1) );
                        % first layer ... nth layer
                        if nlayers > 1
                            for i = 2:nlayers
                                M = M * get_transfer_matrix_at_interface( kz(j,i-1)  , kz(j,i), roughness(i) );
                                M = M * get_transfer_matrix_in_medium(                 kz(j,i), thickness(i) );
                            end
                        end
                        intensity(j) = abs(M(2,1)./M(1,1)).^2;
                    end
                    
                    th=th(:).';
                    intensity=intensity(:).';
                    
                case 2 % recursive Parratt method without transfer matrix (faster)
                    % number of layers
                    nlayers = numel(layer);
                    lambda = get_photon_wavelength(hv);
                    th2kz = @(th,lambda) 4*pi/lambda .* sind(th(:));
                    %----- Wavevector transfer in each layer
                    nths = numel(th); Q = zeros(nlayers,nths);
                    Q(1,:) = th2kz(th,lambda);
                    for j=1:nlayers
                        Q(j+1,:)= sqrt(Q(1,:).^2 + 2*(4*pi/lambda).^2 * (layer(j).xray_refractive_index(:).'-1) * filling(j) );
                    end
                    %----- Reflection coefficients (no multiple scattering)
                    r = zeros(nlayers,nths);
                    for j=1:nlayers
                        r(j,:)=(  (Q(j,:)-Q(j+1,:))./(Q(j,:)+Q(j+1,:)) ) .* exp(-0.5*(Q(j,:).*Q(j+1,:))*roughness(j)^2);
                    end
                    %----- Reflectivity
                    R = zeros(nlayers-1,nths);
                    if nlayers>1
                        R(1,:) =  (r(nlayers-1,:)  + r(nlayers,:) .* exp(1i*Q(nlayers,:)*thickness(nlayers-1)) ) ...
                              ./(1+r(nlayers-1,:) .* r(nlayers,:) .* exp(1i*Q(nlayers,:)*thickness(nlayers-1)) );
                    end
                    if nlayers>2
                    for j=2:nlayers-1
                        R(j,:) =  (r(nlayers-j,:)  + R(j-1,:) .* exp(1i*Q(nlayers-j+1,:)*thickness(nlayers-j)) ) ...
                              ./(1+r(nlayers-j,:) .* R(j-1,:) .* exp(1i*Q(nlayers-j+1,:)*thickness(nlayers-j)) );
                    end
                    end
                    %------ Intensity reflectivity
                    if nlayers==1
                        intensity = abs(r(1,:)).^2;
                    else
                        intensity = abs(R(nlayers-1,:)).^2;
                    end
                    
                    th=th(:).';
                    intensity=intensity(:).';
                    
            end
            
            % add intensity correction due to finite sample width
            if ~isempty(sample_length)
                xray_beam_height = 0.1; % [mm] height of x-ray beam
                th_b = asind(xray_beam_height/sample_length);
                intensity(th<th_b) = sind(th(th<th_b)) ./ sind(th_b) .* intensity(th<th_b);
            end
        end

        function [h]     = plot_bragg_reflections(uc)
            import am_mbe.*
            import am_lib.*
            % generate hkl list
            hv = get_atomic_emission_line_energy(get_atomic_number('Cu'),'kalpha1');
            hkl = permn_([0:4],3).'; hkl=hkl(:,2:end); k = hkl; k_mag = normc_(inv(uc.bas).'*k);
            [~,i] = unique(rnd_(k_mag)); k=k(:,i); k_mag=k_mag(:,i); hkl=hkl(:,i);
            th2 = 2*asind( get_photon_energy(hv) * k_mag / 2 * 10 ); ex_ = abs(imag(th2(:)))<1E-8; 
            k=k(:,ex_); hkl=hkl(:,ex_); th2=th2(ex_);
            % get atomic scattering factors and structure factor
            [Z,~,j] = unique(get_atomic_number({uc.symb{uc.species}}));
            f0 = permute(get_atomic_xray_form_factor(Z,hv,th2),[1,3,2]);
            Fhkl = sum(f0(j,:).*exp(2i*pi*uc.tau.'*k),1); 
            Fhkl2= abs(Fhkl).^2; ex_ = abs(Fhkl)>1E-5;% Fhkl2 = Fhkl2./max(Fhkl2(:));
            % plot Bragg peaks
            h = hggroup;
            plot(th2(ex_),Fhkl2(ex_),'.','markersize',10,'Parent',h); hold on;
            for i = 1:numel(th2); if ex_(i)
                line([th2(i),th2(i)],[0,Fhkl2(i)],'Parent',h);
                text(th2(i),Fhkl2(i),sprintf('  %i%i%i',hkl(:,i)),'Rotation',90,'Parent',h);
            end; end
            hold off;
        end
        
        function [thc]   = get_xray_critical_angle(delta)
            %
            % thc = photoabsorbtion_crosssection(Z,lambda)
            % 
            % thc    [deg]          critical angle
            % delta  [unitless]     real contribution to dielectric function
            %
            % Eq 3.3 in J. Als-Nielsen and Des McMorrow, Elements of Modern
            % X-Ray Physics (John Wiley & Sons, 2011). 
            %
            
            thc = sqrt( 2 * delta ) * 180/pi;
            
        end
        
        function [n]     = get_xray_refractive_index(Z,atomic_density,hv,th)
            %
            % n = get_xray_refractive_index(Z,atomic_density,hv,th)
            % 
            % n              [unitless]     delta contribution to real part of dielectric function (n = 1 - delta + i beta)
            % Z              [Z]            list of atomic numbers
            % atomic_density [atoms/nm^3]   list of atomic density
            % hv             [eV]           photon energy
            %
            % Eq. 3 in B. L. Henke, E. M. Gullikson, and J. C. Davis,
            % Atomic Data and Nuclear Data Tables 54, 181 (1993). 
            %
            
            import am_mbe.*
            
            [f0,f1,f2] = get_atomic_xray_form_factor(Z,hv,th); % [nZs,nhvs,nths]
            n = 1 - am_mbe.r_0 ./ (2*pi) .* get_photon_wavelength(hv).^2 .* sum(( + f0 + f1 - f2*1i ).*atomic_density(:),1);
            n = permute(n,[3,2,1]); % [nhvs,nths]
        end

        function [mu]    = get_xray_linear_absorbtion_coeffcient(Z,atomic_density,hv)
            %
            % mu = get_xray_linear_absorbtion_coeffcient(Z,n,hv)
            % 
            % mu             [1/nm]         atomic linear absorbtion coefficient
            % Z              [Z]            list of atomic numbers
            % atomic_density [atoms/nm^3]   list of atomic number density
            % hv             [eV]           photon energy
            %

            import am_mbe.*
            
            mu = (4*pi) * get_xray_dielectric_function_beta(Z,atomic_density,hv) / get_photon_wavelength(hv);
            
        end

        function [qz]    = get_qz(th,hv)
            
            import am_mbe.*
            
            qz = 4*pi/get_photon_wavelength(hv) * sind(th);
            
        end

        function [th]    = get_th(qz,hv)
            
            import am_mbe.*
            
            th = asind( get_photon_wavelength(hv)*qz/(4*pi) );
            
        end
        
        function [th2]   = get_bragg_angle(hv,a,hkl)
            
            import am_mbe.*
            
            lambda = get_photon_wavelength(hv);
            
            % for cubic materials
            th2 = 2*asind(lambda*sqrt(sum(hkl.^2))/(2*a));
            
        end
        
    end
    
    % materials-related things
    
    methods (Static)
        
        function [layer] = define_material(species,stoichiometry,mass_density)

            import am_mbe.*

            define_layer = @(species,stoichiometry,mass_density) struct(...
                'species',{species},'Z',get_atomic_number(species),...
                'stoichiometry',stoichiometry,'mass_density',mass_density);
            layer = define_layer(species,stoichiometry,mass_density);
            
            % define layer basic properties
            % get molecular weight [amu/f.u.] = [amu/atom * atom/f.u.]
            layer.molecular_weight = sum(get_atomic_mass(layer.Z) .* layer.stoichiometry) ./ sum(layer.stoichiometry);
            % atomic density [atoms/nm^3]: (g/cm^3 * f.u./amu * atoms/f.u.) / (atoms/nm^3) = 602.24
            layer.atomic_density = layer.mass_density / layer.molecular_weight * 602.24;

        end
        
        function [h]     = cubic2hex(a)
            % hexagonal in-plane lattice parameter for 111-oriented of a cubic material
            h = a/sqrt(1-cosd(120));
            
        
        end
    end
    
    % atom-related things
    
    methods (Static)
        
        function [Z]     = get_atomic_number(symb)
            %
            % Z = get_atomic_number(symb)
            % 
            % symb  atomic symbol
            % Z     atomic number
            % 
            symbol_database = {...
                 'h'  ,'he' ,'li' ,'be' ,'b'  ,'c'  ,'n'  ,'o'  ,'f'  ,'ne' ,'na' ,'mg' ,'al' ,'si' ,'p'  ,'s'  , ...
                 'cl' ,'ar' ,'k'  ,'ca' ,'sc' ,'ti' ,'v'  ,'cr' ,'mn' ,'fe' ,'co' ,'ni' ,'cu' ,'zn' ,'ga' ,'ge' , ...
                 'as' ,'se' ,'br' ,'kr' ,'rb' ,'sr' ,'y'  ,'zr' ,'nb' ,'mo' ,'tc' ,'ru' ,'rh' ,'pd' ,'ag' ,'cd' , ...
                 'in' ,'sn' ,'sb' ,'te' ,'i'  ,'xe' ,'cs' ,'ba' ,'la' ,'ce' ,'pr' ,'nd' ,'pm' ,'sm' ,'eu' ,'gd' , ...
                 'tb' ,'dy' ,'ho' ,'er' ,'tm' ,'yb' ,'lu' ,'hf' ,'ta' ,'w'  ,'re' ,'os' ,'ir' ,'pt' ,'au' ,'hg' , ...
                 'tl' ,'pb' ,'bi' ,'po' ,'at' ,'rn' ,'fr' ,'ra' ,'ac' ,'th' ,'pa' ,'u'  ,'np' ,'pu' ,'am' ,'cm' , ...
                 'bk' ,'cf' ,'es' ,'fm' ,'md' ,'no' ,'lr' ,'rf' ,'db' ,'sg' ,'bh' ,'hs' ,'mt' ,'ds' ,'rg' ,'uub', ...
                 'uut','uuq','uup','uuh'}; 
            if iscell(symb)
                nZs = numel(symb); Z = zeros(1,nZs);
                for i = 1:nZs
                    Z(i) = find(strcmp(strtrim(lower(symb{i})),symbol_database));
                end
            else
                Z = find(strcmp(strtrim(lower(symb)),symbol_database));
            end
        end
        
        function [symb]  = get_atomic_symbol(Z)
            %
            % Z = get_atomic_number(symb)
            % 
            % symb  atomic symbol
            % Z     atomic number
            % 
            symbol_database = {...
                 'h'  ,'he' ,'li' ,'be' ,'b'  ,'c'  ,'n'  ,'o'  ,'f'  ,'ne' ,'na' ,'mg' ,'al' ,'si' ,'p'  ,'s'  , ...
                 'cl' ,'ar' ,'k'  ,'ca' ,'sc' ,'ti' ,'v'  ,'cr' ,'mn' ,'fe' ,'co' ,'ni' ,'cu' ,'zn' ,'ga' ,'ge' , ...
                 'as' ,'se' ,'br' ,'kr' ,'rb' ,'sr' ,'y'  ,'zr' ,'nb' ,'mo' ,'tc' ,'ru' ,'rh' ,'pd' ,'ag' ,'cd' , ...
                 'in' ,'sn' ,'sb' ,'te' ,'i'  ,'xe' ,'cs' ,'ba' ,'la' ,'ce' ,'pr' ,'nd' ,'pm' ,'sm' ,'eu' ,'gd' , ...
                 'tb' ,'dy' ,'ho' ,'er' ,'tm' ,'yb' ,'lu' ,'hf' ,'ta' ,'w'  ,'re' ,'os' ,'ir' ,'pt' ,'au' ,'hg' , ...
                 'tl' ,'pb' ,'bi' ,'po' ,'at' ,'rn' ,'fr' ,'ra' ,'ac' ,'th' ,'pa' ,'u'  ,'np' ,'pu' ,'am' ,'cm' , ...
                 'bk' ,'cf' ,'es' ,'fm' ,'md' ,'no' ,'lr' ,'rf' ,'db' ,'sg' ,'bh' ,'hs' ,'mt' ,'ds' ,'rg' ,'uub', ...
                 'uut','uuq','uup','uuh'}; 
            nZs = numel(Z); symb = cell(1,nZs);
            for i = 1:nZs
                symb(i) = symbol_database(Z(i));
            end
        end

        function [m]     = get_atomic_mass(Z)
            %
            % m = get_atomic_mass(Z)
            % 
            % Z     atomic number
            % m     atomic mass
            % 
            mass_database = [...
                    1.007947000,     4.002602000,     6.941200000,     9.012182000,    10.811500000, ...
                   12.011100000,    14.006747000,    15.999430000,    18.998403000,    20.179760000, ...
                   22.989769000,    24.305060000,    26.981540000,    28.085530000,    30.973762000, ...
                   32.066600000,    35.452790000,    39.948100000,    39.098310000,    40.078900000, ...
                   44.955911000,    47.883000000,    50.941510000,    51.996160000,    54.938051000, ...
                   55.847300000,    58.933201000,    58.693400000,    63.546300000,    65.392000000, ...
                   69.723100000,    72.612000000,    74.921592000,    78.963000000,    79.904000000, ...
                   83.801000000,    85.467830000,    87.621000000,    88.905852000,    91.224200000, ...
                   92.906382000,    95.941000000,    98.000000000,   101.072000000,   102.905503000, ...
                  106.421000000,   107.868220000,   112.411800000,   114.821000000,   118.710700000, ...
                  121.757000000,   127.603000000,   126.904473000,   131.292000000,   132.905435000, ...
                  137.327700000,   138.905520000,   140.115400000,   140.907653000,   144.243000000, ...
                  145.000000000,   150.363000000,   151.965900000,   157.253000000,   158.925343000, ...
                  162.503000000,   164.930323000,   167.263000000,   168.934213000,   173.043000000, ...
                  174.967100000,   178.492000000,   180.947910000,   183.853000000,   186.207100000, ...
                  190.210000000,   192.223000000,   195.083000000,   196.966543000,   200.593000000, ...
                  204.383320000,   207.210000000,   208.980373000,   209.000000000,   210.000000000, ...
                  222.000000000,   223.000000000,   226.025000000,   227.028000000,   232.038110000, ...
                  231.035900000,   238.028910000,   237.048000000,   244.000000000,   243.000000000, ...
                  247.000000000,   247.000000000,   251.000000000,   252.000000000,   257.000000000, ...
                  258.000000000,   259.000000000,   262.000000000,   261.000000000,   262.000000000, ...
                  263.000000000,   262.000000000,   265.000000000,   266.000000000]; 
            nZs = numel(Z); m = zeros(1,nZs);
            for i = 1:nZs
                m(i) = mass_database(Z(i));
            end
        end
        
        function [f0,f1,f2] = get_atomic_xray_form_factor(Z,hv,th)
            %
            % [f1,f2] = get_atomic_xray_form_factor(Z,hv,th) [nZs,nhvs,nths]
            % 
            % Z      [Z]            atomic number
            % hv     [eV]           photon energy
            % f0     [e-/atom]      atomic scattering factor
            % f1     [e-/atom]      atomic scattering factor
            % f2     [e-/atom]      atomic scattering factor
            %
            % Conversion factor = 
            %       Plank's constant * speed of light / nm / eV
            %
            % B. L. Henke, E. M. Gullikson, and J. C. Davis, Atomic Data
            % and Nuclear Data Tables 54, 181 (1993). Note: There is a
            % difference in convention between this reference and Warren.
            % Warren defines the anomolous atomic scattering factor as:
            %
            %       f = f0 + f1 + i f2
            %
            % This article clumps f0 and f1 together:
            %
            %       f = f' + i f''
            %
            % f' and f'' are provided in the dataset. To recover Warren's
            % original definition:
            % 
            %       f1 = f' - f0
            %
            % using f0 values computed from Cromer Mann coefficients.
            
            import am_mbe.*
            
            nZs = numel(Z); nths = numel(th); nhvs = numel(hv); f0 = zeros(nZs,nhvs,nths); 
            for i = 1:nZs
                % get f0 analytically from Cromer-Mann coefficients
                CM = get_atomic_cromer_mann_coefficients(Z(i));
                ai = reshape(CM(1:4),1,1,[]);
                bi = reshape(CM(5:8),1,1,[]);
                for j = 1:nhvs
                    f0(i,j,:) = sum( ai .* exp( - bi .* (sind(th)./get_photon_wavelength(hv(j))).^2 ) , 3) + CM(9);
                end
            end
            
            if nargout > 1
                f1 = zeros(nZs,nhvs,nths); f2 = zeros(nZs,nhvs,nths);
                for i = 1:nZs
                    % get f1 and f2 by interpolation
                    symb = get_atomic_symbol(Z);
                    fid = fopen(['data_atomic_scattering_factor/',strtrim(lower(symb{i})),'.nff']);
                        [~] = fgetl(fid); buffer = reshape(fscanf(fid,'%f'),3,[]).';
                    fclose(fid);

                    for j = 1:nhvs
                        f1(i,j,:) = interp1(buffer(:,1),buffer(:,2),hv(j)) - f0(i,j,:);
                        f2(i,j,:) = interp1(buffer(:,1),buffer(:,3),hv(j));
                    end
                end
            end
            
        
            function [CM]    = get_atomic_cromer_mann_coefficients(Z)
                %
                % CM = get_atomic_cromer_mann_coefficients(Z)
                % 
                % Z     atomic number
                % CM    array of Cromer-Mann coefficients [a1, a2, a3, a4, b1, b2, b3, b4, c]
                %
                % Atomic form factors are defined as:
                %
                %           4
                %     f0 = sum [ ai*exp(-bi*s^2) ] + c
                %          i=1
                %
                % for s = sin(theta)/lambda and Cromer-Mann coefficients ai,
                % bi, and c. 
                %
                % S. L. Morelhão, Computer Simulation Tools for X-Ray Analysis
                % (Springer International Publishing, Cham, 2016). 
                %

                CM_database = [ ...
                   0.489918,  0.262003,  0.196767, 0.049879,20.659300, 7.740390, 49.551900,  2.201590,  0.001305;  % H        
                   0.873400,  0.630900,  0.311200, 0.178000, 9.103700, 3.356800, 22.927600,  0.982100,  0.006400;  % He       
                   1.128200,  0.750800,  0.617500, 0.465300, 3.954600, 1.052400, 85.390500,168.261000,  0.037700;  % Li       
                   1.591900,  1.127800,  0.539100, 0.702900,43.642700, 1.862300,103.483000,  0.542000,  0.038500;  % Be       
                   2.054500,  1.332600,  1.097900, 0.706800,23.218500, 1.021000, 60.349800,  0.140300, -0.193200;  % B        
                   2.260690,  1.561650,  1.050750, 0.839259,22.690700, 0.656665,  9.756180, 55.594900,  0.286977;  % C        
                  12.212600,  3.132200,  2.012500, 1.166300, 0.005700, 9.893300, 28.997500,  0.582600,-11.529000;  % N        
                   3.048500,  2.286800,  1.546300, 0.867000,13.277100, 5.701100,  0.323900, 32.908900,  0.250800;  % O        
                   3.539200,  2.641200,  1.517000, 1.024300,10.282500, 4.294400,  0.261500, 26.147600,  0.277600;  % F        
                   3.955300,  3.112500,  1.454600, 1.125100, 8.404200, 3.426200,  0.230600, 21.718400,  0.351500;  % Ne       
                   4.762600,  3.173600,  1.267400, 1.112800, 3.285000, 8.842200,  0.313600,129.424000,  0.676000;  % Na       
                   5.420400,  2.173500,  1.226900, 2.307300, 2.827500,79.261100,  0.380800,  7.193700,  0.858400;  % Mg       
                   6.420200,  1.900200,  1.593600, 1.964600, 3.038700, 0.742600, 31.547200, 85.088600,  1.115100;  % Al       
                   5.662690,  3.071640,  2.624460, 1.393200, 2.665200,38.663400,  0.916946, 93.545800,  1.247070;  % Si       
                   6.434500,  4.179100,  1.780000, 1.490800, 1.906700,27.157000,  0.526000, 68.164500,  1.114900;  % P        
                   6.905300,  5.203400,  1.437900, 1.586300, 1.467900,22.215100,  0.253600, 56.172000,  0.866900;  % S        
                  11.460400,  7.196400,  6.255600, 1.645500, 0.010400, 1.166200, 18.519400, 47.778400, -9.557400;  % Cl       
                   7.484500,  6.772300,  0.653900, 1.644200, 0.907200,14.840700, 43.898300, 33.392900,  1.444500;  % Ar       
                   8.218600,  7.439800,  1.051900, 0.865900,12.794900, 0.774800,213.187000, 41.684100,  1.422800;  % K        
                   8.626600,  7.387300,  1.589900, 1.021100,10.442100, 0.659900, 85.748400,178.437000,  1.375100;  % Ca       
                   9.189000,  7.367900,  1.640900, 1.468000, 9.021300, 0.572900,136.108000, 51.353100,  1.332900;  % Sc       
                   9.759500,  7.355800,  1.699100, 1.902100, 7.850800, 0.500000, 35.633800,116.105000,  1.280700;  % Ti       
                  10.297100,  7.351100,  2.070300, 2.057100, 6.865700, 0.438500, 26.893800,102.478000,  1.219900;  % V        
                  10.640600,  7.353700,  3.324000, 1.492200, 6.103800, 0.392000, 20.262600, 98.739900,  1.183200;  % Cr       
                  11.281900,  7.357300,  3.019300, 2.244100, 5.340900, 0.343200, 17.867400, 83.754300,  1.089600;  % Mn       
                  11.769500,  7.357300,  3.522200, 2.304500, 4.761100, 0.307200, 15.353500, 76.880500,  1.036900;  % Fe       
                  12.284100,  7.340900,  4.003400, 2.348800, 4.279100, 0.278400, 13.535900, 71.169200,  1.011800;  % Co       
                  12.837600,  7.292000,  4.443800, 2.380000, 3.878500, 0.256500, 12.176300, 66.342100,  1.034100;  % Ni       
                  13.338000,  7.167600,  5.615800, 1.673500, 3.582800, 0.247000, 11.396600, 64.812600,  1.191000;  % Cu       
                  14.074300,  7.031800,  5.165200, 2.410000, 3.265500, 0.233300, 10.316300, 58.709700,  1.304100;  % Zn       
                  15.235400,  6.700600,  4.359100, 2.962300, 3.066900, 0.241200, 10.780500, 61.413500,  1.718900;  % Ga       
                  16.081600,  6.374700,  3.706800, 3.683000, 2.850900, 0.251600, 11.446800, 54.762500,  2.131300;  % Ge       
                  16.672300,  6.070100,  3.431300, 4.277900, 2.634500, 0.264700, 12.947900, 47.797200,  2.531000;  % As       
                  17.000600,  5.819600,  3.973100, 4.354300, 2.409800, 0.272600, 15.237200, 43.816300,  2.840900;  % Se       
                  17.178900,  5.235800,  5.637700, 3.985100, 2.172300,16.579600,  0.260900, 41.432800,  2.955700;  % Br       
                  17.355500,  6.728600,  5.549300, 3.537500, 1.938400,16.562300,  0.226100, 39.397200,  2.825000;  % Kr       
                  17.178400,  9.643500,  5.139900, 1.529200, 1.788800,17.315100,  0.274800,164.934000,  3.487300;  % Rb       
                  17.566300,  9.818400,  5.422000, 2.669400, 1.556400,14.098800,  0.166400,132.376000,  2.506400;  % Sr       
                  17.776000, 10.294600,  5.726290, 3.265880, 1.402900,12.800600,  0.125599,104.354000,  1.912130;  % Y        
                  17.876500, 10.948000,  5.417320, 3.657210, 1.276180,11.916000,  0.117622, 87.662700,  2.069290;  % Zr       
                  17.614200, 12.014400,  4.041830, 3.533460, 1.188650,11.766000,  0.204785, 69.795700,  3.755910;  % Nb       
                   3.702500, 17.235600, 12.887600, 3.742900, 0.277200, 1.095800, 11.004000, 61.658400,  4.387500;  % Mo       
                  19.130100, 11.094800,  4.649010, 2.712630, 0.864132, 8.144870, 21.570700, 86.847200,  5.404280;  % Tc       
                  19.267400, 12.918200,  4.863370, 1.567560, 0.808520, 8.434670, 24.799700, 94.292800,  5.378740;  % Ru       
                  19.295700, 14.350100,  4.734250, 1.289180, 0.751536, 8.217580, 25.874900, 98.606200,  5.328000;  % Rh       
                  19.331900, 15.501700,  5.295370, 0.605844, 0.698655, 7.989290, 25.205200, 76.898600,  5.265930;  % Pd       
                  19.280800, 16.688500,  4.804500, 1.046300, 0.644600, 7.472600, 24.660500, 99.815600,  5.179000;  % Ag       
                  19.221400, 17.644400,  4.461000, 1.602900, 0.594600, 6.908900, 24.700800, 87.482500,  5.069400;  % Cd       
                  19.162400, 18.559600,  4.294800, 2.039600, 0.547600, 6.377600, 25.849900, 92.802900,  4.939100;  % In       
                  19.188900, 19.100500,  4.458500, 2.466300, 5.830300, 0.503100, 26.890900, 83.957100,  4.782100;  % Sn       
                  19.641800, 19.045500,  5.037100, 2.682700, 5.303400, 0.460700, 27.907400, 75.282500,  4.590900;  % Sb       
                  19.964400, 19.013800,  6.144870, 2.523900, 4.817420, 0.420885, 28.528400, 70.840300,  4.352000;  % Te       
                  20.147200, 18.994900,  7.513800, 2.273500, 4.347000, 0.381400, 27.766000, 66.877600,  4.071200;  % I        
                  20.293300, 19.029800,  8.976700, 1.990000, 3.928200, 0.344000, 26.465900, 64.265800,  3.711800;  % Xe       
                  20.389200, 19.106200, 10.662000, 1.495300, 3.569000, 0.310700, 24.387900,213.904000,  3.335200;  % Cs       
                  20.336100, 19.297000, 10.888000, 2.695900, 3.216000, 0.275600, 20.207300,167.202000,  2.773100;  % Ba       
                  20.578000, 19.599000, 11.372700, 3.287190, 2.948170, 0.244475, 18.772600,133.124000,  2.146780;  % La       
                  21.167100, 19.769500, 11.851300, 3.330490, 2.812190, 0.226836, 17.608300,127.113000,  1.862640;  % Ce       
                  22.044000, 19.669700, 12.385600, 2.824280, 2.773930, 0.222087, 16.766900,143.644000,  2.058300;  % Pr       
                  22.684500, 19.684700, 12.774000, 2.851370, 2.662480, 0.210628, 15.885000,137.903000,  1.984860;  % Nd       
                  23.340500, 19.609500, 13.123500, 2.875160, 2.562700, 0.202088, 15.100900,132.721000,  2.028760;  % Pm       
                  24.004200, 19.425800, 13.439600, 2.896040, 2.472740, 0.196451, 14.399600,128.007000,  2.209630;  % Sm       
                  24.627400, 19.088600, 13.760300, 2.922700, 2.387900, 0.194200, 13.754600,123.174000,  2.574500;  % Eu       
                  25.070900, 19.079800, 13.851800, 3.545450, 2.253410, 0.181951, 12.933100,101.398000,  2.419600;  % Gd       
                  25.897600, 18.218500, 14.316700, 2.953540, 2.242560, 0.196143, 12.664800,115.362000,  3.583240;  % Tb       
                  26.507000, 17.638300, 14.559600, 2.965770, 2.180200, 0.202172, 12.189900,111.874000,  4.297280;  % Dy       
                  26.904900, 17.294000, 14.558300, 3.638370, 2.070510, 0.197940, 11.440700, 92.656600,  4.567960;  % Ho       
                  27.656300, 16.428500, 14.977900, 2.982330, 2.073560, 0.223545, 11.360400,105.703000,  5.920460;  % Er       
                  28.181900, 15.885100, 15.154200, 2.987060, 2.028590, 0.238849, 10.997500,102.961000,  6.756210;  % Tm       
                  28.664100, 15.434500, 15.308700, 2.989630, 1.988900, 0.257119, 10.664700,100.417000,  7.566720;  % Yb       
                  28.947600, 15.220800, 15.100000, 3.716010, 1.901820, 9.985190,  0.261033, 84.329800,  7.976280;  % Lu       
                  29.144000, 15.172600, 14.758600, 4.300130, 1.832620, 9.599900,  0.275116, 72.029000,  8.581540;  % Hf       
                  29.202400, 15.229300, 14.513500, 4.764920, 1.773330, 9.370460,  0.295977, 63.364400,  9.243540;  % Ta       
                  29.081800, 15.430000, 14.432700, 5.119820, 1.720290, 9.225900,  0.321703, 57.056000,  9.887500;  % W        
                  28.762100, 15.718900, 14.556400, 5.441740, 1.671910, 9.092270,  0.350500, 52.086100, 10.472000;  % Re       
                  28.189400, 16.155000, 14.930500, 5.675890, 1.629030, 8.979480,  0.382661, 48.164700, 11.000500;  % Os       
                  27.304900, 16.729600, 15.611500, 5.833770, 1.592790, 8.865530,  0.417916, 45.001100, 11.472200;  % Ir       
                  27.005900, 17.763900, 15.713100, 5.783700, 1.512930, 8.811740,  0.424593, 38.610300, 11.688300;  % Pt       
                  16.881900, 18.591300, 25.558200, 5.860000, 0.461100, 8.621600,  1.482600, 36.395600, 12.065800;  % Au       
                  20.680900, 19.041700, 21.657500, 5.967600, 0.545000, 8.448400,  1.572900, 38.324600, 12.608900;  % Hg       
                  27.544600, 19.158400, 15.538000, 5.525930, 0.655150, 8.707510,  1.963470, 45.814900, 13.174600;  % Tl       
                  31.061700, 13.063700, 18.442000, 5.969600, 0.690200, 2.357600,  8.618000, 47.257900, 13.411800;  % Pb       
                  33.368900, 12.951000, 16.587700, 6.469200, 0.704000, 2.923800,  8.793700, 48.009300, 13.578200;  % Bi       
                  34.672600, 15.473300, 13.113800, 7.025880, 0.700999, 3.550780,  9.556420, 47.004500, 13.677000;  % Po       
                  35.316300, 19.021100,  9.498870, 7.425180, 0.685870, 3.974580, 11.382400, 45.471500, 13.710800;  % At       
                  35.563100, 21.281600,  8.003700, 7.443300, 0.663100, 4.069100, 14.042200, 44.247300, 13.690500;  % Rn       
                  35.929900, 23.054700, 12.143900, 2.112530, 0.646453, 4.176190, 23.105200,150.645000, 13.724700;  % Fr       
                  35.763000, 22.906400, 12.473900, 3.210970, 0.616341, 3.871350, 19.988700,142.325000, 13.621100;  % Ra       
                  35.659700, 23.103200, 12.597700, 4.086550, 0.589092, 3.651550, 18.599000,117.020000, 13.526600;  % Ac       
                  35.564500, 23.421900, 12.747300, 4.807030, 0.563359, 3.462040, 17.830900, 99.172200, 13.431400;  % Th       
                  35.884700, 23.294800, 14.189100, 4.172870, 0.547751, 3.415190, 16.923500,105.251000, 13.428700;  % Pa       
                  36.022800, 23.412800, 14.949100, 4.188000, 0.529300, 3.325300, 16.092700,100.613000, 13.396600;  % U        
                  36.187400, 23.596400, 15.640200, 4.185500, 0.511929, 3.253960, 15.362200, 97.490800, 13.357300;  % Np       
                  35.510300, 22.578700, 12.776600, 4.921590, 0.498626, 2.966270, 11.948400, 22.750200, 13.211600;  % Pu       
                  36.670600, 24.099200, 17.341500, 3.493310, 0.483629, 3.206470, 14.313600,102.273000, 13.359200;  % Am       
                  36.648800, 24.409600, 17.399000, 4.216650, 0.465154, 3.089970, 13.434600, 88.483400, 13.288700;  % Cm       
                  36.788100, 24.773600, 17.891900, 4.232840, 0.451018, 3.046190, 12.894600, 86.003000, 13.275400;  % Bk       
                  36.918500, 25.199500, 18.331700, 4.243910, 0.437533, 3.007750, 12.404400, 83.788100, 13.267400]; % Cf       
              CM = CM_database(Z,:);
            end
        end

        function [hv]    = get_atomic_emission_line_energy(Z,emission_line)
            %
            % hv = get_xray_energy(Z,emission_line)
            % 
            % hv            [eV]           energy
            % Z             [Z]            atomic number
            % emission_line [string]       Ka1, Ka2, Kb1, La1, La2, Lb1, Lb2, Lg1
            %
            % J. A. Bearden, Reviews of Modern Physics 39, 78 (1967).
            %
            
            xray_database = [ ...
                0      , 0       , 0      , 0      , 0      , 0      , 0      , 0      ; % 1   H   
                0      , 0       , 0      , 0      , 0      , 0      , 0      , 0      ; % 2   He  
                0.0543 , 0       , 0      , 0      , 0      , 0      , 0      , 0      ; % 3   Li  
                0.1085 , 0       , 0      , 0      , 0      , 0      , 0      , 0      ; % 4   Be  
                0.1833 , 0       , 0      , 0      , 0      , 0      , 0      , 0      ; % S   B   
                0.277  , 0       , 0      , 0      , 0      , 0      , 0      , 0      ; % 6   C   
                0.3924 , 0       , 0      , 0      , 0      , 0      , 0      , 0      ; % 7   N   
                0.5249 , 0       , 0      , 0      , 0      , 0      , 0      , 0      ; % 8   O   
                0.6768 , 0       , 0      , 0      , 0      , 0      , 0      , 0      ; % 9   F   
                0.8486 , 0.8486  , 0      , 0      , 0      , 0      , 0      , 0      ; % 10  Ne  
                1.04098, 1.04098 , 1.0711 , 0      , 0      , 0      , 0      , 0      ; % 11  Na  
                1.25360, 1.25360 , 1.3022 , 0      , 0      , 0      , 0      , 0      ; % 12  Mg  
                1.48670, 1.48627 , 1.55745, 0      , 0      , 0      , 0      , 0      ; % 13  Al  
                1.73998, 1.73938 , 1.83594, 0      , 0      , 0      , 0      , 0      ; % 14  Si  
                2.0137 , 2.0127  , 2.1391 , 0      , 0      , 0      , 0      , 0      ; % 15  P   
                2.30784, 2.30664 , 2.46404, 0      , 0      , 0      , 0      , 0      ; % 16  S   
                2.62239, 2.62078 , 2.8156 , 0      , 0      , 0      , 0      , 0      ; % 17  Cl  
                2.95770, 2.95563 , 3.1905 , 0      , 0      , 0      , 0      , 0      ; % 18  Ar  
                3.3138 , 3.3111  , 3.5896 , 0      , 0      , 0      , 0      , 0      ; % 19  K   
                3.69168, 3.68809 , 4.0127 , 0.3413 , 0.3413 , 0.3449 , 0      , 0      ; % 20  Ca  
                4.0906 , 4.0861  , 4.4605 , 0.3954 , 0.3954 , 0.3996 , 0      , 0      ; % 21  Sc  
                4.51084, 4.50486 , 4.93181, 0.4522 , 0.4522 , 0.4584 , 0      , 0      ; % 22  Ti  
                4.95220, 4.94464 , 5.42729, 0.5113 , 0.5113 , 0.5192 , 0      , 0      ; % 23  V   
                5.41472, 5.405509, 5.94671, 0.5728 , 0.5728 , 0.5828 , 0      , 0      ; % 24  Cr  
                5.89875, 5.88765 , 6.49045, 0.6374 , 0.6374 , 0.6488 , 0      , 0      ; % 25  Mn  
                6.40384, 6.39084 , 7.05798, 0.7050 , 0.7050 , 0.7185 , 0      , 0      ; % 26  Fe  
                6.93032, 6.91530 , 7.64943, 0.7762 , 0.7762 , 0.7914 , 0      , 0      ; % 27  Co  
                7.47815, 7.46089 , 8.26466, 0.8515 , 0.8515 , 0.8688 , 0      , 0      ; % 28  Ni  
                8.04778, 8.02783 , 8.90529, 0.9297 , 0.9297 , 0.9498 , 0      , 0      ; % 29  Cu  
                8.63886, 8.61578 , 9.5720 , 1.0117 , 1.0117 , 1.0347 , 0      , 0      ; % 30  Zn  
                9.25174, 9.22482 , 10.2642, 1.09792, 1.09792, 1.1248 , 0      , 0      ; % 31  Ga  
                9.88642, 9.85532 , 10.9821, 1.18800, 1.18800, 1.2185 , 0      , 0      ; % 32  Ge  
                10.5437, 10.50799, 11.7262, 1.2820 , 1.2820 , 1.3170 , 0      , 0      ; % 33  As  
                11.2224, 11.1814 , 12.4959, 1.37910, 1.37910, 1.41923, 0      , 0      ; % 34  Se  
                11.9242, 11.8776 , 13.2914, 1.48043, 1.48043, 1.52590, 0      , 0      ; % 35  Br  
                12.649 , 12.598  , 14.112 , 1.5860 , 1.5860 , 1.6366 , 0      , 0      ; % 36  Kr  
                13.3953, 13.3358 , 14.9613, 1.69413, 1.69256, 1.75217, 0      , 0      ; % 37  Rb  
                14.1650, 14.0979 , 15.8357, 1.80656, 1.80474, 1.87172, 0      , 0      ; % 38  Sr  
                14.9584, 14.8829 , 16.7378, 1.92256, 1.92047, 1.99584, 0      , 0      ; % 39  Y   
                15.7751, 15.6909 , 17.6678, 2.04236, 2.0399 , 2.1244 , 2.2194 , 2.3027 ; % 40  Zr  
                16.6151, 16.5210 , 18.6225, 2.16589, 2.1630 , 2.2574 , 2.3670 , 2.4618 ; % 41  Nb  
                17.4793, 17.3743 , 19.6083, 2.29316, 2.28985, 2.39481, 2.5183 , 2.6235 ; % 42  Mo  
                18.3671, 18.2508 , 20.619 , 2.4240 , 0      , 2.5368 , 0      , 0      ; % 43  Tc  
                19.2792, 19.1504 , 21.6568, 2.55855, 2.55431, 2.68323, 2.8360 , 2.9645 ; % 44  Ru  
                20.2161, 20.0737 , 22.7236, 2.69674, 2.69205, 2.83441, 3.0013 , 3.1438 ; % 45  Rh  
                21.1771, 21.0201 , 23.8187, 2.83861, 2.83325, 2.99022, 3.17179, 3.3287 ; % 46  Pd  
                22.1629, 21.9903 , 24.9424, 2.98431, 2.97821, 3.15094, 3.34781, 3.51959; % 47  Ag  
                23.1736, 22.9841 , 26.0955, 3.13373, 3.12691, 3.31657, 3.52812, 3.71686; % 48  Cd  
                24.2097, 24.0020 , 27.2759, 3.28694, 3.27929, 3.48721, 3.71381, 3.92081; % 49  In  
                25.2713, 25.0440 , 28.4860, 3.44398, 3.43542, 3.66280, 3.90486, 4.13112; % 50  Sn  
                26.3591, 26.1108 , 29.7256, 3.60472, 3.59532, 3.84357, 4.10078, 4.34779; % 51  Sb  
                27.4723, 27.2017 , 30.9957, 3.76933, 3.7588 , 4.02958, 4.3017 , 4.5709 ; % 52  Te  
                28.6120, 28.3172 , 32.2947, 3.93765, 3.92604, 4.22072, 4.5075 , 4.8009 ; % 53  I   
                29.779 , 29.458  , 33.624 , 4.1099 , 0      , 0      , 0      , 0      ; % 54  Xe  
                30.9728, 30.6251 , 34.9869, 4.2865 , 4.2722 , 4.6198 , 4.9359 , 5.2804 ; % 55  Cs  
                32.1936, 31.8171 , 36.3782, 4.46626, 4.45090, 4.82753, 5.1565 , 5.5311 ; % 56  Ba  
                33.4418, 33.0341 , 37.8010, 4.65097, 4.63423, 5.0421 , 5.3835 , 5.7885 ; % 57  La  
                34.7197, 34.2789 , 39.2573, 4.8402 , 4.8230 , 5.2622 , 5.6134 , 6.052  ; % 58  Ce  
                36.0263, 35.5502 , 40.7482, 5.0337 , 5.0135 , 5.4889 , 5.850  , 6.3221 ; % 59  Pr  
                37.3610, 36.8474 , 42.2713, 5.2304 , 5.2077 , 5.7216 , 6.0894 , 6.6021 ; % 60  Nd  
                38.7247, 38.1712 , 43.826 , 5.4325 , 5.4078 , 5.961  , 6.339  , 6.892  ; % 61  Pm  
                40.1181, 39.5224 , 45.413 , 5.6361 , 5.6090 , 6.2051 , 6.586  , 7.178  ; % 62  Sm  
                41.5422, 40.9019 , 47.0379, 5.8457 , 5.8166 , 6.4564 , 6.8432 , 7.4803 ; % 63  Eu  
                42.9962, 42.3089 , 48.697 , 6.0572 , 6.0250 , 6.7132 , 7.1028 , 7.7858 ; % 64  Gd  
                44.4816, 43.7441 , 50.382 , 6.2728 , 6.2380 , 6.978  , 7.3667 , 8.102  ; % 65  Tb  
                45.9984, 45.2078 , 52.119 , 6.4952 , 6.4577 , 7.2477 , 7.6357 , 8.4188 ; % 66  Dy  
                47.5467, 46.6997 , 53.877 , 6.7198 , 6.6795 , 7.5253 , 7.911  , 8.747  ; % 67  Ho  
                49.1277, 48.2211 , 55.681 , 6.9487 , 6.9050 , 7.8109 , 8.1890 , 9.089  ; % 68  Er  
                50.7416, 49.7726 , 57.517 , 7.1799 , 7.1331 , 8.101  , 8.468  , 9.426  ; % 69  Tm  
                52.3889, 51.3540 , 59.37  , 7.4156 , 7.3673 , 8.4018 , 8.7588 , 9.7801 ; % 70  Yb  
                54.0698, 52.9650 , 61.283 , 7.6555 , 7.6049 , 8.7090 , 9.0489 , 10.1434; % 71  Lu  
                55.7902, 54.6114 , 63.234 , 7.8990 , 7.8446 , 9.0227 , 9.3473 , 10.5158; % 72  Hf  
                57.532 , 56.277  , 65.223 , 8.1461 , 8.0879 , 9.3431 , 9.6518 , 10.8952; % 73  Ta  
                59.3182, 57.9817 , 67.2443, 8.3976 , 8.3352 , 9.67235, 9.9615 , 11.2859; % 74  W   
                61.1403, 59.7179 , 69.310 , 8.6525 , 8.5862 , 10.0100, 10.2752, 11.6854; % 75  Re  
                63.0005, 61.4867 , 71.413 , 8.9117 , 8.8410 , 10.3553, 10.5985, 12.0953; % 76  Os  
                64.8956, 63.2867 , 73.5608, 9.1751 , 9.0995 , 10.7083, 10.9203, 12.5126; % 77  Ir  
                66.832 , 65.112  , 75.748 , 9.4423 , 9.3618 , 11.0707, 11.2505, 12.9420; % 78  Pt  
                68.8037, 66.9895 , 77.984 , 9.7133 , 9.6280 , 11.4423, 11.5847, 13.3817; % 79  Au  
                70.819 , 68.895  , 80.253 , 9.9888 , 9.8976 , 11.8226, 11.9241, 13.8301; % 80  Hg  
                72.8715, 70.8319 , 82.576 , 10.2685, 10.1728, 12.2133, 12.2715, 14.2915; % 81  Tl  
                74.9694, 72.8042 , 84.936 , 10.5515, 10.4495, 12.6137, 12.6226, 14.7644; % 82  Pb  
                77.1079, 74.8148 , 87.343 , 10.8388, 10.7309, 13.0235, 12.9799, 15.2477; % 83  Bi  
                79.290 , 76.862  , 89.80  , 11.1308, 11.0158, 13.447 , 13.3404, 15.744 ; % 84  Po  
                81.52  , 78.95   , 92.30  , 11.4268, 11.3048, 13.876 , 0      , 16.251 ; % 8S  At  
                83.78  , 81.07   , 94.87  , 11.7270, 11.5979, 14.316 , 0      , 16.770 ; % 86  Rn  
                86.10  , 83.23   , 97.47  , 12.0313, 11.8950, 14.770 , 14.45  , 17.303 ; % 87  Fr  
                88.47  , 85.43   , 100.13 , 12.3397, 12.1962, 15.2358, 14.8414, 17.849 ; % 88  Ra  
                90.884 , 87.67   , 102.85 , 12.6520, 12.5008, 15.713 , 0      , 18.408 ; % 89  Ac  
                93.350 , 89.953  , 105.609, 12.9687, 12.8096, 16.2022, 15.6237, 18.9825; % 90  Th  
                95.868 , 92.287  , 108.427, 13.2907, 13.1222, 16.702 , 16.024 , 19.568 ; % 91  Pa  
                98.439 , 94.665  , 111.300, 13.6147, 13.4388, 17.2200, 16.4283, 20.16 ]; % 92  U
            
            line_database = {'kalpha1','kalpha2','kbeta1','lalpha1','lalpha2','lbeta1','lbeta2','lgamma1'};
            j = find(strcmp(strtrim(lower(emission_line)),line_database));

            hv = xray_database(Z,j)*1E3;

        end
        
        function [hv]    = get_photon_energy(lambda)
            %
            % E = get_photon_energy(lambda)
            % 
            % E      [eV]           photon energy
            % lambda [nm]           photon wavelength
            %
            % Conversion factor = 
            %       Plank's constant * speed of light / nm / eV
            %
            
            hv = 1239.842 ./ lambda;
            
        end
        
        function [lambda]= get_photon_wavelength(hv)
            %
            % E = get_photon_energy(lambda)
            % 
            % E      [eV]           photon energy
            % lambda [nm]           photon wavelength
            %
            % Conversion factor = 
            %       Plank's constant * speed of light / nm / eV
            %
            
            lambda = 1239.842 ./ hv;
            
        end
        
        
        
    end
    
end














