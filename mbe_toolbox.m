classdef mbe_toolbox

    properties (Constant)
        tiny      = 1E-4; % precision of atomic coordinates
    end
    
    % pressure related things

    methods (Static)
        
        function P = flux2pressure(J,T,m,th)
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
        
        function J = pressure2flux(P,T,m,th)
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
        
        function mfp = pressure2mfp(P,T,sigma)
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
    
    % periodic table related things
    
    methods (Static)
        
        function [Z] = symb2Z(symb)
        s = {'h'  ,'he' ,'li' ,'be' ,'b'  ,'c'  ,'n'  ,'o'  ,'f'  ,'ne' ,'na' ,'mg' ,'al' ,'si' ,'p'  ,'s'  , ...
             'cl' ,'ar' ,'k'  ,'ca' ,'sc' ,'ti' ,'v'  ,'cr' ,'mn' ,'fe' ,'co' ,'ni' ,'cu' ,'zn' ,'ga' ,'ge' , ...
             'as' ,'se' ,'br' ,'kr' ,'rb' ,'sr' ,'y'  ,'zr' ,'nb' ,'mo' ,'tc' ,'ru' ,'rh' ,'pd' ,'ag' ,'cd' , ...
             'in' ,'sn' ,'sb' ,'te' ,'i'  ,'xe' ,'cs' ,'ba' ,'la' ,'ce' ,'pr' ,'nd' ,'pm' ,'sm' ,'eu' ,'gd' , ...
             'tb' ,'dy' ,'ho' ,'er' ,'tm' ,'yb' ,'lu' ,'hf' ,'ta' ,'w'  ,'re' ,'os' ,'ir' ,'pt' ,'au' ,'hg' , ...
             'tl' ,'pb' ,'bi' ,'po' ,'at' ,'rn' ,'fr' ,'ra' ,'ac' ,'th' ,'pa' ,'u'  ,'np' ,'pu' ,'am' ,'cm' , ...
             'bk' ,'cf' ,'es' ,'fm' ,'md' ,'no' ,'lr' ,'rf' ,'db' ,'sg' ,'bh' ,'hs' ,'mt' ,'ds' ,'rg' ,'uub', ...
             'uut','uuq','uup','uuh'}; Z = find(strcmp(strtrim(lower(symb)),s));
        end
        
        function [mass] = Z2mass(Z)
        m = [   1.007947000,     4.002602000,     6.941200000,     9.012182000,    10.811500000, ...
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
              263.000000000,   262.000000000,   265.000000000,   266.000000000]; mass = m(Z);
        end
        
    end
    
end








