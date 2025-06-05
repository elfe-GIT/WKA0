classdef model
    %all properties of system
    %   being computed here ...
    
    properties

        % system matrices
        MO
        MC
        MS
        GO
        GC
        GS
        KK
        KO
        KC
        KS
        % Omega
        Oga
        i
        tEnd
    end
    
    methods
        function obj = model(data)
            %initialize
            %   class
            SIunits = containers.Map;
            SIunits('1'      ) = 1;
            SIunits('mm'     ) = 1/1000;
            SIunits('cm'     ) = 1/100;
            SIunits('m'      ) = 1;
            SIunits('m^2'    ) = 1;
            SIunits('kg'     ) = 1;
            SIunits('km/h'   ) = 1/3.6;
            SIunits('N'      ) = 1;
            SIunits('N*m^2'  ) = 1;
            SIunits('N/m^2'  ) = 1;
            SIunits('kg*m^2' ) = 1;
            SIunits('kg/m^3' ) = 1;
            SIunits('m/s'    ) = 1;
            SIunits('s'      ) = 1;
            SIunits('N*m/rad') = 1;
            SIunits('N/(m/s)' ) = 1;
            SIunits('U/min'   ) = 1/60;
                     
            % ++++++++++++++++ read from database ++++++++++++++++++++++++

            % tower
            l_T = data({'l[T]'},:).number*SIunits(data({'l[T]'},:).unit{1});
            E_T = data({'E[T]'},:).number*SIunits(data({'E[T]'},:).unit{1});
            G_T = data({'G[T]'},:).number*SIunits(data({'G[T]'},:).unit{1});
            rho_T = data({'ρ[T]'},:).number*SIunits(data({'ρ[T]'},:).unit{1});
            D_T = data({'D[T]'},:).number*SIunits(data({'D[T]'},:).unit{1});
            t_T = data({'t[T]'},:).number*SIunits(data({'t[T]'},:).unit{1});

            % nacelle (rigid body)
            m_N = data({'m[N]'},:).number*SIunits(data({'m[N]'},:).unit{1});
            J_N1 = data({'J[N,1]'},:).number*SIunits(data({'J[N,1]'},:).unit{1});
            J_N2 = data({'J[N,2]'},:).number*SIunits(data({'J[N,2]'},:).unit{1});
            J_N3 = data({'J[N,3]'},:).number*SIunits(data({'J[N,3]'},:).unit{1});
            b = data({'b'},:).number*SIunits(data({'b'},:).unit{1});

            % rotor-shaft
            c = data({'c'},:).number*SIunits(data({'c'},:).unit{1});
            D_S = data({'D[S]'},:).number*SIunits(data({'D[S]'},:).unit{1});
            E_S = data({'E[S]'},:).number*SIunits(data({'E[S]'},:).unit{1});

            % rotor
            m_R = data({'m[R]'},:).number*SIunits(data({'m[R]'},:).unit{1});
            R = data({'R'},:).number*SIunits(data({'R'},:).unit{1});
            J_R = data({'J[R]'},:).number*SIunits(data({'J[R]'},:).unit{1});
      
            % blade
            l_B = data({'l[B]'},:).number*SIunits(data({'l[B]'},:).unit{1});
            m_B = data({'m[B]'},:).number*SIunits(data({'m[B]'},:).unit{1});
            EI_B = data({'EI[B]'},:).number*SIunits(data({'EI[B]'},:).unit{1}); % calculated from eigenfreuquency
            A_B = data({'A[B]'},:).number*SIunits(data({'A[B]'},:).unit{1});
            E_B = data({'E[B]'},:).number*SIunits(data({'E[B]'},:).unit{1});

            % parameter study
            N(1) = data({'N[0]'},:).number*SIunits(data({'N[0]'},:).unit{1});
            N(2) = data({'N[1]'},:).number*SIunits(data({'N[1]'},:).unit{1});
            N(3) = data({'ΔN'},:).number*SIunits(data({'ΔN'},:).unit{1});

            % derived properties

            obj.Oga = 2*pi()*linspace(N(1),N(2),ceil((N(2)-N(1))/N(3))+1);
      
            EI_T = E_T*pi()*(D_T/2)^3*t_T;
            i_T  = rho_T * pi()*D_T*t_T * (D_T/2)^2;
            GI_T = G_T* pi()* D_T^3*t_T/16;
            EI_S = E_S*pi()* D_S^4/64;
            %EI_B = ...
            EA_B = E_B*A_B;
            m_T = rho_T*2*pi()*(D_T/2)*t_T*l_B;
                        
            % ................................................
            
            
            
            
            MO = zeros(14,14);
            MC = zeros(14,14);
            MS = zeros(14,14);

            GO = zeros(14,14);
            GC = zeros(14,14);
            GS = zeros(14,14);

            KK = zeros(14,14);
            KO = zeros(14,14);
            KC = zeros(14,14);
            KS = zeros(14,14);

            % start: from matrices.m

            MO( 1 , 1 ) =  m_T/5+m_N+(2*l_B^2*m_B)/l_T^2+(6*R*l_B*m_B)/l_T^2+(12*c^2*m_B)/l_T^2+(24*b*c*m_B)/l_T^2+(12*b^2*m_B)/l_T^2+(6*R^2*m_B)/l_T^2+3*m_B+(4*J_N2)/l_T^2 ; 
            MO( 1 , 6 ) =  m_B/3 ; 
            MO( 1 , 9 ) =  m_B/3 ; 
            MO( 1 , 12 ) =  m_B/3 ; 
            MO( 2 , 2 ) =  m_T/5+m_N+(4*l_B^2*m_B)/l_T^2+(12*R*l_B*m_B)/l_T^2+(12*R^2*m_B)/l_T^2+3*m_B+(4*J_N1)/l_T^2 ; 
            MO( 2 , 3 ) =  3*c*m_B+3*b*m_B ; 
            MO( 2 , 7 ) =  (l_B*m_B)/(2*l_T)+(2*R*m_B)/(3*l_T) ; 
            MO( 2 , 10 ) =  (l_B*m_B)/(2*l_T)+(2*R*m_B)/(3*l_T) ; 
            MO( 2 , 13 ) =  (l_B*m_B)/(2*l_T)+(2*R*m_B)/(3*l_T) ; 
            MO( 3 , 2 ) =  3*c*m_B+3*b*m_B ; 
            MO( 3 , 3 ) =  (l_B^2*m_B)/2+(3*R*l_B*m_B)/2+3*c^2*m_B+6*b*c*m_B+3*b^2*m_B+(3*R^2*m_B)/2+(i_T*l_T)/3+J_N3 ; 
            MO( 4 , 4 ) =  2*m_R+(2*l_B^2*m_B)/l_T^2+(6*R*l_B*m_B)/l_T^2+(6*R^2*m_B)/l_T^2+3*m_B+(4*J_R)/l_T^2 ; 
            MO( 4 , 7 ) =  m_B/3 ; 
            MO( 4 , 9 ) =  (sqrt(3)*l_B*m_B)/(4*l_T)+(R*m_B)/(sqrt(3)*l_T) ; 
            MO( 4 , 10 ) =  -(m_B/6) ; 
            MO( 4 , 11 ) =  -(m_B/(2*sqrt(3))) ; 
            MO( 4 , 12 ) =  -((sqrt(3)*l_B*m_B)/(4*l_T))-(R*m_B)/(sqrt(3)*l_T) ; 
            MO( 4 , 13 ) =  -(m_B/6) ; 
            MO( 4 , 14 ) =  m_B/(2*sqrt(3)) ; 
            MO( 5 , 5 ) =  (2*l_B^2*m_B)/l_T^2+(6*R*l_B*m_B)/l_T^2+(6*R^2*m_B)/l_T^2+3*m_B+(4*J_R)/l_T^2 ; 
            MO( 5 , 6 ) =  -((l_B*m_B)/(2*l_T))-(2*R*m_B)/(3*l_T) ; 
            MO( 5 , 8 ) =  m_B/3 ; 
            MO( 5 , 9 ) =  (l_B*m_B)/(4*l_T)+(R*m_B)/(3*l_T) ; 
            MO( 5 , 10 ) =  m_B/(2*sqrt(3)) ; 
            MO( 5 , 11 ) =  -(m_B/6) ; 
            MO( 5 , 12 ) =  (l_B*m_B)/(4*l_T)+(R*m_B)/(3*l_T) ; 
            MO( 5 , 13 ) =  -(m_B/(2*sqrt(3))) ; 
            MO( 5 , 14 ) =  -(m_B/6) ; 
            MO( 6 , 1 ) =  m_B/3 ; 
            MO( 6 , 5 ) =  -((l_B*m_B)/(2*l_T))-(2*R*m_B)/(3*l_T) ; 
            MO( 6 , 6 ) =  m_B/5 ; 
            MO( 7 , 2 ) =  (l_B*m_B)/(2*l_T)+(2*R*m_B)/(3*l_T) ; 
            MO( 7 , 4 ) =  m_B/3 ; 
            MO( 7 , 7 ) =  m_B/5 ; 
            MO( 8 , 5 ) =  m_B/3 ; 
            MO( 8 , 8 ) =  m_B/5 ; 
            MO( 9 , 1 ) =  m_B/3 ; 
            MO( 9 , 4 ) =  (sqrt(3)*l_B*m_B)/(4*l_T)+(R*m_B)/(sqrt(3)*l_T) ; 
            MO( 9 , 5 ) =  (l_B*m_B)/(4*l_T)+(R*m_B)/(3*l_T) ; 
            MO( 9 , 9 ) =  m_B/5 ; 
            MO( 10 , 2 ) =  (l_B*m_B)/(2*l_T)+(2*R*m_B)/(3*l_T) ; 
            MO( 10 , 4 ) =  -(m_B/6) ; 
            MO( 10 , 5 ) =  m_B/(2*sqrt(3)) ; 
            MO( 10 , 10 ) =  m_B/5 ; 
            MO( 11 , 4 ) =  -(m_B/(2*sqrt(3))) ; 
            MO( 11 , 5 ) =  -(m_B/6) ; 
            MO( 11 , 11 ) =  m_B/5 ; 
            MO( 12 , 1 ) =  m_B/3 ; 
            MO( 12 , 4 ) =  -((sqrt(3)*l_B*m_B)/(4*l_T))-(R*m_B)/(sqrt(3)*l_T) ; 
            MO( 12 , 5 ) =  (l_B*m_B)/(4*l_T)+(R*m_B)/(3*l_T) ; 
            MO( 12 , 12 ) =  m_B/5 ; 
            MO( 13 , 2 ) =  (l_B*m_B)/(2*l_T)+(2*R*m_B)/(3*l_T) ; 
            MO( 13 , 4 ) =  -(m_B/6) ; 
            MO( 13 , 5 ) =  -(m_B/(2*sqrt(3))) ; 
            MO( 13 , 13 ) =  m_B/5; 
            MO( 14 , 4 ) =  m_B/(2*sqrt(3)) ; 
            MO( 14 , 5 ) =  -(m_B/6) ; 
            MO( 14 , 14 ) =  m_B/5 ; 
            MC( 1 , 5 ) =  -((6*c*m_B)/l_T)-(6*b*m_B)/l_T-(2*l_B^2*m_B)/l_T^2-(6*R*l_B*m_B)/l_T^2-(6*R^2*m_B)/l_T^2 ; 
            MC( 1 , 6 ) =  (l_B*m_B)/(2*l_T)+(2*R*m_B)/(3*l_T) ; 
            MC( 1 , 8 ) =  -((2*c*m_B)/(3*l_T))-(2*b*m_B)/(3*l_T) ; 
            MC( 1 , 9 ) =  -((l_B*m_B)/(4*l_T))-(R*m_B)/(3*l_T) ; 
            MC( 1 , 10 ) =  -((c*m_B)/(sqrt(3)*l_T))-(b*m_B)/(sqrt(3)*l_T) ; 
            MC( 1 , 11 ) =  (c*m_B)/(3*l_T)+(b*m_B)/(3*l_T) ; 
            MC( 1 , 12 ) =  -((l_B*m_B)/(4*l_T))-(R*m_B)/(3*l_T) ; 
            MC( 1 , 13 ) =  (c*m_B)/(sqrt(3)*l_T)+(b*m_B)/(sqrt(3)*l_T) ; 
            MC( 1 , 14 ) =  (c*m_B)/(3*l_T)+(b*m_B)/(3*l_T) ; 
            MC( 2 , 4 ) =  3*m_B ; 
            MC( 2 , 7 ) =  m_B/3 ; 
            MC( 2 , 10 ) =  -(m_B/6) ; 
            MC( 2 , 11 ) =  -(m_B/(2*sqrt(3))) ; 
            MC( 2 , 13 ) =  -(m_B/6) ; 
            MC( 2 , 14 ) =  m_B/(2*sqrt(3)) ; 
            MC( 3 , 4 ) =  (l_B^2*m_B)/l_T+(3*R*l_B*m_B)/l_T+(3*R^2*m_B)/l_T+3*c*m_B+3*b*m_B ; 
            MC( 3 , 7 ) =  (c*m_B)/3+(b*m_B)/3 ; 
            MC( 3 , 9 ) =  (sqrt(3)*l_B*m_B)/8+(R*m_B)/(2*sqrt(3)) ; 
            MC( 3 , 10 ) =  -((c*m_B)/6)-(b*m_B)/6 ; 
            MC( 3 , 11 ) =  -((c*m_B)/(2*sqrt(3)))-(b*m_B)/(2*sqrt(3)) ; 
            MC( 3 , 12 ) =  -((sqrt(3)*l_B*m_B)/8)-(R*m_B)/(2*sqrt(3)) ; 
            MC( 3 , 13 ) =  -((c*m_B)/6)-(b*m_B)/6 ; 
            MC( 3 , 14 ) =  (c*m_B)/(2*sqrt(3))+(b*m_B)/(2*sqrt(3)) ; 
            MC( 4 , 2 ) =  3*m_B ; 
            MC( 4 , 3 ) =  (l_B^2*m_B)/l_T+(3*R*l_B*m_B)/l_T+(3*R^2*m_B)/l_T+3*c*m_B+3*b*m_B ; 
            MC( 5 , 1 ) =  -((6*c*m_B)/l_T)-(6*b*m_B)/l_T-(2*l_B^2*m_B)/l_T^2-(6*R*l_B*m_B)/l_T^2-(6*R^2*m_B)/l_T^2 ; 
            MC( 6 , 1 ) =  (l_B*m_B)/(2*l_T)+(2*R*m_B)/(3*l_T) ; 
            MC( 7 , 2 ) =  m_B/3 ; 
            MC( 7 , 3 ) =  (c*m_B)/3+(b*m_B)/3 ; 
            MC( 8 , 1 ) =  -((2*c*m_B)/(3*l_T))-(2*b*m_B)/(3*l_T) ; 
            MC( 9 , 1 ) =  -((l_B*m_B)/(4*l_T))-(R*m_B)/(3*l_T) ; 
            MC( 9 , 3 ) =  (sqrt(3)*l_B*m_B)/8+(R*m_B)/(2*sqrt(3)) ; 
            MC( 10 , 1 ) =  -((c*m_B)/(sqrt(3)*l_T))-(b*m_B)/(sqrt(3)*l_T) ; 
            MC( 10 , 2 ) =  -(m_B/6) ; 
            MC( 10 , 3 ) =  -((c*m_B)/6)-(b*m_B)/6 ; 
            MC( 11 , 1 ) =  (c*m_B)/(3*l_T)+(b*m_B)/(3*l_T) ; 
            MC( 11 , 2 ) =  -(m_B/(2*sqrt(3))) ; 
            MC( 11 , 3 ) =  -((c*m_B)/(2*sqrt(3)))-(b*m_B)/(2*sqrt(3)) ; 
            MC( 12 , 1 ) =  -((l_B*m_B)/(4*l_T))-(R*m_B)/(3*l_T) ; 
            MC( 12 , 3 ) =  -((sqrt(3)*l_B*m_B)/8)-(R*m_B)/(2*sqrt(3)) ; 
            MC( 13 , 1 ) =  (c*m_B)/(sqrt(3)*l_T)+(b*m_B)/(sqrt(3)*l_T) ; 
            MC( 13 , 2 ) =  -(m_B/6) ; 
            MC( 13 , 3 ) =  -((c*m_B)/6)-(b*m_B)/6 ; 
            MC( 14 , 1 ) =  (c*m_B)/(3*l_T)+(b*m_B)/(3*l_T) ; 
            MC( 14 , 2 ) =  m_B/(2*sqrt(3)) ; 
            MC( 14 , 3 ) =  (c*m_B)/(2*sqrt(3))+(b*m_B)/(2*sqrt(3)) ; 
            MS( 1 , 4 ) =  -((6*c*m_B)/l_T)-(6*b*m_B)/l_T-(2*l_B^2*m_B)/l_T^2-(6*R*l_B*m_B)/l_T^2-(6*R^2*m_B)/l_T^2 ; 
            MS( 1 , 7 ) =  -((2*c*m_B)/(3*l_T))-(2*b*m_B)/(3*l_T) ; 
            MS( 1 , 9 ) =  -((sqrt(3)*l_B*m_B)/(4*l_T))-(R*m_B)/(sqrt(3)*l_T) ; 
            MS( 1 , 10 ) =  (c*m_B)/(3*l_T)+(b*m_B)/(3*l_T) ; 
            MS( 1 , 11 ) =  (c*m_B)/(sqrt(3)*l_T)+(b*m_B)/(sqrt(3)*l_T) ; 
            MS( 1 , 12 ) =  (sqrt(3)*l_B*m_B)/(4*l_T)+(R*m_B)/(sqrt(3)*l_T) ; 
            MS( 1 , 13 ) =  (c*m_B)/(3*l_T)+(b*m_B)/(3*l_T) ; 
            MS( 1 , 14 ) =  -((c*m_B)/(sqrt(3)*l_T))-(b*m_B)/(sqrt(3)*l_T) ; 
            MS( 2 , 5 ) =  -(3*m_B) ; 
            MS( 2 , 8 ) =  -(m_B/3) ; 
            MS( 2 , 10 ) =  -(m_B/(2*sqrt(3))) ; 
            MS( 2 , 11 ) =  m_B/6 ; 
            MS( 2 , 13 ) =  m_B/(2*sqrt(3)) ; 
            MS( 2 , 14 ) =  m_B/6 ; 
            MS( 3 , 5 ) =  -((l_B^2*m_B)/l_T)-(3*R*l_B*m_B)/l_T-(3*R^2*m_B)/l_T-3*c*m_B-3*b*m_B ; 
            MS( 3 , 6 ) =  (l_B*m_B)/4+(R*m_B)/3 ; 
            MS( 3 , 8 ) =  -((c*m_B)/3)-(b*m_B)/3 ; 
            MS( 3 , 9 ) =  -((l_B*m_B)/8)-(R*m_B)/6 ; 
            MS( 3 , 10 ) =  -((c*m_B)/(2*sqrt(3)))-(b*m_B)/(2*sqrt(3)) ; 
            MS( 3 , 11 ) =  (c*m_B)/6+(b*m_B)/6 ; 
            MS( 3 , 12 ) =  -((l_B*m_B)/8)-(R*m_B)/6 ; 
            MS( 3 , 13 ) =  (c*m_B)/(2*sqrt(3))+(b*m_B)/(2*sqrt(3)) ; 
            MS( 3 , 14 ) =  (c*m_B)/6+(b*m_B)/6 ; 
            MS( 4 , 1 ) =  -((6*c*m_B)/l_T)-(6*b*m_B)/l_T-(2*l_B^2*m_B)/l_T^2-(6*R*l_B*m_B)/l_T^2-(6*R^2*m_B)/l_T^2 ; 
            MS( 5 , 2 ) =  -(3*m_B) ; 
            MS( 5 , 3 ) =  -((l_B^2*m_B)/l_T)-(3*R*l_B*m_B)/l_T-(3*R^2*m_B)/l_T-3*c*m_B-3*b*m_B ; 
            MS( 6 , 3 ) =  (l_B*m_B)/4+(R*m_B)/3 ; 
            MS( 7 , 1 ) =  -((2*c*m_B)/(3*l_T))-(2*b*m_B)/(3*l_T) ; 
            MS( 8 , 2 ) =  -(m_B/3) ; 
            MS( 8 , 3 ) =  -((c*m_B)/3)-(b*m_B)/3 ; 
            MS( 9 , 1 ) =  -((sqrt(3)*l_B*m_B)/(4*l_T))-(R*m_B)/(sqrt(3)*l_T) ; 
            MS( 9 , 3 ) =  -((l_B*m_B)/8)-(R*m_B)/6 ; 
            MS( 10 , 1 ) =  (c*m_B)/(3*l_T)+(b*m_B)/(3*l_T) ; 
            MS( 10 , 2 ) =  -(m_B/(2*sqrt(3))) ; 
            MS( 10 , 3 ) =  -((c*m_B)/(2*sqrt(3)))-(b*m_B)/(2*sqrt(3)) ; 
            MS( 11 , 1 ) =  (c*m_B)/(sqrt(3)*l_T)+(b*m_B)/(sqrt(3)*l_T) ; 
            MS( 11 , 2 ) =  m_B/6 ; 
            MS( 11 , 3 ) =  (c*m_B)/6+(b*m_B)/6 ; 
            MS( 12 , 1 ) =  (sqrt(3)*l_B*m_B)/(4*l_T)+(R*m_B)/(sqrt(3)*l_T) ; 
            MS( 12 , 3 ) =  -((l_B*m_B)/8)-(R*m_B)/6 ; 
            MS( 13 , 1 ) =  (c*m_B)/(3*l_T)+(b*m_B)/(3*l_T) ; 
            MS( 13 , 2 ) =  m_B/(2*sqrt(3)) ; 
            MS( 13 , 3 ) =  (c*m_B)/(2*sqrt(3))+(b*m_B)/(2*sqrt(3)) ; 
            MS( 14 , 1 ) =  -((c*m_B)/(sqrt(3)*l_T))-(b*m_B)/(sqrt(3)*l_T) ; 
            MS( 14 , 2 ) =  m_B/6 ; 
            MS( 14 , 3 ) =  (c*m_B)/6+(b*m_B)/6 ; 
            GO( 1 , 3 ) =  (2*l_B^2*m_B)/l_T+(6*R*l_B*m_B)/l_T+(6*R^2*m_B)/l_T ; 
            GO( 2 , 8 ) =  -((l_B*m_B)/l_T)-(4*R*m_B)/(3*l_T) ; 
            GO( 2 , 11 ) =  -((l_B*m_B)/l_T)-(4*R*m_B)/(3*l_T) ; 
            GO( 2 , 14 ) =  -((l_B*m_B)/l_T)-(4*R*m_B)/(3*l_T) ; 
            GO( 3 , 1 ) =  -((2*l_B^2*m_B)/l_T)-(6*R*l_B*m_B)/l_T-(6*R^2*m_B)/l_T ; 
            GO( 4 , 5 ) =  -(6*m_B) ; 
            GO( 4 , 8 ) =  -((2*m_B)/3) ; 
            GO( 4 , 10 ) =  -(m_B/sqrt(3)) ; 
            GO( 4 , 11 ) =  m_B/3 ; 
            GO( 4 , 13 ) =  m_B/sqrt(3) ; 
            GO( 4 , 14 ) =  m_B/3 ; 
            GO( 5 , 4 ) =  6*m_B ; 
            GO( 5 , 7 ) =  (2*m_B)/3 ; 
            GO( 5 , 10 ) =  -(m_B/3) ; 
            GO( 5 , 11 ) =  -(m_B/sqrt(3)) ; 
            GO( 5 , 13 ) =  -(m_B/3) ; 
            GO( 5 , 14 ) =  m_B/sqrt(3) ; 
            GO( 7 , 5 ) =  -((2*m_B)/3) ; 
            GO( 7 , 8 ) =  -((2*m_B)/5) ; 
            GO( 8 , 2 ) =  (l_B*m_B)/l_T+(4*R*m_B)/(3*l_T) ; 
            GO( 8 , 4 ) =  (2*m_B)/3 ; 
            GO( 8 , 7 ) =  (2*m_B)/5 ; 
            GO( 10 , 4 ) =  m_B/sqrt(3) ; 
            GO( 10 , 5 ) =  m_B/3 ; 
            GO( 10 , 11 ) =  -((2*m_B)/5) ; 
            GO( 11 , 2 ) =  (l_B*m_B)/l_T+(4*R*m_B)/(3*l_T) ; 
            GO( 11 , 4 ) =  -(m_B/3) ; 
            GO( 11 , 5 ) =  m_B/sqrt(3) ; 
            GO( 11 , 10 ) =  (2*m_B)/5 ; 
            GO( 13 , 4 ) =  -(m_B/sqrt(3)) ; 
            GO( 13 , 5 ) =  m_B/3 ; 
            GO( 13 , 14 ) =  -((2*m_B)/5) ; 
            GO( 14 , 2 ) =  (l_B*m_B)/l_T+(4*R*m_B)/(3*l_T) ; 
            GO( 14 , 4 ) =  -(m_B/3) ; 
            GO( 14 , 5 ) =  -(m_B/sqrt(3)) ; 
            GO( 14 , 13 ) =  (2*m_B)/5 ; 
            GC( 1 , 4 ) =  -((12*c*m_B)/l_T)-(12*b*m_B)/l_T ; 
            GC( 1 , 7 ) =  -((4*c*m_B)/(3*l_T))-(4*b*m_B)/(3*l_T) ; 
            GC( 1 , 10 ) =  (2*c*m_B)/(3*l_T)+(2*b*m_B)/(3*l_T) ; 
            GC( 1 , 11 ) =  (2*c*m_B)/(sqrt(3)*l_T)+(2*b*m_B)/(sqrt(3)*l_T) ; 
            GC( 1 , 13 ) =  (2*c*m_B)/(3*l_T)+(2*b*m_B)/(3*l_T) ; 
            GC( 1 , 14 ) =  -((2*c*m_B)/(sqrt(3)*l_T))-(2*b*m_B)/(sqrt(3)*l_T) ; 
            GC( 2 , 5 ) =  -(6*m_B) ; 
            GC( 2 , 8 ) =  -((2*m_B)/3) ; 
            GC( 2 , 10 ) =  -(m_B/sqrt(3)) ; 
            GC( 2 , 11 ) =  m_B/3 ; 
            GC( 2 , 13 ) =  m_B/sqrt(3) ; 
            GC( 2 , 14 ) =  m_B/3 ; 
            GC( 3 , 5 ) =  -(6*c*m_B)-6*b*m_B ; 
            GC( 3 , 8 ) =  -((2*c*m_B)/3)-(2*b*m_B)/3 ; 
            GC( 3 , 10 ) =  -((c*m_B)/sqrt(3))-(b*m_B)/sqrt(3) ; 
            GC( 3 , 11 ) =  (c*m_B)/3+(b*m_B)/3 ; 
            GC( 3 , 13 ) =  (c*m_B)/sqrt(3)+(b*m_B)/sqrt(3) ; 
            GC( 3 , 14 ) =  (c*m_B)/3+(b*m_B)/3 ; 
            GC( 4 , 1 ) =  -((4*l_B^2*m_B)/l_T^2)-(12*R*l_B*m_B)/l_T^2-(12*R^2*m_B)/l_T^2 ; 
            GC( 5 , 3 ) =  -((2*l_B^2*m_B)/l_T)-(6*R*l_B*m_B)/l_T-(6*R^2*m_B)/l_T ; 
            GC( 6 , 3 ) =  (l_B*m_B)/2+(2*R*m_B)/3 ; 
            GC( 9 , 1 ) =  -((sqrt(3)*l_B*m_B)/(2*l_T))-(2*R*m_B)/(sqrt(3)*l_T) ; 
            GC( 9 , 3 ) =  -((l_B*m_B)/4)-(R*m_B)/3 ; 
            GC( 12 , 1 ) =  (sqrt(3)*l_B*m_B)/(2*l_T)+(2*R*m_B)/(sqrt(3)*l_T) ; 
            GC( 12 , 3 ) =  -((l_B*m_B)/4)-(R*m_B)/3 ; 
            GS( 1 , 5 ) =  (12*c*m_B)/l_T+(12*b*m_B)/l_T ; 
            GS( 1 , 8 ) =  (4*c*m_B)/(3*l_T)+(4*b*m_B)/(3*l_T) ; 
            GS( 1 , 10 ) =  (2*c*m_B)/(sqrt(3)*l_T)+(2*b*m_B)/(sqrt(3)*l_T) ; 
            GS( 1 , 11 ) =  -((2*c*m_B)/(3*l_T))-(2*b*m_B)/(3*l_T) ; 
            GS( 1 , 13 ) =  -((2*c*m_B)/(sqrt(3)*l_T))-(2*b*m_B)/(sqrt(3)*l_T) ; 
            GS( 1 , 14 ) =  -((2*c*m_B)/(3*l_T))-(2*b*m_B)/(3*l_T) ; 
            GS( 2 , 4 ) =  -(6*m_B) ; 
            GS( 2 , 7 ) =  -((2*m_B)/3) ; 
            GS( 2 , 10 ) =  m_B/3 ; 
            GS( 2 , 11 ) =  m_B/sqrt(3) ; 
            GS( 2 , 13 ) =  m_B/3 ; 
            GS( 2 , 14 ) =  -(m_B/sqrt(3)) ; 
            GS( 3 , 4 ) =  -(6*c*m_B)-6*b*m_B ; 
            GS( 3 , 7 ) =  -((2*c*m_B)/3)-(2*b*m_B)/3 ; 
            GS( 3 , 10 ) =  (c*m_B)/3+(b*m_B)/3 ; 
            GS( 3 , 11 ) =  (c*m_B)/sqrt(3)+(b*m_B)/sqrt(3) ; 
            GS( 3 , 13 ) =  (c*m_B)/3+(b*m_B)/3 ; 
            GS( 3 , 14 ) =  -((c*m_B)/sqrt(3))-(b*m_B)/sqrt(3) ; 
            GS( 4 , 3 ) =  -((2*l_B^2*m_B)/l_T)-(6*R*l_B*m_B)/l_T-(6*R^2*m_B)/l_T ; 
            GS( 5 , 1 ) =  (4*l_B^2*m_B)/l_T^2+(12*R*l_B*m_B)/l_T^2+(12*R^2*m_B)/l_T^2 ; 
            GS( 6 , 1 ) =  -((l_B*m_B)/l_T)-(4*R*m_B)/(3*l_T) ; 
            GS( 9 , 1 ) =  (l_B*m_B)/(2*l_T)+(2*R*m_B)/(3*l_T) ; 
            GS( 9 , 3 ) =  -((sqrt(3)*l_B*m_B)/4)-(R*m_B)/sqrt(3) ; 
            GS( 12 , 1 ) =  (l_B*m_B)/(2*l_T)+(2*R*m_B)/(3*l_T) ; 
            GS( 12 , 3 ) =  (sqrt(3)*l_B*m_B)/4+(R*m_B)/sqrt(3) ; 
            KK( 1 , 1 ) =  (4*EI_T)/l_T^3 ; 
            KK( 2 , 2 ) =  (4*EI_T)/l_T^3 ; 
            KK( 3 , 3 ) =  GI_T/l_T ; 
            KK( 4 , 4 ) =  (4*EI_S)/c^3 ; 
            KK( 5 , 5 ) =  (4*EI_S)/c^3 ; 
            KK( 6 , 6 ) =  (4*EI_B)/l_B^3 ; 
            KK( 7 , 7 ) =  (4*EI_B)/l_B^3 ; 
            KK( 8 , 8 ) =  (4*EA_B)/(3*l_B) ; 
            KK( 9 , 9 ) =  (4*EI_B)/l_B^3 ; 
            KK( 10 , 10 ) =  (4*EI_B)/l_B^3 ; 
            KK( 11 , 11 ) =  (4*EA_B)/(3*l_B) ; 
            KK( 12 , 12 ) =  (4*EI_B)/l_B^3 ; 
            KK( 13 , 13 ) =  (4*EI_B)/l_B^3 ; 
            KK( 14 , 14 ) =  (4*EA_B)/(3*l_B) ; 
            KO( 1 , 1 ) =  -((2*l_B^2*m_B)/l_T^2)-(6*R*l_B*m_B)/l_T^2-(6*R^2*m_B)/l_T^2 ; 
            KO( 2 , 2 ) =  -((4*l_B^2*m_B)/l_T^2)-(12*R*l_B*m_B)/l_T^2-(12*R^2*m_B)/l_T^2 ;
            KO( 2 , 7 ) =  -((l_B*m_B)/(2*l_T))-(2*R*m_B)/(3*l_T) ;
            KO( 2 , 10 ) =  -((l_B*m_B)/(2*l_T))-(2*R*m_B)/(3*l_T) ;
            KO( 2 , 13 ) =  -((l_B*m_B)/(2*l_T))-(2*R*m_B)/(3*l_T) ;
            KO( 3 , 3 ) =  -((l_B^2*m_B)/2)-(3*R*l_B*m_B)/2-(3*R^2*m_B)/2 ;
            KO( 4 , 4 ) =  -(3*m_B) ;
            KO( 4 , 7 ) =  -(m_B/3) ;
            KO( 4 , 10 ) =  m_B/6 ;
            KO( 4 , 11 ) =  m_B/(2*sqrt(3)) ;
            KO( 4 , 13 ) =  m_B/6 ;
            KO( 4 , 14 ) =  -(m_B/(2*sqrt(3))) ;
            KO( 5 , 5 ) =  -(3*m_B) ;
            KO( 5 , 8 ) =  -(m_B/3) ;
            KO( 5 , 10 ) =  -(m_B/(2*sqrt(3))) ;
            KO( 5 , 11 ) =  m_B/6 ;
            KO( 5 , 13 ) =  m_B/(2*sqrt(3)) ;
            KO( 5 , 14 ) =  m_B/6 ;
            KO( 6 , 6 ) =  (R*m_B)/(3*l_B)+(4*m_B)/15 ;
            KO( 7 , 2 ) =  -((l_B*m_B)/(2*l_T))-(2*R*m_B)/(3*l_T) ;
            KO( 7 , 4 ) =  -(m_B/3) ;
            KO( 7 , 7 ) =  (R*m_B)/(3*l_B)+m_B/15 ;
            KO( 8 , 5 ) =  -(m_B/3) ;
            KO( 8 , 8 ) =  -(m_B/5) ;
            KO( 9 , 9 ) =  (R*m_B)/(3*l_B)+(4*m_B)/15 ;
            KO( 10 , 2 ) =  -((l_B*m_B)/(2*l_T))-(2*R*m_B)/(3*l_T) ;
            KO( 10 , 4 ) =  m_B/6 ;
            KO( 10 , 5 ) =  -(m_B/(2*sqrt(3))) ;
            KO( 10 , 10 ) =  (R*m_B)/(3*l_B)+m_B/15 ;
            KO( 11 , 4 ) =  m_B/(2*sqrt(3)) ;
            KO( 11 , 5 ) =  m_B/6 ;
            KO( 11 , 11 ) =  -(m_B/5) ;
            KO( 12 , 12 ) =  (R*m_B)/(3*l_B)+(4*m_B)/15 ;
            KO( 13 , 2 ) =  -((l_B*m_B)/(2*l_T))-(2*R*m_B)/(3*l_T) ;
            KO( 13 , 4 ) =  m_B/6 ;
            KO( 13 , 5 ) =  m_B/(2*sqrt(3)) ;
            KO( 13 , 13 ) =  (R*m_B)/(3*l_B)+m_B/15 ;
            KO( 14 , 4 ) =  -(m_B/(2*sqrt(3))) ;
            KO( 14 , 5 ) =  m_B/6 ;
            KO( 14 , 14 ) =  -(m_B/5) ;
            KC( 1 , 5 ) =  (6*c*m_B)/l_T+(6*b*m_B)/l_T ; 
            KC( 1 , 8 ) =  (2*c*m_B)/(3*l_T)+(2*b*m_B)/(3*l_T) ; 
            KC( 1 , 10 ) =  (c*m_B)/(sqrt(3)*l_T)+(b*m_B)/(sqrt(3)*l_T) ; 
            KC( 1 , 11 ) =  -((c*m_B)/(3*l_T))-(b*m_B)/(3*l_T) ; 
            KC( 1 , 13 ) =  -((c*m_B)/(sqrt(3)*l_T))-(b*m_B)/(sqrt(3)*l_T) ; 
            KC( 1 , 14 ) =  -((c*m_B)/(3*l_T))-(b*m_B)/(3*l_T) ; 
            KC( 2 , 4 ) =  -(3*m_B) ; 
            KC( 2 , 7 ) =  -(m_B/3) ; 
            KC( 2 , 10 ) =  m_B/6 ; 
            KC( 2 , 11 ) =  m_B/(2*sqrt(3)) ; 
            KC( 2 , 13 ) =  m_B/6 ; 
            KC( 2 , 14 ) =  -(m_B/(2*sqrt(3))) ; 
            KC( 3 , 4 ) =  -(3*c*m_B)-3*b*m_B ; 
            KC( 3 , 7 ) =  -((c*m_B)/3)-(b*m_B)/3 ; 
            KC( 3 , 10 ) =  (c*m_B)/6+(b*m_B)/6 ; 
            KC( 3 , 11 ) =  (c*m_B)/(2*sqrt(3))+(b*m_B)/(2*sqrt(3)) ; 
            KC( 3 , 13 ) =  (c*m_B)/6+(b*m_B)/6 ; 
            KC( 3 , 14 ) =  -((c*m_B)/(2*sqrt(3)))-(b*m_B)/(2*sqrt(3)) ; 
            KC( 4 , 3 ) =  -((l_B^2*m_B)/l_T)-(3*R*l_B*m_B)/l_T-(3*R^2*m_B)/l_T ; 
            KC( 5 , 1 ) =  (2*l_B^2*m_B)/l_T^2+(6*R*l_B*m_B)/l_T^2+(6*R^2*m_B)/l_T^2 ; 
            KC( 6 , 1 ) =  -((l_B*m_B)/(2*l_T))-(2*R*m_B)/(3*l_T) ; 
            KC( 9 , 1 ) =  (l_B*m_B)/(4*l_T)+(R*m_B)/(3*l_T) ; 
            KC( 9 , 3 ) =  -((sqrt(3)*l_B*m_B)/8)-(R*m_B)/(2*sqrt(3)) ; 
            KC( 12 , 1 ) =  (l_B*m_B)/(4*l_T)+(R*m_B)/(3*l_T) ; 
            KC( 12 , 3 ) =  (sqrt(3)*l_B*m_B)/8+(R*m_B)/(2*sqrt(3)) ; 
            KS( 1 , 4 ) =  (6*c*m_B)/l_T+(6*b*m_B)/l_T ; 
            KS( 1 , 7 ) =  (2*c*m_B)/(3*l_T)+(2*b*m_B)/(3*l_T) ; 
            KS( 1 , 10 ) =  -((c*m_B)/(3*l_T))-(b*m_B)/(3*l_T) ; 
            KS( 1 , 11 ) =  -((c*m_B)/(sqrt(3)*l_T))-(b*m_B)/(sqrt(3)*l_T) ; 
            KS( 1 , 13 ) =  -((c*m_B)/(3*l_T))-(b*m_B)/(3*l_T) ; 
            KS( 1 , 14 ) =  (c*m_B)/(sqrt(3)*l_T)+(b*m_B)/(sqrt(3)*l_T) ; 
            KS( 2 , 5 ) =  3*m_B ; 
            KS( 2 , 8 ) =  m_B/3 ; 
            KS( 2 , 10 ) =  m_B/(2*sqrt(3)) ; 
            KS( 2 , 11 ) =  -(m_B/6) ; 
            KS( 2 , 13 ) =  -(m_B/(2*sqrt(3))) ; 
            KS( 2 , 14 ) =  -(m_B/6) ; 
            KS( 3 , 5 ) =  3*c*m_B+3*b*m_B ; 
            KS( 3 , 8 ) =  (c*m_B)/3+(b*m_B)/3 ; 
            KS( 3 , 10 ) =  (c*m_B)/(2*sqrt(3))+(b*m_B)/(2*sqrt(3)) ; 
            KS( 3 , 11 ) =  -((c*m_B)/6)-(b*m_B)/6 ; 
            KS( 3 , 13 ) =  -((c*m_B)/(2*sqrt(3)))-(b*m_B)/(2*sqrt(3)) ; 
            KS( 3 , 14 ) =  -((c*m_B)/6)-(b*m_B)/6 ; 
            KS( 4 , 1 ) =  (2*l_B^2*m_B)/l_T^2+(6*R*l_B*m_B)/l_T^2+(6*R^2*m_B)/l_T^2 ; 
            KS( 5 , 3 ) =  (l_B^2*m_B)/l_T+(3*R*l_B*m_B)/l_T+(3*R^2*m_B)/l_T ; 
            KS( 6 , 3 ) =  -((l_B*m_B)/4)-(R*m_B)/3 ; 
            KS( 9 , 1 ) =  (sqrt(3)*l_B*m_B)/(4*l_T)+(R*m_B)/(sqrt(3)*l_T) ; 
            KS( 9 , 3 ) =  (l_B*m_B)/8+(R*m_B)/6 ; 
            KS( 12 , 1 ) =  -((sqrt(3)*l_B*m_B)/(4*l_T))-(R*m_B)/(sqrt(3)*l_T) ; 
            KS( 12 , 3 ) =  (l_B*m_B)/8+(R*m_B)/6 ; 

            % end: from matrices.m

            % pivotal coefficients of system matrices are too large -> normalize with MO(i,i)
            %for row = 1:14
                %MO(row,:) = MO(row,:)/MO(row,row);
                %MC(row,:) = MC(row,:)/MO(row,row);
                %MS(row,:) = MS(row,:)/MO(row,row);
                %GO(row,:) = GO(row,:)/MO(row,row);
                %GC(row,:) = GC(row,:)/MO(row,row);
                %GS(row,:) = GS(row,:)/MO(row,row);
                %KK(row,:) = KK(row,:)/MO(row,row);
                %O(row,:) = KO(row,:)/MO(row,row);
                %KC(row,:) = KC(row,:)/MO(row,row);
                %KS(row,:) = KS(row,:)/MO(row,row);
            %end

            obj.MO = MO;
            obj.MC = MC;
            obj.MS = MS;

            obj.GO = GO;
            obj.GC = GC;
            obj.GS = GS;

            obj.KK = KK;
            obj.KO = KO;
            obj.KC = KC;
            obj.KS = KS;

        end

    end
end

