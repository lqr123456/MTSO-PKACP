
%% 10�ֻ���ӳ������ѡ��1-10�ֱ�Ϊ��tent��Logistic��Cubic��chebyshev��Piecewise��sinusoidal��Sine,ICMIC, Circle,Bernoulli
function result = ch(index,N)%��������Ⱥ������ά��

switch index
    case 1
        % tentӳ��
        tent=0.4;  %tent����ϵ��
        Tent(1)=rand;
        for i=1:N-1
            
                if Tent(i)<tent
                    Tent(i+1)=Tent(i)/tent;
                elseif Tent(i)>=tent
                    Tent(i+1)=(1-Tent(i))/(1-tent);
                end
            
        end
        result = Tent;
        
    case 2
        % Logisticӳ��
        miu=4;  %����ϵ��
        Logistic(1)=rand;
        for i=1:N-1
            
                Logistic(i+1)=miu.* Logistic(i).*(1-Logistic(i));
            
        end
        result = Logistic;
        
    case 3
        % Cubicӳ��
        cubic=2.59;%����ϵ��
        Cubic(1)=rand;
        for i=1:N-1
            
                Cubic(i+1)=cubic.*Cubic(i).*(1-Cubic(i).^2);
            
        end
        result = Cubic;
        
    case 4
        %chebyshevӳ��
        chebyshev=5;%����ϵ��
        Chebyshev(1)=rand;
        for i=1:N-1
            
                Chebyshev(i+1)=cos(chebyshev.*acos(Chebyshev(i)));
            
        end
        result = Chebyshev;
        
    case 5
        %Piecewiseӳ��
        p=0.4;
        Piecewise(1)=rand;
        for i=1:N-1
            
                if Piecewise(i)>0&&Piecewise(i)<p
                    Piecewise(i+1)=Piecewise(i)/p;
                elseif Piecewise(i)>=p&&Piecewise(i)<0.5
                    Piecewise(i+1)=(Piecewise(i)-p)/(0.5-p);
                elseif Piecewise(i)>=0.5&&Piecewise(i)<1-p
                    Piecewise(i+1)=(1-p-Piecewise(i))/(0.5-p);
                elseif Piecewise(i)>=1-p&&Piecewise(i)<1
                    Piecewise(i+1)=(1-Piecewise(i))/p;
                
            end
        end
        result = Piecewise;
        
        
    case 6
        %sinusoidalӳ��
        sinusoidal=2.3;
        Sinusoidal(1)=rand;
        for i=1:N-1
            
                Sinusoidal(i+1)=sinusoidal*Sinusoidal(i).^2*(sin(pi*Sinusoidal(i)));
           
        end
        result = Sinusoidal;
        
    case 7
        %Sineӳ��
        sine=4;
        Sine(1)=rand;
        for i=1:N-1
            
                Sine(i+1)=(sine/4)*sin(pi*Sine(i));
            
        end
        result = Sine;
        
        
    case 8
        %         ICMIC ӳ�䣨�����۵���������ӳ�䣩
        icmic=2;
        ICMIC(1)=rand;
        for i=1:N-1
            
                ICMIC(i+1)=sin(icmic/ICMIC(i));
            
        end
        result = ICMIC;
        
        
    case 9
        % Circleӳ��
        a = 0.5; b=2.2;
        Circle(1)=rand;
        for i=1:N-1
            
                Circle(i+1)=mod(Circle(i)+a-b/(2*pi)*sin(2*pi*Circle(i)),1);
            
        end
        result = Circle;
    case 10
        %Bernoulliӳ��
        lammda = 0.4;
        Bernoulli(1)=rand;
        for i=1:N
            
                if Bernoulli(i) <  1-lammda
                    Bernoulli(i+1)= Bernoulli(i)/(1-lammda);
                else
                    Bernoulli(i+1)= (Bernoulli(i)-1+lammda)/lammda;
                end
            
        end
        result = Bernoulli;
        
    case 11
        %Gaussianӳ��
        Gaussian(1)=rand;
        for i=1:N
            
                
                Gaussian(i+1)=1/Gaussian(i)-floor(1/Gaussian(i));
            
        end
        result=Gaussian;
        
    case 12
        %Singerӳ��
        u=1.073;
        Singer(1)=rand;
        for i=1:N
            
                Singer(i+1)=u*(7.86*Singer(i)-23.31*Singer(i)^2+28.75*Singer(i)^3-13.302875*Singer(i)^4);
          
        end
        result=Singer;
                
        
        
end
end

