function [W,failure] = W(Z)

% W function [W,failure] = W(Z)
%
% This function computes the complex scaled complementary error function
% also known as the FADDEEVA function
% W(Z) = exp(-Z^2) * (1-erf(-i*Z)).
% Note that W(i*X) == erfcx(X) for real-valued X. erfcx is a build-in MatLab function.
%
% Inputs and outputs are complex matrices of any size.
% The failure flag indicates numerical overflow.
%
% Source: http://www.netlib.org/toms/680.
% Work published in TRANSACTIONS ON MATHEMATICAL SOFTWARE, VOL. 16, NO. 1, PP. 47.
%
% Translated from F77 to MatLab by Thomas Winiecki 10/05/2006.

% 
% some constants to prevent numerical overflows
RMAX     = 1e308;
RMAXREAL = sqrt(RMAX);
RMAXEXP  = log(RMAX) - log(2);
RMAXGONI = 1e10*pi;
FACTOR   = 2/sqrt(pi);

% split input into real/imag
XI = real(Z);
YI = imag(Z);

% reset outputs
sizein  = size(Z);
failure = repmat(0  ,sizein);
W       = repmat(NaN,sizein);

for idx = 1:prod(sizein)
    
    XABS = abs(XI(idx));
    YABS = abs(YI(idx));
    X    = XABS/6.3;
    Y    = YABS/4.4;
    
    if ((XABS > RMAXREAL)|(YABS > RMAXREAL))
        failure(idx) = 1;
        break;
    else
        QRHO  = X^2 + Y^2;
        XABSQ = XABS^2;
        XQUAD = XABSQ - YABS^2;
        YQUAD = 2*XABS*YABS;
        A     = QRHO < 0.085264;
        if A
            QRHO  = (1-0.85*Y)*sqrt(QRHO);
            N     = round(6 + 72*QRHO);
            J     = 2*N+1;
            XSUM  = 1.0/J;
            YSUM  = 0.0;
            for I = N:-1:1
                J    = J - 2;
                XAUX = (XSUM*XQUAD - YSUM*YQUAD)/I;
                YSUM = (XSUM*YQUAD + YSUM*XQUAD)/I;
                XSUM = XAUX + 1.0/J;
            end
            U1   = -FACTOR*(XSUM*YABS + YSUM*XABS) + 1.0;
            V1   =  FACTOR*(XSUM*XABS - YSUM*YABS);
            DAUX =  exp(-XQUAD);
            U2   =  DAUX*cos(YQUAD);
            V2   = -DAUX*sin(YQUAD);
            U    = U1*U2 - V1*V2;
            V    = U1*V2 + V1*U2;
        else
            if QRHO > 1.0
                H    = 0.0;
                KAPN = 0;
                QRHO = sqrt(QRHO);
                NU   = floor(3 + (1442/(26*QRHO+77)));
            else
                QRHO = (1-Y)*sqrt(1-QRHO);
                H    = 1.88*QRHO;
                H2   = 2*H;
                KAPN = round(7  + 34*QRHO);
                NU   = round(16 + 26*QRHO);
            end
            B = H > 0;
            if B 
                QLAMBDA = H2^KAPN;
            end
            RX = 0.0;
            RY = 0.0;
            SX = 0.0;
            SY = 0.0;
            for N = NU:-1:0
                NP1 = N + 1;
                TX  = YABS + H + NP1*RX;
                TY  = XABS - NP1*RY;
                C   = 0.5/(TX^2 + TY^2);
                RX  = C*TX;
                RY  = C*TY;
                if B & (N <= KAPN)
                    TX = QLAMBDA + SX;
                    SX = RX*TX - RY*SY;
                    SY = RY*TX + RX*SY;
                    QLAMBDA = QLAMBDA/H2;
                end
            end
            if H == 0.0
                U = FACTOR*RX;
                V = FACTOR*RY;
            else
                U = FACTOR*SX;
                V = FACTOR*SY;
            end
            if YABS == 0 
                U = exp(-XABS^2);
            end
        end
        if (YI(idx) < 0.0)
            if A
                U2    = 2*U2;
                V2    = 2*V2;
            else
                XQUAD =  -XQUAD;
                if (YQUAD > RMAXGONI)|(XQUAD > RMAXEXP)
                    failure(idx) = 1;
                    break
                end
                W1 =  2*exp(XQUAD);
                U2 =  W1*cos(YQUAD);
                V2 = -W1*sin(YQUAD);
            end
            U = U2 - U;
            V = V2 - V;
            if XI(idx) > 0 
                V = -V;
            end
        else
            if XI(idx) < 0 
                V = -V;
            end
        end
    end
    
    W(idx) = U + i*V;
    
end

return
