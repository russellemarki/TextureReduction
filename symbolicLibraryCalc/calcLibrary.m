function[almnsp,blmn,cLin] = calcLibrary(maxL, writechk)
if nargin<2
    writechk = 0;
end

if writechk>0
    fprintf('Calculating Qlmn\n')
end
%Calculates Q Coefficients Bunge 14.3
Qlmn = calcQlmn(maxL, writechk);


if writechk>0
    fprintf('Calculating Almnsp\n')
end
%Calculates Almnsp Coefficients Bunge 14.(72,81:84)
almnsp = calcAlmnsp(maxL, Qlmn, writechk);%[];%

if writechk>0
    fprintf('Calculating Blmn\n')
end
Qlmn = double(Qlmn);
%Calculates Blmn through a method devised by Esling
[blmn, cLin] = getCubicBs(maxL, Qlmn);
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function[Q] = calcQlmn(maxL, writechk)
%Copied From Bunge


Q = sym(zeros(alms_ind(maxL, maxL, maxL), 1));


ind = zeros(3, length(Q));
cnt = 0;
for l = 0:2:maxL
    for m = 0:l
        for n = 0:l
            cnt = cnt + 1;
            ind(1, cnt) = l;
            ind(2, cnt) = m;
            ind(3, cnt) = n;
        end
    end
end

cnt - alms_ind(maxL, maxL, maxL)

parfor i = 1:alms_ind(maxL, maxL, maxL)

    l = ind(1, i);
    m = ind(2, i);
    n = ind(3, i);

    x = sym('x');

    fun = ((1-x)^(l-m))*((1+x)^(l+m));
    for j = 1:(l-n)
        fun = diff(fun);
    end

    mult = sym(1i)^(n-m)*sym(1i)^(n+m);
    %Bunge 14.3
    locQ = mult*(-1)^(l-m)/(2^l*symfac(l-m))*...
        sqrt(symfac(l-m)*symfac(l+n)/(symfac(l+m)*symfac(l-n)))...
        *(1-x)^((m-n)/2)*(1+x)^(-(n+m)/2)*fun;

    Q(i) = (simplify(subs(locQ, x, 0)));
end




end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function[almnsp] = calcAlmnsp(maxL, Qlmn, writechk)
%Copied From Bunge

almnsp = (zeros(almnsp_ind(maxL, maxL, maxL, maxL),1));

ind = zeros(4, length(almnsp));
cnt = 0;
for l = 0:2:maxL
    for m = 0:2:l
        for n = 0:l
            for s = 0:l
                cnt = cnt + 1;
                ind(1, cnt) = l;
                ind(2, cnt) = m;
                ind(3, cnt) = n;
                ind(4, cnt) = s;

                %err(cnt) = cnt - almnsp_ind(l, m, n, s);
            end
        end
    end
end

cnt - almnsp_ind(maxL, maxL, maxL, maxL)

parfor i = 1:size(ind, 2)
    
    
    l = ind(1, i);
    m = ind(2, i);
    n = ind(3, i);
    s = ind(4, i);
    %Bunge 14.72
    if mod(m+n, 2)==0
        if s==0
            %Bunge 14.81
            loc_a = Qlmn(alms_ind(l, m, s))*Qlmn(alms_ind(l, n, s));
        else
            %Bunge 14.82
            loc_a = 2*Qlmn(alms_ind(l, m, s))*Qlmn(alms_ind(l, n, s));
        end
    else
        if s==0
            %Bunge 14.83
            loc_a = sym(0);
        else
            %Bunge 14.84
            %We do not multiply by 1i because it is already
            %taken care of.  *1i
            loc_a = 2*Qlmn(alms_ind(l, m, s))*Qlmn(alms_ind(l, n, s));%2*Qlmn(l,m,s)*Qlmn(l,n,s);
        end
    end

    loc_a = simplify(loc_a);

    almnsp(i) = double(loc_a);
end


end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function[blmn, mLin] = getCubicBs(maxL, QLMN)
%Directly Copied From F. WAGNER, C. ESLING AND R.BARO 1977
%Bunge Library
IDIMH = floor(maxL/4 + 1);

DPROJ = (zeros(IDIMH,IDIMH));
DREDUC = (zeros(IDIMH,IDIMH));
blmn = (zeros(2, 2, (maxL-4)/2));

mLin = zeros(1, maxL + 1);

for L = 4:2:maxL
    IDIMH = floor(L/4 + 1);
    LPS1 = L + 1;
    for I = 0:4:LPS1
        II = I/4 + 1;
        for J = 0:4:I
            JJ = J/4 + 1;
            if false
                if J~=0
                    if I==J
                        DPROJ(II,JJ) = (1 + 4*QLMN(alms_ind(L,J,I)))/3;
                    else
                        DPROJ(II,JJ) = (4*QLMN(alms_ind(L,J,I)))/3;
                    end
                else
                    if I==0
                        DPROJ(II,JJ) = (1 + 2*QLMN(alms_ind(L,J,I)))/3;
                    else
                        DPROJ(II,JJ) = (2*sqrt(2)*QLMN(alms_ind(L,J,I)))/3;
                    end
                end
            else
                if J~=0
                    if I==J
                        DPROJ(II,JJ) = (1 + 4*QLMN(alms_ind(L,I,J)))/3;
                    else
                        DPROJ(II,JJ) = (4*QLMN(alms_ind(L,I,J)))/3;
                    end
                else
                    if I==0
                        DPROJ(II,JJ) = (1 + 2*QLMN(alms_ind(L,I,J)))/3;
                    else
                        DPROJ(II,JJ) = (2*sqrt(2)*QLMN(alms_ind(L,I,J)))/3;
                    end
                end
            end
            DPROJ(JJ,II) = DPROJ(II,JJ);
        end
    end
    
%     DPROJ = simplify(DPROJ);
    DDIMS = (0);
    for I = 1:IDIMH
        DDIMS = DDIMS + DPROJ(I,I);
    end
    
    IDIMS = floor(DDIMS + 0.5);
    IMINS = IDIMH + 1 - IDIMS;
    I = IDIMH + 1;
    IMINCHK = true;
    
    while IMINCHK || I-IMINS>0
        
        IMINCHK = false;
        
        I = I - 1;
        
        while I>1 && (DPROJ(I,I)<1e-8)%isAlways
            I = I - 1;
        end
        
        if I>0
            DNORM = sqrt(DPROJ(I,I));
            for J = 1:IDIMH
                DREDUC(I,J) = DPROJ(I,J)/DNORM;
                DPROJ(I,J) = DPROJ(I,J)/DPROJ(I,I);
            end
            
            
            for ILINE = 1:I-1
                for JCOLON = 1:I
                    DPROJ(ILINE,JCOLON) = DPROJ(ILINE,JCOLON) - DPROJ(I,JCOLON)*DPROJ(ILINE,I);
                end
            end
%             DPROJ = simplify(DPROJ);
        else
            I = IMINS;
        end
    end
    
%     DREDUC = simplify(DREDUC);
    for I = IMINS:IDIMH
        II = IDIMH-I+1;
        DREDUC(I,1) = DREDUC(I,1)/sqrt((pi)*2);%sym
        blmn(1, II, (L-4)/2 + 1) = DREDUC(I,1);
        for J = 2:IDIMH
            DREDUC(I,J) = DREDUC(I,J)/sqrt((pi));%sym
            blmn(J,II,(L-4)/2 + 1) = DREDUC(I,J);
        end
        
    end

    mLin(L+1) = IDIMH - IMINS + 1;
end

% blmn = simplify(blmn);
end

function[fac] = symfac(N)
fac = 1;
if N>18
    fac = sym(fac);
end

for mult = 2:N
    fac = fac*mult;
end

fac = sym(fac);
end