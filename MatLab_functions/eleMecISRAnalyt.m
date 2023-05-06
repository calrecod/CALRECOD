function eleMec=eleMecISRAnalyt(c,a,fdpc,h,b,rec,E,t,dUno,dDos)

%------------------------------------------------------------------------
% Syntax:
% eleMec=eleMecISRAnalyt(c,a,fdpc,h,b,rec,E,t,dUno,dDos)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: Any.
%
%------------------------------------------------------------------------
% PURPOSE: To determine the axial load and bending resistance of a 
% rectangular reinforced column cross section with an ISR for a correspon-
% ding neutral axis value c. An analytical mathematical approach is used
% for the ISR analysis.
% 
% OUTPUT: eleMec:             array containing the sum of resistance forces
%                             (axial and being) as [ sum Fs,   sum Ms;
%                                                   sum Fconc, sum Mconc]
%
% INPUT:  b,h:                cross-section dimensions of column (width 
%                             and height)
%
%         t:                  is the column ISR width (considered as 
%                             constant)
%
%         c:                  is the neutral axis depth value
%
%         a:                  is the reduced netural axis depth value,
%                             equal to a=beta1*c
%
%         rec:                is the concrete cover
%
%         fdpc:               is the reduced value of f'c with the factor 
%                             0.85 as prescribed in the ACI 318-19 code:
%                             fdpc=0.85*fc
%
%         dUno:               is the distance from the upper cross-section
%                             border to the centroid of the upper boundary
%                             ISR
%
%         dDos:               is the distance from the upper cross-section
%                             border to the centroid of the lower boundary
%                             ISR
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2023-03-27
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------

%% Main function:
if (c>=-inf && c<=30/51*dUno)
    fa=0.0021*E*t*(b-2*rec);
    fb=0.0021*E*t*(b-2*rec);
    elemAc=casoUnoRec(E,t,dUno,dDos,fa,fb);
elseif(c>30/51*dUno && c<=dUno)
    hTres=21/30*c;
    eps=0.003*(dUno/c-1);
    epsa=0.003*(dUno/c-1);
    if (epsa>0.0021)
        epsa=0.0021;
    end
    epsb=0.0021;
    fa=epsa*E*t*(b-2*rec);
    fb=epsb*E*t*(b-2*rec);
    elemAc=casoDosRec(c,E,t,dUno,dDos,hTres,eps,h,fa,fb);
elseif(c>dUno && c<=30/9*dUno)
    if ((60/102*t)>=(9/51*h-60/51*rec))
        if (c>30/51*dDos && c<=30/9*dUno)
            epsa=(1-dUno/c)*0.003;
            if (epsa>0.0021)
                epsa=0.0021;
            end
            epsb=0.003*(dDos/c-1);
            if (epsb>0.0021)
                epsb=0.0021;
            end
            fa=epsa*E*t*(b-2*rec);
            fb=epsb*E*t*(b-2*rec);
            elemAc=casoEspecial(c,E,t,dUno,dDos,h,fa,fb,epsa,epsb);
        else
            hDos=21/30*c;
            epsa=(1-dUno/c)*0.003;
            if (epsa>0.0021)
                epsa=0.0021;
            end
            epsb=0.0021;
            fa=epsa*E*t*(b-2*rec);
            fb=epsb*E*t*(b-2*rec);
            elemAc=casoTres(c,E,t,dUno,dDos,h,fa,fb,epsa,hDos);
        end
    else
        hDos=21/30*c;
        epsa=(1-dUno/c)*0.003;
        if (epsa>0.0021) 
        epsa=0.0021;
        end
        epsb=0.0021;
        fa=epsa*E*t*(b-2*rec);
        fb=epsb*E*t*(b-2*rec);
        elemAc=casoTres(c,E,t,dUno,dDos,h,fa,fb,epsa,hDos);
    end
    
elseif(c>30/9*dUno && c<=30/51*dDos)
        fa=0.0021*E*t*(b-2*rec);
        fb=0.0021*E*t*(b-2*rec);
        elemAc=casoTresCuatro(c,E,t,dUno,dDos,h,fa,fb);
elseif(c>30/51*dDos && c<=dDos)
        hCuatro=(1-21/30-dUno/c)*c;
        epsa=0.0021;
        epsb=0.003*(dDos/c-1);
        if (epsb>0.0021)
            epsb=0.0021;
        end
        fa=epsa*E*t*(b-2*rec);
        fb=epsb*E*t*(b-2*rec);
        elemAc=casoCuatro(c,E,t,dUno,dDos,h,fa,fb,epsb,hCuatro);
elseif(c>dDos && c<=(30/9*dDos))  
        hCinco=9/30*c;
        epsTres=0.003/c*(c-dDos);
        epsa=0.0021;
        epsb=epsTres;
        if (epsb>0.0021)
            epsb=0.0021;
        end
        fa=epsa*E*t*(b-2*rec);
        fb=epsb*E*t*(b-2*rec);
        elemAc=casoCinco(E,t,dUno,dDos,h,fa,fb,hCinco,epsb);
elseif(c>30/9*dDos)
        fa=0.0021*E*t*(b-2*rec);
        fb=0.0021*E*t*(b-2*rec);
        elemAc=casoSeis(dDos,dUno,E,t,fa,fb);
end

elemConc=casoConcreto(a,fdpc,b,h);

eleMec=[elemAc;
        elemConc];

end

%% Sub-functions:
% (See documentation)
function elemAc=casoUnoRec(E,t,dUno,dDos,fa,fb)

fr=(0.0021*(dDos-dUno)*E*2*t+fa+fb);
mr=0;
elemAc=[fr mr];
end

function elemAc=casoDosRec(c,E,t,dUno,dDos,hTres,eps,h,fa,fb)

funo=2*t*eps*(hTres-(dUno-c))*E;
fdos=2*t*0.5*(0.0021-eps)*(hTres-(dUno-c))*E;
ftres=2*t*(dDos-hTres-c)*0.0021*E;

fr=(funo+fdos+ftres+fa+fb);

brazoUno=dUno+0.5*(hTres-(dUno-c));
brazoDos=dUno+2/3*(hTres-(dUno-c));
brazoTres=c+hTres+0.5*(dDos-c-hTres);

mr=(funo*(brazoUno-0.5*h)+fdos*(brazoDos-0.5*h)+ftres*(brazoTres-0.5*h)+...
    fa*(dUno-0.5*h)+fb*(dDos-0.5*h));

elemAc=[fr mr];
end

function elemAc=casoTres(c,E,t,dUno,dDos,h,fa,fb,epsa,hDos)

funo=-0.5*(c-dUno)*epsa*2*t*E;
fdos=0.5*hDos*0.0021*2*t*E;
ftres=2*t*(dDos-c-hDos)*0.0021*E;

fr=(funo+fdos+ftres-fa+fb);

brazoUno=dUno+1/3*(c-dUno);
brazoDos=c+2/3*hDos;
brazoTres=(c+hDos)+0.5*(dDos-c-hDos);

mr=(funo*(brazoUno-0.5*h)+fdos*(brazoDos-0.5*h)+ftres*(brazoTres-0.5*h)...
    -fa*(dUno-0.5*h)+fb*(dDos-0.5*h));

elemAc=[fr mr];
end

function elemAc=casoTresCuatro(c,E,t,dUno,dDos,h,fa,fb)

sUno=9/30*c-dUno;
sDos=h-51/30*c-dUno;

funo=-sUno*0.0021*2*t*E;
fdos=-0.5*(c-sUno-dUno)*0.0021*2*t*E;
ftres=0.5*2*t*(dDos-c-sDos)*0.0021*E;
fcuatro=sDos*0.0021*2*t*E;

fr=(funo+fdos+ftres+fcuatro-fa+fb);
brazoUno=dUno+0.5*sUno;
brazoDos=dUno+sUno+1/3*(c-dUno-sUno);
brazoTres=(dDos-sDos)-1/3*(dDos-c-sDos);
brazoCuatro=(dDos-0.5*sDos);

mr=(funo*(brazoUno-0.5*h)+fdos*(brazoDos-0.5*h)+ftres*(brazoTres-0.5*h)+...
    fcuatro*(brazoCuatro-0.5*h)-fa*(dUno-0.5*h)+...
        fb*(dDos-0.5*h));

elemAc=[fr mr];
end

function elemAc=casoCuatro(c,E,t,dUno,dDos,h,fa,fb,epsb,hCuatro)

brazoUno=dUno+0.5*hCuatro;
brazoDos=2/3*(dUno+hCuatro+0.5*c);
brazoTres=1/3*c+2/3*dDos;

funo=-2*t*0.0021*hCuatro*E;
fdos=-0.5*2*t*0.0021*(c-dUno-hCuatro)*E;
ftres=0.5*epsb*(dDos-c)*2*t*E;

fr=(funo+fdos+ftres-fa+fb);
mr=(funo*(brazoUno-0.5*h)+fdos*(brazoDos-0.5*h)+ftres*(brazoTres-0.5*h)-...
    fa*(dUno-0.5*h)+fb*(dDos-0.5*h));

elemAc=[fr mr];
end

function elemAc=casoCinco(E,t,dUno,dDos,h,fa,fb,hCinco,epsb)

funo=-2*t*0.0021*(hCinco-dUno)*E;
fdos=-2*t*(epsb)*(dDos-(hCinco))*E;
ftres=-0.5*(0.0021-epsb)*(dDos-(hCinco))*2*t*E;
fr=(funo+fdos+ftres-fa-fb);

brazoUno=dUno+0.5*(hCinco-dUno);
brazoDos=0.5*(hCinco+dDos);
brazoTres=1/3*dDos+2/3*hCinco;

mr=(funo*(brazoUno-0.5*h)+fdos*(brazoDos-0.5*h)+ftres*(brazoTres-0.5*h)-...
    fa*(dUno-0.5*h)-fb*(dDos-0.5*h));

elemAc=[fr mr];
end

function elemAc=casoEspecial(c,E,t,dUno,dDos,h,fa,fb,epsa,epsb)

funo=-0.5*epsa*E*(c-dUno)*2*t;
fdos=0.5*epsb*E*(dDos-c)*2*t;
fr=(funo+fdos-fa+fb);

brazoUno=dUno+1/3*(c-dUno);
brazoDos=dDos-1/3*(dDos-c);

mr=(funo*(brazoUno-0.5*h)+fdos*(brazoDos-0.5*h)-fa*(dUno-0.5*h)+...
    +fb*(dDos-0.5*h));

elemAc=[fr mr];
end

function elemAc=casoSeis(dDos,dUno,E,t,fa,fb)

fr=-(0.0021*(dDos-dUno)*E*2*t+fa+fb);
mr=0;

elemAc=[fr mr];
end