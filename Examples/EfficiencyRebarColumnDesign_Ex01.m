
b=40;
h=50;
E=2.1e6;
npuntos=30;

numberRebars_hdimension=4;
numberRebars_bdimension=3;


nv=2*numberRebars_hdimension+2*numberRebars_bdimension;
dv=6/8*2.54;
concreteCover=[4 4];

[dispositionRebar]=RebarDisposition(b,...
    h,concreteCover,dv,nv,numberRebars_hdimension,numberRebars_bdimension)

load_conditions=[1 0.89 21 14];

fdpc=280*0.85;
beta=0.85;
av=pi/4*(6/8*2.54)^2;
ov=6;
RebarTypeIndex=3;

As=nv*av;

[diagrama,mexef,eficiencia,cxy]=diagramasDisposicion(As,b,h,E,npuntos,...
               fdpc,nv,beta,ov,av,dispositionRebar,load_conditions);
           
bestArrangement=zeros(nv,1)+RebarTypeIndex;

diagramsFinalRebarCols(load_conditions,diagrama,dispositionRebar,...
                h,b,bestArrangement);