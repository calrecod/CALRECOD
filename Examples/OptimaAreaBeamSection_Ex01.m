clc
clear all

b=20;
h=40;
duct=3;
b_rec=4;
h_rec=3;
fc=280;
load_conditions=[1 15];
factor_fc=0.85;
E=2e6;
plotOptimConv=1;
[cbest,bestMr,bestef03,best_Area03,tbest,h]=SGD1tBeamsISR(b,h,duct,...
            b_rec,h_rec,fc,load_conditions,factor_fc,E,plotOptimConv)
                    
