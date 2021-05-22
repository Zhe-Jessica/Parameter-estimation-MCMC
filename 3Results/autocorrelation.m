close all;  clc;
figure
tiledlayout(1,3)
nexttile
 autocorr(record(400:1000,1))           %cut first 400
nexttile
autocorr(record(400:1000,2))
nexttile
 autocorr(record(400:1000,3))
 clear all;
