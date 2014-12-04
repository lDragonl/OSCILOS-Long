function [indexHA,indexLiner,indexDamper] = Fcn_TP_interface_location
global CI
% This function is used to locate the required interfaces
%
% first created: 2014-12-03
% last modified: 2014-12-03
% author: Jingxuan Li (jingxuan.li@imperial.ac.uk)
%
indexHA     = find(CI.CD.index == 10 | CI.CD.index == 11);                  % heat addition
indexLiner  = find(CI.CD.index == 30);                                      % left side of a liner
indexDamper = find(CI.CD.index == 2);                                       % damper
%
% -------------------------------end---------------------------------------