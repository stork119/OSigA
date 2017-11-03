function R = JAKSTAT_rates(x, par, t, stimulus)
R =[
par(1)*(par(13)*x(7)+par(11)*x(8))*x(1)/(par(6)+x(1));
par(2)*par(11)*x(2)*x(2);
par(3)*par(11)*x(3);
par(4)*par(12)*x(4);
par(14)*par(12)*x(5);
par(5)*par(7)*stimulus(1)*par(13)*x(6);
par(8)*par(13)*x(7);
par(9)*par(13)*x(7);
par(10)*par(11)*x(8);
];
 end
