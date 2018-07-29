function mpc = PJM_5Bus
mpc.baseMVA = 100
mpc.BusLen = 5
% Generator Data							
% 	Bus No	Gen_index	a 	b 	c	Pmin	Pmax
mpc.gen = [
	1	1	0.1	12	100	50	800;
	1	2	0.1	20	30	20	800;
	1	3	0.1	8	40	10	800;
	4	1	0.1	10	50	50	800;
	4	2	0.1	25	150	10	800;
	4	3	0.1	19	70	20	800;
	5	1	0.1	17	60	20	800;
	5	2	0.12	18	50	10	800;
	5	3	0.085	20	80	10	800;
];

% Line Data							
% 	fbus	tbus	r	x	b	rateA	
mpc.branch = [
	1	2	0	0.0181	0	100000;	
	1	4	0	0.0304	0	100000;	
	1	5	0	0.0064	0	300;	
	2	3	0	0.0108	0	100000;	
	3	4	0	0.0297	0	100000;	
	4	5	0	0.0297	0	100000;	
];

% Load Bus							
mpc.load = [
	2;
	3;
];