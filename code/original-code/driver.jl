
#parameters
thetar= 426.8693338968694;
k_cm= 0.005990373118888;
gmax= 1260.0;
cl= 0;
thetax= 4.379733394834643;
Kt= 1.0e3;
M= 1.0e8;
we= 4.139172187824451;
Km= 1.0e3;
vm= 5800.0;
nx= 300.0;
Kq= 1.522190403737490e+05;
Kp= 0; #180.1378030928276;
vt= 726.0;
wr= 929.9678874564831;
wq= 948.9349882947897;
wp= 0.0;
nq= 4;
nr= 7549.0;
ns= 10;
dN= 0.02;
parameters= [thetar k_cm gmax cl thetax Kt M we Km vm nx Kq Kp vt wr wq wp nq nr ns dN];

# define rate constants
b= 0;
dm= 0.1;
kb= 1.0;
ku= 1.0;
f= cl*k_cm;
rates= [b dm kb ku f];

# define initial conditions
rmr_0= 0;
em_0= 0;
rmp_0= 0;
rmq_0= 0;
rmt_0= 0;
et_0= 0;
rmm_0= 0;
zmm_0= 0;
zmr_0= 0;
zmp_0= 0;
zmq_0= 0;
zmt_0= 0;
mt_0= 0;
mm_0= 0;
q_0= 0;
p_0= 0;
si_0= 0;
mq_0= 0;
mp_0= 0;
mr_0= 0;
r_0= 10.0;
a_0= 1000;
N_0= 5000;
so_0= 1.0e4;
init= [rmr_0 em_0 rmp_0 rmq_0 rmt_0 et_0 rmm_0 zmm_0 zmr_0 zmp_0 zmq_0 zmt_0 mt_0 mm_0 q_0 p_0 si_0 mq_0 mp_0 mr_0 r_0 a_0 N_0 so_0];

# call solver routine 
t0= 0;
tf= 1e9;
[t,y]= ode15s(@(t,y) cellmodel_odes(t, y, rates, parameters), [t0 tf], init);
rmr= y(:,1);
em= y(:,2);
rmp= y(:,3);
rmq= y(:,4);
rmt= y(:,5);
et= y(:,6);
rmm= y(:,7);
zmm= y(:,8);
zmr= y(:,9);
zmp= y(:,10);
zmq= y(:,11);
zmt= y(:,12);
mt= y(:,13);
mm= y(:,14);
q= y(:,15);
p= y(:,16);
si= y(:,17);
mq= y(:,18);
mp= y(:,19);
mr= y(:,20);
r= y(:,21);
a= y(:,22);
N= y(:,23);
so= y(:,24);

