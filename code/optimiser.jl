##parameter optimisation for AA model##

using JuMP
#~ using GLPK
using Ipopt

#~ modelAA = Model(with_optimizer(GLPK.Optimizer))
modelAA = Model(with_optimizer(Ipopt.Optimizer))

#~ @variable(modelAA, 0<=  )
#the end of the burn in phase
#~ [1.0e11, 4192.83, 333.316, 1624.11, 7.45172, 333.316, 7.45172, 24.5484, 24.5484, 72646.5, 143.077, 5350.34, 4646.88, 1.0342, 1000.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.59708e8, 12757.5, 72474.6, 5337.69, 1620.26]

#~ a = 1.59708e8
#~ AA = 1000.0

rmr_0= 0.0
em_0= 0.0
rmq_0= 0.0
rmt_0= 0.0
et_0= 0.0
rmm_0= 0.0
mt_0= 0.0
mm_0= 0.0
q_0= 0.0
si_0= 0.0
mq_0= 0.0
mr_0= 0.0
r_0= 10.0
a_0= 1000.0
NH4_0 = 1000.0 #num of ammonia particles fixed by nitrogenase
nit_0 = 0.0 #num nitrogenase proteins
nit_mrna_0 = 0.0 #num of nitrogenase coding mRNA molecules
nit_mrna_ribo_0 = 0.0#num of nitrogenase mRNA-ribosome complexes
s_out_0 = 10e10 #external substrate
exported_0= 0.0#total amount of NH4 exported
N_0 = 1.0 #num of bacteria cells to start with
AA_0 = 20000.0#9.6e8 #num of amino acids to start with
AA_prot_0 = 1.0
AA_mrna_0=0.0
AA_mrna_ribo_0=0.0

init= [s_out_0,rmr_0,em_0,rmq_0,rmt_0,et_0,rmm_0,
    mt_0,mm_0,q_0,si_0,mq_0,mr_0,r_0,
    NH4_0,nit_mrna_0,nit_mrna_ribo_0,nit_0,exported_0,N_0,a_0, AA_0, AA_prot_0,AA_mrna_0,AA_mrna_ribo_0]


@variable(modelAA, 0 <= k_ribo_a_var)
@variable(modelAA, 0 <= k_ribo_a_AA_var)
@variable(modelAA, 0 <= k_ribo_AA_a_var)
@variable(modelAA, 0 <= k_a_var)
@variable(modelAA, 0 <= k_cat_AA_var)
@variable(modelAA, 0 <= k_NH4_var)
@variable(modelAA, 0 <= k_a_NH4_var)
@variable(modelAA, 0 <= k_a_AA_var)
@variable(modelAA, 0 <= k_NH4_AA_var)


#~ params = [k_ribo_a_var,k_ribo_a_AA_var,k_ribo_AA_a_var,k_a_var,k_cat_AA_var,k_NH4_var,k_a_AA_var,k_NH4_AA_var]
#~ @variable(modelAA, 0<= params[1:9])

include("optimise-funtion.jl")
JuMP.register(modelAA, :gamma_finder, 9, gamma_finder, autodiff = true)

@NLobjective(modelAA, Max, gamma_finder(k_ribo_a_var,k_ribo_a_AA_var,k_ribo_AA_a_var,
    k_a_var,k_cat_AA_var,k_a_NH4_var,k_NH4_var,k_a_AA_var,k_NH4_AA_var))
    #s_out_0,rmr_0,em_0,rmq_0,rmt_0,et_0,rmm_0,mt_0,mm_0,q_0,si_0,mq_0,mr_0,r_0,
    #NH4_0,nit_mrna_0,nit_mrna_ribo_0,nit_0,exported_0,N_0,a_0, AA_0, AA_prot_0,AA_mrna_0,AA_mrna_ribo_0))

#~ @NLobjective(modelAA, Max, (gmax*a*AA)/(k_ribo_a*k_ribo_a_AA+k_ribo_a_AA*a+k_ribo_AA_a*AA+a*AA))
#~ @NLconstraint(modelAA, con1, k_ribo_a + k_ribo_a_AA + k_ribo_AA_a >= 1)
println("model set up")
#~ @variable(modelAA, 0<= x ) #create a jump variable called x with a lower bound at zero and upper at 2
#~ @variable(model1, 0<= y ) #create a second variable y with its own bounds
#~ @objective(model1, Max, x + y) #define what the object of the optimiser is, here maximising the giving function
#~ @constraint(model1, con, 1x + 2y <= 3) #add a constraint to the model

println("model all set up. optimising now")
optimize!(modelAA)
println("optimising complete")

objective_value(modelAA)
a_end = value(k_ribo_a)
println("optimal k_ribo_a value is $a_end")
a_AA_end = value(k_ribo_a_AA)
println("optimal k_ribo_a_AA value is $a_AA_end")
AA_end = value(k_ribo_AA_a)
println("optimal k_ribo_AA value is $AA_end")

print((gmax*a*AA)/(a_end*a_AA_end+a_AA_end*a+AA_end*AA+a*AA))

