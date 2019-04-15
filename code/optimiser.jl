##parameter optimisation for AA model##

using JuMP
using GLPK

modelAA = Model(with_optimizer(GLPK.Optimizer))


#~ @variable(modelAA, 0<=  )
#the end of the burn in phase
#~ [1.0e11, 4192.83, 333.316, 1624.11, 7.45172, 333.316, 7.45172, 24.5484, 24.5484, 72646.5, 143.077, 5350.34, 4646.88, 1.0342, 1000.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.59708e8, 12757.5, 72474.6, 5337.69, 1620.26]

gmax = 1260.0
a = 1.59708e8
AA = 1000.0


@variable(modelAA, 0 <= k_ribo_a <= 10000)
@variable(modelAA, 0 <= k_ribo_a_AA<= 10000)
@variable(modelAA, 0 <= k_ribo_AA<= 10000)
@NLobjective(modelAA, Max, (gmax*a*AA)/(k_ribo_a*k_ribo_a_AA+k_ribo_a_AA*a+k_ribo_AA_a*AA+a*AA))
@constraint(modelAA, con1, k_ribo_a + k_ribo_a_AA + k_ribo_AA >= 1)

#~ @variable(modelAA, 0<= x ) #create a jump variable called x with a lower bound at zero and upper at 2
#~ @variable(model1, 0<= y ) #create a second variable y with its own bounds
#~ @objective(model1, Max, x + y) #define what the object of the optimiser is, here maximising the giving function
#~ @constraint(model1, con, 1x + 2y <= 3) #add a constraint to the model

println("model all set up. optimising now")
optimize!(modelAA)
println("optimising complete")

objective_value(model1)
a_end = value(k_ribo_a)
println("optimal k_ribo_a value is $a_end")
a_AA_end = value(k_ribo_a_AA)
println("optimal k_ribo_a_AA value is $a_AA_end")
AA_end = value(k_ribo_AA)
println("optimal k_ribo_AA value is $AA_end")



#~ @variable(modelAA, 0<= s_out )
#~ @variable(modelAA, 0<= rmr )
#~ @variable(modelAA, 0<= em )
#~ @variable(modelAA, 0<= rmq )
#~ @variable(modelAA, 0<= rmt )
#~ @variable(modelAA, 0<= et )
#~ @variable(modelAA, 0<= rmm )
#~ @variable(modelAA, 0<= mt )
#~ @variable(modelAA, 0<= mm )
#~ @variable(modelAA, 0<= q )
#~ @variable(modelAA, 0<= si )
#~ @variable(modelAA, 0<= mq )
#~ @variable(modelAA, 0<= mr )
#~ @variable(modelAA, 0<= r )
#~ @variable(modelAA, 0<= NH4 )
#~ @variable(modelAA, 0<= nit_mrna )
#~ @variable(modelAA, 0<= nit_mrna_ribo )
#~ @variable(modelAA, 0<= nit )
#~ @variable(modelAA, 0<= exported )
#~ @variable(modelAA, 0<= N )
#~ @variable(modelAA, 0<= a )
#~ @variable(modelAA, 0<= AA )
#~ @variable(modelAA, 0<= AA_prot )
#~ @variable(modelAA, 0<= AA_mrna )
#~ @variable(modelAA, 0<= AA_mrna_ribo )
