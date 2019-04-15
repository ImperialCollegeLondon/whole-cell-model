println("Precompiling packages. This may take a minute...")
using ODE
using DifferentialEquations
using Pandas
using Plots
using Dates

include("multilevel-functions.jl")

# parameters
thetar= 426.8693338968694 #ribosome transcription threshold, molecules per cell
# s0= 1.0e4 #external substrate level, molecules
gmax= 1260.0 #max translational elongation rate, amino acid per minute molecule
cl= 0 #chloramphenicol level?
thetax= 4.379733394834643 #non ribosomal transcription threshold, molecules per cell
M= 1.0e8 #total mass of the cell (fixed)
Ms=1.0e4 #mid log mass of the cell

vm= 5800.0 #max (metabolic?) enzymatic rate, per minute
vmax_nit = 0.0#650 #max rate of nitrogenase protein., nmol NH4 per minute per mg from Barney, Yurth et al. 2009
vt= 726.0 #max rate of movement of substrate by transporter protein
vmax_AA = 0.0#1000 #max turnover of subtrates and NH4 into amino acids by AA protein


Kq= 152219.0403737490 #housekeeping proteins autoinhibition threshold
#Kp= 0 #180.1378030928276 #net rate of translating gratuitous proteins
Kt= 1.0e3 #the half maximal nutrient import threshold, molecules per cell minute
Km= 1.0e3 #half maximal enzymatic threshold (for michealas menton kinetics)

Km_nit = 0.0#Km #0.18 #half maximal enzymatic threshold for nitrogenase. from Barney, Yurth et al. 2009
#Km_AA = 100.0 #half maximal threshold for the amino acid making protein. need to look this up!


we= 4.139172187824451#max enzyme mRNA transcription rate
wr= 929.9678874564831 #rate of transcription of ribosomal protein mRNAs
wq= 948.9349882947897 #rate of transcription of housekeeping protein mRNAs
v_nit= 0.0#900 #rate of transcription of nitrogenase coding mRNAs
w_AA= 900.0 #rate of transcription of AA making protein
############################################################
############need to look this ^ up!#########################
############################################################

nq= 4 #housekeeping protein autoinhibition Hill coefficient
nr= 7549.0 #ribosome length (num amino acids)
ns= 100.0 #0.5 #substrate use effiency
nx= 300.0 #length of non ribosomal proteins (num amino acids per protein)
# parameters= [thetar,k_cm,gmax,cl,thetax,Kt,M,we,Km,vm,nx,Kq,Kp,vt,wr,wq,wp,nq,nr,ns,v_nit]

# define rate constants
# b= 0 #maybe inhibition efficiency by chloramphenicol?
dm= 0.1 #mRNA degredation rate
ds= 0.01 #rate of chemostat dilution
death_rate = ds #death rate of bacteria
kb= 1.0 #rate of mRNA-ribosome binding
ku= 1.0 #rate of mRNA-ribosome unbinding #was set to one originally

#params for gamma equation
#k_ribo_a = 0.01


#params for amino acid protein equation

# set initial conditions
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
s_out = 10e10 #external substrate
exported_0= 0.0#total amount of NH4 exported
N_0 = 1.0 #num of bacteria cells to start with
AA_0 = 20000.0#9.6e8 #num of amino acids to start with
AA_prot_0 = 1.0
AA_mrna_0=0.0
AA_mrna_ribo_0=0.0

init= [s_out,rmr_0,em_0,rmq_0,rmt_0,et_0,rmm_0,
    mt_0,mm_0,q_0,si_0,mq_0,mr_0,r_0,
    NH4_0,nit_mrna_0,nit_mrna_ribo_0,nit_0,exported_0,N_0,a_0, AA_0, AA_prot_0,AA_mrna_0,AA_mrna_ribo_0]


time1= 5000.0

#~ for i1 in (0.1,10.0,1000.0)
#~     global(k_cat_AA)= i1
#~     for i2 in (0.1,10.0,1000.0)
#~         global(k_a_NH4)=i2
#~             for i3 in (0.1,10.0,1000.0)
#~                 global(k_NH4)= i3
#~                 for i4 in (0.1,10.0,1000.0)
#~                     global(k_a_AA)=i4
#~                     for i5 in (0.1,10.0,1000.0)
#~                         global(k_NH4_AA)=i5
#~ 						for i6 in (0.1,10.0,1000.0)
#~ 							global(k_ribo_a)=i6


k_ribo_a_AA = 10000.0
k_ribo_AA_a = 10000.0

k_a = 10.0
k_cat_AA = i1 = 10.0
k_a_NH4 = i2 = 10.0
k_NH4 = i3 = 10.0
k_a_AA = i4 = 10.0
k_NH4_AA = i5 = 10.0
k_ribo_a = i6 = 10.0


problm = ODEProblem(AA_simple,init,(0.,time1))
#the line above runs the single cell model with AA. no cell growth or nitrogenase
println("solving param-sweep-$i1-$i2-$i3-$i4-$i5-$i6")
#~ start_time = time()

#~ while time() - start_time < 300.0 #stop the model run if its been running for more than 5 mins
#~ 	println("time")
	solved = ODE.solve(problm)
	df1= DataFrame(solved)
#~ 	Pandas.to_csv(df1, "../data/param-sweep-$i1-$i2-$i3-$i4-$i5-$i6.csv")#"testfile.csv"
	println("model burn in successful")
#~ 	break
#~ end

#~ println("Please input how many timesteps you would like the second burn in phase to be (10000 recommended)")
#~ input2= parse(Float64,readline(stdin))
input2 = 5000.0
time2 = time1 + input2
endstate = size(solved,2)
init2 =solved[endstate]
#~ println("setting transporter protein to 1")
#~ init2[6]=1 #this is to replicate andrea's multilevel model. setting the transporters back to 1 induces a lag in population growth
println("running cell growth section")
#global(et) = 1.0 trying to set transporter protein to 1 for the multiscale model. In 
#andrea's paper they set this to 1 to induce a lag in cell growth
problm2 = ODEProblem(AA_popn_growth,init2,(time1,time2))
println("burn in 2 in progress")
solved2 = solve(problm2)
println("burn in 2 complete")

    df1= DataFrame(solved)
    df2= DataFrame(solved2)
    dfs = [df1,df2]
    output= DataFrame(Pandas.concat(dfs))
    Pandas.to_csv(output, "multilevel-output.csv")

 plt= Plots.plot(solved,
        lab= ["S_external 1","ribo mrna comp 2","metab enzyme 3","housekpng mrna comp 4","trans mrna comp 5","transporter prot 6","metab mrna comp 7","trans mrna 8","metab mrna 9","housekepng prot 10","si11","housekpng mrna 12","ribo mrna 13","free ribo 14","NH4 int 15","nit mrna 16","nit mrna comp 17","nitrogenase18","cumulative NH4 19","num cells 20","ATP 21","AA 22","AA prot 23","AA mRNA 24","AA mrna comp 25"],
        size = (1000,500),
        ylims = (0,5e4),
        ylab="Number of Molecules",
        xlab="Timesteps",
        legend=:best,
    #     palette=:heat
    #     seriescolor=:auto
        color = ["#ffe119" "#000000" "#3cb44b" "#0082c8" "#808080" "#46f0f0" "#f032e6" "#d2f53c" "#fabebe" "#008080" "#e6beff" "#aa6e28" "#800000" "#808000" "#ffd8b1" "#000080" "#f58231" "#911eb4" "#e6194b" :grey "#aaffc3" :brown :green :blue :purple]
    )    
    plt =plot!(solved2,
        lab="",
        xlims = (0,time2),
        color = ["#ffe119" "#000000" "#3cb44b" "#0082c8" "#808080" "#46f0f0" "#f032e6" "#d2f53c" "#fabebe" "#008080" "#e6beff" "#aa6e28" "#800000" "#808000" "#ffd8b1" "#000080" "#f58231" "#911eb4" "#e6194b" :grey "#aaffc3" :brown :green :blue :purple]
    )

    savefig(plt, "current-plt.png")



#~ end
#~ end
#~ end
#~ end
#~ end
#~ end
