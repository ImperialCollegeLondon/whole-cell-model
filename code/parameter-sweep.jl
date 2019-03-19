println("Precompiling packages. This may take a minute...")
using ODE
using DifferentialEquations
using Pandas
using Plots
using Dates

include("parameter-sweep-functions.jl")

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

#k_ribo_a = 0.01
#k_ribo_AA= 0.01
#k_ribo_a_AA = 0.01
#k_ribo_AA_a = 1000.0


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


for i1 in (0.1,10.0,1000.0)
    global(k_cat_AA)= i1
    for i2 in (0.1,10.0,1000.0)
        global(k_a_NH4)=i2 #this one can be fixed maybe
        for i3 in (0.1,10.0,1000.0)
            global(k_a)= i3 # k_a set to 1000. maybe try higher in future
            for i4 in (1000.0)
                global(k_NH4)= i4
                for i5 in (0.1,10.0,1000.0)
                    global(k_a_AA)=i5
                    for i6 in (0.1,10.0,1000.0)
                        global(k_NH4_AA)=i6
						for i7 in (0.1,10.0,1000.0)
							global(k_ribo_a)=i7
							for i8 in (1.0)
								global(k_ribo_a_AA)=i8
								for i9 in (1.0)
									global(k_ribo_AA_a)=i9




time1= 25000.0
problm = ODEProblem(AA_simple,init,(0.,time1))
#the line above runs the single cell model with AA. no cell growth or nitrogenase
solved = solve(problm)
df1= DataFrame(solved)
Pandas.to_csv(df1, "../data/$i1/param-sweep-$i1-$i2-$i3-$i4-$i5-$i6-$i7-$i8-$i9.csv")#"testfile.csv"
println("param-sweep-$i1-$i2-$i3-$i4-$i5-$i6-$i7-$i8-$i9 was a success!")
#    open("../data/param-sweep-$i1-$i2-$i3-$i4-$i5-$i6-$i7-$i8-$i9.csv","w") do f
 #   write(f,"AA, ATP, num_cells\n")
 #   end


end
end
end
end
end
end
end
end
end
#end



