##Amino acid model, burn in phases##

println("Precompiling packages. This may take a minute...")
using ODE
using DifferentialEquations
using Pandas
using Plots
using Dates

println("Would you like to run the full model or the reduced?")
println("Please input 'f' or 'r')")
model_type= string(readline(stdin))
if model_type == string('f')
	include("aa-functions.jl")
	println("running full model")
elseif model_type == string('r')
	include("reduced-function.jl")
	println("running reduced model")
else
	println("Please input either 'f' or 'r'")
end
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




################################################################
Km= 1.0e3 #1.0e4 #half maximal enzymatic threshold (for michealas menton kinetics)
#was 1.0e3 in orginal code. why did I change it??
################################################################




Km_nit = 0.0#Km #0.18 #half maximal enzymatic threshold for nitrogenase. from Barney, Yurth et al. 2009
Km_AA = 0.0#1#half maximal threshold for the amino acid making protein. need to look this up!


we= 4.139172187824451#max enzyme mRNA transcription rate
wr= 929.9678874564831 #rate of transcription of ribosomal protein mRNAs
wq= 948.9349882947897 #rate of transcription of housekeeping protein mRNAs
v_nit= 0.0#900 #rate of transcription of nitrogenase coding mRNAs
w_AA= 0#10 #rate of transcription of AA making protein
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
# f= cl*k_cm #chloramphenical level * binding rate = inhibition rate?
# rates= [b,dm,kb,ku,f]

AA_bind_rate = 0.5
atp_bind_rate = 0.5
AA_efficiency = 1/4.6 #how many AA are produced from one ATP so 1/25 means 25 atp required for one aa
# export_rate = 0.1#this much of fixed nitrogen gets exported for use by the plant
############################################################
############this parameter needs to vary #########################
############################################################



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
NH4_0 = 0.0 #num of ammonia particles fixed by nitrogenase
nit_0 = 0.0 #num nitrogenase proteins
nit_mrna_0 = 0.0 #num of nitrogenase coding mRNA molecules
nit_mrna_ribo_0 = 0.0#num of nitrogenase mRNA-ribosome complexes
s_out = 1e5 #external substrate
exported_0= 0.0#total amount of NH4 exported
N_0 = 1.0 #num of bacteria cells to start with
AA_0 = 0.0#9634500.0 #9.6e8 #num of amino acids to start with
AA_prot_0 = 0.0
AA_mrna_0=0.0
AA_mrna_ribo_0=0.0

init= [s_out,rmr_0,em_0,rmq_0,rmt_0,et_0,rmm_0,
    mt_0,mm_0,q_0,si_0,mq_0,mr_0,r_0,
    NH4_0,nit_mrna_0,nit_mrna_ribo_0,nit_0,exported_0,N_0,a_0, AA_0, AA_prot_0,AA_mrna_0,AA_mrna_ribo_0]


println("Please input how many timesteps you would like this initial burn in phase to be (10000 recommended)")
time1= parse(Float64,readline(stdin))
println("initializing burn in 1")
# println(time1)
problm = ODEProblem(model_AA,init,(0.,time1))
println("burn in 1 in progress")
solved = solve(problm)
# solved = solve(model_AA, Rodas4(),reltol=1e-8,abstol=1e-8)
# solved = solve(problm, ode45())
# solved = solve(problm, Tsit5())
# solved = solve(problm, CVODE_BDF())
println("burn in 1 completed")


println("Please input how many timesteps you would like the second burn in phase to be (10000 recommended)")
input2= parse(Float64,readline(stdin))
time2 = time1 + input2
endstate = size(solved,2)
init2 =solved[endstate]
println("initializing burn in 2")
problm2 = ODEProblem(model_AA2,init2,(time1,time2))
println("burn in 2 in progress")
solved2 = solve(problm2)
println("burn in 2 complete")


if model_type == string('r')
    df1= DataFrame(solved)
    df2= DataFrame(solved2)
    dfs = [df1,df2]
    output= DataFrame(Pandas.concat(dfs))
    Pandas.to_csv(output, "reduced-model-output.csv")
    println("output of reduced model saved as reduced-model-output.csv Thanks for running the model!")
exit()
end

println("Please input how many timesteps you would like the amino acid model to run (10000+ recommended)")

input3= parse(Float64,readline(stdin))
time3 = time2 + input3
endstate2 = size(solved2,2)
init3 =solved2[endstate2]
# println("setting problem")
# println("problem being solved")
# solved3 = solve(problm3, CVODE_BDF());
# println("Please input what start value you would like to run the parameter sweep to:")
# input_start1= parse(Float64,readline(stdin))

println("Running parameter sweep. Please hold.")


for i1 in (0.01)#,100.0,10.0,100.0,1000.0,10000.0)
    global(k_cat_AA)= i1
    for i2 in (0.01)#,100.0,10.0,100.0,1000.0,10000.0)
        global(k_a_NH4)=i2
        for i3 in (0.01)#,100.0,10.0,100.0,1000.0,10000.0)
            global(k_a)= i3
            for i4 in (0.01)#,100.0,10.0,100.0,1000.0,10000.0)
                global(k_NH4)= i4
                for i5 in (0.01)#,100.0,10.0,100.0,1000.0,10000.0)
                    global(k_a_AA)=i5
                    for i6 in (0.01)#,100.0,10.0,100.0,1000.0,10000.0)
                        global(k_NH4_AA)=i6
for i7 in (0.01)#,10.0,100.0,1000.0,10000.0)
                        global(k_ribo_a)=i7
for i8 in (0.01)#,10.0,100.0,1000.0,10000.0)
                        global(k_ribo_AA)=i8
for i9 in (0.01)#,10.0,100.0,1000.0,10000.0)
                        global(k_ribo_a_AA)=i9
for i10 in (10000.0)#,0.01,10.0,100.0,1000.0,10000.0)
                        global(k_ribo_AA_a)=i10

#    println("current values are k_cat_AA: $i1 k_a_NH4: $i2")
    problm3 = ODEProblem(model_AA3,init3,(time2,time3))

    open("../data/high-AA-params-file-$i1-$i2-$i3-$i4-$i5-$i6-$i7-$i8-$i9-$i10.csv","w") do f
    write(f,"new_AA, AA_req, AA_prot\n")
    end

    solved3 = solve(problm3)
    endstate3 = size(solved3,2)
    output =solved3[endstate3]
    # println("Results being saved in Data folder with current datetime as filename.")

    df1= DataFrame(solved)
    df2= DataFrame(solved2)
    df3= DataFrame(solved3)

    # append!(df1,df2)
    dfs = [df1,df2,df3]
    output= DataFrame(Pandas.concat(dfs))

    # using Dates
    # filename= "../data/AA-model-results-"* string(now())
    filename = "../data/high-AA-results-k_AA_"*string(i1)*"k_a_NH4"*string(i2)*"k_a"*string(i3)*"k_NH4"* string(i4)*"k_a_AA"*string(i5)*"k_NH4_AA"*string(i6)*"k_ribo_a"*string(i7)*"k_ribo_AA"*string(i8)*"k_ribo_a_AA"*string(i9)*"k_ribo_AA_a"*string(i10)
    filename2= filename* ".csv"
    Pandas.to_csv(output, filename2)




# println("Do you want to plot the results? (y/n)")
# input5= readline(stdin)
# if input5 == "y"
    # println("Plotting. This may take a minute.")
    # using Plots
    plt= Plots.plot(solved,
    #     color= [:black, :orange],
        lab= ["S_external 1","ribo mrna comp 2","metab enzyme 3","housekpng mrna comp 4","trans mrna comp 5","transporter prot 6","metab mrna comp 7","trans mrna 8","metab mrna 9","housekepng prot 10","si11","housekpng mrna 12","ribo mrna 13","free ribo 14","NH4 int 15","nit mrna 16","nit mrna comp 17","nitrogenase18","cumulative NH4 19","num cells 20","ATP 21","AA 22","AA prot 23","AA mRNA 24","AA mrna comp 25"],
        size = (1000,500),
        ylims = (0,5e3),
        ylab="Number of Molecules",
        xlab="Timesteps",
        legend=:best,
    #     palette=:heat
    #     seriescolor=:auto
        color = ["#ffe119" "#000000" "#3cb44b" "#0082c8" "#808080" "#46f0f0" "#f032e6" "#d2f53c" "#fabebe" "#008080" "#e6beff" "#aa6e28" "#800000" "#808000" "#ffd8b1" "#000080" "#f58231" "#911eb4" "#e6194b" :grey "#aaffc3" :brown :green :blue :purple]
    )
    plt =plot!(solved2,
        lab="",
        xlims = (0,time3),
        color = ["#ffe119" "#000000" "#3cb44b" "#0082c8" "#808080" "#46f0f0" "#f032e6" "#d2f53c" "#fabebe" "#008080" "#e6beff" "#aa6e28" "#800000" "#808000" "#ffd8b1" "#000080" "#f58231" "#911eb4" "#e6194b" :grey "#aaffc3" :brown :green :blue :purple]
    )
    plt =plot!(solved3,
        lab="",
        xlims = (0,time3),
        color = ["#ffe119" "#000000" "#3cb44b" "#0082c8" "#808080" "#46f0f0" "#f032e6" "#d2f53c" "#fabebe" "#008080" "#e6beff" "#aa6e28" "#800000" "#808000" "#ffd8b1" "#000080" "#f58231" "#911eb4" "#e6194b" :grey "#aaffc3" :brown :green :blue :purple]
    )
    # gui(plt)
    # display(plt)


# else
# end

# println("Do you want to save the plot? (y/n)")
# input4= readline(stdin)
# if input4 == "y"
    # using Dates
    # using CSV
    # filename= "../data/AA-model-plot-"* string(now())
    filename= "../data/high-AA-plot-k_AA_"*string(i1)*"k_a_NH4"*string(i2)*"k_a"*string(i3)*"k_NH4"*string(i4)*"k_a_AA"*string(i5)*"k_NH4_AA" *string(i6)*"k_ribo_a"*string(i7)*"k_ribo_AA"*string(i8)*"k_ribo_a_AA"*string(i9)*"k_ribo_AA_a"*string(i10)
    savefig(plt, filename)
    println("Plotted and saved model k_cat_AA: ", i1," k_a_NH4: ", i2,"k_a: ",i3,"k_NH4: ",i4,"k_a_AA: ",i5,"k_NH4_AA: ",i6,"k_ribo_a: ",i7,"k_ribo_AA: ",i8,"k_ribo_a_AA: ",i9,"k_ribo_AA_a: ",i10)
    # filename2= filename* ".csv"
    # CSV.write(filename2, output)

end
end
end
end
end
end
end
end
end
end

println("Thanks for running the model!")
