#~Kgamma= gmax/180 #Kgamma is the half maximal translational elongation threshold, gmax is the max rate of translational elongation and Kp is net rate of translation of housekeeping proteins
Kgamma = 3.0e8
#should equal 7 molecules per cell
thetar= 426.8693338968694 #ribosome transcription threshold, molecules per cell
gmax= 1260.0 #max translational elongation rate, amino acid per minute molecule
thetax= 4.379733394834643 #non ribosomal transcription threshold, molecules per cell
M= 1.0e8 #total mass of the cell (fixed)
#Ms=1.0e4 #mid log mass of the cell
vm= 5800.0 #max metabolic enzymatic rate, per minute
vmax_nit = 650 #max rate of nitrogenase protein., nmol NH4 per minute per mg from Barney, Yurth et al. 2009
vt= 726.0 #max rate of movement of substrate by transporter protein
vmax_AA = 1000 #max turnover of subtrates and NH4 into amino acids by AA protein
Kq= 152219.0403737490 #housekeeping proteins autoinhibition threshold
#Kp= 0 #180.1378030928276 #net rate of translating gratuitous proteins
Kt= 1e10 #1.0e3 #the half maximal nutrient import threshold, molecules per cell minute
Km= 1.0e3 #half maximal enzymatic threshold (for michealas menton kinetics)
Km_nit = 0.18 #half maximal enzymatic threshold for nitrogenase. from Barney, Yurth et al. 2009
we= 4.139172187824451#max enzyme mRNA transcription rate
wr= 929.9678874564831 #rate of transcription of ribosomal protein mRNAs
wq= 948.9349882947897 #rate of transcription of housekeeping protein mRNAs
v_nit= 100.0 #900 #rate of transcription of nitrogenase coding mRNAs
w_AA= 10.0 #900.0 #rate of transcription of AA making protein
a_per_AA = 2. #the number of ATP used to make one AA
NH4_per_AA = 2.
############################################################
############need to look this ^ up!#########################
############################################################
nq= 4 #housekeeping protein autoinhibition Hill coefficient
nr= 7549.0 #ribosome length (num amino acids)
ns= 100.0 #0.5 #substrate use effiency
nx= 300.0 #length of non ribosomal proteins (num amino acids per protein)
dm= 0.1 #mRNA degredation rate
ds= 0.01 #rate of chemostat dilution
death_rate = ds #death rate of bacteria
kb= 1.0 #rate of mRNA-ribosome binding
ku= 1.0 #rate of mRNA-ribosome unbinding #was set to one originally
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
NH4_0 = 1e10 #num of ammonia particles
nit_0 = 0.0 #num nitrogenase proteins
nit_mrna_0 = 0.0 #num of nitrogenase coding mRNA molecules
nit_mrna_ribo_0 = 0.0#num of nitrogenase mRNA-ribosome complexes
s_out = 1e10 #external substrate
exported_0= 0.0#total amount of NH4 exported
N_0 = 10.0 #num of bacteria cells to start with
a_0= 1000.0
AA_0 = 1.0 #9.6e8 #num of amino acids to start with
AA_prot_0 = 0.0
AA_mrna_0=0.0
AA_mrna_ribo_0=0.0

init= [s_out,rmr_0,em_0,rmq_0,rmt_0,et_0,rmm_0,
    mt_0,mm_0,q_0,si_0,mq_0,mr_0,r_0,
    NH4_0,nit_mrna_0,nit_mrna_ribo_0,nit_0,exported_0,N_0,a_0, AA_0, AA_prot_0,AA_mrna_0,AA_mrna_ribo_0]

init_original= [s_out,rmr_0,em_0,rmq_0,rmt_0,et_0,rmm_0,
    mt_0,mm_0,q_0,si_0,mq_0,mr_0,r_0,
    0.0,0.0,0.0,0.0,0.0,N_0,a_0, 0.0, 0.0,0.0,0.0]


# k_ribo_a_AA = 1e4 #100000000.0
# k_ribo_AA_a = 1e4 #10000000000.0
# k_ribo_a = 1.0

k_ribo_a_AA = 10000.0
k_ribo_AA_a = 10000.0
k_ribo_a = 10.0


"""The most basic ODE model of Amino Acid metabolism with no cell growth and no changes in external carbon or nitrogen sources.
Nitrogenase equations are still present but not active (inital conditions = 0, dx/dt = 0). Updated gamma equation (rate of translational elongation) that uses both ATP and amino acids as substrates.
    """
function AA_simple(results,y,p,t)
    s_out = y[1] #external level of substrate, a proxy for plant growth state since more plant = more substrate
    rmr= y[2]#num of ribosome/ribosomal mRNA complex molecules
    em= y[3]#metabolic enzyme molecumes
    rmq= y[4]#num of ribosome/housekeeping mRNA complex molecules
    rmt= y[5]#num of ribosome/transporter protein mRNA complex molecules
    et= y[6]#num of transporter enzyme molecules
    rmm= y[7]#num of ribosome/metabolic enzyme mRNA complexes
    mt= y[8]#num of particles of transporter protein mRNA
    mm= y[9]#num of particles of metabolic enzyme mRNA
    q= y[10]#num of housekeeping protein molecules
    si= y[11]#num of substrate molecules inside the cell (internal)
    mq= y[12]#num of mRNA coding for housekeeping proteins
    mr= y[13]#num of mRNA coding for ribosomes
    r= y[14]#num of free ribosomes
    NH4 = y[15] #num of ammonia molecules internally
    nit_mrna = y[16]#num of molecules of nitrogen fixing protein
    nit_mrna_ribo = y[17]#num of molecules of mRNA coding for nitrogen fixing protein
    nit = y[18]#num of molecules of mRNA-ribosome complexes coding for nitrogenase
    exported= y[19]#total amount of NH4 exported
    N = y[20] #num of bacteria cells
    a= y[21]#num of ATP (proxy for level of energy)
    AA = y[22] #num of amino acid molecules
    AA_prot=y[23]#num of protein molecules that make new amino acids
    AA_mrna=y[24]#num of mRNA coding for AA making protein
    AA_mrna_ribo=y[25]#num of mRNA coding for AA making protein bound to ribosomes

    #define more constants

    #gam= gmax*a/(Kgamma + a) #gamma is the current rate of translational elongation? gmax is max rate of translational elongation and a is ATP (energy)level

    gam= (gmax*a*AA)/(k_ribo_a*k_ribo_a_AA+k_ribo_a_AA*a+k_ribo_AA_a*AA+a*AA) #updated gamma equation that uses both ATP and AA

    ttrate= (rmq + rmr + rmt + rmm + nit_mrna_ribo+AA_mrna_ribo)*gam #total translation rate (sum of the mRNA/ribosome complexes times translation rate)
    #lam= ttrate/(M-AA) #lambda, the growth rate/dilution rate, the ratio of total translation rate to total cell mass
    lam= ttrate/M
    #if lam<0.00000001 lam=0.00000001 end
    nucat= em*vm*si/(Km + si) #nucat is the num of internal substrate molecules used by catabolism. em is the num of molecules of metabolic enzyme, vm is the rate of the metabolic enzyme and si is the internal substrate level
    #this is all over the half maximal metabolic rate threshold plus the internal substrate level
    fixation = 0#vmax_nit*a/(Kgamma +a)

    AA_vo = ((k_cat_AA*2*a*NH4)/(k_a_NH4*k_a))/(1+(NH4/k_a_NH4)+(NH4/k_NH4)+2+AA*(k_x))
#     new_AA = AA_vo
    new_AA = AA_vo*AA_prot

    AA_a_use = new_AA*a_per_AA #each AA produced uses x number of ATP molecules

    ################################################
    #export rate needs to vary with AA levels#
    ################################################

    ##################################################################

    results[1]= 0#((-et*vt*s_out/(Kt + s_out))*N)-ds*s_out +(export_rate*NH4*50) +1e7 #change in amount of substrate, a proxy for amount of host plant
    results[2]= kb*r*mr-ku*rmr-(gam/nr)*rmr-lam*rmr#-dp*r
    #num of ribosome-mRNA complexes coding for ribosomes
    results[3]= (gam/nx)*rmm-lam*em#-dp*em
    #num of metabolic proteins
    results[4]= kb*r*mq-ku*rmq-(gam/nx)*rmq-lam*rmq
    #num of ribosome-mRNA complexes coding for the housekeeping proteins
    results[5]= kb*r*mt-ku*rmt-gam/nx*rmt-lam*rmt
    #num of ribosome-mRNA complexes coding for the transporter protein
    results[6]= (gam/nx)*rmt-lam*et#-dp*et
    #num of transporter proteins
    results[7]= kb*r*mm-ku*rmm-gam/nx*rmm-lam*rmm
    #num of ribosome-mRNA complexes coding for the metabolic protein
    results[8]= (we*a/(thetax + a))+ku*rmt+(gam/nx)*rmt-kb*r*mt-dm*mt-lam*mt
    #change in num of mRNAs coding for transporter proteins
    results[9]= (we*a/(thetax + a))+ku*rmm+(gam/nx)*rmm-kb*r*mm-dm*mm-lam*mm
    #change in num of mRNAs coding for metabolic proteins
    results[10]= (gam/nx)*rmq-lam*q#-dp*q
    #num of housekeeping proteins
    results[11]= (et*vt*s_out/(Kt + s_out))-nucat-lam*si#-degredation*si
    results[12]= (wq*a/(thetax + a)/(1 + (q/Kq)^nq))+ku*rmq+(gam/nx)*rmq-kb*r*mq-dm*mq-lam*mq
    #num of free housekeeping mRNAs
    results[13]= (wr*a/(thetar + a))+ku*rmr+(gam/nr)*rmr-kb*r*mr-dm*mr-lam*mr
    #num of free ribosomal mRNAs
    results[14]= ku*rmr +ku*rmt+ ku*rmm+ ku*rmq +(gam/nr)*rmr+
    (gam/nr)*rmr +(gam/nx)*rmt +(gam/nx)*rmm +(gam/nx)*rmq -
    kb*r*mr- kb*r*mt- kb*r*mm -kb*r*mq +ku*nit_mrna_ribo +
    (gam/nx)*nit_mrna_ribo-kb*r*nit_mrna +ku*AA_mrna_ribo +
    (gam/nx)*AA_mrna_ribo-kb*r*AA_mrna -lam*r
    #num of free ribosomes

    #new nitrogen fixing equations. also opposite term in ATP equation
    ##########################################
    results[15]= 0 #fixation*nit-export_rate*NH4 -lam*NH4
    #change in num of NH4 molecules in the cell from fixation rate - export rate
    results[16]= 0#(v_nit*a/(thetax + a))+(ku*nit_mrna_ribo)+(gam/nx*nit_mrna_ribo)-(kb*r*nit_mrna)-(dm*nit_mrna)-(lam*nit_mrna)
    #num of free nitrogenase coding mRNA
    results[17]= 0#(kb*r*nit_mrna)-(ku*nit_mrna_ribo)-(gam/nx*nit_mrna_ribo)-(lam*nit_mrna_ribo)
    #num of nit mRNA-ribosome complexes
    results[18]= 0 #(gam/nx)*nit_mrna_ribo-lam*nit
    #num of nitrogenase proteins
    results[19]= 0#(export_rate*NH4)#total num of NH4 molecules exported this timestep
    results[20]=0#lam*N-death_rate*N #the change in num of bacterial cells, instead of just loss to dilution
    results[21]= ns*nucat-ttrate-fixation*nit-AA_a_use-lam*a
    ##########################################
    ## AAA equations##
    ##########################################
    results[22]= new_AA -ttrate-lam*AA #num of AA
    results[23]=(gam/nx)*AA_mrna_ribo-lam*AA_prot #num of AA making proteins
    results[24]=((w_AA*a/(thetax + a)))+(ku*AA_mrna_ribo)+
    (gam/nx*AA_mrna_ribo)-(kb*r*AA_mrna)-(dm*AA_mrna)-(lam*AA_mrna)#num of free AA mRNAs
    results[25]=(kb*r*AA_mrna)-(ku*AA_mrna_ribo)-(gam/nx*AA_mrna_ribo)-(lam*AA_mrna_ribo)#num of AA mRNAs bound to ribosomes


end




"""The most basic ODE model of Amino Acid metabolism with no cell growth and no changes in external carbon or nitrogen sources.
Nitrogenase equations are still present but not active (inital conditions = 0, dx/dt = 0). Updated gamma equation (rate of translational elongation) that uses both ATP and amino acids as substrates.
    """
function AA_simple_with_popn_growth(results,y,p,t)
    s_out = y[1] #external level of substrate, a proxy for plant growth state since more plant = more substrate
    rmr= y[2]#num of ribosome/ribosomal mRNA complex molecules
    em= y[3]#metabolic enzyme molecumes
    rmq= y[4]#num of ribosome/housekeeping mRNA complex molecules
    rmt= y[5]#num of ribosome/transporter protein mRNA complex molecules
    et= y[6]#num of transporter enzyme molecules
    rmm= y[7]#num of ribosome/metabolic enzyme mRNA complexes
    mt= y[8]#num of particles of transporter protein mRNA
    mm= y[9]#num of particles of metabolic enzyme mRNA
    q= y[10]#num of housekeeping protein molecules
    si= y[11]#num of substrate molecules inside the cell (internal)
    mq= y[12]#num of mRNA coding for housekeeping proteins
    mr= y[13]#num of mRNA coding for ribosomes
    r= y[14]#num of free ribosomes
    NH4 = y[15] #num of ammonia molecules fixed
    nit_mrna = y[16]#num of molecules of nitrogen fixing protein
    nit_mrna_ribo = y[17]#num of molecules of mRNA coding for nitrogen fixing protein
    nit = y[18]#num of molecules of mRNA-ribosome complexes coding for nitrogenase
    exported= y[19]#total amount of NH4 exported
    N = y[20] #num of bacteria cells
    a= y[21]#num of ATP (proxy for level of energy)
    AA = y[22] #num of amino acid molecules
    AA_prot=y[23]#num of protein molecules that make new amino acids
    AA_mrna=y[24]#num of mRNA coding for AA making protein
    AA_mrna_ribo=y[25]#num of mRNA coding for AA making protein bound to ribosomes
    #define more constants

    #gam= gmax*a/(Kgamma + a) #gamma is the current rate of translational elongation? gmax is max rate of translational elongation and a is ATP (energy)level
    gam= (gmax*a*AA)/(k_ribo_a*k_ribo_a_AA+k_ribo_a_AA*a+k_ribo_AA_a*AA+a*AA) #two substrate gamma equation that uses both ATP and AA
    ttrate= (rmq + rmr + rmt + rmm + nit_mrna_ribo+AA_mrna_ribo)*gam #total translation rate (sum of the mRNA/ribosome complexes times translation rate)
    #lam= ttrate/(M-AA) #lambda, the growth rate/dilution rate, the ratio of total translation rate to total cell mass
    lam= ttrate/M
    nucat= em*vm*si/(Km + si) #nucat is the num of internal substrate molecules used by catabolism. em is the num of molecules of metabolic enzyme, vm is the rate of the metabolic enzyme and si is the internal substrate level
    #this is all over the half maximal metabolic rate threshold plus the internal substrate level
    fixation = 0#vmax_nit*a/(Kgamma +a)

    AA_vo = ((k_cat_AA*2*a*NH4)/(k_a_NH4*k_a))/(1+(NH4/k_a_NH4)+(NH4/k_NH4)+2+AA*(k_x))
    new_AA = AA_vo*AA_prot

    AA_a_use = new_AA*a_per_AA #each AA produced uses x number of ATP molecules
    AA_NH4_use = new_AA*NH4_per_AA

    ################################################
    #export rate needs to vary with AA levels#
    ################################################

    ##################################################################

    results[1]= ((-et*vt*s_out/(Kt + s_out))*N)-ds*s_out +sub_inflow #+(export_rate*NH4*50)  #change in amount of substrate, a proxy for amount of host plant
    results[2]= kb*r*mr-ku*rmr-(gam/nr)*rmr-lam*rmr#-dp*r
    #num of ribosome-mRNA complexes coding for ribosomes
    results[3]= (gam/nx)*rmm-lam*em#-dp*em
    #num of metabolic proteins
    results[4]= kb*r*mq-ku*rmq-(gam/nx)*rmq-lam*rmq
    #num of ribosome-mRNA complexes coding for the housekeeping proteins
    results[5]= kb*r*mt-ku*rmt-gam/nx*rmt-lam*rmt
    #num of ribosome-mRNA complexes coding for the transporter protein
    results[6]= (gam/nx)*rmt-lam*et#-dp*et
    #num of transporter proteins
    results[7]= kb*r*mm-ku*rmm-gam/nx*rmm-lam*rmm
    #num of ribosome-mRNA complexes coding for the metabolic protein
    results[8]= (we*a/(thetax + a))+ku*rmt+(gam/nx)*rmt-kb*r*mt-dm*mt-lam*mt
    #change in num of mRNAs coding for transporter proteins
    results[9]= (we*a/(thetax + a))+ku*rmm+(gam/nx)*rmm-kb*r*mm-dm*mm-lam*mm
    #change in num of mRNAs coding for metabolic proteins
    results[10]= (gam/nx)*rmq-lam*q#-dp*q
    #num of housekeeping proteins
    results[11]= (et*vt*s_out/(Kt + s_out))-nucat-lam*si#-degredation*si
    results[12]= (wq*a/(thetax + a)/(1 + (q/Kq)^nq))+ku*rmq+(gam/nx)*rmq-kb*r*mq-dm*mq-lam*mq
    #num of free housekeeping mRNAs
    results[13]= (wr*a/(thetar + a))+ku*rmr+(gam/nr)*rmr-kb*r*mr-dm*mr-lam*mr
    #num of free ribosomal mRNAs
    results[14]= ku*rmr +ku*rmt+ ku*rmm+ ku*rmq +(gam/nr)*rmr+
    (gam/nr)*rmr +(gam/nx)*rmt +(gam/nx)*rmm +(gam/nx)*rmq -
    kb*r*mr- kb*r*mt- kb*r*mm -kb*r*mq +ku*nit_mrna_ribo +
    (gam/nx)*nit_mrna_ribo-kb*r*nit_mrna +ku*AA_mrna_ribo +
    (gam/nx)*AA_mrna_ribo-kb*r*AA_mrna -lam*r

    #num of free ribosomes

    #new nitrogen fixing equations. also opposite term in ATP equation
    ##########################################
    results[15]= NH4_inflow -AA_NH4_use -lam*NH4 -ds*NH4 #fixation*nit-export_rate*NH4 -lam*NH4
    #change in num of NH4 molecules in the cell from fixation rate - export rate
    results[16]= 0#(v_nit*a/(thetax + a))+(ku*nit_mrna_ribo)+(gam/nx*nit_mrna_ribo)-(kb*r*nit_mrna)-(dm*nit_mrna)-(lam*nit_mrna)
    #num of free nitrogenase coding mRNA
    results[17]= 0#(kb*r*nit_mrna)-(ku*nit_mrna_ribo)-(gam/nx*nit_mrna_ribo)-(lam*nit_mrna_ribo)
    #num of nit mRNA-ribosome complexes
    results[18]= 0 #(gam/nx)*nit_mrna_ribo-lam*nit
    #num of nitrogenase proteins
    results[19]= 0#(export_rate*NH4)#total num of NH4 molecules exported this timestep
    results[20]=lam*N -death_rate*N #the change in num of bacterial cells, instead of just loss to dilution
    results[21]= ns*nucat-ttrate-AA_a_use-lam*a #ATP equation
    ##########################################
    ## AAA equations##
    ##########################################
    results[22]= new_AA -ttrate-lam*AA #num of AA
    results[23]=(gam/nx)*AA_mrna_ribo-lam*AA_prot #num of AA making proteins
    results[24]=((w_AA*a/(thetax + a)))+(ku*AA_mrna_ribo)+
    (gam/nx*AA_mrna_ribo)-(kb*r*AA_mrna)-(dm*AA_mrna)-(lam*AA_mrna)#num of free AA mRNAs
    results[25]=(kb*r*AA_mrna)-(ku*AA_mrna_ribo)-(gam/nx*AA_mrna_ribo)-(lam*AA_mrna_ribo)#num of AA mRNAs bound to ribosomes


end

"""The most basic ODE model of Amino Acid metabolism with no cell growth and no changes in external carbon or nitrogen sources.
Nitrogenase equations are still present but not active (inital conditions = 0, dx/dt = 0). Updated gamma equation (rate of translational elongation) that uses both ATP and amino acids as substrates.
    """
function AA_simple_with_popn_growth_no_replenish(results,y,p,t)
    s_out = y[1] #external level of substrate, a proxy for plant growth state since more plant = more substrate
    rmr= y[2]#num of ribosome/ribosomal mRNA complex molecules
    em= y[3]#metabolic enzyme molecumes
    rmq= y[4]#num of ribosome/housekeeping mRNA complex molecules
    rmt= y[5]#num of ribosome/transporter protein mRNA complex molecules
    et= y[6]#num of transporter enzyme molecules
    rmm= y[7]#num of ribosome/metabolic enzyme mRNA complexes
    mt= y[8]#num of particles of transporter protein mRNA
    mm= y[9]#num of particles of metabolic enzyme mRNA
    q= y[10]#num of housekeeping protein molecules
    si= y[11]#num of substrate molecules inside the cell (internal)
    mq= y[12]#num of mRNA coding for housekeeping proteins
    mr= y[13]#num of mRNA coding for ribosomes
    r= y[14]#num of free ribosomes
    NH4 = y[15] #num of ammonia molecules fixed
    nit_mrna = y[16]#num of molecules of nitrogen fixing protein
    nit_mrna_ribo = y[17]#num of molecules of mRNA coding for nitrogen fixing protein
    nit = y[18]#num of molecules of mRNA-ribosome complexes coding for nitrogenase
    exported= y[19]#total amount of NH4 exported
    N = y[20] #num of bacteria cells
    a= y[21]#num of ATP (proxy for level of energy)
    AA = y[22] #num of amino acid molecules
    AA_prot=y[23]#num of protein molecules that make new amino acids
    AA_mrna=y[24]#num of mRNA coding for AA making protein
    AA_mrna_ribo=y[25]#num of mRNA coding for AA making protein bound to ribosomes

    #define more constants

    #gam= gmax*a/(Kgamma + a) #gamma is the current rate of translational elongation? gmax is max rate of translational elongation and a is ATP (energy)level
    gam= (gmax*a*AA)/(k_ribo_a*k_ribo_a_AA+k_ribo_a_AA*a+k_ribo_AA_a*AA+a*AA) #two substrate gamma equation that uses both ATP and AA
    ttrate= (rmq + rmr + rmt + rmm + nit_mrna_ribo+AA_mrna_ribo)*gam #total translation rate (sum of the mRNA/ribosome complexes times translation rate)
    #lam= ttrate/(M-AA) #lambda, the growth rate/dilution rate, the ratio of total translation rate to total cell mass
    lam= ttrate/M
    nucat= em*vm*si/(Km + si) #nucat is the num of internal substrate molecules used by catabolism. em is the num of molecules of metabolic enzyme, vm is the rate of the metabolic enzyme and si is the internal substrate level
    #this is all over the half maximal metabolic rate threshold plus the internal substrate level
    fixation = 0#vmax_nit*a/(Kgamma +a)

    AA_vo = ((k_cat_AA*2*a*NH4)/(k_a_NH4*k_a))/(1+(NH4/k_a_NH4)+(NH4/k_NH4)+2+AA*(k_x))
    new_AA = AA_vo*AA_prot

    AA_a_use = new_AA*a_per_AA #each AA produced uses x number of ATP molecules
    AA_NH4_use = new_AA*NH4_per_AA
    ################################################
    #export rate needs to vary with AA levels#
    ################################################

    ##################################################################

    results[1]= ((-et*vt*s_out/(Kt + s_out))*N)# +sub_inflow #+(export_rate*NH4*50)  #change in amount of substrate, a proxy for amount of host plant
    results[2]= kb*r*mr-ku*rmr-(gam/nr)*rmr-lam*rmr#-dp*r
    #num of ribosome-mRNA complexes coding for ribosomes
    results[3]= (gam/nx)*rmm-lam*em#-dp*em
    #num of metabolic proteins
    results[4]= kb*r*mq-ku*rmq-(gam/nx)*rmq-lam*rmq
    #num of ribosome-mRNA complexes coding for the housekeeping proteins
    results[5]= kb*r*mt-ku*rmt-gam/nx*rmt-lam*rmt
    #num of ribosome-mRNA complexes coding for the transporter protein
    results[6]= (gam/nx)*rmt-lam*et#-dp*et
    #num of transporter proteins
    results[7]= kb*r*mm-ku*rmm-gam/nx*rmm-lam*rmm
    #num of ribosome-mRNA complexes coding for the metabolic protein
    results[8]= (we*a/(thetax + a))+ku*rmt+(gam/nx)*rmt-kb*r*mt-dm*mt-lam*mt
    #change in num of mRNAs coding for transporter proteins
    results[9]= (we*a/(thetax + a))+ku*rmm+(gam/nx)*rmm-kb*r*mm-dm*mm-lam*mm
    #change in num of mRNAs coding for metabolic proteins
    results[10]= (gam/nx)*rmq-lam*q#-dp*q
    #num of housekeeping proteins
    results[11]= (et*vt*s_out/(Kt + s_out))-nucat-lam*si#-degredation*si
    results[12]= (wq*a/(thetax + a)/(1 + (q/Kq)^nq))+ku*rmq+(gam/nx)*rmq-kb*r*mq-dm*mq-lam*mq
    #num of free housekeeping mRNAs
    results[13]= (wr*a/(thetar + a))+ku*rmr+(gam/nr)*rmr-kb*r*mr-dm*mr-lam*mr
    #num of free ribosomal mRNAs
    results[14]= ku*rmr +ku*rmt+ ku*rmm+ ku*rmq +(gam/nr)*rmr+
    (gam/nr)*rmr +(gam/nx)*rmt +(gam/nx)*rmm +(gam/nx)*rmq -
    kb*r*mr- kb*r*mt- kb*r*mm -kb*r*mq +ku*nit_mrna_ribo +
    (gam/nx)*nit_mrna_ribo-kb*r*nit_mrna +ku*AA_mrna_ribo +
    (gam/nx)*AA_mrna_ribo-kb*r*AA_mrna -lam*r

    #num of free ribosomes

    #new nitrogen fixing equations. also opposite term in ATP equation
    ##########################################
    results[15]= -AA_NH4_use -lam*NH4 #fixation*nit-export_rate*NH4 -lam*NH4
    #change in num of NH4 molecules in the cell from fixation rate - export rate
    results[16]= 0#(v_nit*a/(thetax + a))+(ku*nit_mrna_ribo)+(gam/nx*nit_mrna_ribo)-(kb*r*nit_mrna)-(dm*nit_mrna)-(lam*nit_mrna)
    #num of free nitrogenase coding mRNA
    results[17]= 0#(kb*r*nit_mrna)-(ku*nit_mrna_ribo)-(gam/nx*nit_mrna_ribo)-(lam*nit_mrna_ribo)
    #num of nit mRNA-ribosome complexes
    results[18]= 0 #(gam/nx)*nit_mrna_ribo-lam*nit
    #num of nitrogenase proteins
    results[19]= 0#(export_rate*NH4)#total num of NH4 molecules exported this timestep
    results[20]=lam*N -death_rate*N #the change in num of bacterial cells, instead of just loss to dilution
    results[21]= ns*nucat-ttrate-AA_a_use-lam*a #ATP equation
    ##########################################
    ## AAA equations##
    ##########################################
    results[22]= new_AA -ttrate-lam*AA #num of AA
    results[23]=(gam/nx)*AA_mrna_ribo-lam*AA_prot #num of AA making proteins
    results[24]=((w_AA*a/(thetax + a)))+(ku*AA_mrna_ribo)+
    (gam/nx*AA_mrna_ribo)-(kb*r*AA_mrna)-(dm*AA_mrna)-(lam*AA_mrna)#num of free AA mRNAs
    results[25]=(kb*r*AA_mrna)-(ku*AA_mrna_ribo)-(gam/nx*AA_mrna_ribo)-(lam*AA_mrna_ribo)#num of AA mRNAs bound to ribosomes


end



"""ODE model of with no AA and no cell growth as well as no changes in external carbon. No nitrogen source provided.
Nitrogenase equations are removed. Original gamma equation (rate of translational elongation) that uses only ATP.
    """
function original_model(results,y,p,t)
    s_out = y[1] #external level of substrate, a proxy for plant growth state since more plant = more substrate
    rmr= y[2]#num of ribosome/ribosomal mRNA complex molecules
    em= y[3]#metabolic enzyme molecumes
    rmq= y[4]#num of ribosome/housekeeping mRNA complex molecules
    rmt= y[5]#num of ribosome/transporter protein mRNA complex molecules
    et= y[6]#num of transporter enzyme molecules
    rmm= y[7]#num of ribosome/metabolic enzyme mRNA complexes
    mt= y[8]#num of particles of transporter protein mRNA
    mm= y[9]#num of particles of metabolic enzyme mRNA
    q= y[10]#num of housekeeping protein molecules
    si= y[11]#num of substrate molecules inside the cell (internal)
    mq= y[12]#num of mRNA coding for housekeeping proteins
    mr= y[13]#num of mRNA coding for ribosomes
    r= y[14]#num of free ribosomes
    NH4 = y[15] #num of ammonia molecules fixed
    nit_mrna = y[16]#num of molecules of nitrogen fixing protein
    nit_mrna_ribo = y[17]#num of molecules of mRNA coding for nitrogen fixing protein
    nit = y[18]#num of molecules of mRNA-ribosome complexes coding for nitrogenase
    exported= y[19]#total amount of NH4 exported
    N = y[20] #num of bacteria cells
    a= y[21]#num of ATP (proxy for level of energy)
    AA = y[22] #num of amino acid molecules
    AA_prot=y[23]#num of protein molecules that make new amino acids
    AA_mrna=y[24]#num of mRNA coding for AA making protein
    AA_mrna_ribo=y[25]#num of mRNA coding for AA making protein bound to ribosomes
 
    #define more constants

    gam= gmax*a/(Kgamma + a) #gamma is the current rate of translational elongation? gmax is max rate of translational elongation and a is ATP (energy)level
    ttrate = (rmq + rmr + rmt + rmm)*gam #total translation rate (sum of the mRNA/ribosome complexes times translation rate) 
    lam= ttrate/M
    nucat= em*vm*si/(Km + si) #nucat is the num of internal substrate molecules used by catabolism. em is the num of molecules of metabolic enzyme, vm is the rate of the metabolic enzyme and si is the internal substrate level
    ##################################################################

    results[1]= 0 #((-et*vt*s_out/(Kt + s_out))*N)-ds*s_out +1e7 #change in amount of substrate, a proxy for amount of host plant
    results[2]= kb*r*mr-ku*rmr-(gam/nr)*rmr-lam*rmr #num of ribosome-mRNA complexes coding for ribosomes
    results[3]= (gam/nx)*rmm-lam*em #num of metabolic proteins
    results[4]= kb*r*mq-ku*rmq-(gam/nx)*rmq-lam*rmq #num of ribosome-mRNA complexes coding for the housekeeping proteins
    results[5]= kb*r*mt-ku*rmt-gam/nx*rmt-lam*rmt #num of ribosome-mRNA complexes coding for the transporter protein
    results[6]= (gam/nx)*rmt-lam*et #num of transporter proteins
    results[7]= kb*r*mm-ku*rmm-gam/nx*rmm-lam*rmm #num of ribosome-mRNA complexes coding for the metabolic protein
    results[8]= (we*a/(thetax + a))+ku*rmt+(gam/nx)*rmt-kb*r*mt-dm*mt-lam*mt #change in num of mRNAs coding for transporter proteins
    results[9]= (we*a/(thetax + a))+ku*rmm+(gam/nx)*rmm-kb*r*mm-dm*mm-lam*mm #change in num of mRNAs coding for metabolic proteins
    results[10]= (gam/nx)*rmq-lam*q #num of housekeeping proteins
    results[11]= (et*vt*s_out/(Kt + s_out))-nucat-lam*si #internal substrate
    results[12]= (wq*a/(thetax + a)/(1 + (q/Kq)^nq))+ku*rmq+(gam/nx)*rmq-kb*r*mq-dm*mq-lam*mq #num of free housekeeping mRNAs
    results[13]= (wr*a/(thetar + a))+ku*rmr+(gam/nr)*rmr-kb*r*mr-dm*mr-lam*mr #num of free ribosomal mRNAs
    results[14]= ku*rmr +ku*rmt+ ku*rmm+ ku*rmq +(gam/nr)*rmr+(gam/nr)*rmr +(gam/nx)*rmt +(gam/nx)*rmm +(gam/nx)*rmq -kb*r*mr- kb*r*mt- kb*r*mm -kb*r*mq -lam*r #num of free ribosomes
    results[15]= 0 #NH4 molecules
    results[16]= 0 #num of free nitrogenase coding mRNA
    results[17]= 0 #num of nit mRNA-ribosome complexes
    results[18]= 0 #num of nitrogenase proteins
    results[19]= 0 #exported NH4
    results[20]= 0 #lam*N-death_rate*N #num of bacterial cells
    results[21]= ns*nucat-ttrate-lam*a
    results[22]= 0 #num of AA
    results[23]= 0 #num of AA making proteins
    results[24]= 0 #num of free AA mRNAs
    results[25]= 0 #num of AA mRNAs bound to ribosomes
end


"""ODE model of with no AA and no cell growth as well as no changes in external carbon. No nitrogen source provided.
Nitrogenase equations are removed. Original gamma equation (rate of translational elongation) that uses only ATP.
    """
function original_model_popn_growth(results,y,p,t)
    s_out = y[1] #external level of substrate, a proxy for plant growth state since more plant = more substrate
    rmr= y[2]#num of ribosome/ribosomal mRNA complex molecules
    em= y[3]#metabolic enzyme molecumes
    rmq= y[4]#num of ribosome/housekeeping mRNA complex molecules
    rmt= y[5]#num of ribosome/transporter protein mRNA complex molecules
    et= y[6]#num of transporter enzyme molecules
    rmm= y[7]#num of ribosome/metabolic enzyme mRNA complexes
    mt= y[8]#num of particles of transporter protein mRNA
    mm= y[9]#num of particles of metabolic enzyme mRNA
    q= y[10]#num of housekeeping protein molecules
    si= y[11]#num of substrate molecules inside the cell (internal)
    mq= y[12]#num of mRNA coding for housekeeping proteins
    mr= y[13]#num of mRNA coding for ribosomes
    r= y[14]#num of free ribosomes
    NH4 = y[15] #num of ammonia molecules fixed
    nit_mrna = y[16]#num of molecules of nitrogen fixing protein
    nit_mrna_ribo = y[17]#num of molecules of mRNA coding for nitrogen fixing protein
    nit = y[18]#num of molecules of mRNA-ribosome complexes coding for nitrogenase
    exported= y[19]#total amount of NH4 exported
    N = y[20] #num of bacteria cells
    a= y[21]#num of ATP (proxy for level of energy)
    AA = y[22] #num of amino acid molecules
    AA_prot=y[23]#num of protein molecules that make new amino acids
    AA_mrna=y[24]#num of mRNA coding for AA making protein
    AA_mrna_ribo=y[25]#num of mRNA coding for AA making protein bound to ribosomes
    
    #define more constants

    gam= gmax*a/(Kgamma + a) #gamma is the current rate of translational elongation? gmax is max rate of translational elongation and a is ATP (energy)level
    ttrate = (rmq + rmr + rmt + rmm)*gam #total translation rate (sum of the mRNA/ribosome complexes times translation rate) 
    lam= ttrate/M
    nucat= em*vm*si/(Km + si) #nucat is the num of internal substrate molecules used by catabolism. em is the num of molecules of metabolic enzyme, vm is the rate of the metabolic enzyme and si is the internal substrate level
    ##################################################################

    results[1]= sub_inflow -((et*vt*s_out/(Kt + s_out))*N) -ds*s_out #change in amount of substrate, a proxy for amount of host plant
    results[2]= kb*r*mr-ku*rmr-(gam/nr)*rmr-lam*rmr #num of ribosome-mRNA complexes coding for ribosomes
    results[3]= (gam/nx)*rmm-lam*em #num of metabolic proteins
    results[4]= kb*r*mq-ku*rmq-(gam/nx)*rmq-lam*rmq #num of ribosome-mRNA complexes coding for the housekeeping proteins
    results[5]= kb*r*mt-ku*rmt-gam/nx*rmt-lam*rmt #num of ribosome-mRNA complexes coding for the transporter protein
    results[6]= (gam/nx)*rmt-lam*et #num of transporter proteins
    results[7]= kb*r*mm-ku*rmm-gam/nx*rmm-lam*rmm #num of ribosome-mRNA complexes coding for the metabolic protein
    results[8]= (we*a/(thetax + a))+ku*rmt+(gam/nx)*rmt-kb*r*mt-dm*mt-lam*mt #change in num of mRNAs coding for transporter proteins
    results[9]= (we*a/(thetax + a))+ku*rmm+(gam/nx)*rmm-kb*r*mm-dm*mm-lam*mm #change in num of mRNAs coding for metabolic proteins
    results[10]= (gam/nx)*rmq-lam*q #num of housekeeping proteins
    results[11]= (et*vt*s_out/(Kt + s_out))-nucat-lam*si #internal substrate
    results[12]= (wq*a/(thetax + a)/(1 + (q/Kq)^nq))+ku*rmq+(gam/nx)*rmq-kb*r*mq-dm*mq-lam*mq #num of free housekeeping mRNAs
    results[13]= (wr*a/(thetar + a))+ku*rmr+(gam/nr)*rmr-kb*r*mr-dm*mr-lam*mr #num of free ribosomal mRNAs
    results[14]= ku*rmr +ku*rmt+ ku*rmm+ ku*rmq +(gam/nr)*rmr+(gam/nr)*rmr +(gam/nx)*rmt +(gam/nx)*rmm +(gam/nx)*rmq -kb*r*mr- kb*r*mt- kb*r*mm -kb*r*mq -lam*r #num of free ribosomes
    results[15]= 0 #NH4 molecules
    results[16]= 0 #num of free nitrogenase coding mRNA
    results[17]= 0 #num of nit mRNA-ribosome complexes
    results[18]= 0 #num of nitrogenase proteins
    results[19]= 0 #exported NH4
    results[20]= lam*N-death_rate*N #num of bacterial cells
    results[21]= ns*nucat-ttrate-lam*a
    results[22]= 0 #num of AA
    results[23]= 0 #num of AA making proteins
    results[24]= 0 #num of free AA mRNAs
    results[25]= 0 #num of AA mRNAs bound to ribosomes
end




"""ODE model of with no AA and no cell growth as well as no changes in external carbon. No nitrogen source provided.
Nitrogenase equations are removed. Original gamma equation (rate of translational elongation) that uses only ATP.
    """
function original_model_popn_growth_no_replenish(results,y,p,t)
    s_out = y[1] #external level of substrate, a proxy for plant growth state since more plant = more substrate
    rmr= y[2]#num of ribosome/ribosomal mRNA complex molecules
    em= y[3]#metabolic enzyme molecumes
    rmq= y[4]#num of ribosome/housekeeping mRNA complex molecules
    rmt= y[5]#num of ribosome/transporter protein mRNA complex molecules
    et= y[6]#num of transporter enzyme molecules
    rmm= y[7]#num of ribosome/metabolic enzyme mRNA complexes
    mt= y[8]#num of particles of transporter protein mRNA
    mm= y[9]#num of particles of metabolic enzyme mRNA
    q= y[10]#num of housekeeping protein molecules
    si= y[11]#num of substrate molecules inside the cell (internal)
    mq= y[12]#num of mRNA coding for housekeeping proteins
    mr= y[13]#num of mRNA coding for ribosomes
    r= y[14]#num of free ribosomes
    NH4 = y[15] #num of ammonia molecules fixed
    nit_mrna = y[16]#num of molecules of nitrogen fixing protein
    nit_mrna_ribo = y[17]#num of molecules of mRNA coding for nitrogen fixing protein
    nit = y[18]#num of molecules of mRNA-ribosome complexes coding for nitrogenase
    exported= y[19]#total amount of NH4 exported
    N = y[20] #num of bacteria cells
    a= y[21]#num of ATP (proxy for level of energy)
    AA = y[22] #num of amino acid molecules
    AA_prot=y[23]#num of protein molecules that make new amino acids
    AA_mrna=y[24]#num of mRNA coding for AA making protein
    AA_mrna_ribo=y[25]#num of mRNA coding for AA making protein bound to ribosomes
    
    #define more constants

    gam= gmax*a/(Kgamma + a) #gamma is the current rate of translational elongation? gmax is max rate of translational elongation and a is ATP (energy)level
    ttrate = (rmq + rmr + rmt + rmm)*gam #total translation rate (sum of the mRNA/ribosome complexes times translation rate) 
    lam= ttrate/M
    nucat= em*vm*si/(Km + si) #nucat is the num of internal substrate molecules used by catabolism. em is the num of molecules of metabolic enzyme, vm is the rate of the metabolic enzyme and si is the internal substrate level
    ##################################################################

    results[1]= ((-et*vt*s_out/(Kt + s_out))*N) #change in amount of substrate, a proxy for amount of host plant
    results[2]= kb*r*mr-ku*rmr-(gam/nr)*rmr-lam*rmr #num of ribosome-mRNA complexes coding for ribosomes
    results[3]= (gam/nx)*rmm-lam*em #num of metabolic proteins
    results[4]= kb*r*mq-ku*rmq-(gam/nx)*rmq-lam*rmq #num of ribosome-mRNA complexes coding for the housekeeping proteins
    results[5]= kb*r*mt-ku*rmt-gam/nx*rmt-lam*rmt #num of ribosome-mRNA complexes coding for the transporter protein
    results[6]= (gam/nx)*rmt-lam*et #num of transporter proteins
    results[7]= kb*r*mm-ku*rmm-gam/nx*rmm-lam*rmm #num of ribosome-mRNA complexes coding for the metabolic protein
    results[8]= (we*a/(thetax + a))+ku*rmt+(gam/nx)*rmt-kb*r*mt-dm*mt-lam*mt #change in num of mRNAs coding for transporter proteins
    results[9]= (we*a/(thetax + a))+ku*rmm+(gam/nx)*rmm-kb*r*mm-dm*mm-lam*mm #change in num of mRNAs coding for metabolic proteins
    results[10]= (gam/nx)*rmq-lam*q #num of housekeeping proteins
    results[11]= (et*vt*s_out/(Kt + s_out))-nucat-lam*si #internal substrate
    results[12]= (wq*a/(thetax + a)/(1 + (q/Kq)^nq))+ku*rmq+(gam/nx)*rmq-kb*r*mq-dm*mq-lam*mq #num of free housekeeping mRNAs
    results[13]= (wr*a/(thetar + a))+ku*rmr+(gam/nr)*rmr-kb*r*mr-dm*mr-lam*mr #num of free ribosomal mRNAs
    results[14]= ku*rmr +ku*rmt+ ku*rmm+ ku*rmq +(gam/nr)*rmr+(gam/nr)*rmr +(gam/nx)*rmt +(gam/nx)*rmm +(gam/nx)*rmq -kb*r*mr- kb*r*mt- kb*r*mm -kb*r*mq -lam*r #num of free ribosomes
    results[15]= 0 #NH4 molecules
    results[16]= 0 #num of free nitrogenase coding mRNA
    results[17]= 0 #num of nit mRNA-ribosome complexes
    results[18]= 0 #num of nitrogenase proteins
    results[19]= 0 #exported NH4
    results[20]= lam*N #num of bacterial cells
    
    results[21]= ns*nucat-ttrate-lam*a
    results[22]= 0 #num of AA
    results[23]= 0 #num of AA making proteins
    results[24]= 0 #num of free AA mRNAs
    results[25]= 0 #num of AA mRNAs bound to ribosomes
end

"""The most basic ODE model of Nitrogen fixing no amino acids, no cell growth, and no changes in external carbon or nitrogen sources.
Nitrogenase equations are still present but not active (inital conditions = 0, dx/dt = 0). Updated gamma equation (rate of translational elongation) that uses both ATP and amino acids as substrates.
    """

function nitrogenase_simple_no_AA(results,y,p,t)
    s_out = y[1] #external level of substrate, a proxy for plant growth state since more plant = more substrate
    rmr= y[2]#num of ribosome/ribosomal mRNA complex molecules
    em= y[3]#metabolic enzyme molecumes
    rmq= y[4]#num of ribosome/housekeeping mRNA complex molecules
    rmt= y[5]#num of ribosome/transporter protein mRNA complex molecules
    et= y[6]#num of transporter enzyme molecules
    rmm= y[7]#num of ribosome/metabolic enzyme mRNA complexes
    mt= y[8]#num of particles of transporter protein mRNA
    mm= y[9]#num of particles of metabolic enzyme mRNA
    q= y[10]#num of housekeeping protein molecules
    si= y[11]#num of substrate molecules inside the cell (internal)
    mq= y[12]#num of mRNA coding for housekeeping proteins
    mr= y[13]#num of mRNA coding for ribosomes
    r= y[14]#num of free ribosomes
    NH4 = y[15] #num of ammonia molecules fixed
    nit_mrna = y[16]#num of molecules of nitrogen fixing protein
    nit_mrna_ribo = y[17]#num of molecules of mRNA coding for nitrogen fixing protein
    nit = y[18]#num of molecules of mRNA-ribosome complexes coding for nitrogenase
    exported= y[19]#total amount of NH4 exported
    N = y[20] #num of bacteria cells
    a= y[21]#num of ATP (proxy for level of energy)
    AA = y[22] #num of amino acid molecules
    AA_prot=y[23]#num of protein molecules that make new amino acids
    AA_mrna=y[24]#num of mRNA coding for AA making protein
    AA_mrna_ribo=y[25]#num of mRNA coding for AA making protein bound to ribosomes
 
    #define more constants
    gam= gmax*a/(Kgamma + a)
#     gam= (gmax*a*AA)/(k_ribo_a*k_ribo_a_AA+k_ribo_a_AA*a+k_ribo_AA_a*AA+a*AA) #two substrate gamma equation that uses both ATP and AA
    ttrate= (rmq + rmr + rmt + rmm + nit_mrna_ribo+AA_mrna_ribo)*gam #total translation rate (sum of the mRNA/ribosome complexes times translation rate)
    lam= ttrate/M
    nucat= em*vm*si/(Km + si) #nucat is the num of internal substrate molecules used by catabolism. em is the num of molecules of metabolic enzyme, vm is the rate of the metabolic enzyme and si is the internal substrate level
    fixation = vmax_nit*a/(Kgamma +a)
    AA_vo = ((k_cat_AA*2*a*NH4)/(k_a_NH4*k_a))/(1+(NH4/k_a_NH4)+(NH4/k_NH4)+2+AA*(k_x))
    new_AA = AA_vo*AA_prot

    AA_a_use = new_AA*a_per_AA #each AA produced uses x number of ATP molecules
    AA_NH4_use = new_AA*NH4_per_AA
    ##################################################################

    results[1]= 0#((-et*vt*s_out/(Kt + s_out))*N)-ds*s_out +(export_rate*NH4*exchange_rate) +sub_inflow #change in amount of substrate, a proxy for amount of host plant
    results[2]= kb*r*mr-ku*rmr-(gam/nr)*rmr-lam*rmr#-dp*r
    #num of ribosome-mRNA complexes coding for ribosomes
    results[3]= (gam/nx)*rmm-lam*em#-dp*em
    #num of metabolic proteins
    results[4]= kb*r*mq-ku*rmq-(gam/nx)*rmq-lam*rmq
    #num of ribosome-mRNA complexes coding for the housekeeping proteins
    results[5]= kb*r*mt-ku*rmt-gam/nx*rmt-lam*rmt
    #num of ribosome-mRNA complexes coding for the transporter protein
    results[6]= (gam/nx)*rmt-lam*et#-dp*et
    #num of transporter proteins
    results[7]= kb*r*mm-ku*rmm-gam/nx*rmm-lam*rmm
    #num of ribosome-mRNA complexes coding for the metabolic protein
    results[8]= (we*a/(thetax + a))+ku*rmt+(gam/nx)*rmt-kb*r*mt-dm*mt-lam*mt
    #change in num of mRNAs coding for transporter proteins
    results[9]= (we*a/(thetax + a))+ku*rmm+(gam/nx)*rmm-kb*r*mm-dm*mm-lam*mm
    #change in num of mRNAs coding for metabolic proteins
    results[10]= (gam/nx)*rmq-lam*q#-dp*q
    #num of housekeeping proteins
    results[11]= (et*vt*s_out/(Kt + s_out))-nucat-lam*si#-degredation*si
    results[12]= (wq*a/(thetax + a)/(1 + (q/Kq)^nq))+ku*rmq+(gam/nx)*rmq-kb*r*mq-dm*mq-lam*mq
    #num of free housekeeping mRNAs
    results[13]= (wr*a/(thetar + a))+ku*rmr+(gam/nr)*rmr-kb*r*mr-dm*mr-lam*mr
    #num of free ribosomal mRNAs
    results[14]= ku*rmr +ku*rmt+ ku*rmm+ ku*rmq +(gam/nr)*rmr+
    (gam/nr)*rmr +(gam/nx)*rmt +(gam/nx)*rmm +(gam/nx)*rmq -
    kb*r*mr- kb*r*mt- kb*r*mm -kb*r*mq +ku*nit_mrna_ribo +
    (gam/nx)*nit_mrna_ribo-kb*r*nit_mrna +ku*AA_mrna_ribo +
    (gam/nx)*AA_mrna_ribo-kb*r*AA_mrna -lam*r
    #num of free ribosomes
    results[15]= 0 #NH4_inflow +fixation*nit -export_rate*NH4 -lam*NH4
    #change in num of NH4 molecules in the cell from fixation rate - export rate
    results[16]= (v_nit*a/(thetax + a))+(ku*nit_mrna_ribo)+(gam/nx*nit_mrna_ribo)-(kb*r*nit_mrna)-(dm*nit_mrna)-(lam*nit_mrna)
    #num of free nitrogenase coding mRNA
    results[17]= (kb*r*nit_mrna)-(ku*nit_mrna_ribo)-(gam/nx*nit_mrna_ribo)-(lam*nit_mrna_ribo)
    #num of nit mRNA-ribosome complexes
    results[18]= 0#(gam/nx)*nit_mrna_ribo-lam*nit
    #num of nitrogenase proteins
    results[19]= 0#(export_rate*NH4)#total num of NH4 molecules exported this timestep
    results[20]= 0#lam*N-death_rate*N #the change in num of bacterial cells, instead of just loss to dilution
    results[21]= ns*nucat-ttrate-fixation*nit-lam*a
    results[22]= 0 #num of AA
    results[23]= 0 #num of AA making proteins
    results[24]= 0 #num of free AA mRNAs
    results[25]= 0 #num of AA mRNAs bound to ribosomes
end


"""ODE model of nitrogen fixing no Amino Acid metabolism with cell growth and changes allowed in external carbon or nitrogen sources.
Nitrogenase equations are still present but not active (inital conditions = 0, dx/dt = 0). Updated gamma equation (rate of translational elongation) that uses both ATP and amino acids as substrates.
    """
function nitrogenase_no_AA_second_burn_in(results,y,p,t)
    s_out = y[1] #external level of substrate, a proxy for plant growth state since more plant = more substrate
    rmr= y[2]#num of ribosome/ribosomal mRNA complex molecules
    em= y[3]#metabolic enzyme molecumes
    rmq= y[4]#num of ribosome/housekeeping mRNA complex molecules
    rmt= y[5]#num of ribosome/transporter protein mRNA complex molecules
    et= y[6]#num of transporter enzyme molecules
    rmm= y[7]#num of ribosome/metabolic enzyme mRNA complexes
    mt= y[8]#num of particles of transporter protein mRNA
    mm= y[9]#num of particles of metabolic enzyme mRNA
    q= y[10]#num of housekeeping protein molecules
    si= y[11]#num of substrate molecules inside the cell (internal)
    mq= y[12]#num of mRNA coding for housekeeping proteins
    mr= y[13]#num of mRNA coding for ribosomes
    r= y[14]#num of free ribosomes
    NH4 = y[15] #num of ammonia molecules fixed
    nit_mrna = y[16]#num of molecules of nitrogen fixing protein
    nit_mrna_ribo = y[17]#num of molecules of mRNA coding for nitrogen fixing protein
    nit = y[18]#num of molecules of mRNA-ribosome complexes coding for nitrogenase
    exported= y[19]#total amount of NH4 exported
    N = y[20] #num of bacteria cells
    a= y[21]#num of ATP (proxy for level of energy)
    AA = y[22] #num of amino acid molecules
    AA_prot=y[23]#num of protein molecules that make new amino acids
    AA_mrna=y[24]#num of mRNA coding for AA making protein
    AA_mrna_ribo=y[25]#num of mRNA coding for AA making protein bound to ribosomes
 
    #define more constants
    gam= gmax*a/(Kgamma + a)
#     gam= (gmax*a*AA)/(k_ribo_a*k_ribo_a_AA+k_ribo_a_AA*a+k_ribo_AA_a*AA+a*AA) #two substrate gamma equation that uses both ATP and AA
    ttrate= (rmq + rmr + rmt + rmm + nit_mrna_ribo+AA_mrna_ribo)*gam #total translation rate (sum of the mRNA/ribosome complexes times translation rate)
    lam= ttrate/M
    nucat= em*vm*si/(Km + si) #nucat is the num of internal substrate molecules used by catabolism. em is the num of molecules of metabolic enzyme, vm is the rate of the metabolic enzyme and si is the internal substrate level
    fixation = vmax_nit*a/(Kgamma +a)
    AA_vo = ((k_cat_AA*2*a*NH4)/(k_a_NH4*k_a))/(1+(NH4/k_a_NH4)+(NH4/k_NH4)+2+AA*(k_x))
    new_AA = AA_vo*AA_prot

    AA_a_use = new_AA*a_per_AA #each AA produced uses x number of ATP molecules
    AA_NH4_use = new_AA*NH4_per_AA
    
    ##################################################################

    results[1]= 0#((-et*vt*s_out/(Kt + s_out))*N)-ds*s_out +(export_rate*NH4*exchange_rate) +sub_inflow #change in amount of substrate, a proxy for amount of host plant
    results[2]= kb*r*mr-ku*rmr-(gam/nr)*rmr-lam*rmr#-dp*r
    #num of ribosome-mRNA complexes coding for ribosomes
    results[3]= (gam/nx)*rmm-lam*em#-dp*em
    #num of metabolic proteins
    results[4]= kb*r*mq-ku*rmq-(gam/nx)*rmq-lam*rmq
    #num of ribosome-mRNA complexes coding for the housekeeping proteins
    results[5]= kb*r*mt-ku*rmt-gam/nx*rmt-lam*rmt
    #num of ribosome-mRNA complexes coding for the transporter protein
    results[6]= (gam/nx)*rmt-lam*et#-dp*et
    #num of transporter proteins
    results[7]= kb*r*mm-ku*rmm-gam/nx*rmm-lam*rmm
    #num of ribosome-mRNA complexes coding for the metabolic protein
    results[8]= (we*a/(thetax + a))+ku*rmt+(gam/nx)*rmt-kb*r*mt-dm*mt-lam*mt
    #change in num of mRNAs coding for transporter proteins
    results[9]= (we*a/(thetax + a))+ku*rmm+(gam/nx)*rmm-kb*r*mm-dm*mm-lam*mm
    #change in num of mRNAs coding for metabolic proteins
    results[10]= (gam/nx)*rmq-lam*q#-dp*q
    #num of housekeeping proteins
    results[11]= (et*vt*s_out/(Kt + s_out))-nucat-lam*si#-degredation*si
    results[12]= (wq*a/(thetax + a)/(1 + (q/Kq)^nq))+ku*rmq+(gam/nx)*rmq-kb*r*mq-dm*mq-lam*mq
    #num of free housekeeping mRNAs
    results[13]= (wr*a/(thetar + a))+ku*rmr+(gam/nr)*rmr-kb*r*mr-dm*mr-lam*mr
    #num of free ribosomal mRNAs
    results[14]= ku*rmr +ku*rmt+ ku*rmm+ ku*rmq +(gam/nr)*rmr+
    (gam/nr)*rmr +(gam/nx)*rmt +(gam/nx)*rmm +(gam/nx)*rmq -
    kb*r*mr- kb*r*mt- kb*r*mm -kb*r*mq +ku*nit_mrna_ribo +
    (gam/nx)*nit_mrna_ribo-kb*r*nit_mrna +ku*AA_mrna_ribo +
    (gam/nx)*AA_mrna_ribo-kb*r*AA_mrna -lam*r
    #num of free ribosomes
    results[15]= 0 #NH4_inflow +fixation*nit -export_rate*NH4 -lam*NH4
    #change in num of NH4 molecules in the cell from fixation rate - export rate
    results[16]= (v_nit*a/(thetax + a))+(ku*nit_mrna_ribo)+(gam/nx*nit_mrna_ribo)-(kb*r*nit_mrna)-(dm*nit_mrna)-(lam*nit_mrna)
    #num of free nitrogenase coding mRNA
    results[17]= (kb*r*nit_mrna)-(ku*nit_mrna_ribo)-(gam/nx*nit_mrna_ribo)-(lam*nit_mrna_ribo)
    #num of nit mRNA-ribosome complexes
    results[18]= (gam/nx)*nit_mrna_ribo-lam*nit
    #num of nitrogenase proteins
    results[19]= 0#(export_rate*NH4)#total num of NH4 molecules exported this timestep
    results[20]= 0 #lam*N-death_rate*N #the change in num of bacterial cells, instead of just loss to dilution
    results[21]= ns*nucat-ttrate-fixation*nit-AA_a_use-lam*a
    results[22]= 0 #num of AA
    results[23]= 0 #num of AA making proteins
    results[24]= 0 #num of free AA mRNAs
    results[25]= 0 #num of AA mRNAs bound to ribosomes
end

"""ODE model of nitrogen fixing no Amino Acid metabolism with cell growth and changes allowed in external carbon or nitrogen sources.
Nitrogenase equations are still present but not active (inital conditions = 0, dx/dt = 0). Updated gamma equation (rate of translational elongation) that uses both ATP and amino acids as substrates.
    """
function nitrogenase_no_AA_with_popn_growth(results,y,p,t)
    s_out = y[1] #external level of substrate, a proxy for plant growth state since more plant = more substrate
    rmr= y[2]#num of ribosome/ribosomal mRNA complex molecules
    em= y[3]#metabolic enzyme molecumes
    rmq= y[4]#num of ribosome/housekeeping mRNA complex molecules
    rmt= y[5]#num of ribosome/transporter protein mRNA complex molecules
    et= y[6]#num of transporter enzyme molecules
    rmm= y[7]#num of ribosome/metabolic enzyme mRNA complexes
    mt= y[8]#num of particles of transporter protein mRNA
    mm= y[9]#num of particles of metabolic enzyme mRNA
    q= y[10]#num of housekeeping protein molecules
    si= y[11]#num of substrate molecules inside the cell (internal)
    mq= y[12]#num of mRNA coding for housekeeping proteins
    mr= y[13]#num of mRNA coding for ribosomes
    r= y[14]#num of free ribosomes
    NH4 = y[15] #num of ammonia molecules fixed
    nit_mrna = y[16]#num of molecules of nitrogen fixing protein
    nit_mrna_ribo = y[17]#num of molecules of mRNA coding for nitrogen fixing protein
    nit = y[18]#num of molecules of mRNA-ribosome complexes coding for nitrogenase
    exported= y[19]#total amount of NH4 exported
    N = y[20] #num of bacteria cells
    a= y[21]#num of ATP (proxy for level of energy)
    AA = y[22] #num of amino acid molecules
    AA_prot=y[23]#num of protein molecules that make new amino acids
    AA_mrna=y[24]#num of mRNA coding for AA making protein
    AA_mrna_ribo=y[25]#num of mRNA coding for AA making protein bound to ribosomes
 
    #define more constants
    gam= gmax*a/(Kgamma + a)
#     gam= (gmax*a*AA)/(k_ribo_a*k_ribo_a_AA+k_ribo_a_AA*a+k_ribo_AA_a*AA+a*AA) #two substrate gamma equation that uses both ATP and AA
    ttrate= (rmq + rmr + rmt + rmm + nit_mrna_ribo+AA_mrna_ribo)*gam #total translation rate (sum of the mRNA/ribosome complexes times translation rate)
    lam= ttrate/M
    nucat= em*vm*si/(Km + si) #nucat is the num of internal substrate molecules used by catabolism. em is the num of molecules of metabolic enzyme, vm is the rate of the metabolic enzyme and si is the internal substrate level
    fixation = vmax_nit*a/(Kgamma +a)
    AA_vo = ((k_cat_AA*2*a*NH4)/(k_a_NH4*k_a))/(1+(NH4/k_a_NH4)+(NH4/k_NH4)+2+AA*(k_x))
    new_AA = AA_vo*AA_prot

    AA_a_use = new_AA*a_per_AA #each AA produced uses x number of ATP molecules
    AA_NH4_use = new_AA*NH4_per_AA
    ##################################################################

    results[1]= ((-et*vt*s_out/(Kt + s_out))*N)-ds*s_out +(export_rate*NH4*exchange_rate)+sub_inflow #change in amount of substrate, a proxy for amount of host plant
    results[2]= kb*r*mr-ku*rmr-(gam/nr)*rmr-lam*rmr#-dp*r
    #num of ribosome-mRNA complexes coding for ribosomes
    results[3]= (gam/nx)*rmm-lam*em#-dp*em
    #num of metabolic proteins
    results[4]= kb*r*mq-ku*rmq-(gam/nx)*rmq-lam*rmq
    #num of ribosome-mRNA complexes coding for the housekeeping proteins
    results[5]= kb*r*mt-ku*rmt-gam/nx*rmt-lam*rmt
    #num of ribosome-mRNA complexes coding for the transporter protein
    results[6]= (gam/nx)*rmt-lam*et#-dp*et
    #num of transporter proteins
    results[7]= kb*r*mm-ku*rmm-gam/nx*rmm-lam*rmm
    #num of ribosome-mRNA complexes coding for the metabolic protein
    results[8]= (we*a/(thetax + a))+ku*rmt+(gam/nx)*rmt-kb*r*mt-dm*mt-lam*mt
    #change in num of mRNAs coding for transporter proteins
    results[9]= (we*a/(thetax + a))+ku*rmm+(gam/nx)*rmm-kb*r*mm-dm*mm-lam*mm
    #change in num of mRNAs coding for metabolic proteins
    results[10]= (gam/nx)*rmq-lam*q#-dp*q
    #num of housekeeping proteins
    results[11]= (et*vt*s_out/(Kt + s_out))-nucat-lam*si#-degredation*si
    results[12]= (wq*a/(thetax + a)/(1 + (q/Kq)^nq))+ku*rmq+(gam/nx)*rmq-kb*r*mq-dm*mq-lam*mq
    #num of free housekeeping mRNAs
    results[13]= (wr*a/(thetar + a))+ku*rmr+(gam/nr)*rmr-kb*r*mr-dm*mr-lam*mr
    #num of free ribosomal mRNAs
    results[14]= ku*rmr +ku*rmt+ ku*rmm+ ku*rmq +(gam/nr)*rmr+
    (gam/nr)*rmr +(gam/nx)*rmt +(gam/nx)*rmm +(gam/nx)*rmq -
    kb*r*mr- kb*r*mt- kb*r*mm -kb*r*mq +ku*nit_mrna_ribo +
    (gam/nx)*nit_mrna_ribo-kb*r*nit_mrna +ku*AA_mrna_ribo +
    (gam/nx)*AA_mrna_ribo-kb*r*AA_mrna -lam*r
    #num of free ribosomes
    results[15]= NH4_inflow +fixation*nit -export_rate*NH4 -lam*NH4
    #change in num of NH4 molecules in the cell from fixation rate - export rate
    results[16]= (v_nit*a/(thetax + a))+(ku*nit_mrna_ribo)+(gam/nx*nit_mrna_ribo)-(kb*r*nit_mrna)-(dm*nit_mrna)-(lam*nit_mrna)
    #num of free nitrogenase coding mRNA
    results[17]= (kb*r*nit_mrna)-(ku*nit_mrna_ribo)-(gam/nx*nit_mrna_ribo)-(lam*nit_mrna_ribo)
    #num of nit mRNA-ribosome complexes
    results[18]= (gam/nx)*nit_mrna_ribo-lam*nit
    #num of nitrogenase proteins
    results[19]= (export_rate*NH4)#total num of NH4 molecules exported this timestep
    results[20]=lam*N-death_rate*N #the change in num of bacterial cells, instead of just loss to dilution
    results[21]= ns*nucat-ttrate-fixation*nit-lam*a
    results[22]= 0 #num of AA
    results[23]= 0 #num of AA making proteins
    results[24]= 0 #num of free AA mRNAs
    results[25]= 0 #num of AA mRNAs bound to ribosomes
end


"""ODE model of nitrogen fixing no Amino Acid metabolism with cell growth and changes allowed in external carbon or nitrogen sources.
Nitrogenase equations are still present but not active (inital conditions = 0, dx/dt = 0). Updated gamma equation (rate of translational elongation) that uses both ATP and amino acids as substrates.
    """
function nitrogenase_no_AA_with_popn_growth_run_out(results,y,p,t)
    s_out = y[1] #external level of substrate, a proxy for plant growth state since more plant = more substrate
    rmr= y[2]#num of ribosome/ribosomal mRNA complex molecules
    em= y[3]#metabolic enzyme molecumes
    rmq= y[4]#num of ribosome/housekeeping mRNA complex molecules
    rmt= y[5]#num of ribosome/transporter protein mRNA complex molecules
    et= y[6]#num of transporter enzyme molecules
    rmm= y[7]#num of ribosome/metabolic enzyme mRNA complexes
    mt= y[8]#num of particles of transporter protein mRNA
    mm= y[9]#num of particles of metabolic enzyme mRNA
    q= y[10]#num of housekeeping protein molecules
    si= y[11]#num of substrate molecules inside the cell (internal)
    mq= y[12]#num of mRNA coding for housekeeping proteins
    mr= y[13]#num of mRNA coding for ribosomes
    r= y[14]#num of free ribosomes
    NH4 = y[15] #num of ammonia molecules fixed
    nit_mrna = y[16]#num of molecules of nitrogen fixing protein
    nit_mrna_ribo = y[17]#num of molecules of mRNA coding for nitrogen fixing protein
    nit = y[18]#num of molecules of mRNA-ribosome complexes coding for nitrogenase
    exported= y[19]#total amount of NH4 exported
    N = y[20] #num of bacteria cells
    a= y[21]#num of ATP (proxy for level of energy)
    AA = y[22] #num of amino acid molecules
    AA_prot=y[23]#num of protein molecules that make new amino acids
    AA_mrna=y[24]#num of mRNA coding for AA making protein
    AA_mrna_ribo=y[25]#num of mRNA coding for AA making protein bound to ribosomes
 
    #define more constants
    gam= gmax*a/(Kgamma + a)
#     gam= (gmax*a*AA)/(k_ribo_a*k_ribo_a_AA+k_ribo_a_AA*a+k_ribo_AA_a*AA+a*AA) #two substrate gamma equation that uses both ATP and AA
    ttrate= (rmq + rmr + rmt + rmm + nit_mrna_ribo+AA_mrna_ribo)*gam #total translation rate (sum of the mRNA/ribosome complexes times translation rate)
    lam= ttrate/M
    nucat= em*vm*si/(Km + si) #nucat is the num of internal substrate molecules used by catabolism. em is the num of molecules of metabolic enzyme, vm is the rate of the metabolic enzyme and si is the internal substrate level
    fixation = vmax_nit*a/(Kgamma +a)
    AA_vo = ((k_cat_AA*2*a*NH4)/(k_a_NH4*k_a))/(1+(NH4/k_a_NH4)+(NH4/k_NH4)+2+AA*(k_x))
    new_AA = AA_vo*AA_prot

    AA_a_use = new_AA*a_per_AA #each AA produced uses x number of ATP molecules
    AA_NH4_use = new_AA*NH4_per_AA
    ##################################################################

    results[1]= ((-et*vt*s_out/(Kt + s_out))*N)-ds*s_out +(export_rate*NH4*exchange_rate) #+sub_inflow #change in amount of substrate, a proxy for amount of host plant
    results[2]= kb*r*mr-ku*rmr-(gam/nr)*rmr-lam*rmr#-dp*r
    #num of ribosome-mRNA complexes coding for ribosomes
    results[3]= (gam/nx)*rmm-lam*em#-dp*em
    #num of metabolic proteins
    results[4]= kb*r*mq-ku*rmq-(gam/nx)*rmq-lam*rmq
    #num of ribosome-mRNA complexes coding for the housekeeping proteins
    results[5]= kb*r*mt-ku*rmt-gam/nx*rmt-lam*rmt
    #num of ribosome-mRNA complexes coding for the transporter protein
    results[6]= (gam/nx)*rmt-lam*et#-dp*et
    #num of transporter proteins
    results[7]= kb*r*mm-ku*rmm-gam/nx*rmm-lam*rmm
    #num of ribosome-mRNA complexes coding for the metabolic protein
    results[8]= (we*a/(thetax + a))+ku*rmt+(gam/nx)*rmt-kb*r*mt-dm*mt-lam*mt
    #change in num of mRNAs coding for transporter proteins
    results[9]= (we*a/(thetax + a))+ku*rmm+(gam/nx)*rmm-kb*r*mm-dm*mm-lam*mm
    #change in num of mRNAs coding for metabolic proteins
    results[10]= (gam/nx)*rmq-lam*q#-dp*q
    #num of housekeeping proteins
    results[11]= (et*vt*s_out/(Kt + s_out))-nucat-lam*si#-degredation*si
    results[12]= (wq*a/(thetax + a)/(1 + (q/Kq)^nq))+ku*rmq+(gam/nx)*rmq-kb*r*mq-dm*mq-lam*mq
    #num of free housekeeping mRNAs
    results[13]= (wr*a/(thetar + a))+ku*rmr+(gam/nr)*rmr-kb*r*mr-dm*mr-lam*mr
    #num of free ribosomal mRNAs
    results[14]= ku*rmr +ku*rmt+ ku*rmm+ ku*rmq +(gam/nr)*rmr+
    (gam/nr)*rmr +(gam/nx)*rmt +(gam/nx)*rmm +(gam/nx)*rmq -
    kb*r*mr- kb*r*mt- kb*r*mm -kb*r*mq +ku*nit_mrna_ribo +
    (gam/nx)*nit_mrna_ribo-kb*r*nit_mrna +ku*AA_mrna_ribo +
    (gam/nx)*AA_mrna_ribo-kb*r*AA_mrna -lam*r
    #num of free ribosomes
    results[15]= fixation*nit -export_rate*NH4 -lam*NH4
    #change in num of NH4 molecules in the cell from fixation rate - export rate
    results[16]= (v_nit*a/(thetax + a))+(ku*nit_mrna_ribo)+(gam/nx*nit_mrna_ribo)-(kb*r*nit_mrna)-(dm*nit_mrna)-(lam*nit_mrna)
    #num of free nitrogenase coding mRNA
    results[17]= (kb*r*nit_mrna)-(ku*nit_mrna_ribo)-(gam/nx*nit_mrna_ribo)-(lam*nit_mrna_ribo)
    #num of nit mRNA-ribosome complexes
    results[18]= (gam/nx)*nit_mrna_ribo-lam*nit
    #num of nitrogenase proteins
    results[19]= (export_rate*NH4)#total num of NH4 molecules exported this timestep
    results[20]=lam*N-death_rate*N #the change in num of bacterial cells, instead of just loss to dilution
    results[21]= ns*nucat-ttrate-fixation*nit-lam*a
    results[22]= 0 #num of AA
    results[23]= 0 #num of AA making proteins
    results[24]= 0 #num of free AA mRNAs
    results[25]= 0 #num of AA mRNAs bound to ribosomes
end


"""The most basic ODE model of Nitrogen fixing with amino acids with no cell growth and no changes in external carbon or nitrogen sources.
    """
function nitrogenase_simple_with_AA(results,y,p,t)
    s_out = y[1] #external level of substrate, a proxy for plant growth state since more plant = more substrate
    rmr= y[2]#num of ribosome/ribosomal mRNA complex molecules
    em= y[3]#metabolic enzyme molecumes
    rmq= y[4]#num of ribosome/housekeeping mRNA complex molecules
    rmt= y[5]#num of ribosome/transporter protein mRNA complex molecules
    et= y[6]#num of transporter enzyme molecules
    rmm= y[7]#num of ribosome/metabolic enzyme mRNA complexes
    mt= y[8]#num of particles of transporter protein mRNA
    mm= y[9]#num of particles of metabolic enzyme mRNA
    q= y[10]#num of housekeeping protein molecules
    si= y[11]#num of substrate molecules inside the cell (internal)
    mq= y[12]#num of mRNA coding for housekeeping proteins
    mr= y[13]#num of mRNA coding for ribosomes
    r= y[14]#num of free ribosomes
    NH4 = y[15] #num of ammonia molecules internally
    nit_mrna = y[16]#num of molecules of mRNA coding for nitrogen fixing protein
    nit_mrna_ribo = y[17]#num of molecules of mRNA-ribosome complexes coding for nitrogenase
    nit = y[18] #num of molecules of nitrogen fixing protein
    exported= y[19]#total amount of NH4 exported
    N = y[20] #num of bacteria cells
    a= y[21]#num of ATP (proxy for level of energy)
    AA = y[22] #num of amino acid molecules
    AA_prot=y[23]#num of protein molecules that make new amino acids
    AA_mrna=y[24]#num of mRNA coding for AA making protein
    AA_mrna_ribo=y[25]#num of mRNA coding for AA making protein bound to ribosomes

    #define more constants

    #gam= gmax*a/(Kgamma + a) #gamma is the current rate of translational elongation? gmax is max rate of translational elongation and a is ATP (energy)level

    gam= (gmax*a*AA)/(k_ribo_a*k_ribo_a_AA+k_ribo_a_AA*a+k_ribo_AA_a*AA+a*AA) #updated gamma equation that uses both ATP and AA

    ttrate= (rmq + rmr + rmt + rmm + nit_mrna_ribo+AA_mrna_ribo)*gam #total translation rate (sum of the mRNA/ribosome complexes times translation rate)
    #lam= ttrate/(M-AA) #lambda, the growth rate/dilution rate, the ratio of total translation rate to total cell mass
    lam= ttrate/M
    #if lam<0.00000001 lam=0.00000001 end
    nucat= em*vm*si/(Km + si) #nucat is the num of internal substrate molecules used by catabolism. em is the num of molecules of metabolic enzyme, vm is the rate of the metabolic enzyme and si is the internal substrate level
    #this is all over the half maximal metabolic rate threshold plus the internal substrate level
    fixation = vmax_nit*a/(Kgamma +a)

    AA_vo = ((k_cat_AA*2*a*NH4)/(k_a_NH4*k_a))/(1+(NH4/k_a_NH4)+(NH4/k_NH4)+2+AA*(k_x))
#     new_AA = AA_vo
    new_AA = AA_vo*AA_prot

    AA_a_use = new_AA*a_per_AA #each AA produced uses x number of ATP molecules
    AA_NH4_use = new_AA*NH4_per_AA

    ##################################################################

    results[1]= 0#((-et*vt*s_out/(Kt + s_out))*N)-ds*s_out +(export_rate*NH4*exchange_rate) +sub_inflow #change in amount of substrate, a proxy for amount of host plant
    results[2]= kb*r*mr-ku*rmr-(gam/nr)*rmr-lam*rmr#-dp*r
    #num of ribosome-mRNA complexes coding for ribosomes
    results[3]= (gam/nx)*rmm-lam*em#-dp*em
    #num of metabolic proteins
    results[4]= kb*r*mq-ku*rmq-(gam/nx)*rmq-lam*rmq
    #num of ribosome-mRNA complexes coding for the housekeeping proteins
    results[5]= kb*r*mt-ku*rmt-gam/nx*rmt-lam*rmt
    #num of ribosome-mRNA complexes coding for the transporter protein
    results[6]= (gam/nx)*rmt-lam*et#-dp*et
    #num of transporter proteins
    results[7]= kb*r*mm-ku*rmm-gam/nx*rmm-lam*rmm
    #num of ribosome-mRNA complexes coding for the metabolic protein
    results[8]= (we*a/(thetax + a))+ku*rmt+(gam/nx)*rmt-kb*r*mt-dm*mt-lam*mt
    #change in num of mRNAs coding for transporter proteins
    results[9]= (we*a/(thetax + a))+ku*rmm+(gam/nx)*rmm-kb*r*mm-dm*mm-lam*mm
    #change in num of mRNAs coding for metabolic proteins
    results[10]= (gam/nx)*rmq-lam*q#-dp*q
    #num of housekeeping proteins
    results[11]= (et*vt*s_out/(Kt + s_out))-nucat-lam*si#-degredation*si
    results[12]= (wq*a/(thetax + a)/(1 + (q/Kq)^nq))+ku*rmq+(gam/nx)*rmq-kb*r*mq-dm*mq-lam*mq
    #num of free housekeeping mRNAs
    results[13]= (wr*a/(thetar + a))+ku*rmr+(gam/nr)*rmr-kb*r*mr-dm*mr-lam*mr
    #num of free ribosomal mRNAs
    results[14]= ku*rmr +ku*rmt+ ku*rmm+ ku*rmq +(gam/nr)*rmr+
    (gam/nr)*rmr +(gam/nx)*rmt +(gam/nx)*rmm +(gam/nx)*rmq -
    kb*r*mr- kb*r*mt- kb*r*mm -kb*r*mq +ku*nit_mrna_ribo +
    (gam/nx)*nit_mrna_ribo-kb*r*nit_mrna +ku*AA_mrna_ribo +
    (gam/nx)*AA_mrna_ribo-kb*r*AA_mrna -lam*r
    #num of free ribosomes

    #new nitrogen fixing equations. also opposite term in ATP equation
    ##########################################
    results[15]= 0 #NH4_inflow +fixation*nit -export_rate*NH4 -AA_NH4_use -lam*NH4
    #change in num of NH4 molecules in the cell from inflow + fixation - export
    results[16]= (v_nit*a/(thetax + a))+(ku*nit_mrna_ribo)+(gam/nx*nit_mrna_ribo)-(kb*r*nit_mrna)-(dm*nit_mrna)-(lam*nit_mrna)
    #num of free nitrogenase coding mRNA
    results[17]= (kb*r*nit_mrna)-(ku*nit_mrna_ribo)-(gam/nx*nit_mrna_ribo)-(lam*nit_mrna_ribo)
    #num of nit mRNA-ribosome complexes
    results[18]= 0 #(gam/nx)*nit_mrna_ribo-lam*nit
    #num of nitrogenase proteins
    results[19]= 0#(export_rate*NH4)#total num of NH4 molecules exported this timestep
    results[20]= 0 #lam*N-death_rate*N #the change in num of bacterial cells, instead of just loss to dilution
    results[21]= ns*nucat-ttrate-fixation*nit-AA_a_use-lam*a
    ##########################################
    ## AAA equations##
    ##########################################
    results[22]= 0 #new_AA -ttrate-lam*AA #num of AA
    results[23]=(gam/nx)*AA_mrna_ribo-lam*AA_prot #num of AA making proteins
    results[24]=((w_AA*a/(thetax + a)))+(ku*AA_mrna_ribo)+
    (gam/nx*AA_mrna_ribo)-(kb*r*AA_mrna)-(dm*AA_mrna)-(lam*AA_mrna)#num of free AA mRNAs
    results[25]=(kb*r*AA_mrna)-(ku*AA_mrna_ribo)-(gam/nx*AA_mrna_ribo)-(lam*AA_mrna_ribo)#num of AA mRNAs bound to ribosomes


end


"""ODE model of nitrogen fixing with Amino Acid metabolism with cell growth and changes allowed in external carbon or nitrogen sources.
Nitrogenase equations are still present but not active (inital conditions = 0, dx/dt = 0). Updated gamma equation (rate of translational elongation) that uses both ATP and amino acids as substrates.
    """
function nitrogenase_with_AA_second_burn_in(results,y,p,t)
    s_out = y[1] #external level of substrate, a proxy for plant growth state since more plant = more substrate
    rmr= y[2]#num of ribosome/ribosomal mRNA complex molecules
    em= y[3]#metabolic enzyme molecumes
    rmq= y[4]#num of ribosome/housekeeping mRNA complex molecules
    rmt= y[5]#num of ribosome/transporter protein mRNA complex molecules
    et= y[6]#num of transporter enzyme molecules
    rmm= y[7]#num of ribosome/metabolic enzyme mRNA complexes
    mt= y[8]#num of particles of transporter protein mRNA
    mm= y[9]#num of particles of metabolic enzyme mRNA
    q= y[10]#num of housekeeping protein molecules
    si= y[11]#num of substrate molecules inside the cell (internal)
    mq= y[12]#num of mRNA coding for housekeeping proteins
    mr= y[13]#num of mRNA coding for ribosomes
    r= y[14]#num of free ribosomes
    NH4 = y[15] #num of ammonia molecules fixed
    nit_mrna = y[16]#num of molecules of nitrogen fixing protein
    nit_mrna_ribo = y[17]#num of molecules of mRNA coding for nitrogen fixing protein
    nit = y[18]#num of molecules of mRNA-ribosome complexes coding for nitrogenase
    exported= y[19]#total amount of NH4 exported
    N = y[20] #num of bacteria cells
    a= y[21]#num of ATP (proxy for level of energy)
    AA = y[22] #num of amino acid molecules
    AA_prot=y[23]#num of protein molecules that make new amino acids
    AA_mrna=y[24]#num of mRNA coding for AA making protein
    AA_mrna_ribo=y[25]#num of mRNA coding for AA making protein bound to ribosomes
    #define more constants

    #gam= gmax*a/(Kgamma + a) #gamma is the current rate of translational elongation? gmax is max rate of translational elongation and a is ATP (energy)level
    gam= (gmax*a*AA)/(k_ribo_a*k_ribo_a_AA+k_ribo_a_AA*a+k_ribo_AA_a*AA+a*AA) #two substrate gamma equation that uses both ATP and AA
    ttrate= (rmq + rmr + rmt + rmm + nit_mrna_ribo+AA_mrna_ribo)*gam #total translation rate (sum of the mRNA/ribosome complexes times translation rate)
    #lam= ttrate/(M-AA) #lambda, the growth rate/dilution rate, the ratio of total translation rate to total cell mass
    lam= ttrate/M
    nucat= em*vm*si/(Km + si) #nucat is the num of internal substrate molecules used by catabolism. em is the num of molecules of metabolic enzyme, vm is the rate of the metabolic enzyme and si is the internal substrate level
    #this is all over the half maximal metabolic rate threshold plus the internal substrate level
    fixation = vmax_nit*a/(Kgamma +a)

    AA_vo = ((k_cat_AA*2*a*NH4)/(k_a_NH4*k_a))/(1+(NH4/k_a_NH4)+(NH4/k_NH4)+2+AA*(k_x))
    new_AA = AA_vo*AA_prot

    AA_a_use = new_AA*a_per_AA #each AA produced uses x number of ATP molecules
    AA_NH4_use = new_AA*NH4_per_AA

    ################################################
    #export rate needs to vary with AA levels#
    ################################################

    ##################################################################

    results[1]= 0#((-et*vt*s_out/(Kt + s_out))*N)-ds*s_out +sub_inflow #+(export_rate*NH4*exchange_rate)  #change in amount of substrate, a proxy for amount of host plant
    results[2]= kb*r*mr-ku*rmr-(gam/nr)*rmr-lam*rmr#-dp*r
    #num of ribosome-mRNA complexes coding for ribosomes
    results[3]= (gam/nx)*rmm-lam*em#-dp*em
    #num of metabolic proteins
    results[4]= kb*r*mq-ku*rmq-(gam/nx)*rmq-lam*rmq
    #num of ribosome-mRNA complexes coding for the housekeeping proteins
    results[5]= kb*r*mt-ku*rmt-gam/nx*rmt-lam*rmt
    #num of ribosome-mRNA complexes coding for the transporter protein
    results[6]= (gam/nx)*rmt-lam*et#-dp*et
    #num of transporter proteins
    results[7]= kb*r*mm-ku*rmm-gam/nx*rmm-lam*rmm
    #num of ribosome-mRNA complexes coding for the metabolic protein
    results[8]= (we*a/(thetax + a))+ku*rmt+(gam/nx)*rmt-kb*r*mt-dm*mt-lam*mt
    #change in num of mRNAs coding for transporter proteins
    results[9]= (we*a/(thetax + a))+ku*rmm+(gam/nx)*rmm-kb*r*mm-dm*mm-lam*mm
    #change in num of mRNAs coding for metabolic proteins
    results[10]= (gam/nx)*rmq-lam*q#-dp*q
    #num of housekeeping proteins
    results[11]= (et*vt*s_out/(Kt + s_out))-nucat-lam*si#-degredation*si
    results[12]= (wq*a/(thetax + a)/(1 + (q/Kq)^nq))+ku*rmq+(gam/nx)*rmq-kb*r*mq-dm*mq-lam*mq
    #num of free housekeeping mRNAs
    results[13]= (wr*a/(thetar + a))+ku*rmr+(gam/nr)*rmr-kb*r*mr-dm*mr-lam*mr
    #num of free ribosomal mRNAs
    results[14]= ku*rmr +ku*rmt+ ku*rmm+ ku*rmq +(gam/nr)*rmr+
    (gam/nr)*rmr +(gam/nx)*rmt +(gam/nx)*rmm +(gam/nx)*rmq -
    kb*r*mr- kb*r*mt- kb*r*mm -kb*r*mq +ku*nit_mrna_ribo +
    (gam/nx)*nit_mrna_ribo-kb*r*nit_mrna +ku*AA_mrna_ribo +
    (gam/nx)*AA_mrna_ribo-kb*r*AA_mrna -lam*r

    #num of free ribosomes

    #new nitrogen fixing equations. also opposite term in ATP equation
    ##########################################
    results[15]= 0 #NH4_inflow +fixation*nit -export_rate*NH4 -AA_NH4_use -lam*NH4#fixation*nit-export_rate*NH4 -lam*NH4
    #change in num of NH4 molecules in the cell from fixation rate - export rate
    results[16]= (v_nit*a/(thetax + a))+(ku*nit_mrna_ribo)+(gam/nx*nit_mrna_ribo)-(kb*r*nit_mrna)-(dm*nit_mrna)-(lam*nit_mrna)
    #num of free nitrogenase coding mRNA
    results[17]= (kb*r*nit_mrna)-(ku*nit_mrna_ribo)-(gam/nx*nit_mrna_ribo)-(lam*nit_mrna_ribo)
    #num of nit mRNA-ribosome complexes
    results[18]= (gam/nx)*nit_mrna_ribo-lam*nit
    #num of nitrogenase proteins
    results[19]= 0#(export_rate*NH4)#total num of NH4 molecules exported this timestep
    results[20]= 0.0 #lam*N -death_rate*N #the change in num of bacterial cells, instead of just loss to dilution
    results[21]= ns*nucat-ttrate-AA_a_use-lam*a #ATP equation
    ##########################################
    ## AAA equations##
    ##########################################
    results[22]= new_AA -ttrate-lam*AA #num of AA
    results[23]=(gam/nx)*AA_mrna_ribo-lam*AA_prot #num of AA making proteins
    results[24]=((w_AA*a/(thetax + a)))+(ku*AA_mrna_ribo)+
    (gam/nx*AA_mrna_ribo)-(kb*r*AA_mrna)-(dm*AA_mrna)-(lam*AA_mrna)#num of free AA mRNAs
    results[25]=(kb*r*AA_mrna)-(ku*AA_mrna_ribo)-(gam/nx*AA_mrna_ribo)-(lam*AA_mrna_ribo)#num of AA mRNAs bound to ribosomes

end

function nitrogenase_with_AA_with_popn_growth(results,y,p,t)
    s_out = y[1] #external level of substrate, a proxy for plant growth state since more plant = more substrate
    rmr= y[2]#num of ribosome/ribosomal mRNA complex molecules
    em= y[3]#metabolic enzyme molecumes
    rmq= y[4]#num of ribosome/housekeeping mRNA complex molecules
    rmt= y[5]#num of ribosome/transporter protein mRNA complex molecules
    et= y[6]#num of transporter enzyme molecules
    rmm= y[7]#num of ribosome/metabolic enzyme mRNA complexes
    mt= y[8]#num of particles of transporter protein mRNA
    mm= y[9]#num of particles of metabolic enzyme mRNA
    q= y[10]#num of housekeeping protein molecules
    si= y[11]#num of substrate molecules inside the cell (internal)
    mq= y[12]#num of mRNA coding for housekeeping proteins
    mr= y[13]#num of mRNA coding for ribosomes
    r= y[14]#num of free ribosomes
    NH4 = y[15] #num of ammonia molecules fixed
    nit_mrna = y[16]#num of molecules of nitrogen fixing protein
    nit_mrna_ribo = y[17]#num of molecules of mRNA coding for nitrogen fixing protein
    nit = y[18]#num of molecules of mRNA-ribosome complexes coding for nitrogenase
    exported= y[19]#total amount of NH4 exported
    N = y[20] #num of bacteria cells
    a= y[21]#num of ATP (proxy for level of energy)
    AA = y[22] #num of amino acid molecules
    AA_prot=y[23]#num of protein molecules that make new amino acids
    AA_mrna=y[24]#num of mRNA coding for AA making protein
    AA_mrna_ribo=y[25]#num of mRNA coding for AA making protein bound to ribosomes
    #define more constants

    #gam= gmax*a/(Kgamma + a) #gamma is the current rate of translational elongation? gmax is max rate of translational elongation and a is ATP (energy)level
    gam= (gmax*a*AA)/(k_ribo_a*k_ribo_a_AA+k_ribo_a_AA*a+k_ribo_AA_a*AA+a*AA) #two substrate gamma equation that uses both ATP and AA
    ttrate= (rmq + rmr + rmt + rmm + nit_mrna_ribo+AA_mrna_ribo)*gam #total translation rate (sum of the mRNA/ribosome complexes times translation rate)
    #lam= ttrate/(M-AA) #lambda, the growth rate/dilution rate, the ratio of total translation rate to total cell mass
    lam= ttrate/M
    nucat= em*vm*si/(Km + si) #nucat is the num of internal substrate molecules used by catabolism. em is the num of molecules of metabolic enzyme, vm is the rate of the metabolic enzyme and si is the internal substrate level
    #this is all over the half maximal metabolic rate threshold plus the internal substrate level
    fixation = vmax_nit*a/(Kgamma +a)

    AA_vo = ((k_cat_AA*2*a*NH4)/(k_a_NH4*k_a))/(1+(NH4/k_a_NH4)+(NH4/k_NH4)+2+AA*(k_x))
    new_AA = AA_vo*AA_prot

    AA_a_use = new_AA*a_per_AA #each AA produced uses x number of ATP molecules
    AA_NH4_use = new_AA*NH4_per_AA
    ################################################
    #export rate needs to vary with AA levels#
    ################################################

    ##################################################################

    results[1]= ((-et*vt*s_out/(Kt + s_out))*N)-ds*s_out +sub_inflow +(export_rate*NH4*exchange_rate)  #change in amount of substrate, a proxy for amount of host plant
    results[2]= kb*r*mr-ku*rmr-(gam/nr)*rmr-lam*rmr#-dp*r
    #num of ribosome-mRNA complexes coding for ribosomes
    results[3]= (gam/nx)*rmm-lam*em#-dp*em
    #num of metabolic proteins
    results[4]= kb*r*mq-ku*rmq-(gam/nx)*rmq-lam*rmq
    #num of ribosome-mRNA complexes coding for the housekeeping proteins
    results[5]= kb*r*mt-ku*rmt-gam/nx*rmt-lam*rmt
    #num of ribosome-mRNA complexes coding for the transporter protein
    results[6]= (gam/nx)*rmt-lam*et#-dp*et
    #num of transporter proteins
    results[7]= kb*r*mm-ku*rmm-gam/nx*rmm-lam*rmm
    #num of ribosome-mRNA complexes coding for the metabolic protein
    results[8]= (we*a/(thetax + a))+ku*rmt+(gam/nx)*rmt-kb*r*mt-dm*mt-lam*mt
    #change in num of mRNAs coding for transporter proteins
    results[9]= (we*a/(thetax + a))+ku*rmm+(gam/nx)*rmm-kb*r*mm-dm*mm-lam*mm
    #change in num of mRNAs coding for metabolic proteins
    results[10]= (gam/nx)*rmq-lam*q#-dp*q
    #num of housekeeping proteins
    results[11]= (et*vt*s_out/(Kt + s_out))-nucat-lam*si#-degredation*si
    results[12]= (wq*a/(thetax + a)/(1 + (q/Kq)^nq))+ku*rmq+(gam/nx)*rmq-kb*r*mq-dm*mq-lam*mq
    #num of free housekeeping mRNAs
    results[13]= (wr*a/(thetar + a))+ku*rmr+(gam/nr)*rmr-kb*r*mr-dm*mr-lam*mr
    #num of free ribosomal mRNAs
    results[14]= ku*rmr +ku*rmt+ ku*rmm+ ku*rmq +(gam/nr)*rmr+
    (gam/nr)*rmr +(gam/nx)*rmt +(gam/nx)*rmm +(gam/nx)*rmq -
    kb*r*mr- kb*r*mt- kb*r*mm -kb*r*mq +ku*nit_mrna_ribo +
    (gam/nx)*nit_mrna_ribo-kb*r*nit_mrna +ku*AA_mrna_ribo +
    (gam/nx)*AA_mrna_ribo-kb*r*AA_mrna -lam*r

    #num of free ribosomes

    #new nitrogen fixing equations. also opposite term in ATP equation
    ##########################################
    results[15]= NH4_inflow +fixation*nit -export_rate*NH4 -AA_NH4_use -lam*NH4 #change in num of NH4 molecules in the cell from fixation rate - export rate
    results[16]= (v_nit*a/(thetax + a))+(ku*nit_mrna_ribo)+(gam/nx*nit_mrna_ribo)-(kb*r*nit_mrna)-(dm*nit_mrna)-(lam*nit_mrna)
    #num of free nitrogenase coding mRNA
    results[17]= (kb*r*nit_mrna)-(ku*nit_mrna_ribo)-(gam/nx*nit_mrna_ribo)-(lam*nit_mrna_ribo)
    #num of nit mRNA-ribosome complexes
    results[18]= (gam/nx)*nit_mrna_ribo-lam*nit
    #num of nitrogenase proteins
    results[19]= (export_rate*NH4)#total num of NH4 molecules exported this timestep
    results[20]=lam*N -death_rate*N #the change in num of bacterial cells, instead of just loss to dilution
    results[21]= ns*nucat-ttrate-AA_a_use-lam*a #ATP equation
    ##########################################
    ## AAA equations##
    ##########################################
    results[22]= new_AA -ttrate-lam*AA #num of AA
    results[23]=(gam/nx)*AA_mrna_ribo-lam*AA_prot #num of AA making proteins
    results[24]=((w_AA*a/(thetax + a)))+(ku*AA_mrna_ribo)+
    (gam/nx*AA_mrna_ribo)-(kb*r*AA_mrna)-(dm*AA_mrna)-(lam*AA_mrna)#num of free AA mRNAs
    results[25]=(kb*r*AA_mrna)-(ku*AA_mrna_ribo)-(gam/nx*AA_mrna_ribo)-(lam*AA_mrna_ribo)#num of AA mRNAs bound to ribosomes

end

function nitrogenase_with_AA_with_popn_growth_run_out(results,y,p,t)
    s_out = y[1] #external level of substrate, a proxy for plant growth state since more plant = more substrate
    rmr= y[2]#num of ribosome/ribosomal mRNA complex molecules
    em= y[3]#metabolic enzyme molecumes
    rmq= y[4]#num of ribosome/housekeeping mRNA complex molecules
    rmt= y[5]#num of ribosome/transporter protein mRNA complex molecules
    et= y[6]#num of transporter enzyme molecules
    rmm= y[7]#num of ribosome/metabolic enzyme mRNA complexes
    mt= y[8]#num of particles of transporter protein mRNA
    mm= y[9]#num of particles of metabolic enzyme mRNA
    q= y[10]#num of housekeeping protein molecules
    si= y[11]#num of substrate molecules inside the cell (internal)
    mq= y[12]#num of mRNA coding for housekeeping proteins
    mr= y[13]#num of mRNA coding for ribosomes
    r= y[14]#num of free ribosomes
    NH4 = y[15] #num of ammonia molecules fixed
    nit_mrna = y[16]#num of molecules of nitrogen fixing protein
    nit_mrna_ribo = y[17]#num of molecules of mRNA coding for nitrogen fixing protein
    nit = y[18]#num of molecules of mRNA-ribosome complexes coding for nitrogenase
    exported= y[19]#total amount of NH4 exported
    N = y[20] #num of bacteria cells
    a= y[21]#num of ATP (proxy for level of energy)
    AA = y[22] #num of amino acid molecules
    AA_prot=y[23]#num of protein molecules that make new amino acids
    AA_mrna=y[24]#num of mRNA coding for AA making protein
    AA_mrna_ribo=y[25]#num of mRNA coding for AA making protein bound to ribosomes
    #define more constants

    #gam= gmax*a/(Kgamma + a) #gamma is the current rate of translational elongation? gmax is max rate of translational elongation and a is ATP (energy)level
    gam= (gmax*a*AA)/(k_ribo_a*k_ribo_a_AA+k_ribo_a_AA*a+k_ribo_AA_a*AA+a*AA) #two substrate gamma equation that uses both ATP and AA
    ttrate= (rmq + rmr + rmt + rmm + nit_mrna_ribo+AA_mrna_ribo)*gam #total translation rate (sum of the mRNA/ribosome complexes times translation rate)
    #lam= ttrate/(M-AA) #lambda, the growth rate/dilution rate, the ratio of total translation rate to total cell mass
    lam= ttrate/M
    nucat= em*vm*si/(Km + si) #nucat is the num of internal substrate molecules used by catabolism. em is the num of molecules of metabolic enzyme, vm is the rate of the metabolic enzyme and si is the internal substrate level
    #this is all over the half maximal metabolic rate threshold plus the internal substrate level
    fixation = vmax_nit*a/(Kgamma +a)

    AA_vo = ((k_cat_AA*2*a*NH4)/(k_a_NH4*k_a))/(1+(NH4/k_a_NH4)+(NH4/k_NH4)+2+AA*(k_x))
    new_AA = AA_vo*AA_prot

    AA_a_use = new_AA*a_per_AA #each AA produced uses x number of ATP molecules
    AA_NH4_use = new_AA*NH4_per_AA
    ################################################
    #export rate needs to vary with AA levels#
    ################################################

    ##################################################################

    results[1]= ((-et*vt*s_out/(Kt + s_out))*N)-ds*s_out +(export_rate*NH4*exchange_rate) #+sub_inflow  #change in amount of substrate, a proxy for amount of host plant
    results[2]= kb*r*mr-ku*rmr-(gam/nr)*rmr-lam*rmr#-dp*r
    #num of ribosome-mRNA complexes coding for ribosomes
    results[3]= (gam/nx)*rmm-lam*em#-dp*em
    #num of metabolic proteins
    results[4]= kb*r*mq-ku*rmq-(gam/nx)*rmq-lam*rmq
    #num of ribosome-mRNA complexes coding for the housekeeping proteins
    results[5]= kb*r*mt-ku*rmt-gam/nx*rmt-lam*rmt
    #num of ribosome-mRNA complexes coding for the transporter protein
    results[6]= (gam/nx)*rmt-lam*et#-dp*et
    #num of transporter proteins
    results[7]= kb*r*mm-ku*rmm-gam/nx*rmm-lam*rmm
    #num of ribosome-mRNA complexes coding for the metabolic protein
    results[8]= (we*a/(thetax + a))+ku*rmt+(gam/nx)*rmt-kb*r*mt-dm*mt-lam*mt
    #change in num of mRNAs coding for transporter proteins
    results[9]= (we*a/(thetax + a))+ku*rmm+(gam/nx)*rmm-kb*r*mm-dm*mm-lam*mm
    #change in num of mRNAs coding for metabolic proteins
    results[10]= (gam/nx)*rmq-lam*q#-dp*q
    #num of housekeeping proteins
    results[11]= (et*vt*s_out/(Kt + s_out))-nucat-lam*si#-degredation*si
    results[12]= (wq*a/(thetax + a)/(1 + (q/Kq)^nq))+ku*rmq+(gam/nx)*rmq-kb*r*mq-dm*mq-lam*mq
    #num of free housekeeping mRNAs
    results[13]= (wr*a/(thetar + a))+ku*rmr+(gam/nr)*rmr-kb*r*mr-dm*mr-lam*mr
    #num of free ribosomal mRNAs
    results[14]= ku*rmr +ku*rmt+ ku*rmm+ ku*rmq +(gam/nr)*rmr+
    (gam/nr)*rmr +(gam/nx)*rmt +(gam/nx)*rmm +(gam/nx)*rmq -
    kb*r*mr- kb*r*mt- kb*r*mm -kb*r*mq +ku*nit_mrna_ribo +
    (gam/nx)*nit_mrna_ribo-kb*r*nit_mrna +ku*AA_mrna_ribo +
    (gam/nx)*AA_mrna_ribo-kb*r*AA_mrna -lam*r

    #num of free ribosomes

    #new nitrogen fixing equations. also opposite term in ATP equation
    ##########################################
    results[15]= fixation*nit -export_rate*NH4 -AA_NH4_use -lam*NH4 # +NH4_inflow #change in num of NH4 molecules in the cell 
    results[16]= (v_nit*a/(thetax + a))+(ku*nit_mrna_ribo)+(gam/nx*nit_mrna_ribo)-(kb*r*nit_mrna)-(dm*nit_mrna)-(lam*nit_mrna)
    #num of free nitrogenase coding mRNA
    results[17]= (kb*r*nit_mrna)-(ku*nit_mrna_ribo)-(gam/nx*nit_mrna_ribo)-(lam*nit_mrna_ribo)
    #num of nit mRNA-ribosome complexes
    results[18]= (gam/nx)*nit_mrna_ribo-lam*nit
    #num of nitrogenase proteins
    results[19]= (export_rate*NH4)#total num of NH4 molecules exported this timestep
    results[20]=lam*N -death_rate*N #the change in num of bacterial cells, instead of just loss to dilution
    results[21]= ns*nucat-ttrate-AA_a_use-lam*a #ATP equation
    ##########################################
    ## AAA equations##
    ##########################################
    results[22]= new_AA -ttrate-lam*AA #num of AA
    results[23]=(gam/nx)*AA_mrna_ribo-lam*AA_prot #num of AA making proteins
    results[24]=((w_AA*a/(thetax + a)))+(ku*AA_mrna_ribo)+
    (gam/nx*AA_mrna_ribo)-(kb*r*AA_mrna)-(dm*AA_mrna)-(lam*AA_mrna)#num of free AA mRNAs
    results[25]=(kb*r*AA_mrna)-(ku*AA_mrna_ribo)-(gam/nx*AA_mrna_ribo)-(lam*AA_mrna_ribo)#num of AA mRNAs bound to ribosomes

end

"""
given the values of metabolic enzyme and internal subtrate, calculates the number of new ATP molecules produced by metabolism of said substrate
"""
find_nucat(em,si)=begin
    output= (em*vm*si/(Km + si))*ns #nucat is the num of internal substrate molecules used by catabolism. em is the num of molecules of metabolic enzyme, vm is the rate of the metabolic enzyme and si is the internal substrate level
end

"""
Calculates the total transcription rate for given numbers of mRNAs bound to ribosomes and current ATP level.
"""
find_ttrate_no_AA(rmq,rmr,rmt,rmm,nit_mrna_ribo,AA_mrna_ribo,a)=begin
    output = ((rmq + rmr + rmt + rmm + nit_mrna_ribo+AA_mrna_ribo)*(gmax*a/(Kgamma + a)))
end

    
"""
Calculates the total transcription rate for given numbers of mRNAs bound to ribosomes and current ATP level, and current AA level.
"""
find_ttrate_with_AA(rmq,rmr,rmt,rmm,nit_mrna_ribo,AA_mrna_ribo,a,AA)=begin
    ((rmq + rmr + rmt + rmm + nit_mrna_ribo+AA_mrna_ribo)*((gmax*a*AA)/(k_ribo_a*k_ribo_a_AA+k_ribo_a_AA*a+k_ribo_AA_a*AA+a*AA)))
end

find_gamma_no_AA(a)=begin
    (gmax*a/(Kgamma + a))
end

find_gamma_with_AA(a,AA)=begin
    (gmax*a*AA)/(k_ribo_a*k_ribo_a_AA+k_ribo_a_AA*a+k_ribo_AA_a*AA+a*AA)
end
            
find_new_AA(a,NH4,AA,AA_prot)=begin
    output = (((k_cat_AA*2*a*NH4)/(k_a_NH4*k_a))/(1+(NH4/k_a_NH4)+(NH4/k_NH4)+2+AA*(k_x)))*AA_prot
end
################
#try to make equation such that if there is very little ATP it shouldn't make AA even if there is an excess of NH4

find_transport(et,s_out,N)= begin
    ((et*vt*s_out/(Kt + s_out))*N)
end

find_sub_dilution(s_out)= begin
    ds*s_out
end 
#     ((-et*vt*s_out/(Kt + s_out))*N)-ds*s_out +sub_inflow

