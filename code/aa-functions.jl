##Functions for AA model##

function model_AA(results,y,p,t)
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
#     Kgamma= gmax/180 #Kgamma is the half maximal translational elongation threshold, gmax is the max rate of translational elongation and Kp is net rate of translation of housekeeping proteins
    Kgamma = 3.0e8
    #should equal 7 molecules per cell
    gam= gmax*a/(Kgamma + a) #gamma is the current rate of translational elongation? gmax is max rate of translational elongation and a is ATP (energy)level


    #######
    #need to make gam limited like we did with the carrying capacity of a population. gamma should plateau as it approches the number of amino acids
    #######


    ttrate= (rmq + rmr + rmt + rmm + nit_mrna_ribo+AA_mrna_ribo)*gam #total translation rate (sum of the mRNA/ribosome complexes times translation rate)
    lam= ttrate/Ms #lambda, the growth rate/dilution rate, the ratio of total translation rate to total cell mass
    nucat= em*vm*si/(Km + abs(si)) #nucat is the num of internal substrate molecules used by catabolism. em is the num of molecules of metabolic enzyme, vm is the rate of the metabolic enzyme and si is the internal substrate level
    #this is all over the half maximal metabolic rate threshold plus the internal substrate level
    fixation = vmax_nit*a/(Kgamma +a)

    #AA equations
    AA_req = ttrate #the number of aa required at this time time for all translation to go forward
    AA_diff = AA-AA_req #aa difference is the required num of AA minus the num of AA in the cell

#     println("amino acid requirment is $AA_req")
#     println("amino acids available is $AA")
    if AA_diff > 0 #if there are enough or more than enough AAs to fufill all translation then dont make new AAs
        new_AA = 0
        export_rate = 0.9
    else #AA_diff < 0 #if there are not enough AAs for all translation then make new AAs
        AA_a_use = vmax_AA*a/(Kgamma+a)
        new_AA = AA_a_use * AA_efficiency
        export_rate = 0.1
        gam=0 #gamma is the current rate of translational elongation? gmax is max rate of translational elongation and a is ATP (energy)level
        ttrate= (rmq + rmr + rmt + rmm + nit_mrna_ribo+AA_mrna_ribo)*gam #total translation rate (sum of the mRNA/ribosome complexes times translation rate)
        AA_req = ttrate #the number of aa required at this time time for all translation to go forward
    end
#     println("num of AA produced is $new_AA")

    ##################################################################
    results[1]= 0 #((-et*vt*s_out/(Kt + abs(s_out)))*N)-ds*s_out #+(export_rate*NH4*50) #change in amount of substrate, a proxy for amount of host plant
#     var1 = 10+(-et*vt*s_out/(Kt + abs(s_out)))*N-ds*s_out
#     println("net change in external substrate is $var1")

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
    ######################################
    degredation= 0.0
    ###### setting degredation of internal substrate here to make the cells go extinct eventually
    results[11]= (et*vt*s_out/(Kt + abs(s_out)))-nucat-lam*si#-degredation*si
    #num of internal substrate molecules
#     var2 = (et*vt*s_out/(Kt + abs(s_out)))
#     println("import rate is $var2")

    results[12]= (wq*a/(thetax + a)/(1 + (q/Kq)^nq))+ku*rmq+(gam/nx)*rmq-kb*r*mq-dm*mq-lam*mq
    #num of free housekeeping mRNAs
    results[13]= (wr*a/(thetar + a))+ku*rmr+(gam/nr)*rmr-kb*r*mr-dm*mr-lam*mr
    #num of free ribosomal mRNAs
    results[14]= ku*rmr +ku*rmt+ ku*rmm+ ku*rmq +(gam/nr)*rmr +(gam/nr)*rmr +(gam/nx)*rmt +(gam/nx)*rmm +(gam/nx)*rmq -kb*r*mr- kb*r*mt- kb*r*mm -kb*r*mq +ku*nit_mrna_ribo +(gam/nx)*nit_mrna_ribo-kb*r*nit_mrna +ku*AA_mrna_ribo +(gam/nx)*AA_mrna_ribo-kb*r*AA_mrna  -lam*r
    #num of free ribosomes

    #new nitrogen fixing equations. also opposite term in ATP equation
    ##########################################
    results[15]= fixation*nit-export_rate*NH4 -lam*NH4
    #change in num of NH4 molecules in the cell from fixation rate - export rate
    results[16]= (v_nit*a/(thetax + a))+(ku*nit_mrna_ribo)+(gam/nx*nit_mrna_ribo)-(kb*r*nit_mrna)-(dm*nit_mrna)-(lam*nit_mrna)
    #num of free nitrogenase coding mRNA
    results[17]= (kb*r*nit_mrna)-(ku*nit_mrna_ribo)-(gam/nx*nit_mrna_ribo)-(lam*nit_mrna_ribo)
    #num of nit mRNA-ribosome complexes
    results[18]= gam/nx*nit_mrna_ribo-lam*nit
    #num of nitrogenase proteins
    results[19]= (export_rate*NH4)#total num of NH4 molecules exported this timestep
    results[20]= 0#lam*N-death_rate*N #the change in num of bacterial cells, instead of just loss to dilution
    results[21]= ns*nucat-ttrate-fixation*nit-AA_efficiency-lam*a
    ##########################################
    ## AAA equations##
    ##########################################


    results[22]= 0#new_AA-AA_req-lam*AA #num of AA
    results[23]=gam/nx*AA_mrna_ribo-lam*AA_prot #num of AA making proteins
    results[24]=(w_AA*a/(thetax + a))+(ku*AA_mrna_ribo)+(gam/nx*AA_mrna_ribo)-(kb*r*AA_mrna)-(dm*AA_mrna)-(lam*AA_mrna)#num of free AA mRNAs
    results[25]=(kb*r*AA_mrna)-(ku*AA_mrna_ribo)-(gam/nx*AA_mrna_ribo)-(lam*AA_mrna_ribo)#num of AA mRNAs bound to ribosomes



end



function model_AA2(results,y,p,t)
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
#     Kgamma= gmax/180 #Kgamma is the half maximal translational elongation threshold, gmax is the max rate of translational elongation and Kp is net rate of translation of housekeeping proteins
    Kgamma = 3.0e8
    #should equal 7 molecules per cell
    gam= gmax*a/(Kgamma + a) #gamma is the current rate of translational elongation? gmax is max rate of translational elongation and a is ATP (energy)level

    #######
    #need to make gam limited like we did with the carrying capacity of a population. gamma should plateau as it approches the number of amino acids
    #######

    ttrate= (rmq + rmr + rmt + rmm + nit_mrna_ribo+AA_mrna_ribo)*gam #total translation rate (sum of the mRNA/ribosome complexes times translation rate)
    lam= ttrate/Ms #lambda, the growth rate/dilution rate, the ratio of total translation rate to total cell mass
    nucat= em*vm*si/(Km + abs(si)) #nucat is the num of internal substrate molecules used by catabolism. em is the num of molecules of metabolic enzyme, vm is the rate of the metabolic enzyme and si is the internal substrate level
    #this is all over the half maximal metabolic rate threshold plus the internal substrate level
    fixation = vmax_nit*a/(Kgamma +a)

    #AA equations
    AA_req = ttrate #the number of aa required at this time time for all translation to go forward
    AA_diff = AA-AA_req #aa difference is the required num of AA minus the num of AA in the cell

    if AA_diff > 0 #if there are enough or more than enough AAs to fufill all translation then dont make new AAs
        new_AA = 0
        export_rate = 0.9
        AA_a_use = 0
    else #AA_diff < 0 #if there are not enough AAs for all translation then make new AAs
        AA_a_use = vmax_AA*a/(Kgamma+ a)
        new_AA = AA_a_use * AA_efficiency*AA_prot
        export_rate = 0.1
        gam=0
        ttrate= (rmq + rmr + rmt + rmm + nit_mrna_ribo+AA_mrna_ribo)*gam #total translation rate (sum of the mRNA/ribosome complexes times translation rate)
        AA_req = ttrate #the number of aa required at this time time for all translation to go forward
    end

    ##################################################################
    results[1]= ((-et*vt*s_out/(Kt + abs(s_out)))*N)-ds*s_out +(export_rate*NH4*50) #change in amount of substrate, a proxy for amount of host plant
#     var1 = 10+(-et*vt*s_out/(Kt + abs(s_out)))*N-ds*s_out
#     println("net change in external substrate is $var1")

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
    ######################################
    degredation= 0.0
    ###### setting degredation of internal substrate here to make the cells go extinct eventually
    results[11]= (et*vt*s_out/(Kt + abs(s_out)))-nucat-lam*si#-degredation*si
    results[12]= (wq*a/(thetax + a)/(1 + (q/Kq)^nq))+ku*rmq+(gam/nx)*rmq-kb*r*mq-dm*mq-lam*mq
    #num of free housekeeping mRNAs
    results[13]= (wr*a/(thetar + a))+ku*rmr+(gam/nr)*rmr-kb*r*mr-dm*mr-lam*mr
    #num of free ribosomal mRNAs
    results[14]= ku*rmr +ku*rmt+ ku*rmm+ ku*rmq +(gam/nr)*rmr +(gam/nr)*rmr +(gam/nx)*rmt +(gam/nx)*rmm +(gam/nx)*rmq -kb*r*mr- kb*r*mt- kb*r*mm -kb*r*mq +ku*nit_mrna_ribo +(gam/nx)*nit_mrna_ribo-kb*r*nit_mrna +ku*AA_mrna_ribo +(gam/nx)*AA_mrna_ribo-kb*r*AA_mrna  -lam*r
    #num of free ribosomes

    #new nitrogen fixing equations. also opposite term in ATP equation
    ##########################################
    results[15]= fixation*nit-export_rate*NH4 -lam*NH4
    #change in num of NH4 molecules in the cell from fixation rate - export rate
    results[16]= (v_nit*a/(thetax + a))+(ku*nit_mrna_ribo)+(gam/nx*nit_mrna_ribo)-(kb*r*nit_mrna)-(dm*nit_mrna)-(lam*nit_mrna)
    #num of free nitrogenase coding mRNA
    results[17]= (kb*r*nit_mrna)-(ku*nit_mrna_ribo)-(gam/nx*nit_mrna_ribo)-(lam*nit_mrna_ribo)
    #num of nit mRNA-ribosome complexes
    results[18]= gam/nx*nit_mrna_ribo-lam*nit
    #num of nitrogenase proteins
    results[19]= (export_rate*NH4)#total num of NH4 molecules exported this timestep
    results[20]=lam*N-death_rate*N #the change in num of bacterial cells, instead of just loss to dilution
    results[21]= ns*nucat-ttrate-fixation*nit-AA_a_use-lam*a
    ##########################################
    ## AAA equations##
    ##########################################


    results[22]= 0#new_AA-AA_req-lam*AA #num of AA
    results[23]=gam/nx*AA_mrna_ribo-lam*AA_prot #num of AA making proteins
    results[24]=(w_AA*a/(thetax + a))+(ku*AA_mrna_ribo)+(gam/nx*AA_mrna_ribo)-(kb*r*AA_mrna)-(dm*AA_mrna)-(lam*AA_mrna)#num of free AA mRNAs
    results[25]=(kb*r*AA_mrna)-(ku*AA_mrna_ribo)-(gam/nx*AA_mrna_ribo)-(lam*AA_mrna_ribo)#num of AA mRNAs bound to ribosomes



end


function model_AA3(results,y,p,t)
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
#     Kgamma= gmax/180 #Kgamma is the half maximal translational elongation threshold, gmax is the max rate of translational elongation and Kp is net rate of translation of housekeeping proteins
    Kgamma = 3.0e8
    #should equal 7 molecules per cell
    #gam= gmax*a/(Kgamma + a) #gamma is the current rate of translational elongation? gmax is max rate of translational elongation and a is ATP (energy)level

    #k_ribo_a = 1
    #k_ribo_AA = 1
    #k_ribo_a_AA = 1
    #k_ribo_AA_a = 1

    gam= (gmax*a*AA)/(k_ribo_a*k_ribo_a_AA+k_ribo_a_AA*a+k_ribo_AA_a*AA+a*AA) #updated gamma equation that uses both ATP and AA


    ttrate= (rmq + rmr + rmt + rmm + nit_mrna_ribo+AA_mrna_ribo)*gam #total translation rate (sum of the mRNA/ribosome complexes times translation rate)
    lam= ttrate/(Ms-AA) #lambda, the growth rate/dilution rate, the ratio of total translation rate to total cell mass
    if lam<0.00000001 lam=0.00000001 end
    nucat= em*vm*si/(Km + abs(si)) #nucat is the num of internal substrate molecules used by catabolism. em is the num of molecules of metabolic enzyme, vm is the rate of the metabolic enzyme and si is the internal substrate level
    #this is all over the half maximal metabolic rate threshold plus the internal substrate level
    fixation = vmax_nit*a/(Kgamma +a)

    #AA equations
    AA_req = ttrate #the number of aa required at this time time for all translation to go forward

    # k_AA= 50.0 #the constant that will determine how many AAs the cell is aiming to have
    # k_cat_AA= 1000.0 #rate of production of product of the AA enzyme when bound to both substrates
    # println("current k_cat_AA value inside ODE is ",k_cat_AA)
    # k_a_NH4 = 10.0 #rate of binding of NH4 to atp-enzyme complex
    # k_a = 10.0 #rate of binding of atp to enzyme
    # k_NH4 = 10.0 #rate of binding of NH4 to enzyme
    # k_a_AA = 5.0 #rate of binding of the end product inhibitor to the enzyme-atp complex
    # k_NH4_AA= 5.0 #rate of binding of the end product acting as an inhibitor to the enzyme-NH4 complex

    AA_vo = ((k_cat_AA*a*NH4)/k_a_NH4*k_a)/(1+((1+(AA/k_a_AA)+(NH4/k_a_NH4)))*(a/k_a)+(1+(AA/k_NH4_AA)*(NH4/k_NH4)))

    new_AA = AA_vo*AA_prot

#     Kgamma = 3.0e8
#     k_AA= 100.0 #this constant...
#     a= 1000.0 #num of ATP
#     cons1=1.0 #... and this constant
    #the product of these two constants determines the num of AA the cell is aiming to have.
#     var1 = vmax_AA*a/(Kgamma*(cons1+(cur_AA/k_AA)) + a)

    # new_AA = max((AA_a_use * AA_efficiency*AA_prot),AA_req)
#     export_rate = max(0,1/AA_diff)
#     println("export rate is $export_rate")
    export_rate=0.3

#     println("")

    ################################################
    #export rate needs to vary with AA levels#
    ################################################

#     gam=0.00000001
#     ttrate= (rmq + rmr + rmt + rmm + nit_mrna_ribo+AA_mrna_ribo)*gam #total translation rate (sum of the mRNA/ribosome complexes times translation rate)
#     AA_req = ttrate #the number of aa required at this time time for all translation to go forward
#     lam= ttrate/(Ms-AA) #lambda, the growth rate/dilution rate, the ratio of total translation rate to total cell mass


    ##################################################################
    results[1]= ((-et*vt*s_out/(Kt + abs(s_out)))*N)-ds*s_out +(export_rate*NH4*50) #change in amount of substrate, a proxy for amount of host plant
    if s_out <0.00000001 s_out=0.00000001 end
    results[2]= kb*r*mr-ku*rmr-(gam/nr)*rmr-lam*rmr#-dp*r
    if rmr<0.00000001 rmr = 0.00000001 end
    #num of ribosome-mRNA complexes coding for ribosomes
    results[3]= (gam/nx)*rmm-lam*em#-dp*em
    if em<0.00000001 em=0.00000001 end
    #num of metabolic proteins
    results[4]= kb*r*mq-ku*rmq-(gam/nx)*rmq-lam*rmq
    if rmq<0.00000001 rmq=0.00000001 end
    #num of ribosome-mRNA complexes coding for the housekeeping proteins
    results[5]= kb*r*mt-ku*rmt-gam/nx*rmt-lam*rmt
    if rmt<0.00000001 rmt=0.00000001 end
    #num of ribosome-mRNA complexes coding for the transporter protein
    results[6]= (gam/nx)*rmt-lam*et#-dp*et
    if et<0.00000001 et=0.00000001 end
    #num of transporter proteins
    results[7]= kb*r*mm-ku*rmm-gam/nx*rmm-lam*rmm
    if rmm<0.00000001 rmm=0.00000001 end
    #num of ribosome-mRNA complexes coding for the metabolic protein
    results[8]= (we*a/(thetax + a))+ku*rmt+(gam/nx)*rmt-kb*r*mt-dm*mt-lam*mt
    if mt<0.00000001 mt=0.00000001 end
    #change in num of mRNAs coding for transporter proteins
    results[9]= (we*a/(thetax + a))+ku*rmm+(gam/nx)*rmm-kb*r*mm-dm*mm-lam*mm
    if mm<0.00000001 mm=0.00000001 end
    #change in num of mRNAs coding for metabolic proteins
    results[10]= (gam/nx)*rmq-lam*q#-dp*q
    if q<0.00000001 q=0.00000001 end
    #num of housekeeping proteins
    ######################################
#     degredation= 0.0
    ###### setting degredation of internal substrate here to make the cells go extinct eventually
    results[11]= (et*vt*s_out/(Kt + abs(s_out)))-nucat-lam*si#-degredation*si
    if si<0.00000001 si=0.00000001 end
    results[12]= (wq*a/(thetax + a)/(1 + (q/Kq)^nq))+ku*rmq+(gam/nx)*rmq-kb*r*mq-dm*mq-lam*mq
    if mq<0.00000001 mq=0.00000001 end
    #num of free housekeeping mRNAs
    results[13]= (wr*a/(thetar + a))+ku*rmr+(gam/nr)*rmr-kb*r*mr-dm*mr-lam*mr
    if mr<0.00000001 mr=0.00000001 end
    #num of free ribosomal mRNAs
    results[14]= ku*rmr +ku*rmt+ ku*rmm+ ku*rmq +(gam/nr)*rmr+
    (gam/nr)*rmr +(gam/nx)*rmt +(gam/nx)*rmm +(gam/nx)*rmq -
    kb*r*mr- kb*r*mt- kb*r*mm -kb*r*mq +ku*nit_mrna_ribo +
    (gam/nx)*nit_mrna_ribo-kb*r*nit_mrna +ku*AA_mrna_ribo +
    (gam/nx)*AA_mrna_ribo-kb*r*AA_mrna -lam*r
    if r < 0.00000001 r = 0.00000001 end
    #num of free ribosomes

    #new nitrogen fixing equations. also opposite term in ATP equation
    ##########################################
    results[15]= fixation*nit-export_rate*NH4 -lam*NH4
    if NH4<0.00000001 NH4=0.00000001 end
    #change in num of NH4 molecules in the cell from fixation rate - export rate
    results[16]= (v_nit*a/(thetax + a))+(ku*nit_mrna_ribo)+(gam/nx*nit_mrna_ribo)-
    (kb*r*nit_mrna)-(dm*nit_mrna)-(lam*nit_mrna)
    if nit_mrna<0.00000001 nit_mrna=0.00000001 end
    #num of free nitrogenase coding mRNA
    results[17]= (kb*r*nit_mrna)-(ku*nit_mrna_ribo)-(gam/nx*nit_mrna_ribo)-(lam*nit_mrna_ribo)
    if nit_mrna_ribo<0.00000001 nit_mrna_ribo=0.00000001 end
    #num of nit mRNA-ribosome complexes
    results[18]= (gam/nx)*nit_mrna_ribo-lam*nit
    if nit<0.00000001 nit=0.00000001 end
    #num of nitrogenase proteins
    results[19]= (export_rate*NH4)#total num of NH4 molecules exported this timestep
    results[20]=lam*N-death_rate*N #the change in num of bacterial cells, instead of just loss to dilution
    if N<0.00000001 N=0.00000001 end
    #########################################
    #set how many ATP it takes to make one AA
    AA_a_use= new_AA*2
    #########################################
    results[21]= ns*nucat-ttrate-fixation*nit-AA_a_use-lam*a
    if a<0.00000001 a=0.00000001 end
    ##########################################
    ## AAA equations##
    ##########################################

    results[22]= new_AA-AA_req-lam*AA #num of AA
	####save the final values of these three variable to file###
    if AA<0.00000001 AA=0.00000001 end

    results[23]=(gam/nx)*AA_mrna_ribo-lam*AA_prot #num of AA making proteins
    if AA_prot<0.00000001 AA_prot=0.00000001 end

    AA_dependent_rate= (10.0+AA_req)/AA #AA_req/AA #this is the first attempt at AA control
    if AA_dependent_rate < 0.00000001 AA_dependent_rate =0.0000001 end

    results[24]=((w_AA*a/(thetax + a))*AA_dependent_rate)+(ku*AA_mrna_ribo)+
    (gam/nx*AA_mrna_ribo)-(kb*r*AA_mrna)-(dm*AA_mrna)-(lam*AA_mrna)#num of free AA mRNAs
    if AA_mrna<0.00000001 AA_mrna=0.0000001 end

    results[25]=(kb*r*AA_mrna)-(ku*AA_mrna_ribo)-(gam/nx*AA_mrna_ribo)-(lam*AA_mrna_ribo)#num of AA mRNAs bound to ribosomes
    if AA_mrna_ribo<0.00000001 AA_mrna_ribo=0.0000001 end

end
