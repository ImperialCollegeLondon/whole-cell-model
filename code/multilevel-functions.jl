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

    if s_out <0.000001 s_out=0 end
    if rmr<0.000001 rmr= 0 end
    if em<0.000001 em=0 end
    if rmq<0.000001 rmq=0 end
    if rmt<0.000001 rmt=0 end
    if et<0.000001 et=0 end
    if rmm<0.000001 rmm=0 end
    if mt<0.000001 mt=0 end
    if mm<0.000001 mm=0 end
    if q<0.000001 q=0 end
    if si<0.000001 si=0 end
    if mq<0.000001 mp=0 end
    if mr<0.000001 mr=0 end
    if r < 0.000001 r = 0 end    
    if NH4<0.000001 NH4=0 end
    if nit_mrna<0.000001 nit_mrna=0 end
    if nit_mrna_ribo<0.000001 nit_mrna_ribo=0 end
    if nit<0.000001 nit=0 end
    if N<0.000001 N=0 end
    if a<0.000001 a=0 end
    if AA<0.000001 AA=0 end
    if AA_prot<0.000001 AA_prot=0 end
    if AA_mrna<0.0000001 AA_mrna=0 end    
    if AA_mrna_ribo<0.000001 AA_mrna_ribo=0 end

    #define more constants
#~      Kgamma= gmax/180 #Kgamma is the half maximal translational elongation threshold, gmax is the max rate of translational elongation and Kp is net rate of translation of housekeeping proteins
    Kgamma = 3.0e8
    #should equal 7 molecules per cell
    #gam= gmax*a/(Kgamma + a) #gamma is the current rate of translational elongation? gmax is max rate of translational elongation and a is ATP (energy)level

    #k_ribo_a = 1
    #k_ribo_AA = 1
    #k_ribo_a_AA = 1
    #k_ribo_AA_a = 1

    gam= (gmax*a*AA)/(k_ribo_a*k_ribo_a_AA+k_ribo_a_AA*a+k_ribo_AA_a*AA+a*AA) #updated gamma equation that uses both ATP and AA

    ttrate= (rmq + rmr + rmt + rmm + nit_mrna_ribo+AA_mrna_ribo)*gam #total translation rate (sum of the mRNA/ribosome complexes times translation rate)
    #lam= ttrate/(M-AA) #lambda, the growth rate/dilution rate, the ratio of total translation rate to total cell mass
    lam= ttrate/M
    #if lam<0.00000001 lam=0.00000001 end
    nucat= em*vm*si/(Km + abs(si)) #nucat is the num of internal substrate molecules used by catabolism. em is the num of molecules of metabolic enzyme, vm is the rate of the metabolic enzyme and si is the internal substrate level
    #this is all over the half maximal metabolic rate threshold plus the internal substrate level
    fixation = 0#vmax_nit*a/(Kgamma +a)

    #AA equations
     #k_AA= 1000.0 #the constant that will determine how many AAs the cell is aiming to have
     #k_cat_AA= 1000.0 #rate of production of product of the AA enzyme when bound to both substrates
    # println("current k_cat_AA value inside ODE is ",k_cat_AA)
     #k_a_NH4 = 1000.0 #rate of binding of NH4 to atp-enzyme complex
     #k_a = 1000.0 #rate of binding of atp to enzyme
     #k_NH4 = 1000.0 #rate of binding of NH4 to enzyme

     #k_a_AA = 0.01 #rate of binding of the end product inhibitor to the enzyme-atp complex
     #k_NH4_AA= 0.01 #rate of binding of the end product acting as an inhibitor to the enzyme-NH4 complex

    AA_vo = ((k_cat_AA*a*NH4)/k_a_NH4*k_a)/(1+((1+(AA/k_a_AA)+(NH4/k_a_NH4)))*(a/k_a)+(1+(AA/k_NH4_AA)*(NH4/k_NH4)))

    new_AA = AA_vo*AA_prot

    #open("../data/param-sweep-$k_cat_AA-$k_a_NH4-$k_a-$k_NH4-$k_a_AA-$k_NH4_AA-$k_ribo_a-$k_ribo_AA-$k_ribo_a_AA-$k_ribo_AA_a.csv","a") do f
    #write(f, "\n$AA, $a, $N")
    #end

    AA_a_use = new_AA*2 #each AA produced uses x number of ATP molecules


    export_rate=0#0.3
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
    results[15]= 0#fixation*nit-export_rate*NH4 -lam*NH4
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
	#println("new AA is $new_AA")
	#AA_used = ttrate+lam*AA
	#println("AA used is $AA_used")
	#tempvar= results[22]
	#println("results 22 is $tempvar")
	####save the final values of these three variable to file###
    
    #open("../data/high-AA-params-file-$k_cat_AA-$k_a_NH4-$k_a-$k_NH4-$k_a_AA-$k_NH4_AA-$k_ribo_a-$k_ribo_AA-$k_ribo_a_AA-$k_ribo_AA_a.csv","a") do f
    #write(f, "\n$new_AA, $ttrate, $AA_prot")
    #end


    results[23]=(gam/nx)*AA_mrna_ribo-lam*AA_prot #num of AA making proteins


    results[24]=((w_AA*a/(thetax + a)))+(ku*AA_mrna_ribo)+
    (gam/nx*AA_mrna_ribo)-(kb*r*AA_mrna)-(dm*AA_mrna)-(lam*AA_mrna)#num of free AA mRNAs



    results[25]=(kb*r*AA_mrna)-(ku*AA_mrna_ribo)-(gam/nx*AA_mrna_ribo)-(lam*AA_mrna_ribo)#num of AA mRNAs bound to ribosomes


end


function AA_popn_growth(results,y,p,t)
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

    if s_out <0.000001 s_out=0 end
    if rmr<0.000001 rmr= 0 end
    if em<0.000001 em=0 end
    if rmq<0.000001 rmq=0 end
    if rmt<0.000001 rmt=0 end
    if et<0.000001 et=0 end
    if rmm<0.000001 rmm=0 end
    if mt<0.000001 mt=0 end
    if mm<0.000001 mm=0 end
    if q<0.000001 q=0 end
    if si<0.000001 si=0 end
    if mq<0.000001 mp=0 end
    if mr<0.000001 mr=0 end
    if r < 0.000001 r = 0 end    
    if NH4<0.000001 NH4=0 end
    if nit_mrna<0.000001 nit_mrna=0 end
    if nit_mrna_ribo<0.000001 nit_mrna_ribo=0 end
    if nit<0.000001 nit=0 end
    if N<0.000001 N=0 end
    if a<0.000001 a=0 end
    if AA<0.000001 AA=0 end
    if AA_prot<0.000001 AA_prot=0 end
    if AA_mrna<0.0000001 AA_mrna=0 end    
    if AA_mrna_ribo<0.000001 AA_mrna_ribo=0 end

    #define more constants
#~      Kgamma= gmax/180 #Kgamma is the half maximal translational elongation threshold, gmax is the max rate of translational elongation and Kp is net rate of translation of housekeeping proteins
    Kgamma = 3.0e8
    #should equal 7 molecules per cell
    #gam= gmax*a/(Kgamma + a) #gamma is the current rate of translational elongation? gmax is max rate of translational elongation and a is ATP (energy)level

    #k_ribo_a = 1
    #k_ribo_AA = 1
    #k_ribo_a_AA = 1
    #k_ribo_AA_a = 1

    gam= (gmax*a*AA)/(k_ribo_a*k_ribo_a_AA+k_ribo_a_AA*a+k_ribo_AA_a*AA+a*AA) #updated gamma equation that uses both ATP and AA

    ttrate= (rmq + rmr + rmt + rmm + nit_mrna_ribo+AA_mrna_ribo)*gam #total translation rate (sum of the mRNA/ribosome complexes times translation rate)
    #lam= ttrate/(M-AA) #lambda, the growth rate/dilution rate, the ratio of total translation rate to total cell mass
    lam= ttrate/M
    #if lam<0.00000001 lam=0.00000001 end
    nucat= em*vm*si/(Km + abs(si)) #nucat is the num of internal substrate molecules used by catabolism. em is the num of molecules of metabolic enzyme, vm is the rate of the metabolic enzyme and si is the internal substrate level
    #this is all over the half maximal metabolic rate threshold plus the internal substrate level
    fixation = 0#vmax_nit*a/(Kgamma +a)

    #AA equations
     #k_AA= 1000.0 #the constant that will determine how many AAs the cell is aiming to have
     #k_cat_AA= 1000.0 #rate of production of product of the AA enzyme when bound to both substrates
    # println("current k_cat_AA value inside ODE is ",k_cat_AA)
     #k_a_NH4 = 1000.0 #rate of binding of NH4 to atp-enzyme complex
     #k_a = 1000.0 #rate of binding of atp to enzyme
     #k_NH4 = 1000.0 #rate of binding of NH4 to enzyme

     #k_a_AA = 0.01 #rate of binding of the end product inhibitor to the enzyme-atp complex
     #k_NH4_AA= 0.01 #rate of binding of the end product acting as an inhibitor to the enzyme-NH4 complex

    AA_vo = ((k_cat_AA*a*NH4)/k_a_NH4*k_a)/(1+((1+(AA/k_a_AA)+(NH4/k_a_NH4)))*(a/k_a)+(1+(AA/k_NH4_AA)*(NH4/k_NH4)))

    new_AA = AA_vo*AA_prot

    #open("../data/param-sweep-$k_cat_AA-$k_a_NH4-$k_a-$k_NH4-$k_a_AA-$k_NH4_AA-$k_ribo_a-$k_ribo_AA-$k_ribo_a_AA-$k_ribo_AA_a.csv","a") do f
    #write(f, "\n$AA, $a, $N")
    #end

    AA_a_use = new_AA*2 #each AA produced uses x number of ATP molecules
    AA_NH4_use = new_AA*2 #each AA produced uses 2 NH4 molecules

    export_rate=0#0.3
    ################################################
    #export rate needs to vary with AA levels#
    ################################################

    ##################################################################

    results[1]= ((-et*vt*s_out/(Kt + s_out))*N)-ds*s_out +(export_rate*NH4*50) +1e7 #change in amount of substrate, a proxy for amount of host plant

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
    #####################################Rplot#####
    results[15]= fixation*nit-export_rate*NH4 -lam*NH4 -AA_NH4_use + 5e6
    #change in num of NH4 molecules in the cell from fixation rate - export rate
    results[16]= 0#(v_nit*a/(thetax + a))+(ku*nit_mrna_ribo)+(gam/nx*nit_mrna_ribo)-(kb*r*nit_mrna)-(dm*nit_mrna)-(lam*nit_mrna)
    #num of free nitrogenase coding mRNA
    results[17]= 0#(kb*r*nit_mrna)-(ku*nit_mrna_ribo)-(gam/nx*nit_mrna_ribo)-(lam*nit_mrna_ribo)
    #num of nit mRNA-ribosome complexes
    results[18]= 0 #(gam/nx)*nit_mrna_ribo-lam*nit
    #num of nitrogenase proteins
    results[19]= 0#(export_rate*NH4)#total num of NH4 molecules exported this timestep
    results[20]=lam*N-death_rate*N #the change in num of bacterial cells, instead of just loss to dilution
    results[21]= ns*nucat-ttrate-fixation*nit-AA_a_use-lam*a
    ##########################################
    ## AAA equations##
    ##########################################
    results[22]= new_AA -ttrate-lam*AA #num of AA
	#println("new AA is $new_AA")
	#AA_used = ttrate+lam*AA
	#println("AA used is $AA_used")
	#tempvar= results[22]
	#println("results 22 is $tempvar")
	####save the final values of these three variable to file###
    
    #open("../data/high-AA-params-file-$k_cat_AA-$k_a_NH4-$k_a-$k_NH4-$k_a_AA-$k_NH4_AA-$k_ribo_a-$k_ribo_AA-$k_ribo_a_AA-$k_ribo_AA_a.csv","a") do f
    #write(f, "\n$new_AA, $ttrate, $AA_prot")
    #end


    results[23]=(gam/nx)*AA_mrna_ribo-lam*AA_prot #num of AA making proteins


    results[24]=((w_AA*a/(thetax + a)))+(ku*AA_mrna_ribo)+
    (gam/nx*AA_mrna_ribo)-(kb*r*AA_mrna)-(dm*AA_mrna)-(lam*AA_mrna)#num of free AA mRNAs



    results[25]=(kb*r*AA_mrna)-(ku*AA_mrna_ribo)-(gam/nx*AA_mrna_ribo)-(lam*AA_mrna_ribo)#num of AA mRNAs bound to ribosomes


end
