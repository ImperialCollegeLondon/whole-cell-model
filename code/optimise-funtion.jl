#write function where input is the 8 parameters 
#the function should take the parameters, solve the 25 ODEs 
#output of function should be value of gamma (proxy for growth rate)
#this function should be used in the optimiser to find maximum growth rate



function gamma_finder(one,two,three,four,five,six,seven,eight,nine)
    #s_out_0,rmr_0,em_0,rmq_0,rmt_0,et_0,rmm_0,
    #mt_0,mm_0,q_0,si_0,mq_0,mr_0,r_0,
    #NH4_0,nit_mrna_0,nit_mrna_ribo_0,nit_0,exported_0,N_0,a_0, AA_0, 
    #AA_prot_0,AA_mrna_0,AA_mrna_ribo_0)
    #input should be three arrays, the first with the 8 parameter values that 
#we are varying and the second array should be the current values of the 25 molecules
    


    k_ribo_a= one#params[1]
    k_ribo_a_AA= two#params[2]
    k_ribo_AA_a= three#params[3]
    k_a= four#params[4]
    k_cat_AA= five#params[5]
    k_a_NH4= six#params[6]
    k_NH4= seven#params[7]
    k_a_AA= eight#params[8]
    k_NH4_AA= nine#params[9]
    #unpack the 9 enzyme kinetic parameters that i want to optimise
    
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
    
    s_out = s_out_0 #external level of substrate, a proxy for plant growth state since more plant = more substrate
    rmr= rmr_0#num of ribosome/ribosomal mRNA complex molecules
    em= em_0#metabolic enzyme molecumes
    rmq= rmq_0#num of ribosome/housekeeping mRNA complex molecules
    rmt= rmt_0#num of ribosome/transporter protein mRNA complex molecules
    et= et_0#num of transporter enzyme molecules
    rmm= rmm_0#num of ribosome/metabolic enzyme mRNA complexes
    mt= mt_0#num of particles of transporter protein mRNA
    mm= mm_0#num of particles of metabolic enzyme mRNA
    q= q_0#num of housekeeping protein molecules
    si= si_0#num of substrate molecules inside the cell (internal)
    mq= mq_0#num of mRNA coding for housekeeping proteins
    mr= mr_0#num of mRNA coding for ribosomes
    r= r_0#num of free ribosomes
    NH4 = NH4_0 #num of ammonia molecules fixed
    nit_mrna = nit_mrna_0#num of molecules of nitrogen fixing protein
    nit_mrna_ribo = nit_mrna_ribo_0#num of molecules of mRNA coding for nitrogen fixing protein
    nit = nit_0#num of molecules of mRNA-ribosome complexes coding for nitrogenase
    exported= exported_0#total amount of NH4 exported
    N = N_0 #num of bacteria cells
    a= a_0#num of ATP (proxy for level of energy)
    AA = AA_0 #num of amino acid molecules
    AA_prot= AA_prot_0#num of protein molecules that make new amino acids
    AA_mrna= AA_mrna_0#num of mRNA coding for AA making protein
    AA_mrna_ribo= AA_mrna_ribo_0#num of mRNA coding for AA making protein bound to ribosomes
    
    results = [s_out_0,rmr_0,em_0,rmq_0,rmt_0,et_0,rmm_0,
    mt_0,mm_0,q_0,si_0,mq_0,mr_0,r_0,
    NH4_0,nit_mrna_0,nit_mrna_ribo_0,nit_0,exported_0,N_0,a_0, AA_0, AA_prot_0,AA_mrna_0,AA_mrna_ribo_0]
    #making an array to record the results in
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

    
    
    gmax = 1260.0
    gam= (gmax*a*AA)/(k_ribo_a*k_ribo_a_AA+k_ribo_a_AA*a+k_ribo_AA_a*AA+a*AA)
    ttrate= (rmq + rmr + rmt + rmm + nit_mrna_ribo+AA_mrna_ribo)*gam
    lam= ttrate/M
    nucat= em*vm*si/(Km + abs(si))
    fixation = 0

    AA_vo = ((k_cat_AA*a*NH4)/k_a_NH4*k_a)/(1+((1+(AA/k_a_AA)+(NH4/k_a_NH4)))*(a/k_a)+(1+(AA/k_NH4_AA)*(NH4/k_NH4)))
    new_AA = AA_vo*AA_prot
    AA_a_use = new_AA*2
    AA_NH4_use = new_AA*2
    export_rate=0


    results[1]= ((-et*vt*s_out/(Kt + s_out))*N)-ds*s_out +(export_rate*NH4*50) +1e7 #change in amount of substrate, a proxy for amount of host plant
    results[2]= kb*r*mr-ku*rmr-(gam/nr)*rmr-lam*rmr#-dp*r
    results[3]= (gam/nx)*rmm-lam*em#-dp*em
    results[4]= kb*r*mq-ku*rmq-(gam/nx)*rmq-lam*rmq
    results[5]= kb*r*mt-ku*rmt-gam/nx*rmt-lam*rmt
    results[6]= (gam/nx)*rmt-lam*et#-dp*et
    results[7]= kb*r*mm-ku*rmm-gam/nx*rmm-lam*rmm
    results[8]= (we*a/(thetax + a))+ku*rmt+(gam/nx)*rmt-kb*r*mt-dm*mt-lam*mt
    results[9]= (we*a/(thetax + a))+ku*rmm+(gam/nx)*rmm-kb*r*mm-dm*mm-lam*mm
    results[10]= (gam/nx)*rmq-lam*q#-dp*q
    results[11]= (et*vt*s_out/(Kt + s_out))-nucat-lam*si#-degredation*si
    results[12]= (wq*a/(thetax + a)/(1 + (q/Kq)^nq))+ku*rmq+(gam/nx)*rmq-kb*r*mq-dm*mq-lam*mq
    results[13]= (wr*a/(thetar + a))+ku*rmr+(gam/nr)*rmr-kb*r*mr-dm*mr-lam*mr
    results[14]= ku*rmr +ku*rmt+ ku*rmm+ ku*rmq +(gam/nr)*rmr+(gam/nr)*rmr +(gam/nx)*rmt +(gam/nx)*rmm +(gam/nx)*rmq -kb*r*mr- kb*r*mt- kb*r*mm -kb*r*mq +ku*nit_mrna_ribo +(gam/nx)*nit_mrna_ribo-kb*r*nit_mrna +ku*AA_mrna_ribo +(gam/nx)*AA_mrna_ribo-kb*r*AA_mrna -lam*r
    results[15]= fixation*nit-export_rate*NH4 -lam*NH4 -AA_NH4_use + 3e6
    results[16]= 0#(v_nit*a/(thetax + a))+(ku*nit_mrna_ribo)+(gam/nx*nit_mrna_ribo)-(kb*r*nit_mrna)-(dm*nit_mrna)-(lam*nit_mrna)
    results[17]= 0#(kb*r*nit_mrna)-(ku*nit_mrna_ribo)-(gam/nx*nit_mrna_ribo)-(lam*nit_mrna_ribo)
    results[18]= 0 #(gam/nx)*nit_mrna_ribo-lam*nit
    results[19]= 0#(export_rate*NH4)#total num of NH4 molecules exported this timestep
    results[20]=lam*N-death_rate*N #the change in num of bacterial cells, instead of just loss to dilution
    results[21]= ns*nucat-ttrate-fixation*nit-AA_a_use-lam*a
    results[22]= new_AA -ttrate-lam*AA #num of AA
    results[23]=(gam/nx)*AA_mrna_ribo-lam*AA_prot #num of AA making proteins
    results[24]=((w_AA*a/(thetax + a)))+(ku*AA_mrna_ribo)+(gam/nx*AA_mrna_ribo)-(kb*r*AA_mrna)-(dm*AA_mrna)-(lam*AA_mrna)#num of free AA mRNAs
    results[25]=(kb*r*AA_mrna)-(ku*AA_mrna_ribo)-(gam/nx*AA_mrna_ribo)-(lam*AA_mrna_ribo)#num of AA mRNAs bound to ribosomes

    #need to add the change back to the original value to get the actual value after that timestep
    #which can then be used to calculate the growth rate after that step
    a = a + results[21]
    AA= AA + results[22]
    rmq = rmq + results[4]
    rmr = rmr + results[2]
    rmt = rmt + results[5]
    rmm = rmm + results[7]
    nit_mrna_ribo = nit_mrna_ribo + results[17]
    AA_mrna_ribo = AA_mrna_ribo + results[25]

    
    end_gam= (gmax*a*AA)/(k_ribo_a*k_ribo_a_AA+k_ribo_a_AA*a+k_ribo_AA_a*AA+a*AA)
    end_ttrate= (rmq + rmr + rmt + rmm + nit_mrna_ribo+AA_mrna_ribo)*end_gam
    end_lam= end_ttrate/M
#~     println("new gamma is $end_gam")
#~     println("new ttrate is $end_ttrate")
#~     println("new lambda is $end_lam")

#~     output = [end_lam, results]
#~     return(end_lam)


    end  
    
