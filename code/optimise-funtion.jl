#write function where input is the 8 parameters 
#the function should take the parameters, solve the 25 ODEs 
#output of function should be value of gamma (proxy for growth rate)
#this function should be used in the optimiser to find maximum growth rate



function gamma_finder(one,two,three,four,five,six,seven,eight,nine,
    s_out,rmr_0,em_0,rmq_0,rmt_0,et_0,rmm_0,
    mt_0,mm_0,q_0,si_0,mq_0,mr_0,r_0,
    NH4_0,nit_mrna_0,nit_mrna_ribo_0,nit_0,exported_0,N_0,a_0, AA_0, 
    AA_prot_0,AA_mrna_0,AA_mrna_ribo_0)
    #input should be three arrays, the first with the 8 parameter values that 
#we are varying and the second array should be the current values of the 25 molecules
    
    
    results = [s_out,rmr_0,em_0,rmq_0,rmt_0,et_0,rmm_0,
    mt_0,mm_0,q_0,si_0,mq_0,mr_0,r_0,
    NH4_0,nit_mrna_0,nit_mrna_ribo_0,nit_0,exported_0,N_0,a_0, AA_0, AA_prot_0,AA_mrna_0,AA_mrna_ribo_0]
    #making an array to record the results in

    k_ribo_a= one#params[1]
    k_ribo_a_AA= two#params[2]
    k_ribo_AA= three#params[3]
    k_a= four#params[4]
    k_cat_AA= five#params[5]
    k_a_NH4= six#params[6]
    k_NH4= seven#params[7]
    k_a_AA= eight#params[8]
    k_NH4_AA= nine#params[9]
    #unpack the 9 enzyme kinetic parameters that i want to optimise
    
#~     s_out = y[1] #external level of substrate, a proxy for plant growth state since more plant = more substrate
#~     rmr= y[2]#num of ribosome/ribosomal mRNA complex molecules
#~     em= y[3]#metabolic enzyme molecumes
#~     rmq= y[4]#num of ribosome/housekeeping mRNA complex molecules
#~     rmt= y[5]#num of ribosome/transporter protein mRNA complex molecules
#~     et= y[6]#num of transporter enzyme molecules
#~     rmm= y[7]#num of ribosome/metabolic enzyme mRNA complexes
#~     mt= y[8]#num of particles of transporter protein mRNA
#~     mm= y[9]#num of particles of metabolic enzyme mRNA
#~     q= y[10]#num of housekeeping protein molecules
#~     si= y[11]#num of substrate molecules inside the cell (internal)
#~     mq= y[12]#num of mRNA coding for housekeeping proteins
#~     mr= y[13]#num of mRNA coding for ribosomes
#~     r= y[14]#num of free ribosomes
#~     NH4 = y[15] #num of ammonia molecules fixed
#~     nit_mrna = y[16]#num of molecules of nitrogen fixing protein
#~     nit_mrna_ribo = y[17]#num of molecules of mRNA coding for nitrogen fixing protein
#~     nit = y[18]#num of molecules of mRNA-ribosome complexes coding for nitrogenase
#~     exported= y[19]#total amount of NH4 exported
#~     N = y[20] #num of bacteria cells
#~     a= y[21]#num of ATP (proxy for level of energy)
#~     AA = y[22] #num of amino acid molecules
#~     AA_prot=y[23]#num of protein molecules that make new amino acids
#~     AA_mrna=y[24]#num of mRNA coding for AA making protein
#~     AA_mrna_ribo=y[25]#num of mRNA coding for AA making protein bound to ribosomes
    
    
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
    rmm = rmm + resutls[7]
    nit_mrna_ribo = nit_mrna_ribo + results[17]
    AA_mrna_ribo = AA_mrna_ribo + results[25]


    end_gam= (gmax*a*AA)/(k_ribo_a*k_ribo_a_AA+k_ribo_a_AA*a+k_ribo_AA_a*AA+a*AA)
    end_ttrate= (rmq + rmr + rmt + rmm + nit_mrna_ribo+AA_mrna_ribo)*end_gam
    end_lam= end_ttrate/M


#~     output = [end_lam, results]
    return(end_lam)


    end  
    
