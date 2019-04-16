#write function where input is the 8 parameters 
#the function should take the parameters, solve the 25 ODEs 
#output of function should be value of gamma (proxy for growth rate)
#this function should be used in the optimiser to find maximum growth rate


function gamma_finder(params)

    k_ribo_a= params[1]
    k_ribo_a_AA= params[2]
    k_ribo_AA= params[3]
    k_a= params[4]
    k_cat_AA= params[5]
    k_a_NH4= params[6]
    k_NH4= params[7]
    k_a_AA= params[8]
    k_NH4_AA= params[9]
    #unpack the 9 enzyme kinetic parameters that i want to optimise
    
    gam= (gmax*a*AA)/(k_ribo_a*k_ribo_a_AA+k_ribo_a_AA*a+k_ribo_AA_a*AA+a*AA)
    ttrate= (rmq + rmr + rmt + rmm + nit_mrna_ribo+AA_mrna_ribo)*gam
    lam= ttrate/M
    nucat= em*vm*si/(Km + abs(si))
    fixation = 0

    AA_vo = ((k_cat_AA*a*NH4)/k_a_NH4*k_a)/(1+((1+(AA/k_a_AA)+(NH4/k_a_NH4)))*(a/k_a)+(1+(AA/k_NH4_AA)*(NH4/k_NH4)))
    new_AA = AA_vo*AA_prot
    AA_a_use = new_AA*2 #each AA produced uses x number of ATP molecules
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

    end_gam= (gmax*a*AA)/(k_ribo_a*k_ribo_a_AA+k_ribo_a_AA*a+k_ribo_AA_a*AA+a*AA)

    
    
    
