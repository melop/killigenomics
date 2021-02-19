GENERAL-INFO-START

		seq-file		../PLP.gphocs.1000.50000.3.in
		trace-file		mcmc-trace.out
		num-loci		13509
		random-seed		23333
		burn-in		0
		mcmc-iterations		300000
		mcmc-sample-skip		0
		start-mig		0
		iterations-per-log		100
		logs-per-line		100

		tau-theta-print		10000
		tau-theta-alpha		1
		tau-theta-beta		10000

		mig-rate-print		0.001
		mig-rate-alpha		0.002
		mig-rate-beta		0.00001

		locus-mut-rate		VAR		1.0

		find-finetunes		TRUE
		find-finetunes-num-steps		100
		find-finetunes-samples-per-step		100


GENERAL-INFO-END


CURRENT-POPS-START

		POP-START
				name		B
				samples		PLP_B_PLP-B3_D504_D704 d 
		POP-END

		POP-START
				name		C
				samples		PLP_C_PLP-C1_D506_D706 d PLP_C_PLP-C2_D507_D707 d PLP_C_PLP-C3_D508_D708 d PLP_C_PLP-C4_D502_D701 d 
		POP-END

		POP-START
				name		A
				samples		PLP_A_PLP-A3_D501_D701 d PLP_A_PLP-A4_D502_D702 d 
		POP-END

		POP-START
				name		F
				samples		PLP_F_PLP-F1_D503_D701 d PLP_F_PLP-F2_D504_D702 d PLP_F_PLP-F3_D505_D703 d 
		POP-END

		POP-START
				name		D
				samples		PLP_D_PLP-D0_ref d PLP_D_PLP-D1_D503_D702 d 
		POP-END

		POP-START
				name		E
				samples		PLP_E_PLP-E1_D505_D704 d PLP_E_PLP-E2_D506_D705 d PLP_E_PLP-E3_D507_D706 d PLP_E_PLP-E4_D508_D707 d PLP_E_PLP-E5_D501_D708 d 
		POP-END

		POP-START
				name		G
				samples		PLP_G_PLP-G1_D506_D704 d PLP_G_PLP-G3_D508_D706 d PLP_G_PLP-G5_D502_D708 d 
		POP-END

CURRENT-POPS-END


ANCESTRAL-POPS-START

		POP-START
				name		Pras_N
				children		B		C
				tau-initial		0.001
		POP-END

		POP-START
				name		Praslin
				children		Pras_N		A
				tau-initial		0.003
		POP-END

		POP-START
				name		Praslin_Curieuse
				children		Praslin		F
				tau-initial		0.005
		POP-END

		POP-START
				name		PCD
				children		Praslin_Curieuse		D
				tau-initial		0.007
		POP-END

		POP-START
				name		Mahe
				children		E		G
				tau-initial		0.007
		POP-END

		POP-START
				name		root
				children		PCD		Mahe
				tau-initial		0.01
		POP-END

ANCESTRAL-POPS-END


MIG-BANDS-START

		BAND-START
				source		D
				target		Pras_N
		BAND-END

MIG-BANDS-END



