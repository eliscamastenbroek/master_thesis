## True proportions in simulated data
true_proportions = c(0.6061153, 0.2577143, 0.1361704)
names(true_proportions) = c("Permanent","Other","Flexible")

#Create ME matrices in global environment
ME_matrix1 = create_ME_matrix(-log(18),-log(18),log(324),log(324),log(18),log(18)) #10% ME
ME_matrix2 = create_ME_matrix(-3*log(2),-3*log(2),6*log(2),6*log(2),3*log(2),3*log(2)) #20% ME
ME_matrix3 = create_ME_matrix(-1.54045,-1.54045,3.0809,3.0809,1.54045,1.54045) #30% ME
ME_matrix4a = create_ME_matrix(-4.4917,-4.94368,7.64123,5.6678,3.09876,4.55064) #Realistic 7% ME
ME_matrix4b = create_ME_matrix(-6.14311,-2.63157,11.9482,5.03275,1.72427,1.53296) #Realistic 7% ME
ME_matrix4 = list(ME_matrix4a,ME_matrix4b) #Realistic 7% ME
