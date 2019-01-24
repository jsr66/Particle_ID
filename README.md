# Particle-ID
The primary purpose of this project is to identify different species of particle emerging from collisions in a particle accelerator. The project consists of three files: ParticleID.cpp, LikelihoodFcn.cpp, MonteCarloHitGenerator.cpp.

The main file is ParticleID.cpp. This code takes as input data from an unknown track in a particle detector and outputs the probability that the particle is an electron, pion, muon, kaon, or proton. More technical details are given in the file comments. 

The file LikelihoodFcn.cpp provides a means of aggregrating measurements of ionization from a single track in a particle detector (which indirectly enable us to determine the particle's velocity through what is known as the Universal Curve). Specifically, it returns the MPV of the probability distribution (known as the Landau distribution) that most likely produced the multiple ionization measurements associated with the track through Maximum Likelihood fitting. This MPV enables us to meaningfully assign a single value for the ionization to the particle track, and, ultimately, to determine its velocity through the Universal Curve, which relates ionization deposits left by a track to particle velocity (more precisely, the relativistic factor betagamma). More technical details are given in the file comments. 

The file MonteCarloHitGenerator.cpp simulates particle track measurements using Monte Carlo Methods. Specifically, it draws particle track measurements randomly from a distribution whose MPV has been specified by the user. This provides a useful means of testing LikelihoodFcn.cpp and ParticleID.cpp: after inputting simulated measurements drawn from a known distribution into LikelihoodFcn.cpp, we should expect to recover as output an MPV that closely approximates the MPV used generate those measurements. More technical details are given in the file comments. 


