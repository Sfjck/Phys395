Assignment 5
Ensure the following empty folders exist in the directory. They should be part of the archive already.
Q1Out/
Q2Out/
Q3Out/
Q4Out/
Q5Out/
They are used to store the (many) output files from integrations. They must be present or the code will crash.

Then use this command:
bash runme.sh

For outputs, look at the following files:
output
prob1.pdf
prob2.pdf
prob3.pdf
prob4.pdf
prob5.pdf

Some explanations:
For harmonic oscillator, eigenvalues should be (n+1/2)hw, n=0,1,2...
For anharmonic oscillator, I couldn't find any literature on expected eigenvalues in higher energies with h defined as 1. But the ground state is:
E0 = 1.060 * (lambda/4)^(1/3) * (h^2 / 2m)^(2/3)
Plugging in lambda = h = m = 1 gives ~0.42
(Reference: A Modern Approach to Quantum Mechanics 2e, Townsend, p449)
Of course, all eigenvalues from bisection can be compared to results from Rayleigh.
