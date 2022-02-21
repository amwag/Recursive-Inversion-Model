# Recursive Inversion Model
 
The provided code gives all the tools needed to create, sample from, fit, and find marginal distributions for recursive inversion models (Link to paper to come when publication is finalized). The main method provides some tools for producing basic results that can be fine tuned for simple testing. These methods also provide examples of the different functions and variables needed to produce most necessary values for anyone interested in more involved testing.

The following three tasks can be easily completed from the command line, and should provide most the necessary resources for most basic Recursive Inversion Model functions, fitting to data, reproducing known synthetic models, or finding marginal distributions of rank or inversion matricies. 

- Reproducability testing
Example commandline arguments:
`-synthtest n=12 thetalow=.5 thetahigh=2 Ntrain=100 iters=150 temp=.02`
* -synthtest will produce a random tree parameterized by random values, sample from that tree, and attempt to reproduce it. (Note that in rare instances this method can return a non-canonical, but valid tree)
 * n - (Integer, Optional, default n=15) the number of items to be ordered in the model.
 * thetalow (Positive double, Optional, default thetalow=.5) The lowest value of randomized θ (sampled uniformly) to be used in the model.
 * thetahigh (Positive double, Optional, default thetahigh=2) The highest value of randomized θ (sampled uniformly) to be used in the model.
 * Ntrain (Positive interger, Optional, default Ntrain=5000) The number of rankings to sample from the synthetic model for refitting.
 * iters (Positive integer, Optional, default iters=1000) The number of steps through the model fit process, i.e. the number of unique rank orderings tested.
 * temp (Positive double [0,1], Optional, default temp=.02) The 'temperature' used for simualted annealing process.



This 
