# Recursive Inversion Model
 
The provided code gives all the tools needed to create, sample from, fit, and find marginal distributions for recursive inversion models (Link to paper to come when publication is finalized). The main method provides some tools for producing basic results that can be fine tuned for simple testing. These methods also provide examples of the different functions and variables needed to produce most necessary values for anyone interested in more involved testing.

The following three tasks can be easily completed from the command line, and should provide most the necessary resources for most basic Recursive Inversion Model functions, fitting to data, reproducing known synthetic models, or finding marginal distributions of rank or inversion matricies. 

- Reproducability testing
Example commandline arguments:
`-synthtest n=12 thetalow=.5 thetahigh=2 Ntrain=100 iters=150 temp=.02`
* -synthtest will produce a random tree parameterized by random values, sample from that tree, and attempt to reproduce it. (Note that in rare instances this method can return a non-canonical, but valid tree)
  * n (Integer, Optional, default n=15) The number of items to be ordered in the model.
  * thetalow (Positive double, Optional, default thetalow=.5) The lowest value of randomized θ (sampled uniformly) to be used in the model.
  * thetahigh (Positive double, Optional, default thetahigh=2) The highest value of randomized θ (sampled uniformly) to be used in the model.
  * Ntrain (Positive interger, Optional, default Ntrain=5000) The number of rankings to sample from the synthetic model for refitting.
  * iters (Positive integer, Optional, default iters=1000) The number of steps through the model fit process, i.e. the number of unique rank orderings tested.
  * temp (Positive double [0,1], Optional, default temp=.02) The 'temperature' used for simualted annealing process.

`-permtest file=data/sushi/SushiClean.txt ptest=.3 pvalidation=.2 iters=1000 temp=.03 dataSeed=400 runtimeSeed=400`
* -permtest will import data from a CSV containing rankings of integers, split the data into (up to) 3 parts for training, testing, and validation, then fit the best model possible to the training data set and provide the resultant log-likelihood for all three subsets of the data.
  * file (String, REQUIRED) The file location for a CSV of rankings of integers. (Please see note below concerning the Sushi dataset)
  * ptest (Positive double [0,1], Optional, default ptest=0) The proportion of data to set aside for testing.
  * pvalidation (Positive double [0,1], Optional, default pvalidation=0) The proportion of data to set aside for validation. Remaining data not used for testing or validation is used for training, thus, pvalidation+ptest must be less than 1.
  * iters (Positive integer, Optional, default iters=1000) The number of steps through the model fit process, i.e. the number of unique rank orderings tested.
  * temp (Positive double [0,1], Optional, default temp=.02) The 'temperature' used for simualted annealing process.
  * dataSeed (Integer, Optional, randomized by default) The seed used when splitting the data into training/testing/validation, for reproducing results. Note that if ptest=pvalidation=0, dataSeed serves no purpose.
  * runtimeSeed (Integer, Optional, randomized by default) The seed used to choose new orderings to fit the model to at runtime, for reproducing results.

`-importrim rim=([0.7],0,([1.1],([0.9],1,([0.2],([1.0],2,3),4)),([0.3],([0.9],([0.8],5,6),7),([1.0],8,9))))`
* -importrim will import the provided model and produce the marginal rank and inversion matricies. Intended to mostly serve as an example of how to import a model.
  * rim (String, REQUIRED) The model to imported, in the format ([theta],node,node) where a node is of the form ([theta],node,node) or a unique integer for a leaf node.


A note on the Sushi dataset: The authors of these methods would like to thank the creators of the Sushi rank dataset (found at https://www.kamishima.net/sushi/) which was used in testing and studying these models. Per the request of the data set authors, the data is not distributed alongside this code and must be separately downloaded and preprocessed using the ReformatSushi.R script found within the data/sushi folder. Consult the RIM-README-Sushi.txt file found therein for more details.

A note on future developments: Anyone digging into the code further may discover that it contains functionality not accessed by the examples here. The code authors have developed Recursive Inversion Models further than what is available here, and plan to introduce more methods for different functionality to the public respository over time.
