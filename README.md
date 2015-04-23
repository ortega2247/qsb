#Project 4 (Data from Lopez, Sorger et al, 2013

##Questions (copied from the project descriptions in the course web page)

* Use several mathematical functions (e.g polynomial, exponential) to fit the data. How do they perform
* Assume you will approximate this system with a few enzymatic reactions. Could you write a set of ODEs and integrate them to encode a system with similar behavior? Can you create a "best" model based on goodness-of-fit?
* How do model parameters affect the signal output? Could you make predictions that you could test experimentally based on a model?
* Is your model prognostic, or diagnostic? How would you test it?
* Use a machine-learning approach or data-driven method to derive a minimal model topology and parameters that explain the observed behavior. What could you learn?
* Encode different reaction network motifs and generate the ODEs that behave with a similar behavior as the data. What can you learn?

##I propose that we start putting information here about the status of the questions, and add more questions if necessary

##Proposed "experiments"
Here, I will start adding a few ideas

* Randomly perturb (use parallel MCMC, maybe) the rates at random sites on the system of ODEs and compare the perturbed system with the original one (or with the raw data?) to examine critical parameters

* I think question 1 is already answered. I think we could do a hill equation fit and that could tell us somethings about cooperativity in the reactions.

* Also, about the perturbations, i think we should perturbe the system and compare the new outputs with the output of the model that fits best the experimentl values. a
