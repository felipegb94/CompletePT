# PermTestingToolbox
CPU and GPU implementations of Permutation Testing using the libraries Armadillo and Arrayfire. The permutation testing schemes implemented herein are based on two-sample and one-sample t-test. 

## Setup 

## Background OLD
### General Linear Model
GLM models observed data (dependent variable) as a linear combination
of predictor variables (independent variables, covariates etc). `Y (observed data matrix), X (design matrix, predictor variables), b (parameter estimate vector), and e (error vector).`

```
Y = X * b + e
```

Many statistical techniques are special cases of the general linear model. For example:

* ANOVA asks whether different experimental conditions (X1, X2, etc)
are associated with different scores.

* Multiple regression asks whether scores are related to predictor
variables (X1, X2, etc) 

* T-test is a special case of ANOVA where there are only two groups X1 and X2

The three of them are asking the same question. Is there a relationship between  a dependent variable (Yi) and one or more independent variables (Xi).

### T-tests

__1. Paired Sample t-test:__ Used to compare two population means in the case of two samples that are correlated. Commonly used in 'before and after' studies, case-control studies.The typical workflow for two-sample t-test hypothesis testing is:

__2. One-Sample t-test:__ The one sample t test is an appropriate analysis when the research looks to compare the mean of a sample with a hypothesized mean to assess if differences occur. The assumptions of the one sample t test include: the data must be normally distributed within the population and the data should be independent; scores of one participant are not dependent upon scores of another.

__3. Two-Sample t-test:__ Used to determine if two population means are equal.

### Hypothesis Testing
Hypothesis testing is a group of techniques in Statistics that are often used in functional magnetic resonance imaging (fMRI) data to identify areas in the brain that display statistical significant activity. So how do we classify a voxel as statistically significant?

  *1. Select univariate test-statistic:* The job of this test statistic is to act as the mapping from data to a detection threshold.
   
  *2. Hypothesis setup:* Setup two hypothesis. The null hypothesis (H0) says that there is no mean difference between both samples. The alternate hypothesis says that they are different.   
  
  *3. Select Significance Level:* Choose the significance level. Usually 5% in most studies and 1% in medical studies. This is usually denoted as alpha, and it tells us the probability of making a Type I error; that is the probability of deciding erroneously on the alternative when, in fact, the null hypothesis is true.   
  
  *4. Calculate the Parameter:* To calculate the parameter using the fo 
  
  *5. Testing of hypothesis:* Compare result to table value.



### fMRI
In placebo-controlled clinical trials, the goal is to assess the improvements of the treatment, if any. When evaluation treatments for neurodegenerative disorders, fMRI data can be used to aid the evaluation of the treatment. Since we are trying to detect if there is any difference between the placebo and the control fMRI data a two-sample t-test would be an appropriate test.


## Setup

### OSX

Download homebrew

### Ubuntu 

Download the C++ linear algebra library, Armadillo. Do NOT use `apt-get` to install Armadillo, this will install a very old version (4.2). This software has been tested against Armadillo 5.*, 6.1 and 6.2. So go to http://arma.sourceforge.net/download.html download the latest stable version and follow the `README.txt` instruction on how to install Armadillo. It is VERY important that you install all the requirements (section 3 of the `README.txt`) before you proceed to build and install Armadillo (section 4 of the `README.txt`).
