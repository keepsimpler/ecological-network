Code and Data accompanying the paper:
-------- 
Wenfeng Feng & Kazuhiro Takemoto, Heterogeneity in ecological mutualistic networks dominantly determines community stability, 2014.

###License:
The code is made available under the GNU GPL license:
http://www.gnu.org/licenses/gpl.html

###Requests:
Please send any request/bug/issue/question to
Wenfeng Feng fengwenfeng@gmail.com

###Citation:
If you use the code for your work, please reference the paper:

Wenfeng Feng & Kazuhiro Takemoto, Heterogeneity in ecological mutualistic networks dominantly determines community stability, 2014.


###Structure:
All the code and data files are for the statistical software R (http://cran.us.r-project.org/).

**plustwins.functions.R**

    the functions needed to estimate dominant eigenvalues of bipartite networks, and to measure structural features of null models.

**plustwins.test.R**

    the test scripts for the figures of paper.

**datasets.RData**

    40 empirical mutualistic networks from literatures.

They require the following packages:
>igraph, bipartite, plyr, reshape2, doMC, MASS, ggplot2, akima, etc.


###Usage (in R):
```
load('datasets.RData')
source('plustwins.functions.R')
source('plustwins.test.R')
```

For details on the use, load the files in R, or browse the code.

**Note:** The scripts in *plustwins.test.R* are very time consuming, you may need to evaluate the scripts step by step.

