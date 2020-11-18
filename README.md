# MOPSO
An implementation of multi objective particle swarm optimization technique for a minimization problem

[![View MOPSO on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://se.mathworks.com/matlabcentral/fileexchange/70021-mopso)

In this project there are two ways that I have implemented nonlinear constraints:
1. In MOPSO1 constraints are computed with objectives in one file and a zero or positive infeasabilty value is assinged to a particle where zero means it is feasable.
1. In MOPSO2 constraints are in other function and if a particle is infeasable the onjective function will not be evaluated in order to save computation time, and objective values of that particle will be NaN.
