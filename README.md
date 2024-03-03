# Dual_Simplex
It is a bit modified implementation of dual simplex method. The amazing thing is that it works in all cases albeit in worst case a bit slower than Dual simplex method.

I did this as a part of my MTL103(course at IITD) assignment.

The main function is the simplex_algo(filename) which when runs automatically reads the input in the filename passed and returns a dictionary containing:
{'solution_status', 'optimal_value', 'optimal_solution', 'initial_tableau', 'final_tableau'}
In case of infeasible/unbounded case:
It only reports it.<br />
Assumption: variables are non-negative.<br />
