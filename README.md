# Dual_Simplex
It is a bit modified implementation of dual simplex method. The amazing thing is that it works in all cases albeit in worst case a bit slower than Dual simplex method.

I did this as a part of my MTL103(course at IITD) assignment.

The main function is the simplex_algo() which when runs automatically reads the input in the "input.txt" file and returns a dictionary containing keys:
{'solution_status', 'optimal_value', 'optimal_solution', 'initial_tableau', 'final_tableau'}


In case of 'solution_status':infeasible/unbounded case:<br />
It only reports it.<br />
Assumption: variables are non-negative.
