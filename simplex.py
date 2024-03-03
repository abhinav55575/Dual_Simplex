import numpy as np
from numpy.linalg import inv,matrix_rank  # Matrix inverse
rank=matrix_rank

epsilon = 10**(-10)


def read_input(filename):
    #Assumed constraints have at max 4 decimal places.
    #Please change if required for accuracy.
    maxdec=4
    objective = None
    A = []
    b = []
    constraint_types = []
    c = []
    with open(filename, 'r') as file:
      data = file.read().replace('\n', '')
    #if no decimal is present we can do maxdec=0.
    if "=" in data:  
      ind=data.find("=")
    
      if "." not in data[:ind]:
        maxdec=0

    with open(filename, 'r') as file:
        lines = file.readlines()
        
        idx = 0
        cnt = 0
        while idx < len(lines):
            line = lines[idx].strip()
            if "[" in line:
                cnt += 1
                idx += 1
                curr = []
                while idx < len(lines):
                    line = lines[idx].strip()
                    if "[" in line:
                        break
                    line = lines[idx].strip()
                    if line != "":
                        curr += [line]
                    idx += 1
                if cnt == 1:
                    objective = curr[0].strip()
                elif cnt >= 2 and cnt <= 2:
                    for u in curr:
                        v = u.strip()
                        A.append([float(h)*10**maxdec for h in v.split(",")])
                elif cnt >= 3 and cnt <= 3:
                    for u in curr:
                        v = u.strip()
                        b.append(float(v)*10**maxdec)
                elif cnt == 4:
                    for u in curr:
                        v = u.strip()
                        constraint_types.append(v)
                else:
                    c = [float(h) for h in curr[0].strip().split(",")]
            else:
                idx += 1
    nA=[]
    nb=[]
    
    cnt=0
    st=len(A[0])
    for u in constraint_types:
        if u!="=":
            cnt+=1
            
    for i in range(len(constraint_types)):
        if constraint_types[i]=="=":
            nA.append(A[i]+[0]*cnt)
            nb.append(b[i])
            
    
    for i in range(len(constraint_types)):
        if constraint_types[i]=="<=":
            nA.append(A[i]+[0]*cnt)
            nA[-1][st]=1
            st+=1
            nb.append(b[i])
            c.append(0)
    for i in range(len(constraint_types)):
        if constraint_types[i]==">=":
            nA.append(A[i]+[0]*cnt)
            nA[-1][st]=-1
            st+=1
            nb.append(b[i])
            c.append(0)
    flag=0        
    if objective=="maximize":
        flag=1
        c=[-u for u in c]

    return nA,nb,c,flag,cnt     
            
    
            
            
    

def simplex(A, b, c, flag,cnt):
    """
    Outer "wrapper" for executing the simplex method: phase I and phase II.

    :param A: constraint matrix
    :param b: independent terms in constraints
    :param c: costs vector
    

    This function prints the outcome of each step to stdout.
    
    """
    c=np.array(c)
    
    A=np.array(A)
    b=np.array(b)
    

    m, n = A.shape[0], A.shape[1]
    b1=list(b)
    

    
    """setup"""
    answ=dict()
    answ["solution_status"]=None
    answ["optimal_value"]=None
    answ["optimal_solution"]=None
    answ["initial_tableau"]=None
    answ["final_tableau"]=None
    
    A_I = np.concatenate((A, np.identity(m)), axis=1)  #  constraint matrix
    x_I = np.concatenate((np.zeros(n), b))  #  variable vector
    #changed the cost to c+[1]*m,m=no of aux variables
    c_I = np.concatenate((c, np.ones(m)))  #  c_j vector
    basic_I = set(range(n, n + m))
    

    """execution"""
    
    
    ext, x, basic, z, d, it_II = simplex_core(A_I, c_I, x_I, basic_I,answ,b1,1)
    infeasible_check=trunc(sum(x[n:n+m]))
    if infeasible_check>0:
        answ["solution_status"]="infeasible"
        return answ

    

    if ext == 0:
        answ["solution_status"]="optimal"
        if flag==1:
            answ["optimal_value"]=-z
        else:    
            answ["optimal_value"]=z
        answ["optimal_solution"]=list(x)[:n-cnt]
        
        return answ
    elif ext == 1:
        
        answ["solution_status"]="unbounded"
        return answ
        

    

    return ext, x, z, d


def simplex_core(A, c, x, basic, answ=None,b1=None,flu=0):
    """
    This function executes the simplex algorithm iteratively until it
    terminates. It is the core function of this project.

    :param A: constraint matrix
    :param c: costs vector
    :param x: initial BFS
    :param basic: initial basic index set
    
    :return: a tuple consisting of the exit code, the value of x, basic index set,
    optimal cost (if optimum has been found), and BFD corresponding to
    feasible ray (if unlimited problem)
    """

    m, n = A.shape[0], A.shape[1]  # no. of rows, columns of A, respectively

    B, N = list(basic), set(range(n))  # Basic /nonbasic index lists
    for u in B:
        N.remove(u)
    del basic  
    B_inv = inv(A[:, B])  # Calculate inverse of basic matrix (`A[:, B]`)

    z = np.dot(c, x)  # Value of obj. function
    prices = c[B] @ B_inv
    


    it = 1  # Iteration number
    if flu==1:
            a2=[(c[q],q) for q in N]
            na2=[0]*n
            for u in a2:
                na2[u[1]]=u[0]
                
            a1=-prices@b1
            a3=B_inv@b1
            a4=B_inv@A
            las_tab=[]
            las_tab.append([a1]+na2)
            
            
            for u in range(len(a3)):
                las_tab.append([a3[u]]+list(a4[u]))
            answ["initial_tableau"]=las_tab   
                
    
                
    while it <= 500:
        r_q, q, p, theta, d = None, None, None, None, None


        """Optimality test"""
        prices = c[B] @ B_inv  # Store product for efficiency

        
        optimum = True
            
        for q in N:  # Read in lexicographical index order
                r_q = c[q] - prices @ A[:, q]
                if r_q < 0:
                    optimum = False
                    break    
        if flu==1:
            #print("space")
            a2=[(c[q] - prices @ A[:, q], q) for q in N]
            na2=[0]*n
            for u in a2:
                na2[u[1]]=u[0]
                
            a1=-prices@b1
            a3=B_inv@b1
            a4=B_inv@A
            las_tab=[]
            las_tab.append([a1]+na2)
            
            
            for u in range(len(a3)):
                las_tab.append([a3[u]]+list(a4[u]))
            answ["final_tableau"]=las_tab    
                

        if optimum:
            #print("\tfound optimum")
            return 0, x, set(B), z, None, it  # Found optimal solution
    


        """Feasible basic direction"""
        d = np.zeros(n)
        for i in range(m):
            d[B[i]] = trunc(-B_inv[i, :] @ A[:, q])
        d[q] = 1
        #print(prices @ A[:, q])


        """Maximum step length"""
        # List of tuples of "candidate" theta an corresponding index in basic list:
        neg = [(-x[B[i]] / d[B[i]], i) for i in range(m) if d[B[i]] < 0]

        if len(neg) == 0:
            #print("\tidentified unlimited problem")
            return 1, x, set(B),  None, d, it  # Flag problem as unlimited and return ray

        # Get theta and index (in basis) of exiting basic variable:
        theta, p = min(neg, key=(lambda tup: tup[0]))
        #print(p)


        """Variable updates"""
        x = np.array([trunc(var) for var in (x + theta * d)])  # Update all variables
        assert x[B[p]] == 0

        z = trunc(z + theta * r_q)  # Update obj. function value

        # Update inverse:
        for i in set(range(m)) - {p}:
            B_inv[i, :] -= d[B[i]]/d[B[p]] * B_inv[p, :]
        B_inv[p, :] /= -d[B[p]]

        N = N - {q} | {B[p]}  # Update nonbasic index set
        B[p] = q  # Update basic index list

        



def trunc(x):
    """
    Returns 0 if x is smaller (in absolute value) than a certain global constant.
    """
    return x if abs(x) >= epsilon else 0
def simplex_algo():
    A,b,c,flag,cnt=read_input("input.txt")
    return simplex(A,b,c,flag,cnt)

