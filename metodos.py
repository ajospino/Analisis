import numpy as np
import numpy.matlib
from sympy import sympify, Abs, Symbol
from sympy.core.evalf import evalf


def fixedPoint(fx, gx, x0, tol, numMax):
    results = list()
    results.append(["iter", '{:^15.7}'.format("xi"), '{:^15.7}'.format("g(xi)"), '{:^15.7}'.format("f(xi)"), '{:^15.7}'.format("E")])
    
    x = Symbol('x')
    i = 1
    cond = tol
    error = 1.000
    
    ex = sympify(fx)
    rx = sympify(gx)

    xP = x0
    xA = 0.0

    ea = ex.subs(x, xP)
    ea = ea.evalf()

    ra = rx.subs(x, xP)
    ra = ra.evalf()

    while((error > cond) and (i < numMax)):
    
        ra = rx.subs(x,xP)
        xA = ra.evalf()

        ea = ex.subs(x, xP)
        ea = ea.evalf()

        error = Abs(xA - xP)

        xP = xA

        results.append([i, '{:^15.7f}'.format(float(xA)), '{:^15.7f}'.format(float(ra)), '{:^15.7E}'.format(float(ea)), '{:^15.7E}'.format(float(error))])
        
        i += 1
        
    
    for i in results:
        print(i)

def bisec(a, b, fx, Error,numMax):
    results = list()
    results.append(["iter", '{:^15.7}'.format("a"), '{:^15.7}'.format("xm"), '{:^15.7}'.format("b"), '{:^15.7}'.format("f(xm)"), '{:^15.7}'.format("E")])
        
    x = Symbol('x')
    i = 1
    cond = Error
    error = 1.0000000

    ex = sympify(fx)
    
    ea = ex.subs(x, a)
    ea = ea.evalf()

    xm0 = 0.0
    ex_3 = 0

    xm = (a + b)/2

    ex_3 = ex.subs(x, xm)
    ex_3 = ex_3.evalf()

    while(error > cond) and (i < numMax):
        
        if (ea*ex_3 < 0):
            b = xm
        else:
            a = xm

        xm0 = xm
        xm = (a+b)/2

        ex_3 = ex.subs(x, xm)
        ex_3 = ex_3.evalf()

        error = Abs(xm-xm0)

        results.append([i, '{:^15.7f}'.format(float(a)), '{:^15.7E}'.format(float(xm)), '{:^15.7E}'.format(float(b)), '{:^15.7E}'.format(float(ex_3)), '{:^15.7E}'.format(float(error))])

                        
        i += 1

    for i in results:
        print(i)

def regulaFalsi(a, b, fx, Error,numMax):
    results = list()
    results.append(["iter", '{:^15.7}'.format("a"), '{:^15.7}'.format("xm"), '{:^15.7}'.format("b"), '{:^15.7}'.format("f(xm)"), '{:^15.7}'.format("E")])
    

    x = Symbol('x')
    i = 1
    cond = Error
    error = 1.0000000

    ex = sympify(fx)

    xm = 0
    xm0 = 0
    ex_2 = 0
    ex_3 = 0
    ex_a = 0
    ex_b = 0

    while (error > cond) and (i < numMax):
        if i == 1:
            ex_2 = ex.subs(x, a)
            ex_2 = ex_2.evalf()
            ex_a = ex_2

            ex_2 = ex.subs(x, b)
            ex_2 = ex_2.evalf()
            ex_b = ex_2

            xm = (ex_b*a - ex_a*b)/(ex_b-ex_a)
            ex_3 = ex.subs(x, xm)
            ex_3 = ex_3.evalf()
            results.append([i, '{:^15.7f}'.format(a), '{:^15.7f}'.format(xm), '{:^15.7f}'.format(b), '{:^15.7E}'.format(ex_3)])
        else:

            if (ex_a*ex_3 < 0):
                b = xm
            else:
                a = xm

            xm0 = xm
            ex_2 = ex.subs(x, a)
            ex_2 = ex_2.evalf()
            ex_a = ex_2

            ex_2 = ex.subs(x, b)
            ex_2 = ex_2.evalf()
            ex_b = ex_2

            xm = (ex_b*a - ex_a*b)/(ex_b-ex_a)

            ex_3 = ex.subs(x, xm)
            ex_3 = ex_3.evalf()

            error = Abs(xm-xm0)
            er = sympify(error)
            error = er.evalf()
            results.append([i, '{:^15.7f}'.format(a), '{:^15.7f}'.format(xm), '{:^15.7f}'.format(b), '{:^15.7E}'.format(ex_3), '{:^15.7E}'.format(error)])
        i += 1
    
    for i in results:
        print(i)

def secan(x0, x1, fx, tol, numMax):
    results = list()
    results.append(["iter", '{:^15.7}'.format("xi"), '{:^15.7}'.format("f(xi)"), '{:^15.7}'.format("E")])
    
    x = Symbol('x')
    i = 0
    cond = tol
    error = 1.0000000

    ex = sympify(fx)

    y = x0
    ex_0 = ex
    ex_1 = ex

    while((error > cond) and (i < numMax)):
        if i == 0:
            ex_0 = ex.subs(x, x0)
            ex_0 = ex_0.evalf()
            results.append([i, '{:^15.7f}'.format(float(x0)), '{:^15.7E}'.format(float(ex_0))])
        elif i == 1:
            ex_1 = ex.subs(x, x1)
            ex_1 = ex_1.evalf()
            results.append([i, '{:^15.7f}'.format(float(x1)), '{:^15.7E}'.format(float(ex_1))])
        else:
            y = x1
            x1 = x1 - (ex_1*(x1 - x0)/(ex_1 - ex_0))
            x0 = y

            ex_0 = ex.subs(x, x0)
            ex_0 = ex_1.evalf()

            ex_1 = ex.subs(x, x1)
            ex_1 = ex_1.evalf()

            error = Abs(x1 - x0)
            
            results.append([i, '{:^15.7f}'.format(float(x1)), '{:^15.7E}'.format(float(ex_1)), '{:^15.7E}'.format(float(error))])
        i += 1

    for i in results:
        print(i)

def gaussSimple(Ma, b):
    # Getting matrix dimention
    n = len(Ma)

    # Adding the the vector B at the end of the matrix
    matrixMa = np.matrix(Ma)
    vectorB = np.array(b)
    M = np.column_stack((matrixMa, vectorB))

    steps = {'Step 0': np.copy(M)}

    # Matrix reduction
    for i in range(n-1):
        for j in range(i+1, n):
            if (M[j, i] != 0):
                M[j, i:n+1] = M[j, i:n+1]-(M[j, i]/M[i, i])*M[i, i:n+1]
        steps[f'Step {i+1}'] = np.copy(M)

    x = backSubst(M)

    for i, j in steps.items():
        print(i)
        print(j)
        print("\n")
    print("X")
    print(x)

def vanderMon(Vx,Vy):
    output = list()
    X = np.array(Vx)
    s1 = X.size

    Y = np.array(Vy)

    A = np.zeros((s1,s1))

    i = 0
    for i in range(i,s1):
        A[:,i]=X**(s1-i-1)
    
    Coef = (np.linalg.solve(A,Y)).conj().transpose()
    
    output.append("A")
    output.append(A)
    output.append("Coef")
    output.append(Coef)
    
    for i in output:
        print(i)
        print("\n")
    
def difdivid(Vx,Vy):
    output=list()

    X = np.array(Vx)
    n = X.size

    Y = np.array(Vy)

    D = np.zeros((n,n))

    D[:,0]=Y.T
    for i in range(1,n):
        aux0 = D[i-1:n,i-1]
        aux = np.diff(aux0)
        aux2 = X[i:n] - X[0:n-i]
        D[i:n,i] = aux/aux2.T  

    Coef = np.diag(D)
    
    output.append("D")
    output.append(D)
    output.append("Coef")
    output.append(Coef)
    
    for i in output:
        print(i)
        print("\n")

def LUSimple(Ma, b):
    matrixMa = np.array(Ma).T
    n = matrixMa.shape[0]
    L = np.eye(n)
    U = np.zeros((n,n))
    M = matrixMa
    
    steps = {'Step 0': [np.copy(M)]}
    
    for i in range(n-1):
        for j in range(i+1, n):
            if not (M[j,i] == 0):
                L[j,i]=M[j,i]/M[i,i]
                M[j,i:n]=M[j,i:n]-(M[j,i]/M[i,i])*M[i,i:n]
        U[i, i:n]=M[i,i:n]
        U[i+1,i+1:n]=M[i+1,i+1:n]
        steps[f"Step {i+1}"] = np.copy(M)
        steps[f"Step {i+1}.1"] = {"L:":np.copy(L)}
        steps[f"Step {i+1}.2"] = {"U:":np.copy(U)}
    U[n-1,n-1]=M[n-1,n-1]
    
    z=forSubst(np.column_stack((L,b)))
    x=backSubst(np.column_stack((U,z)))

    for i, j in steps.items():
        print(i)
        print(j)
        print("\n")
    print("X")
    print(x)

def backSubst(M):
    n = M.shape[0]
    x = np.matlib.zeros((n, 1), dtype=complex)
    x[n-1] = M[n-1, n]/M[n-1, n-1]
    for i in range(n-2, -1, -1):
        aux1 = np.hstack((1, np.asarray(x[i+1:n]).reshape(-1)))
        aux2 = np.hstack((M[i, n], np.asarray(-M[i, i+1:n]).reshape(-1)))
        x[i] = np.dot(aux1, aux2)/M[i, i]
    return x

def forSubst(M):
    n = M.shape[0]
    x = np.matlib.zeros((n, 1), dtype=complex)
    x[0] = M[0, n]/M[0, 0]
    for i in range(1, n, 1):
        aux1 = np.hstack((1, np.asarray(x[0:i]).reshape(-1)))
        aux2 = np.hstack((M[i, n], np.asarray(-M[i, 0:i]).reshape(-1)))
        x[i] = np.dot(aux1, aux2)/M[i, i]
    return x

def steffenSen(x0,f,tol,nmax):
    x = Symbol('x')
    
    ex = sympify(f)
    ex1 = ex.subs(x,x0)
    ex1 = ex1.evalf()
    
    x1 = ex1
    
    ex2 = ex.subs(x,x1)
    ex2 = ex2.evalf()
    
    x2 = ex2 
    
    xn = x0 - pow((x1-x0),2)/(x2-2*x1+x0)
    error = abs(x2-xn)
    i = 1
    
    while(error>=tol and i < nmax):
        x0 = xn
        
        ex = sympify(f)
        ex1 = ex.subs(x,x0)
        ex1 = ex1.evalf()
        
        x1 = ex1
        
        ex2 = ex.subs(x,x1)
        ex2 = ex2.evalf()
        
        x2 = ex2 
        xn = x0 - pow((x1-x0),2)/(x2-2*x1+x0)
        error = abs(x0-xn)
        
        i += 1
        
        print(f"El error en la iteración {i} : " + str(error))
        print(f"El valor de x en la iteración {i} : " + str(xn))








#fixedPoint( "ln(sin(x)**2 + 1) - (0.5) -x ", "ln(sin(x)**2 + 1) -(0.5)", -0.5, 0.0000001, 100)
#bisec(0,1,"ln(sin(x)**2 + 1) -(0.5)",0.0000001,100)
#regulaFalsi(0, 1, "ln(sin(x)**2 + 1) -(0.5)", 0.0000001, 100)
#secan(0.5, 1, "ln(sin(x)**2 + 1) -(0.5)", 0.0000001, 100)
A = [
     [2, -1, 0, 3],
     [1, 0.5, 3, 8],
     [0, 13, -2, 11],
     [14, 5, -2, 3]
     ]
b = [1,1,1,1]
#gaussSimple(A,b)
#vanderMon([-1, 0, 3, 4],[15.5, 3, 8, 1])
#difdivid([-1, 0, 3, 4],[15.5, 3, 8, 1])
"""LUSimple([
    [4,-1,0,3],
    [1,15.5,3,8],
    [0,-1.3,-4,1.1],
    [14,5,-2,30]
], [1,1,1,1])
"""
#steffenSen(1.5,'(10/(x+4))^(1/2)',0.00001,10)