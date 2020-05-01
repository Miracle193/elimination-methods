import time

# p - pivot row
# i - column index
# n - number of unknowns
# a - augmented matrix
# finds the pivot which is the smallest integer
# such that i <= p <= n
def findPivotBackward(p, i, n, a):
    if p == n + 1:
        return -1
    if a[p][i] == 0:
        p = p + 1
        return findPivotBackward(p, i, n, a)
    return p

# n - number of unknowns
# matrix - augmented matrix
def gaussianBackwardSubstitution(n, matrix):
    for i in range(0, n): # Elimination process
        p = findPivotBackward(i, i, n, matrix) # find the pivot element
        if p == -1: # if matrix[p][i] == 0 there is no unique solution
            print("No unique solution exits")
            return
        if p != i:
            matrix[p], matrix[i] = matrix[i], matrix[p]

        # simulated row exchange
        for j in range(i+1, n+1):
            m = matrix[j][i] / matrix[i][i]
            multiple = [(elem * m) for elem in matrix[i]]
            matrix[j] = [(elem - elemMultiple) for elem, 
            elemMultiple in zip(matrix[j], multiple)]
    if matrix[n-1][n-1] == 0:
        print("No unique solution exits")
        return

    x = [0] * (n + 1)
    # start backward substitution
    x[n] = matrix[n][n+1]/matrix[n][n] 
    for i in range(n-1, -1, -1):
        for j in range(i+1, n+1):
            x[i] = x[i] + (matrix[i][j] * x[j])
        x[i] = (matrix[i][n+1] - x[i]) / matrix[i][i]
    return x

# i - column index
# n - number of unknowns
# a - augmented matrix
# nrow - row pointer
# finds the pivot which is the smallest integer
# such that i <= p <= n and |a[nrow[p]][i]| = max(|a[nrow[j]][i]|)
# where i <= j <= n
def findPivotPartial(i, n, a, nrow):
    maxnum = 0
    p = 0
    for j in range(i, n+1):
        if abs(a[nrow[j]][i]) > maxnum:
            maxnum = abs(a[nrow[j]][i])
            p = j
    return p

# n - number of unknowns
# matrix - augmented matrix
def gaussianPartialPivoting(n, matrix):
    NROW = [None] * (n + 1) 
    # Initialize row pointer
    for i in range(0, n+1):
        NROW[i] = i

    for i in range(0, n): # Elimination process
        p = findPivotPartial(i, n, matrix, NROW) # find pivot element
        if matrix[NROW[p]][i] == 0:
            print("No unique solution exits")
            return
        if NROW[i] != NROW[p]:
            NCOPY = NROW[i]
            NROW[i] = NROW[p]
            NROW[p] = NCOPY

        # simulated row exchange
        for j in range(i+1, n+1): 
            m = matrix[NROW[j]][i]/matrix[NROW[i]][i]
            multiple = [(m * elem) for elem in matrix[NROW[i]]]
            matrix[NROW[j]] = [(elem - elemMultiple) for 
            elem, elemMultiple in zip(matrix[NROW[j]], multiple)]

    if matrix[NROW[n]][n] == 0:
        print("No unique solution exits")
        return

    x = [0] * (n+1) 
    # start backward substitution
    x[n] = matrix[NROW[n]][n+1]/matrix[NROW[n]][n]
    for i in range(n-1, -1, -1):
        for j in range(i+1, n+1):
            x[i] = x[i] + (matrix[NROW[i]][j] * x[j])
        x[i] = (matrix[NROW[i]][n+1] - x[i])/matrix[NROW[i]][i]
    return x

# i - column index
# n - number of unknowns
# a - augmented matrix
# nrow - row pointer
# s - list of scale factor for each row
# finds the pivot which is the smallest integer
# such that i <= p <= n and 
# |a[nrow[p]][i]|/s[nrow[p]] = max(|a[nrow[j]][i]|/s[nrow[j]])
# where i <= j <= n
def findPivotScaled(i, n, a, nrow, s):
    p = 0
    maxnum = 0
    for j in range(i, n+1):
        div = abs(a[nrow[j]][i])/s[j]
        if div > maxnum:
            max = div
            p = j
    return p

# n - number of unknowns
# matrix - augmented matrix
def gaussianScaledPivoting(n, matrix):
    s = [None] * (n + 1)
    NROW = [None] * (n + 1)
    # Initialize scale factor and row pointer
    for i in range(0, n+1): # set up scale f
        absolute = [abs(elem) for elem in matrix[i]]
        s[i] = max(absolute)
        if s[i] == 0:
            print("No unique solution exits")
            return
        else:
            NROW[i] = i

    for i in range(0, n): # Elimination process
        p = findPivotScaled(i, n, matrix, NROW, s) # find pivot element
        if matrix[NROW[p]][i] == 0:
            print("No unique solution exits")
            return
        if NROW[i] != NROW[p]:
            NCOPY = NROW[i]
            NROW[i] = NROW[p]
            NROW[p] = NCOPY

        # simulated row exchange
        for j in range(i+1, n+1):
            m = matrix[NROW[j]][i]/matrix[NROW[i]][i]
            multiple = [(m * elem) for elem in matrix[NROW[i]]]
            matrix[NROW[j]] = [(elem - elemMultiple) for 
            elem, elemMultiple in zip(matrix[NROW[j]], multiple)]
    if matrix[NROW[n]][n] == 0:
        print("No unique solution exits")
        return

    x = [0] * (n+1)
    # start backward substitution
    x[n] = matrix[NROW[n]][n+1]/matrix[NROW[n]][n] 
    x[n] = round(x[n], 0)
    for i in range(n-1, -1, -1):
        for j in range(i+1, n+1):
            x[i] = x[i] + (matrix[NROW[i]][j] * x[j])
        x[i] = (matrix[NROW[i]][n+1] - x[i])/matrix[NROW[i]][i]
        x[i] = round(x[i], 0)
    return x

# n - number of unknowns
# matrix - augmented matrix
def gaussJordan(n, matrix):
    for i in range(0, n+1): # Elimination process
        p = findPivotBackward(i, i, n, matrix) # find pivot element
        if p == -1:
            print("No unique solution exits")
            return -1
        if p != i:
            matrix[p], matrix[i] = matrix[i], matrix[p]
        loop = list(range(0, n+1))
        loop.remove(i)

        # simulated row exchange
        for j in loop:
            m = matrix[j][i]/matrix[i][i]
            multiple = [(elem * m) for elem in matrix[i]]
            matrix[j] = [(elem - elemMultiple) for 
            elem, elemMultiple in zip(matrix[j], multiple)]
    if matrix[n][n] == 0:
        print("No unique solution exits")
        return -1

    x = [0] * (n+1)
    # start backward substitution
    for i in range(n, -1, -1):
        x[i] = matrix[i][n+1]/matrix[i][i]
    return x

matrix = [[2, -5, 15], 
[3, 1, 31]]

n = 1
start_time = time.time()
print(gaussianBackwardSubstitution(n, matrix))
print(time.time()-start_time)

start_time = time.time()
print(gaussianPartialPivoting(n, matrix))
print(time.time()-start_time)

start_time = time.time()
print(gaussianScaledPivoting(n, matrix))
print(time.time()-start_time)

start_time = time.time()
print(gaussJordan(n, matrix))
print(time.time()-start_time)