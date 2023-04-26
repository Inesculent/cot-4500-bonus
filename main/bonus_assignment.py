#This is the bonus assignment for COT-4500
#The majority of the code has been recycled from previous assignments and modified slightly in order to fit the reqs of this 
#assignment


import numpy as np


def make_diagonally_dominant(matrix, b_vector):
    n = len(matrix)

    for i in range(n):
        pivot: float = matrix[i][i]
        sum_of_other_elements = sum(abs(matrix[i][i+1:]))

        # we can guarantee this pivot is the largest in the row
        if abs(pivot) > abs(sum_of_other_elements):
            continue

        # if we reach this point, this means we need to swap AT LEAST ONCE
        max_value_of_row = 0
        max_index_in_row = 0
        for j in range(n):
            current_value_in_row: float = abs(matrix[i][j])
            if current_value_in_row > max_value_of_row:
                max_value_of_row = current_value_in_row
                max_index_in_row = j

        # now that we have a new "pivot", we swap cur_row with the expected index
        matrix[[i, max_index_in_row]] = matrix[[max_index_in_row, i]]
        b_vector[[i, max_index_in_row]] = b_vector[[max_index_in_row, i]]
        
    return matrix, b_vector

def jacobi(matrix, vector):

    prevx1 = 0
    prevx2 = 0
    prevx3 = 0

    x1 = (vector[0] - (prevx2 * matrix[0][1]) - (prevx3 * matrix[0][2]))/matrix[0][0]
    x2 = (vector[1] - (prevx1 * matrix[1][0]) - (prevx3 * matrix[1][2]))/matrix[1][1]
    x3 = (vector[2] - (prevx1 * matrix[2][0]) - (prevx2 * matrix[2][1]))/matrix[2][2]

    tab = (x1 - prevx1)

    prevx1 = x1
    prevx2 = x2
    prevx3 = x3

    #loop until all xs in array
    #loop until finished
    counter = 2
    while (abs(tab) > 10**-6 and counter <=50):

        x1 = (vector[0] - (prevx2 * matrix[0][1]) - (prevx3 * matrix[0][2]))/matrix[0][0]
        x2 = (vector[1] - (prevx1 * matrix[1][0]) - (prevx3 * matrix[1][2]))/matrix[1][1]
        x3 = (vector[2] - (prevx1 * matrix[2][0]) - (prevx2 * matrix[2][1]))/matrix[2][2]

        tab = (x1 - prevx1)
        prevx1 = x1
        prevx2 = x2
        prevx3 = x3
        counter += 1

    return counter


def gauss_seidel(matrix, vector):

    prevx1 = 0
    prevx2 = 0
    prevx3 = 0

    x1 = (vector[0] - (prevx2 * matrix[0][1]) - (prevx3 * matrix[0][2]))/matrix[0][0]
    tab = (x1 - prevx1)
    prevx1 = x1
    x2 = (vector[1] - (prevx1 * matrix[1][0]) - (prevx3 * matrix[1][2]))/matrix[1][1]
    prevx2 = x2
    x3 = (vector[2] - (prevx1 * matrix[2][0]) - (prevx2 * matrix[2][1]))/matrix[2][2]
    prevx3 = x3


    #loop until all xs in array
    #loop until finished
    counter = 1
    while (abs(tab) > 10**-6 and counter <=50):

        x1 = (vector[0] - (prevx2 * matrix[0][1]) - (prevx3 * matrix[0][2]))/matrix[0][0]
        tab = (x1 - prevx1)
        prevx1 = x1
        x2 = (vector[1] - (prevx1 * matrix[1][0]) - (prevx3 * matrix[1][2]))/matrix[1][1]
        prevx2 = x2
        x3 = (vector[2] - (prevx1 * matrix[2][0]) - (prevx2 * matrix[2][1]))/matrix[2][2]
        prevx3 = x3
        counter += 1

    return counter

#Is simply the derivative of ("x^3 + 4x^2 -10")
def custom_derivative(value):
    return (3 * value* value) - (2 * value)


#Function to utilize newton_raphson method
def newton_raphson(initial: float, tolerance: float, sequence):
    # Track iterations
    counter = 0
    # Set x to the initial approximation
    x = initial
    #Find the derivative f'
    f = eval(sequence)  

    fPrime = custom_derivative(initial) 
    
    #We get our approximation using the intial divided by the derivative
    approximation: float = f / fPrime
    #As long as it's above the tolerance
    while(abs(approximation) >= tolerance):
        #find function f
        x = initial
        f = eval(sequence)
        #find the derivative
        fPrime = custom_derivative(initial)
        
        #Find the approximation
        approximation = f / fPrime
        
        #Subtract the approximation from initial and continue
        initial -= approximation
        counter += 1

    #Return the counter because that's what we're looking for
    return counter
    

def function1(t: float, w: float):
    return w - (t**3)


def do_work(t, w, h):
    basic_function_call = function1(t, w)

    incremented_t = t + h
    incremented_w = w + (h * basic_function_call)
    incremented_function_call = function1(incremented_t, incremented_w)

    return basic_function_call + incremented_function_call

def modified_eulers():
    original_w = .5
    start_of_t, end_of_t = (0, 3)
    num_of_iterations = 100

    # set up h
    h = (end_of_t - start_of_t) / num_of_iterations

    next_w = 0

    for cur_iteration in range(0, num_of_iterations):
        # do we have all values ready?
        t = start_of_t
        w = original_w
        h = h

        # create a function for the inner work
        inner_math = do_work(t, w, h)

        # this gets the next approximation
        next_w = w + ( (h / 2) * inner_math )


        # we need to set the just solved "w" to be the original w
        # and not only that, we need to change t as well
        start_of_t = t + h
        original_w = next_w
        
    return next_w



def apply_div_dif(matrix: np.array):
    size = len(matrix)
    for i in range(2, size):
        for j in range(2, i+2):
            # skip if value is prefilled or if exceeds bounds
            if j >= len(matrix[i]) or matrix[i][j] != 0:
                continue
            
            #Get left
            left: float = matrix[i][j-1]
            #Get right
            diagonal_left: float = matrix[i-1][j-1]
            #Numerator is left - diagonal left
            numerator: float = left - diagonal_left
            #Denominator is as follows...
            denominator = (matrix[i][0] - matrix[i-j+1][0])
            #Divide numerator by denominator and save into the matrix
            operation = numerator / denominator
            matrix[i][j] = operation
    
    return matrix

def hermite_interpolation():
    #Declare inputs
    x_points = [0,1,2]
    y_points = [1,2,4]
    slopes = [1.06,1.23,1.55]
   
    num_pts = len(x_points)
    matrix = np.zeros((2*num_pts, 2*num_pts))
    #Fill x values
    for x in range(2 * num_pts - 1):
        matrix[x][0] = x_points[int(x/2)]
        matrix[x + 1][0] = x_points[int(x/2)]
        x += 1
    
    #Fill y values
    for x in range(2 * num_pts - 1):
        matrix[x][1] = y_points[int(x/2)]
        matrix[x + 1][1] = y_points[int(x/2)]
        x += 1


    #Fill derivatives
    for x in range(num_pts):
        matrix[(x*2) + 1][2] = slopes[int(x)]

    #Now utilize the divided difference method 
    filled_matrix = apply_div_dif(matrix)

    print(filled_matrix, "\n")



if __name__ == "__main__":
    np.set_printoptions(precision=7, suppress=True, linewidth=100)

    matrix = np.array([[3,1,1],
                       [1,4,1],
                       [2,3,7]])
    b_vector = np.array([1,3,0])

    d_matrix, new_b = make_diagonally_dominant(matrix, b_vector)

    print(gauss_seidel(d_matrix, new_b),"\n")
    
    print(jacobi(d_matrix, new_b),"\n")

    function = "x**3 - x**2 + 2"
    accuracy: float = 10**-6
    
    print(newton_raphson(.5, accuracy, function), "\n")

    hermite_interpolation()


    print("%.5f" % modified_eulers(), "\n")
