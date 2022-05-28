import sys
import numpy
import pandas as pd
from ctypes import *


# helper methods

def invalidInput():
    print("Invalid Input!")
    sys.exit()


def errorOccurred():
    print("An Error Has Occurred")
    sys.exit()


def euclidean_norm_squared(vector1, vector2):
    # this function determines the euclidean distance of two vectors, WITHOUT taking a sqrt at the end.
    # CONSIDER: using numpy's built-in?

    dimension = len(vector1)  # supposed to be aligned with dimension outside
    updating_sum = 0

    for coordinate in range(0, dimension):
        updating_sum = updating_sum + numpy.square((vector1[coordinate] - vector2[coordinate]))

    return updating_sum


# input validations

n = len(sys.argv)  # num of cmd line arguments, helps us determine where our different inputs are located

if (n != 5) and (n != 6):  # invalid number of cmd line args
    invalidInput()

k = sys.argv[1]

if not k.isdigit():  # k is not an integer
    invalidInput()

else:
    k = int(k)
    if k <= 0:
        invalidInput()

if n == 6:  # max_iter provided
    max_iter = sys.argv[2]
    epsilon_index = 3
else:  # n == 5, max_iter not provided
    max_iter = 300
    epsilon_index = 2

if (type(max_iter) is str) and (not max_iter.isdigit()):  # max_iter is not an integer
    invalidInput()

else:
    max_iter = int(max_iter)
    if max_iter <= 0:
        invalidInput()

EPSILON = sys.argv[epsilon_index]

try:
    EPSILON = float(EPSILON)
    if EPSILON < 0:
        invalidInput()

except:
    invalidInput()

# no validations required for input files
f1 = sys.argv[epsilon_index + 1]
f2 = sys.argv[epsilon_index + 2]

# after receiving and validating cmd line args, preforming inner join using pandas:

try:
    f1_df = pd.read_csv(f1, sep=',', engine='python', header=None)
    f1_df.columns = [("col " + str(i)) for i in
                     range(len(f1_df.columns))]  # naming the columns, so we have direct access to them by name
    f2_df = pd.read_csv(f2, sep=',', engine='python', header=None)
    f2_df.columns = [("col " + str(i)) for i in
                     range(len(f2_df.columns))]  # naming the columns, so we have direct access to them by name

    # convert file_1 && file_2 to dataframe in pandas, make first column their index column and sort by index
    vector_table = pd.merge(f1_df, f2_df, how='inner', on='col 0')  # merge by first column
    vector_table = vector_table.set_index('col 0')  # setting the first column as index
    vector_table = vector_table.sort_index()

    # implementation algorithm of kmeans_++

    # converting dataframe to numpy 2D-array, in order to use seed and numpy.random.choice() for random choosing
    vector_table = vector_table.to_numpy()
    list_of_all_vectors = vector_table.tolist()  # creating a list of lists to send to c module

    if len(list_of_all_vectors) != 0:
        dimension = len(list_of_all_vectors[0])  # dimension of first vector in the list (assuming input is valid, so all vectors have the same dimension)
    else:  # there are no vectors in the list!
        errorOccurred()

    num_of_vectors = len(vector_table)  # counting the number of vectors provided

    index_list = []  # index list: will be returned to the user at the end of the run

    numpy.random.seed(0)
    curr_miu_index = numpy.random.randint(0, num_of_vectors)  # CONSIDER: instead of choosing a number with randint,
    # maybe choosing directly from the index column in table? note: choosing a random vector does not follow the desired
    # result! (e.g, 26 instead of 44 for input 1)

    index_list.append(curr_miu_index)

    miu_list = [vector_table[curr_miu_index]]  # miu_list now holds the first miu that has been chosen randomly
    curr_miu = miu_list[0]
    prob_list = [0 for i in range(0, num_of_vectors)]  # CONSIDER: change the initial value? for safety reasons
    d_list = [numpy.inf for i in range(0, num_of_vectors)]  # for x_l, d_list[l] is D_l

    # len(miu_list) = i, and it is 1 right now

    while len(miu_list) < k:
        for i in range(0, num_of_vectors):  # for every vector in the list:
            curr_vector = vector_table[i]

            for miu in miu_list:  # for every miu in miu list, we check the min. euclidean distance from curr_vector,
                # and update d_list accordingly
                curr_norm = euclidean_norm_squared(curr_vector, miu)
                if curr_norm < d_list[i]:
                    d_list[i] = curr_norm

        d_sum = sum(d_list)  # sum of all D's (sigma (i from 0 to num_of_vectors) : D_i)

        for i in range(0, num_of_vectors):
            prob_list[i] = d_list[i] / d_sum  # for x_l :  D_l/(sigma (i from 0 to num_of_vectors) : D_i)

        curr_miu_index = numpy.random.choice(vector_table.shape[0], size=1, replace=False, p=prob_list)[
            0]  # choosing a vector randomly, with probabilities (prob_list) taken to account
        index_list.append(curr_miu_index)
        # index_list.append(full_index_list[curr_miu_index])
        curr_miu = vector_table[curr_miu_index]
        miu_list = numpy.vstack([miu_list, [curr_miu]])

    k_first_vectors = miu_list.tolist()

    # CALLING CLUSTERING METHOD FROM HW1 (CONSIDER: need to fix memory freeing)
    c_module = CDLL("./mykmeanssp.so")
    c_module.connect()
    cent_list = c_module.fit(k_first_vectors, list_of_all_vectors, k, num_of_vectors, max_iter, dimension, EPSILON)

    # retrieving original indices of chosen k vectors

    for i in range(0, len(index_list)):  # printing indices (with commas)
        if i == (len(index_list) - 1):
            print(index_list[i])
        else:
            print(index_list[i], ",", sep='', end='')

    for vector in cent_list:  # printing centroids (by coordinates with commas)
        for coordinate in vector:
            if (coordinate == len(vector) - 1):
                print(vector[coordinate])
            else:
                print(vector[coordinate], ",", sep='', end='')

except:
    errorOccurred()
