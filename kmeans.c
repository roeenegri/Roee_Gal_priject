#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdio.h> 
#include <stdlib.h> /*memo alloc*/
#include <math.h> /*math ops*/
#include <string.h> /*string ops*/
#include <ctype.h>


/* function declaration */

double** extract_k_first_vectors(char* , int, int, int*);
int check_centroid_diff(double**, double**, int, int, double);
void calc_centroids(double**, int *, double **, int, int);
int matching_cluster(double*, double **, int, int);
double calc_norm_form(double*, double *, int, int);
int extract_vectors (char*, double**, int, int);
int terminate_with_error ();
int terminate_invalid_input ();
void sum_two_vectors (double **, double *, int, int);
int find_dimension(char*);
int find_num_of_vectors(char*);
int isNumber(char*);



/* constants */
char INVALID [] = "Invalid Input!"; /*print and return 1 if input is invalid*/
char ERROR [] = "An Error Has Occurred"; /*print and return 1 if any other error happend*/





static PyMethodDef _capiMethods[] = {
    {"fit", (PyCFunction) mainC_Py, METH_VARARGS, PyDoc_STR("calc the k meams")},
    {NULL, NULL, 0, NULL}
};


static struct PyModuleDef _moduledef = {
    PyModuleDef_HEAD_INIT,
    "mykmeanssp",
    NULL,
    -1,
    _capiMethods
};


PyMODINIT_FUNC PyInit_capi_kmeans(void) {
    PyObject *m;
    m = PyModule_Create(&_moduledef);
    if(!m){
        return NULL;
    }
    return m;
}

/* main function */

static PyObject *
mainC_Py(PyObject *self, PyObject *args){
    PyObject *file_namePy;
    PyObject *kFirstVectorsIndexPy;
    PyObject *ret;
    char* file_name;
    int* k_first_vectors_index;
    double **data;
    double **kFirstVectors;
    int k, n, max_iter, dimension;
    float EPSILON;
    int cnt ; /*iteration counter of the while-loop*/
    double **cents; /*array of centroids*/
    double **prev_cents; /*array of previos centroids (for later calculations)*/
    double **data_list;
    int i;
    int val;
    int *num_of_vectors;
    double **sum_of_vectors;
    int j;

    if (!PyArg_ParseTuple(args, "sOiiiid", &file_namePy, &kFirstVectorsIndexPy, &k, &n, &max_iter, &dimension, &EPSILON)){
        //print erorr??????????????
        return NULL;
    }

    k_first_vectors_index = calloc(k, sizeof(int));

    for(i=0; i < k; i++){
        k_first_vectors_index[i] = PyLong_AsLong(PyList_GET_ITEM(kFirstVectorsIndexPy, i));
    }

    file_name = file_namePy; 

    cnt = 0;

    data_list = calloc(n, sizeof(double*));
    for(i = 0; i < n; ++i){
        data_list[i] = calloc(dimension, sizeof(double));
    }

    prev_cents = calloc(k, sizeof(double*));
    for(i = 0; i < k; ++i){
        prev_cents[i] = calloc(dimension, sizeof(double));
    }
    
    val = extract_vectors(file_name, data_list, n, dimension); //put all the data from the file to the datalist
    if(val == -2){
        terminate_with_error();
    }

     /*algorithm*/

    /*initialization of centroid list: with first k vectors*/

    cents = extract_k_first_vectors(file_name, k, dimension, k_first_vectors_index); 

    do {
      
        int current_xi_position = -1;
        int clust_num = -1;
      
        sum_of_vectors = calloc(k, sizeof(double*));
        for(i = 0; i < k; ++i){
            sum_of_vectors[i] = calloc(dimension, sizeof(double));
        } 

        num_of_vectors = calloc(k, sizeof(int));

        /* prev_cents = cents; we save the current cents to compare the with new ones later: make prev point to cents,,,, could be a pointer problem here*/ 

        for (i=0; i<k; i++) {
            for (j=0; j<dimension; j++) {
                prev_cents[i][j] = cents[i][j];
            }
        }


        for (current_xi_position = 0; current_xi_position < n; current_xi_position ++) {
            
            clust_num = matching_cluster(data_list[current_xi_position], cents, k, dimension); 
            sum_two_vectors(sum_of_vectors, data_list[current_xi_position], dimension, clust_num);
            num_of_vectors[clust_num] += 1;
        }

        calc_centroids(sum_of_vectors, num_of_vectors, cents, k, dimension);

        cnt = cnt+1;

    } while ((cnt < max_iter) && (check_centroid_diff(prev_cents, cents, k, dimension, EPSILON) == 0));

    ret = PyList_New(k);
    for (i = 0;i < k; i++){
        PyObject *vector_py = PyList_New(dimension);
        for(j = 0;j < dimension; j++){
            PyList_SetItem(vector_py, j, PyFloat_FromDouble(cents[i][j]));
        }
        PyList_SetItem(ret,i, vector_py);
    }




    return ret;
}


/* implementation of helper functions */


int extract_vectors (char *filename, double **data_list, int n, int dimension) {
        
    /*input: filename, a pointer to the data (vector) list
    output: number of vectors
    function will extract all the vectors from txt file to a list of vectors, and return their number*/

    int i = 0;
    double number = 0;
    int j = 0;

    FILE* f = fopen(filename, "r");
    if (f == NULL) {
        return -2; /*error will be printed outside the function call*/
    
    }

    for(i = 0; i < n; i++){
        while(fscanf(f, "%lf", &number) != EOF && (j < dimension)) {
                if ((number != ',') && (number != '\n')){
                data_list[i][j] = number;
                j++;
                number = fgetc(f);
                }
              
        }
        j = 0;
    }

    if (data_list == NULL) {
        return -2;
    } 

    fclose(f);
    return 0;

     /* i (divided by) dimension =  number of vectors, this is what needs to be returned*/
    /* MISSING: we need to also return the dimension outside somehow */
}


int check_centroid_diff(double **prv_cents, double **cents, int k, int dimension, double epsilon) { 
    
    int i = 0;
    double diff;

    for(i = 0; i < k; i++){
        diff = calc_norm_form(prv_cents[i], cents[i], dimension, 1);

        if(diff > epsilon){
            return 0;
        }
    }  

    return 1;
   
    /*check_centroid_diff (prev_cents, curr_cents):
    input: 2 lists of previous centroids and current centroids
    output: 1 if ALL of them hasn't changed more than epsilon, else 0*/
}


void calc_centroids(double **sum_of_vectors,int *num_of_vectors,double **cents, int k, int dimension){

    int i;
    int j;
    for(i = 0; i < k; ++i){
        for(j = 0; j < dimension; ++j){

            cents[i][j] = sum_of_vectors[i][j] / num_of_vectors[i];

        }
    }

    return;
    
    /*input: two arrays: sum_of_vectors, num_of_vectors
    output: void (just updates the centroid value in cents list)*/
}

int matching_cluster(double *current_xi,double **cents, int k, int dimension) {

    int i;
    double curr;
    double min_diff;
    int min_centro;
    
    i = 0;

    min_diff = calc_norm_form(current_xi, cents[0], dimension, 0);
    min_centro = 0;

    for(i = 1; i < k; i++){
        curr = calc_norm_form(current_xi, cents[i], dimension, 0);
        if(curr < min_diff){
            min_diff = curr;
            min_centro = i;
        }
    }

    return min_centro;

    /*input: current vector xi, and the centroid list (data_list)
    output: (int) the number of the new cluster to assign the vector to: the closest centroid to the vector*/
}


double calc_norm_form(double *vec1, double *vec2, int dimension, int with_sqrt) {

    /* calculate the norm of the difference vector */

    double sum;
    int i;
    double num1;
    double num2;

    sum = 0;
    i = 0;

    for(i = 0; i < dimension; i++){
        num1 = vec1[i];
        num2 = vec2[i];
        sum += pow((num1-num2), 2.0);

    }
    
    if (with_sqrt == 1) {
        
        sum = pow(sum, 0.5);

    }

    return sum;
}



int terminate_with_error () {
        
    /* terminates the program: prints the general error message and returns 1 */

    printf("%s", ERROR);
    return 1;
}

int terminate_invalid_input () {
    
    /* terminates the program: prints the invalid input message and returns 1 */

    printf("%s", INVALID);
    return 1;
}

void sum_two_vectors (double **sum_of_vectors, double *new_vector, int dimension, int cluster_number) {

    int i;
    i = 0;
    /*input: a list of vectors, a vector, the dimension and the index of the list to the specific sum-vector we want to add the vector to.
    output: void
    the function takes a vector in the list that represents a sum, and addes a new vector to the sum (the change is done in place)*/

    for (i = 0; i < dimension; i++) {

        sum_of_vectors[cluster_number][i] += new_vector[i];
    }

}

double** extract_k_first_vectors(char* file_name, int k, int dimension, int* index_list) {
    /*extract_first_k_vectors(data_list *, arr *)):
    input: data_list *
    output: will update the double-array arr with first k vectors of the input file*/
    int i;
    int j = 0;
    double **centroids;
    double number;
    int index;
    int counter = 0;
    centroids = calloc(k, sizeof(double *));
    for(i = 0; i < k; ++i){
        centroids[i] = calloc(dimension, sizeof(double));
    }

    FILE* f = fopen(file_name, "r");
    for(i = 0; i < k; i++){
        index = index_list[i];
        while(fscanf(f, "%lf", &number) != EOF && (j < dimension)) {
            if(index == counter){
                if ((number != ',') && (number != '\n')){
                centroids[i][j] = number;
                j++;
                number = fgetc(f);
                }
            }
            else{
                while(number != '\n'){
                    number = fgetc(f);
                }  
            }
            j = 0;
            counter++;
            }
        }
        return centroids;
    }    
   
    








