#include <stdio.h> 
#include <stdlib.h> /*memo alloc*/
#include <math.h> /*math ops*/
#include <string.h> /*string ops*/
#include <ctype.h>


/* function declaration */

double** extract_k_first_vectors(double**, int, int);
int check_centroid_diff(double**, double**, int, int);
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
int EPSILON = 0.001;


/* main function */

int main (int argc, char **argv) {
    
    /*input varibles declaration*/

    int k;
    int max_iter;
    int input_file_index;
    int output_file_index;
    int n; /*default number of vectors in input file is set to 0, will be set according to number of vectors soon*/
    int cnt ; /*iteration counter of the while-loop*/
    int dimension; /*dimension of the vectors*/
    double **cents; /*array of centroids*/
    double **prev_cents; /*array of previos centroids (for later calculations)*/
    double **data_list;
    int i;
    int val;
    int *num_of_vectors;
    double **sum_of_vectors;
    FILE* output;
    int j;

    cnt = 0;

    if(isNumber(argv[1]) == 0){
        return terminate_invalid_input();
    }
    /* MISSING: casting? check that arguments are from the right types (int, etc.) */
    sscanf(argv[1], "%d", &k);

    if (k <= 1) { /* we assume argv[1] is k. validation of smaller than n will be checked later: row 74*/
        return terminate_invalid_input();
    }

    else if (argc == 4) { /*no max_iter provided*/
        
        max_iter = 200; /*default value*/
        input_file_index = 2;
        output_file_index = 3;
    }

    else if (argc == 5) { /*max_iter provided, we assume it's in argv[2]*/
        if(isNumber(argv[2]) == 0){
            return terminate_invalid_input();
        }
        sscanf(argv[2], "%d", &max_iter);
        if (max_iter <= 0) {
            return terminate_invalid_input();
        }
        input_file_index = 3;
        output_file_index = 4;
    }

    else {

        return terminate_invalid_input();
    }

    /*varible decleration*/


    n = find_num_of_vectors(argv[input_file_index]);

    if(k > n){
        return terminate_invalid_input();
    }


    dimension = find_dimension(argv[input_file_index]);
    
    data_list = calloc(n, sizeof(double*));
    for(i = 0; i < n; ++i){
        data_list[i] = calloc(dimension, sizeof(double));
    }

    prev_cents = calloc(k, sizeof(double*));
    for(i = 0; i < k; ++i){
        prev_cents[i] = calloc(dimension, sizeof(double));
    }

    val = extract_vectors(argv[input_file_index], data_list, n, dimension);
    if (val == -1){
        return terminate_invalid_input();
    }

    else if (val == -2) {
        return terminate_with_error();
    }

    if (k >= n) {
        return terminate_invalid_input();
    }

    /*algorithm*/

    /*initialization of centroid list: with first k vectors*/
    cents = extract_k_first_vectors(data_list, k, dimension); 


    /* checked with input_1 until here ! */

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

    } while ((cnt < max_iter) && (check_centroid_diff(prev_cents, cents, k, dimension) == 0));


    /* output file generation */

    output = fopen (argv[output_file_index], "w");
    if (output == NULL) {
        return terminate_invalid_input();
    }

    for (i = 0; i < k; i ++) {
        for (j = 0; j < dimension; j ++){
            if ((j == 0) && (i != 0)) {
                
                fprintf (output, "\n");
            }

            if (j == (dimension-1)) {
                
                fprintf (output,"%.4f", cents[i][j]); 
            }

            else {
                
                fprintf (output,"%.4f,", cents[i][j]);
            }
        }
    }

    fclose(output);

    /* MISSING : other errors? we haven't used it: try/catch? assert? */

    free (data_list);
    free (prev_cents);
    free (cents);

    return 0;
}


/* implementation of helper functions */

double** extract_k_first_vectors(double **data_list, int k, int dimention) {
    int i;
    int j;
    double **centroids;
    centroids = calloc(k, sizeof(double *));
    for(i = 0; i < k; ++i){
        centroids[i] = calloc(dimention, sizeof(double));
    }
    for(i = 0; i < k; ++i){
        for(j = 0; j < dimention; ++j){
            centroids[i][j] = data_list[i][j];
        }
    }

    return centroids;
    
    /*extract_first_k_vectors(data_list *, arr *)):
    input: data_list *
    output: will update the double-array arr with first k vectors of the input file*/
}

int check_centroid_diff(double **prv_cents, double **cents, int k, int dimension) { 
    
    int i = 0;
    double diff;

    for(i = 0; i < k; i++){
        diff = calc_norm_form(prv_cents[i], cents[i], dimension, 1);

        if(diff > 0.001){
            return 0;
        }
    }  

    return 1;
   
    /*check_centroid_diff (prev_cents, curr_cents):
    input: 2 lists of previous centroids and current centroids
    output: 1 if ALL of them hasn't changed more than epsilon, else 0*/
}

int find_dimension(char* input_file){
    char c;
    int dimension = 0;
    FILE *f = fopen(input_file, "r");
    if(f == NULL){
        return terminate_invalid_input();
    }

    c = fgetc(f);

    while ((c != '\n') && (c != EOF)) {
        if(c == ','){
            dimension += 1;
        }
        c = fgetc(f);
    }

    fclose(f);

    dimension += 1;

    return dimension;
}


int find_num_of_vectors(char* input_file){
            
    FILE *f;
    int count = 0; 
    char c; 
 
    f = fopen(input_file, "r");
  
    if (f == NULL){
        return -1;
    }
  
    c = fgetc(f);

    while (c != EOF) {
        if (c == '\n') {
            count = count + 1;
        }
    
        c = fgetc(f);
    }

    fclose(f);

    return count;
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

int extract_vectors (char *filename, double **data_list, int n, int dimension) {
        
    /*input: filename, a pointer to the data (vector) list
    output: number of vectors
    function will extract all the vectors from txt file to a list of vectors, and return their number*/

    int i = 0;
    double number = 0;
    int j = 0;

    FILE* f = fopen(filename, "r");
    if (f == NULL) {
        return -1; /*invalid input will be printed outside the function call*/
    
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


int terminate_with_error () {
        
    /* terminates the program: prints the general error message and returns 1 */

    printf(ERROR);
    return 1;
}

int terminate_invalid_input () {
    
    /* terminates the program: prints the invalid input message and returns 1 */

    printf(INVALID);
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

    int isNumber(char number[]){
    int i = 0;

    if (number[0] == '-')
        return 0;
    for (; number[i] != 0; i++)
    {
        if (!isdigit(number[i]))
            return 0;
    }
    return 1;
}




