#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <assert.h>
#include <math.h>

#define A1 1.0
#define B1 3.0
#define A2 -1.5
#define B2 1.5

struct parameters {
    int M, N;
    double h1, h2;
};

int is_inside_area(double x, double y) {
    return (x * x - 4.0 * y * y > 1.0) && (x > 1.0) && (x < 3.0);
}

double get_xi_by_index(struct parameters *params, int i) {
    double res = A1 + i*(params->h1);
    return res;
}

double get_yj_by_index(struct parameters *params, int j) {
    double res = A2 + j*(params->h2);
    return res;
}

double **initialize_to_zeros(struct parameters *params) {
    int i,j;
    double **w;
    int M = params->M, N = params->N;
    w = malloc(sizeof(double*)*(M+1));
    for(i=0; i<=M; i++) {
        w[i] = malloc(sizeof(double)*(N+1));
        for(j=0; j<=N; j++) {
            w[i][j] = 0;
        }
    }
    return w;
}

void matrix_linear_combination(struct parameters *params, double alpha, double betta, double **M1, double **M2, double **res) {
    int i,j;
    int M = params->M, N = params->N;
    for(i=0; i<=M; i++){
        for(j=0; j<=N; j++){
            res[i][j] = alpha*M1[i][j] + betta*M2[i][j];
        }
    }
}

double **initialize_right_part(struct parameters *params) {
    int i,j;
    int M = params->M, N = params->N;
    double xi, yj;
    double **res;
    res = malloc((sizeof(double*))*(M+1));
    for(i=0; i<=M; i++){
        res[i] = malloc(sizeof(double)*(M+1));
        for(j=0; j<=N; j++){
            xi = get_xi_by_index(params, i);
            yj = get_yj_by_index(params, j);
            if(is_inside_area(xi,yj)){
                res[i][j] = 1;
            } else {
                res[i][j] = 0;
            }
        }
    }
    return res;
}

double scalar_product(struct parameters *params, double **u, double **v) {
    double h1 = params->h1;
    double h2 = params->h2;
    double res = 0.0;
    int i,j;
    for(i=1; i<(params->M); i++){
        for(j=1; j<(params->N); j++){
            res += h1*h2*u[i][j]*v[i][j];
        }
    }
    return res;
}

double euclidean_norm(struct parameters *params, double **u) {
    double res = sqrt(scalar_product(params, u, u));
    return res;
}

double left_derivative_x(struct parameters *params, double **w, int i, int j) {
    double h1 = params->h1;
    return (w[i][j] - w[i-1][j])/h1;
}

double right_derivative_x(struct parameters *params, double **w, int i, int j) {
    double h1 = params->h1;
    return (w[i+1][j] - w[i][j])/h1;
}

double left_derivative_y(struct parameters *params, double **w, int i, int j) {
    double h2 = params->h2;
    return (w[i][j] - w[i][j-1])/h2;
}

double right_derivative_y(struct parameters *params, double **w, int i, int j) {
    double h2 = params->h2;
    return (w[i][j+1] - w[i][j])/h2;
}

double express_x_through_y(double y) {
    return sqrt(1+4*y*y);
}

double express_y_through_x(double x) {
    return sqrt((x*x-1)/4);
}

double max(double a, double b) {
    if(a > b) {
        return a;
    }
    return b;
}

double min(double a, double b) {
    if(a < b) {
        return a;
    }
    return b;
}

double calculate_lij_y(double xi, double y1, double y2) {
    double first_intersection = express_y_through_x(xi);
    double second_intersection = -first_intersection;
    double y_max, y_min;
    if(xi<=A1 || xi >=B1) {
        return 0.0;
    }
    if((y1 > first_intersection) || (y2 < second_intersection)) {
        return 0.0;
    }
    y_max = min(first_intersection, y2);
    y_min = max(second_intersection, y1);
    return y_max-y_min;
}

double calculate_lij_x(double yi, double x1, double x2) {
    double intersection = express_x_through_y(yi);
    double x_max = B1;
    if(intersection >= x_max) {
        return 0.0;
    }
    if(x2 < intersection) {
        return 0.0;
    }
    if(x1 > intersection) {
        return x2-x1;
    } else {
        return x2-intersection;
    }
}

double calculate_aij(struct parameters *params, int i, int j) {
    double h1 = params->h1, h2 = params->h2;
    double h = max(h1,h2);
    double eps = h*h;
    double xi = get_xi_by_index(params, i)-0.5*h1;
    double yj = get_yj_by_index(params, j);
    double y1 = yj - 0.5*h2, y2 = yj + 0.5*h2;
    double lij = calculate_lij_y(xi, y1, y2);
    return (1/h2)*lij+(1/eps)*(1-(1/h2)*lij);
}

double calculate_bij(struct parameters *params, int i, int j) {
    double h1 = params->h1, h2 = params->h2;
    double h = max(h1,h2);
    double eps = h*h;
    double yj = get_yj_by_index(params, j)-0.5*h2;
    double xi = get_xi_by_index(params, i);
    double x1 = xi - 0.5*h1, x2 = xi + 0.5*h1;
    double lij = calculate_lij_x(yj, x1, x2);
    return (1/h1)*lij+(1/eps)*(1-(1/h1)*lij);
}

double calculate_Awij(struct parameters *params, double **w, int i, int j) {
    double aij = calculate_aij(params, i, j);
    double bij = calculate_bij(params, i, j);
    double aij_1 = calculate_aij(params, i+1, j);
    double bij_1 = calculate_bij(params, i, j+1);
    double h1 = params->h1, h2 = params->h2;
    double Awij = (-1/h1)*(aij_1*right_derivative_x(params,w,i,j) - aij*left_derivative_x(params,w,i,j))+
                 (-1/h2)*(bij_1*right_derivative_y(params,w,i,j) - bij*left_derivative_y(params,w,i,j));
    return Awij;
}


void operator_A(struct parameters *params, double **w, double **res) {
    int i,j;
    int M = params->M, N = params->N;
    for(i=1; i<=M-1; i++) {
        for(j=1; j<=N-1; j++) {
            res[i][j] = calculate_Awij(params, w, i, j);
        }
    }
}

double **solution_method(struct parameters *params) {
    double **r, **tmp, **Ar, **w_next, **w, **B, **Aw;
    double Arr, Ar_norm, tau_next, norm;
    r = initialize_to_zeros(params);
    B = initialize_right_part(params);
    w = initialize_to_zeros(params);
    Aw = initialize_to_zeros(params);
    Ar = initialize_to_zeros(params);
    w_next = initialize_to_zeros(params);
    size_t count=0;
    double delta = 0.1;
    do{
        operator_A(params, w, Aw);
        matrix_linear_combination(params, 1.0, -1.0, Aw, B, r);
        operator_A(params,r, Ar);
        Arr = scalar_product(params, Ar, r);
        Ar_norm = euclidean_norm(params, Ar);
        tau_next = Arr/(Ar_norm*Ar_norm);
        matrix_linear_combination(params, 1, -tau_next, w, r, w_next);
        tmp = w;
        w = w_next;
        w_next = tmp;
        norm = euclidean_norm(params, r);
       if(count%1000 == 0) {
        printf("%lf\n",norm);
       }
       count++;

    }
    while(norm>delta);
    return w;
}

void get_params(struct parameters *params, int argc, char **argv) {
    int M = atoi(argv[1]);
    int N = atoi(argv[2]);
    params->M = M;
    params->N = N;
    params->h1 = (B1-A1)/(double)M;
    params->h2 = (B2-A2)/(double)N;
}

void print_result(struct parameters *params, double **res)
{
    int i,j;
    int M = params->M, N=params->N;
    double xi,yj;
    for (i = 0; i <= M; ++i) {
        for (j = 0; j <= N; ++j) {
            xi=get_xi_by_index(params, i);
            yj=get_yj_by_index(params,j);
            fprintf(stderr, "%lf\t%lf\t%lf\n", xi, yj, res[i][j]);
        }
    }
}

int main(int argc, char **argv) {
    if(argc < 3){
        fprintf(stderr, "you need to enter M and N");
        exit(1);
    }
    struct parameters *params = malloc(sizeof(*params));
    get_params(params, argc, argv);
    double **res = solution_method(params);
    print_result(params, res);
}



