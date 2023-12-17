#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <mpi.h>

#define A1 1.0
#define B1 3.0
#define A2 -1.5
#define B2 1.5

int proc_id, proc_amount;
int N;
int M;
double delta_stop;
const char *file_name;
int inside_domain_N;
int inside_domain_M;
int domain_perim;
int sum_perim;
int vertical_domains_count;
int domain_size;
int grid_size;
double step_x;
double step_y;
double epsilon;

inline int inside_domain_i(int i) {
    return i % inside_domain_M;
}

inline int inside_domain_j(int j) {
    return j % inside_domain_N;
}

inline int domain_id(int i, int j) {
    return (i / inside_domain_M) * vertical_domains_count + j / inside_domain_N;
}

inline int global_i(int i, int domain_id) {
    return (domain_id / vertical_domains_count) * inside_domain_M + i;
}

inline int global_j(int j, int domain_id) {
    return (domain_id % vertical_domains_count) * inside_domain_N + j;
}

inline int index_inside_domain(int i, int j) {
    return i * inside_domain_N + j;
}    

inline int global_index(int i, int j){
    return i * N + j;
}

inline int domain_index(int i, int j) {
    return domain_id(i,j) * domain_size + index_inside_domain(inside_domain_i(i), inside_domain_j(j));
}

inline double right_derivative_x(double *w, double *b, int g_i, int g_j, int i, int j) {
    double tmp;
    if(domain_id(g_i+1, g_j) == proc_id) {
        tmp = w[index_inside_domain(i+1,j)];
    } else {
        tmp = b[domain_id(g_i+1, g_j)*domain_perim + 2*inside_domain_M + j];
    }
    return (tmp - w[index_inside_domain(i,j)])/step_x;
}

inline double left_derivative_x(double *w, double *b, int g_i, int g_j, int i, int j) {
    double tmp;
    if(domain_id(g_i-1, g_j) == proc_id) {
        tmp = w[index_inside_domain(i-1,j)];
    } else {
        tmp = b[domain_id(g_i-1, g_j)*domain_perim + 2*inside_domain_M + inside_domain_N + j];
    }
    return (w[index_inside_domain(i,j)]-tmp)/step_x;
}

inline double right_derivative_y(double *w, double *b, int g_i, int g_j, int i, int j) {
    double tmp;
    if(domain_id(g_i, g_j+1) == proc_id) {
        tmp = w[index_inside_domain(i,j+1)];
    } else {
        tmp = b[domain_id(g_i, g_j+1)*domain_perim + i];
    }
    return (tmp - w[index_inside_domain(i,j)])/step_y;
}

inline double left_derivative_y(double *w, double *b, int g_i, int g_j, int i, int j) {
    double tmp;
    if(domain_id(g_i, g_j-1) == proc_id) {
        tmp = w[index_inside_domain(i,j-1)];
    } else {
        tmp = b[domain_id(g_i, g_j-1)*domain_perim + inside_domain_M + i];
    }
    return (w[index_inside_domain(i,j)]-tmp)/step_y;
}

double get_xi_by_index(int i) {
    double res = A1 + i*step_x;
    return res;
}

double get_yj_by_index(int j) {
    double res = A2 + j*step_y;
    return res;
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

void calculate_domain_grid_sizes(void)
{
    int is_horizontal = 0, proc_count = proc_amount / 2;
    inside_domain_N = N;
    inside_domain_M = M;
    while (proc_count) {
        if (is_horizontal)
            inside_domain_N /= 2;
        else
            inside_domain_M /= 2;
        is_horizontal = !is_horizontal;
        proc_count /= 2;
    }
    domain_perim = 2 * (inside_domain_N + inside_domain_M);
    sum_perim = domain_perim * proc_amount;
    vertical_domains_count = N / inside_domain_N;
    domain_size = inside_domain_M * inside_domain_N;
}

int is_inside_area(double x, double y) {
    return (x * x - 4.0 * y * y > 1.0) && (x > 1.0) && (x < 3.0);
}

double express_y_through_x(double x) {
    return sqrt((x*x-1)/4);
}

double express_x_through_y(double y) {
    return sqrt(1+4*y*y);
}

double convert_right_part(double x, double y)
{
    return is_inside_area(x, y) ? 1.0 : 0.0;
}

double scalar_product(double *u, double *v)
{
    double result = 0.0;
    double koeff = step_x * step_y;
    int i, j;
    for (i = 0; i < inside_domain_M; ++i) {
        for (j = 0; j < inside_domain_N; ++j) {
            result += u[index_inside_domain(i, j)] * v[index_inside_domain(i, j)];
        }
    }
    return result * koeff;
}

double sqr_euclidean_norm(double *u)
{
    return scalar_product(u, u);
}

void matrix_linear_combination(double *result, double alpha, double betta,
                 double *M1, double *M2)
{
    int i, j;
    for (i = 0; i < inside_domain_M; ++i) {
        for (j = 0; j < inside_domain_N; ++j) {
            result[index_inside_domain(i, j)] = alpha * M1[index_inside_domain(i, j)] + betta * M2[index_inside_domain(i, j)];
        }
    }
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

void get_aij(double *result, double eps)
{
    double inv_h2 = 1 / step_y;
    int i, j;
    for (i = 1; i < M ; ++i) {
        for (j = 1; j < N; ++j) {
            double x = get_xi_by_index(i), y = get_yj_by_index(j);
            double l = calculate_lij_y(x - 0.5 * step_x, y - 0.5 * step_y, y + 0.5 * step_y);
            result[global_index(i, j)] = inv_h2 * l + (1.0 - inv_h2 * l) / eps;
        }
    }
}

void get_bij(double *result, double eps)
{
    int i, j;
    for (i = 1; i < M; ++i) {
        for (j = 1; j < N; ++j) {
            double x = get_xi_by_index(i), y = get_yj_by_index(j);
            double l = calculate_lij_x(y - 0.5 * step_y, x - 0.5 * step_x, x + 0.5 * step_x);
            result[global_index(i, j)] = (1/step_x) * l + (1.0 - (1/step_x) * l) / eps;
        }
    }
}

double *initialize_to_zeros(size_t size)
{
    size_t i;
    double *grid = malloc(size * sizeof(double));
    for (i = 0; i < size; ++i)
        grid[i] = 0.0;
    return grid;
}

void initialize_right_part(double *grid, double (*grinit)(double, double))
{
    int i, j;
    for (i = 0; i < M; ++i) {
        for (j = 0; j < N; ++j) {
            double x = get_xi_by_index(i), y = get_yj_by_index(j);
            grid[domain_index(i, j)] = grinit(x, y);
        }
    }
}

void operator_A(double *Aw, double *w, double *border_w,
                         double *aij, double *bij)
{
    int i, j;
    for (i = 0; i < inside_domain_M; ++i) {
        for (j = 0; j < inside_domain_N; ++j) {
            double a1, a2, b1, b2;
            double right_der_x, left_der_x, right_der_y, left_der_y;
            int g_i = global_i(i, proc_id), g_j = global_j(j, proc_id);
            if (g_i == 0 || g_i == M - 1 || g_j == 0 || g_j == N - 1)
                continue;
            a1 = aij[global_index(g_i, g_j)];
            a2 = aij[global_index(g_i + 1, g_j)];
            b1 = bij[global_index(g_i, g_j)];
            b2 = bij[global_index(g_i, g_j + 1)];
            right_der_x = right_derivative_x(w, border_w, g_i, g_j, i, j);
            left_der_x = left_derivative_x(w, border_w, g_i, g_j, i, j);
            right_der_y = right_derivative_y(w, border_w, g_i, g_j, i, j);
            left_der_y = left_derivative_y(w, border_w, g_i, g_j, i, j);
            Aw[index_inside_domain(i, j)] = -(a2 * right_der_x - a1 * left_der_x) / step_x - (b2 * right_der_y - b1 * left_der_y) / step_y;
        }
    }
}

void get_result(double *result, const char *file)
{
    int i, j;
    FILE *fp = fopen(file, "w");
    if (!fp) {
        perror("file doesn't open");
        return;
    }
    for (i = 0; i < M; ++i) {
        for (j = 0; j < N; ++j) {
            double x = get_xi_by_index(i), y = get_yj_by_index(j);
            fprintf(fp, "%lf\t%lf\t%lf\n", x, y, result[domain_index(i, j)]);
        }
    }
    fclose(fp);
}

void get_border(double *border, const double *grid)
{
    int i, j;
    for (i = 0; i < inside_domain_M; ++i) {
        border[i] = grid[index_inside_domain(i, 0)];
        border[inside_domain_M + i] = grid[index_inside_domain(i, inside_domain_N - 1)];
    }
    for (j = 0; j < inside_domain_N; ++j) {
        border[2*inside_domain_M + j] = grid[index_inside_domain(0, j)];
        border[2*inside_domain_M + inside_domain_N + j] = grid[index_inside_domain(inside_domain_M - 1, j)];
    }
}

void solution_method(void)
{
    double delta, tau;
    double *tmp, local_vars[3], global_vars[3];
    unsigned long iteration_num = 0;
    double *w, *Aw, *w_next;
    double *r, *Ar, *b;
    double *aij, *bij;
    double *global_b, *global_w, *border_w, *border_r;
    double *global_border_w, *global_border_r;
    
    w = initialize_to_zeros(domain_size);
    Aw = initialize_to_zeros(domain_size);
    w_next = initialize_to_zeros(domain_size);
    border_r = initialize_to_zeros(domain_perim);
    border_w = initialize_to_zeros(domain_perim);
    global_border_r = initialize_to_zeros(sum_perim);
    global_border_w = initialize_to_zeros(sum_perim);
    r = initialize_to_zeros(domain_size);
    Ar = initialize_to_zeros(domain_size);

    aij = initialize_to_zeros(grid_size);
    bij = initialize_to_zeros(grid_size);
    b = initialize_to_zeros(domain_size);

    global_b = NULL;
    global_w = NULL;
    if (proc_id == 0) {
        global_b = initialize_to_zeros(grid_size);
        global_w = initialize_to_zeros(grid_size);
        initialize_right_part(global_b, convert_right_part);
        get_aij(aij, epsilon);
        get_bij(bij, epsilon);
    }

    MPI_Scatter(global_b, domain_size, MPI_DOUBLE, b, domain_size,
                MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Bcast(aij, grid_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(bij, grid_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    do {

        get_border(border_w, w);

        MPI_Allgather(border_w, domain_perim, MPI_DOUBLE,
                      global_border_w, domain_perim, MPI_DOUBLE,
                      MPI_COMM_WORLD);

        operator_A(Aw, w, global_border_w, aij, bij);
        matrix_linear_combination(r, 1.0, -1.0, Aw, b);

        get_border(border_r, r);

        MPI_Allgather(border_r, domain_perim, MPI_DOUBLE,
                      global_border_r, domain_perim, MPI_DOUBLE,
                      MPI_COMM_WORLD);

        operator_A(Ar, r, global_border_r, aij, bij);

        local_vars[0] = scalar_product(Ar, r);
        local_vars[1] = sqr_euclidean_norm(Ar);
        local_vars[2] = sqr_euclidean_norm(r);

        MPI_Allreduce(local_vars, global_vars, 3,
                      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        delta = sqrt(global_vars[2]);

        tau = global_vars[0] / (global_vars[1] + 0.000001);

        matrix_linear_combination(w_next, 1.0, -tau, w, r);

        tmp = w;
        w = w_next;
        w_next = tmp;

        if (proc_id == 0) {
            if (iteration_num % 10000 == 0)
                fprintf(stderr, "delta = %lf\n", delta);
            ++iteration_num;
        }
    } while (delta > delta_stop);

    MPI_Gather(w, domain_size, MPI_DOUBLE, global_w,
               domain_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (proc_id == 0) {
        if (file_name)
            get_result(global_w, file_name);
    }
}

int is_power_of_two(int number) {
    return ((number) && !((number) & ((number) - 1)));
}

void get_params(int argc, char **argv) {
    M = atoi(argv[1]);
    N = atoi(argv[2]);
    delta_stop = atof(argv[3]);
    step_x = (B1-A1)/(M-1);
    step_y = (B2-A2)/(N-1);
    if(argc > 4) {
        file_name = argv[4];
    }
    double h_max = max(step_x, step_y);
    grid_size = N*M;
    epsilon = h_max*h_max;
    calculate_domain_grid_sizes();
}

int main(int argc, char **argv)
{
    double start_time = 0.0;
    if(argc < 4){
        fprintf(stderr, "you need to enter M and N");
        exit(1);
    }
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &proc_amount);
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
    get_params(argc, argv);
    if (proc_id == 0) {
        fprintf(stderr, "The amount of processes: %d\n", proc_amount);
        fprintf(stderr, "Grid sizes = (%d, %d)\n", M, N);
        start_time = MPI_Wtime();
    }
    solution_method();
    MPI_Barrier(MPI_COMM_WORLD);
    if (proc_id == 0) {
        fprintf(stderr,"Time: %lf\n", MPI_Wtime() - start_time);
    }
    MPI_Finalize();
    return 0;
}

