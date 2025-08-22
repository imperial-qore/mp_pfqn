#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>

/* Default ranges */
#define DEFAULT_M_MIN 2
#define DEFAULT_M_MAX 10
#define DEFAULT_R_MIN 2
#define DEFAULT_R_MAX 5
#define DEFAULT_N_MIN 1
#define DEFAULT_N_MAX 15
#define DEFAULT_MI_MIN 1
#define DEFAULT_MI_MAX 2
#define DEFAULT_L_MIN 1
#define DEFAULT_L_MAX 10000
#define DEFAULT_Z_MIN 0
#define DEFAULT_Z_MAX 10000

typedef struct {
    int M_min, M_max;      /* stations range */
    int R_min, R_max;      /* classes range */
    int N_min, N_max;      /* jobs per class range */
    int mi_min, mi_max;    /* multiplicities range */
    int L_min, L_max;      /* demands range */
    int Z_min, Z_max;      /* think times range */
    int M, R, N;           /* specific values (for backward compatibility) */
    int seed;
    bool use_specific_values;
} rndmodel_params;

int randi(int min, int max)
{
    if (min >= max) return min;
    return min + (rand() % (max - min + 1));
}

void print_usage(const char* progname)
{
    printf("USAGE: %s [options] [M R N seed]\n\n", progname);
    printf("Generate random queueing network models.\n\n");
    printf("Backward compatibility mode:\n");
    printf("  %s M R N seed    Generate model with M stations, R classes, N total jobs\n\n", progname);
    printf("Options mode:\n");
    printf("  -h, --help              : Show this help\n");
    printf("  -s, --seed SEED         : Random seed (default: current time)\n");
    printf("  -M, --stations MIN MAX  : Range for number of stations (default: %d-%d)\n", DEFAULT_M_MIN, DEFAULT_M_MAX);
    printf("  -R, --classes MIN MAX   : Range for number of classes (default: %d-%d)\n", DEFAULT_R_MIN, DEFAULT_R_MAX);
    printf("  -N, --jobs MIN MAX      : Range for jobs per class (default: %d-%d)\n", DEFAULT_N_MIN, DEFAULT_N_MAX);
    printf("  -m, --multi MIN MAX     : Range for station multiplicities (default: %d-%d)\n", DEFAULT_MI_MIN, DEFAULT_MI_MAX);
    printf("  -L, --demands MIN MAX   : Range for service demands (default: %d-%d)\n", DEFAULT_L_MIN, DEFAULT_L_MAX);
    printf("  -Z, --think MIN MAX     : Range for think times (default: %d-%d)\n\n", DEFAULT_Z_MIN, DEFAULT_Z_MAX);
    printf("Examples:\n");
    printf("  %s 5 3 20 42                    # 5 stations, 3 classes, 20 total jobs, seed 42\n", progname);
    printf("  %s -M 3 8 -R 2 4 -s 123         # 3-8 stations, 2-4 classes, seed 123\n", progname);
    printf("  %s -N 5 15 -L 100 5000          # 5-15 jobs per class, demands 100-5000\n", progname);
}

int parse_arguments(int argc, char **argv, rndmodel_params *params)
{
    /* Initialize default values */
    params->M_min = DEFAULT_M_MIN; params->M_max = DEFAULT_M_MAX;
    params->R_min = DEFAULT_R_MIN; params->R_max = DEFAULT_R_MAX;
    params->N_min = DEFAULT_N_MIN; params->N_max = DEFAULT_N_MAX;
    params->mi_min = DEFAULT_MI_MIN; params->mi_max = DEFAULT_MI_MAX;
    params->L_min = DEFAULT_L_MIN; params->L_max = DEFAULT_L_MAX;
    params->Z_min = DEFAULT_Z_MIN; params->Z_max = DEFAULT_Z_MAX;
    params->seed = time(NULL);
    params->use_specific_values = false;

    /* Check for backward compatibility mode (4 arguments: M R N seed) */
    if (argc == 5 && argv[1][0] != '-') {
        params->M = atoi(argv[1]);
        params->R = atoi(argv[2]);
        params->N = atoi(argv[3]);
        params->seed = atoi(argv[4]);
        params->use_specific_values = true;
        return 0;
    }

    /* Parse options */
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            print_usage(argv[0]);
            exit(0);
        } else if (strcmp(argv[i], "-s") == 0 || strcmp(argv[i], "--seed") == 0) {
            if (i + 1 >= argc) {
                fprintf(stderr, "Error: %s requires an argument\n", argv[i]);
                return -1;
            }
            params->seed = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-M") == 0 || strcmp(argv[i], "--stations") == 0) {
            if (i + 2 >= argc) {
                fprintf(stderr, "Error: %s requires two arguments (min max)\n", argv[i]);
                return -1;
            }
            params->M_min = atoi(argv[++i]);
            params->M_max = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-R") == 0 || strcmp(argv[i], "--classes") == 0) {
            if (i + 2 >= argc) {
                fprintf(stderr, "Error: %s requires two arguments (min max)\n", argv[i]);
                return -1;
            }
            params->R_min = atoi(argv[++i]);
            params->R_max = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-N") == 0 || strcmp(argv[i], "--jobs") == 0) {
            if (i + 2 >= argc) {
                fprintf(stderr, "Error: %s requires two arguments (min max)\n", argv[i]);
                return -1;
            }
            params->N_min = atoi(argv[++i]);
            params->N_max = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-m") == 0 || strcmp(argv[i], "--multi") == 0) {
            if (i + 2 >= argc) {
                fprintf(stderr, "Error: %s requires two arguments (min max)\n", argv[i]);
                return -1;
            }
            params->mi_min = atoi(argv[++i]);
            params->mi_max = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-L") == 0 || strcmp(argv[i], "--demands") == 0) {
            if (i + 2 >= argc) {
                fprintf(stderr, "Error: %s requires two arguments (min max)\n", argv[i]);
                return -1;
            }
            params->L_min = atoi(argv[++i]);
            params->L_max = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-Z") == 0 || strcmp(argv[i], "--think") == 0) {
            if (i + 2 >= argc) {
                fprintf(stderr, "Error: %s requires two arguments (min max)\n", argv[i]);
                return -1;
            }
            params->Z_min = atoi(argv[++i]);
            params->Z_max = atoi(argv[++i]);
        } else {
            fprintf(stderr, "Error: Unknown option %s\n", argv[i]);
            return -1;
        }
    }

    /* Validate ranges */
    if (params->M_min > params->M_max || params->M_min < 1) {
        fprintf(stderr, "Error: Invalid station range [%d, %d]\n", params->M_min, params->M_max);
        return -1;
    }
    if (params->R_min > params->R_max || params->R_min < 1) {
        fprintf(stderr, "Error: Invalid class range [%d, %d]\n", params->R_min, params->R_max);
        return -1;
    }
    if (params->N_min > params->N_max || params->N_min < 1) {
        fprintf(stderr, "Error: Invalid jobs range [%d, %d]\n", params->N_min, params->N_max);
        return -1;
    }
    if (params->mi_min > params->mi_max || params->mi_min < 1) {
        fprintf(stderr, "Error: Invalid multiplicity range [%d, %d]\n", params->mi_min, params->mi_max);
        return -1;
    }
    if (params->L_min > params->L_max || params->L_min < 1) {
        fprintf(stderr, "Error: Invalid demand range [%d, %d]\n", params->L_min, params->L_max);
        return -1;
    }
    if (params->Z_min > params->Z_max || params->Z_min < 0) {
        fprintf(stderr, "Error: Invalid think time range [%d, %d]\n", params->Z_min, params->Z_max);
        return -1;
    }

    return 0;
}

int main(int argc, char **argv)
{
    rndmodel_params params;
    FILE *f = stdout;
    int L[1000][100];
    int M, R, r, i;

    if (argc < 2) {
        print_usage(argv[0]);
        return 1;
    }

    if (parse_arguments(argc, argv, &params) != 0) {
        return 1;
    }

    srand(params.seed);

    /* Generate or use specific values */
    if (params.use_specific_values) {
        /* Backward compatibility mode */
        M = params.M;
        R = params.R;
        fprintf(f, "%d\n", R);
        for (r = 1; r <= R; r++) {
            fprintf(f, "%d ", params.N / R);  /* Distribute total jobs equally */
        }
        fprintf(f, "\n");
        for (r = 1; r <= R; r++) {
            fprintf(f, "%d ", 0);  /* Think times = 0 for backward compatibility */
        }
        fprintf(f, "\n");
        fprintf(f, "%d\n", M);
        
        /* Generate demands using old range */
        for (r = 1; r <= R; r++) {
            for (i = 1; i <= M; i++) {
                L[i][r] = randi(1, 10000);
            }
        }
        
        for (i = 1; i <= M; i++) {
            fprintf(f, "%d ", 1);  /* Multiplicities = 1 for backward compatibility */
            for (r = 1; r <= R; r++) {
                fprintf(f, "%d ", L[i][r]);
            }
            fprintf(f, "\n");
        }
    } else {
        /* New options mode */
        R = randi(params.R_min, params.R_max);
        M = randi(params.M_min, params.M_max);
        
        fprintf(f, "%d\n", R);
        
        /* Generate job populations for each class */
        for (r = 1; r <= R; r++) {
            fprintf(f, "%d ", randi(params.N_min, params.N_max));
        }
        fprintf(f, "\n");
        
        /* Generate think times for each class */
        for (r = 1; r <= R; r++) {
            fprintf(f, "%d ", randi(params.Z_min, params.Z_max));
        }
        fprintf(f, "\n");
        
        fprintf(f, "%d\n", M);
        
        /* Generate demands */
        for (r = 1; r <= R; r++) {
            for (i = 1; i <= M; i++) {
                L[i][r] = randi(params.L_min, params.L_max);
            }
        }
        
        /* Generate stations with multiplicities and demands */
        for (i = 1; i <= M; i++) {
            fprintf(f, "%d ", randi(params.mi_min, params.mi_max));
            for (r = 1; r <= R; r++) {
                fprintf(f, "%d ", L[i][r]);
            }
            fprintf(f, "\n");
        }
    }

    return 0;
}
