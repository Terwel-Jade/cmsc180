// lab05 - Distributed Z-Score Normalization  (parallel send + receive)
//
// Extends lab04 to:
//   (1) Slave computes Z-Score Normalization column-wise on its submatrix of X,
//       producing its part of matrix T.
//   (2) Slave sends its part of T back to master.
//   (3) Master times from before distributing X to after rebuilding full T.
//   (4) Slave times ONLY the Z-Score computation (after rx, before tx).
//
// Parallel improvement (vs original lab05):
//   - Master pre-splits all X submatrices into flat buffers BEFORE timing.
//   - One pthread per slave handles: connect → send X → receive T → store.
//   - Master clock covers ONLY the parallel network phase.
//
// Build :  gcc -O2 -o lab05 lab05.c -lm -lpthread
// Config:  config.txt
//            master  <ip>   <port>
//            slave   <ip>   <port>
//            ...
// Usage :  ./lab05 <n> <p> <s>
//            n  matrix size (n x n)
//            p  port number this process listens / connects on
//            s  0 = master   1 = slave

#define _POSIX_C_SOURCE 200112L
#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <math.h>
#include <errno.h>
#include <arpa/inet.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <pthread.h>
#include <sched.h>

/* ------------------------------------------------------------------ */
/*  Peer info                                                          */
/* ------------------------------------------------------------------ */
typedef struct {
    char role[16];
    char ip[64];
    int  port;
} PeerInfo;

/* ------------------------------------------------------------------ */
/*  Per-slave task (one per pthread)                                   */
/*                                                                     */
/*  send_buf  : pre-split flat X submatrix  (input,  owned by master) */
/*  recv_buf  : flat T submatrix buffer     (output, allocated here)  */
/* ------------------------------------------------------------------ */
typedef struct {
    /* identity */
    int  slave_idx;
    char ip[64];
    int  port;

    /* dimensions */
    int  sub_rows;
    int  cols;

    /* data */
    double *send_buf;   /* X submatrix — filled before thread spawned  */
    double *recv_buf;   /* T submatrix — malloc'd by thread, read after join */

    /* result: 0 = ok, -1 = error */
    int result;
} SlaveTask;

/* ------------------------------------------------------------------ */
/*  Matrix helpers                                                     */
/* ------------------------------------------------------------------ */
static double **alloc_matrix(int rows, int cols)
{
    double **M = malloc(rows * sizeof(double *));
    if (!M) { perror("malloc"); exit(1); }
    for (int i = 0; i < rows; i++) {
        M[i] = malloc(cols * sizeof(double));
        if (!M[i]) { perror("malloc"); exit(1); }
    }
    return M;
}

static void free_matrix(double **M, int rows)
{
    for (int i = 0; i < rows; i++) free(M[i]);
    free(M);
}

static void gen_rand_matrix(double **M, int rows, int cols)
{
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            M[i][j] = (double)(rand() % 100 + 1);
}

static void print_matrix(double **M, int rows, int cols, const char *label)
{
    printf("\n--- %s (%d x %d) ---\n", label, rows, cols);
    int rlim = (rows < 8) ? rows : 8;
    int clim = (cols < 8) ? cols : 8;
    for (int i = 0; i < rlim; i++) {
        for (int j = 0; j < clim; j++)
            printf("%8.4f ", M[i][j]);
        if (cols > clim) printf("...");
        printf("\n");
    }
    if (rows > rlim) printf("  ... (%d more rows)\n", rows - rlim);
}

/* ------------------------------------------------------------------ */
/*  Z-Score Normalization (column-wise on a submatrix)                */
/* ------------------------------------------------------------------ */
static void zscore_normalize(double **X, double **T, int rows, int cols)
{
    for (int j = 0; j < cols; j++) {
        double mean = 0.0;
        for (int i = 0; i < rows; i++)
            mean += X[i][j];
        mean /= rows;

        double variance = 0.0;
        for (int i = 0; i < rows; i++) {
            double diff = X[i][j] - mean;
            variance += diff * diff;
        }
        double std = sqrt(variance / rows);

        for (int i = 0; i < rows; i++)
            T[i][j] = (std < 1e-12) ? 0.0 : (X[i][j] - mean) / std;
    }
}

/* ------------------------------------------------------------------ */
/*  Config reader                                                      */
/* ------------------------------------------------------------------ */
static int read_config(PeerInfo peers[], int max)
{
    FILE *f = fopen("config.txt", "r");
    if (!f) { perror("fopen config.txt"); exit(1); }

    int count = 0;
    char line[256];
    while (count < max && fgets(line, sizeof(line), f)) {
        char *p = line;
        while (*p == ' ' || *p == '\t') p++;
        if (*p == '#' || *p == '\n' || *p == '\r' || *p == '\0')
            continue;
        if (sscanf(p, "%15s %63s %d",
                   peers[count].role,
                   peers[count].ip,
                   &peers[count].port) == 3)
            count++;
    }
    fclose(f);
    return count;
}

/* ------------------------------------------------------------------ */
/*  Socket helpers                                                     */
/* ------------------------------------------------------------------ */
static int send_all(int fd, const void *buf, size_t len)
{
    const char *p = buf;
    size_t sent = 0;
    while (sent < len) {
        ssize_t n = send(fd, p + sent, len - sent, 0);
        if (n <= 0) return -1;
        sent += n;
    }
    return 0;
}

static int recv_all(int fd, void *buf, size_t len)
{
    char *p = buf;
    size_t got = 0;
    while (got < len) {
        ssize_t n = recv(fd, p + got, len - got, 0);
        if (n <= 0) return -1;
        got += n;
    }
    return 0;
}

static int make_server_socket(int port)
{
    int sfd = socket(AF_INET, SOCK_STREAM, 0);
    if (sfd < 0) { perror("socket"); exit(1); }

    int opt = 1;
    setsockopt(sfd, SOL_SOCKET, SO_REUSEADDR, &opt, sizeof(opt));

    struct sockaddr_in addr = {0};
    addr.sin_family      = AF_INET;
    addr.sin_addr.s_addr = INADDR_ANY;
    addr.sin_port        = htons(port);

    if (bind(sfd, (struct sockaddr *)&addr, sizeof(addr)) < 0) {
        perror("bind"); exit(1);
    }
    if (listen(sfd, 64) < 0) { perror("listen"); exit(1); }
    return sfd;
}

static int connect_to(const char *ip, int port)
{
    int fd = socket(AF_INET, SOCK_STREAM, 0);
    if (fd < 0) { perror("socket"); exit(1); }

    struct sockaddr_in addr = {0};
    addr.sin_family = AF_INET;
    addr.sin_port   = htons(port);
    if (inet_pton(AF_INET, ip, &addr.sin_addr) <= 0) {
        fprintf(stderr, "Invalid IP: %s\n", ip); exit(1);
    }

    for (int attempt = 0; attempt < 10; attempt++) {
        if (connect(fd, (struct sockaddr *)&addr, sizeof(addr)) == 0)
            return fd;
        fprintf(stderr, "[MASTER] Retrying connection to %s:%d (%d/10)...\n",
                ip, port, attempt + 1);
        sleep(1);
    }
    fprintf(stderr, "Could not connect to %s:%d\n", ip, port);
    exit(1);
}

/* ------------------------------------------------------------------ */
/*  Core-affinity helper (for slaves)                                 */
/* ------------------------------------------------------------------ */
static void pin_process_to_core(int core_id)
{
    int num_cores = (int)sysconf(_SC_NPROCESSORS_ONLN);
    if (num_cores <= 1) return;

    core_id = core_id % (num_cores - 1);

    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(core_id, &cpuset);

    if (sched_setaffinity(0, sizeof(cpu_set_t), &cpuset) != 0)
        perror("sched_setaffinity (non-fatal)");
    else
        printf("[SLAVE] Pinned to CPU core %d (of %d)\n", core_id, num_cores);
}

/* ================================================================== */
/*  Thread worker                                                      */
/*                                                                     */
/*  For one slave:                                                     */
/*    1. connect                                                       */
/*    2. send header + flat X submatrix buffer                        */
/*    3. receive header + flat T submatrix into recv_buf              */
/*    4. validate dimensions, store result                            */
/* ================================================================== */
static void *slave_exchange_thread(void *arg)
{
    SlaveTask *task = (SlaveTask *)arg;
    task->result    = -1;

    int fd = connect_to(task->ip, task->port);
    printf("[MASTER] [thread slave=%d] Connected to %s:%d\n",
           task->slave_idx, task->ip, task->port);

    /* ---- send X submatrix ---- */
    int32_t hdr_out[2] = { htonl(task->sub_rows), htonl(task->cols) };
    if (send_all(fd, hdr_out, sizeof(hdr_out)) < 0) {
        fprintf(stderr, "[thread slave=%d] send header failed\n", task->slave_idx);
        close(fd); return NULL;
    }

    size_t data_bytes = (size_t)task->sub_rows * task->cols * sizeof(double);
    if (send_all(fd, task->send_buf, data_bytes) < 0) {
        fprintf(stderr, "[thread slave=%d] send X data failed\n", task->slave_idx);
        close(fd); return NULL;
    }
    printf("[MASTER] [thread slave=%d] X submatrix sent (%d x %d)\n",
           task->slave_idx, task->sub_rows, task->cols);

    /* ---- receive T submatrix ---- */
    int32_t hdr_in[2];
    if (recv_all(fd, hdr_in, sizeof(hdr_in)) < 0) {
        fprintf(stderr, "[thread slave=%d] recv T header failed\n", task->slave_idx);
        close(fd); return NULL;
    }

    int rx_rows = ntohl(hdr_in[0]);
    int rx_cols = ntohl(hdr_in[1]);

    if (rx_rows != task->sub_rows || rx_cols != task->cols) {
        fprintf(stderr,
                "[thread slave=%d] dimension mismatch: expected %dx%d got %dx%d\n",
                task->slave_idx, task->sub_rows, task->cols, rx_rows, rx_cols);
        close(fd); return NULL;
    }

    /* Allocate recv_buf — master reads it after pthread_join */
    task->recv_buf = malloc(data_bytes);
    if (!task->recv_buf) { perror("malloc recv_buf"); close(fd); return NULL; }

    if (recv_all(fd, task->recv_buf, data_bytes) < 0) {
        fprintf(stderr, "[thread slave=%d] recv T data failed\n", task->slave_idx);
        free(task->recv_buf);
        task->recv_buf = NULL;
        close(fd); return NULL;
    }

    printf("[MASTER] [thread slave=%d] T submatrix received (%d x %d)\n",
           task->slave_idx, rx_rows, rx_cols);

    close(fd);
    task->result = 0;
    return NULL;
}

/* ================================================================== */
/*  MASTER                                                             */
/*                                                                     */
/*  Phase 1 (pre-timing): generate X, pre-split into flat buffers.   */
/*  Phase 2 (timed):      spawn threads → join → copy T into place.  */
/* ================================================================== */
static void run_master(int n)
{
    /* ---- read config ---- */
    PeerInfo peers[64];
    int total = read_config(peers, 64);

    PeerInfo slaves[64];
    int t = 0;
    for (int i = 0; i < total; i++)
        if (strcmp(peers[i].role, "slave") == 0)
            slaves[t++] = peers[i];

    if (t == 0) { fprintf(stderr, "No slaves found in config.txt\n"); exit(1); }
    printf("[MASTER] n=%d | slaves=%d\n", n, t);

    /* ---- generate random matrix X ---- */
    srand((unsigned)time(NULL));
    double **X = alloc_matrix(n, n);
    gen_rand_matrix(X, n, n);
    print_matrix(X, n, n, "Matrix X (master, original)");

    /* ---- allocate full T to be rebuilt ---- */
    double **T = alloc_matrix(n, n);

    /* ---- compute row-partition boundaries ---- */
    int base_rows  = n / t;
    int row_starts[64];
    int row_counts[64];
    {
        int rs = 0;
        for (int s = 0; s < t; s++) {
            row_starts[s] = rs;
            row_counts[s] = (s == t - 1) ? (n - rs) : base_rows;
            rs += row_counts[s];
        }
    }

    /* ==============================================================
     * PHASE 1 — Pre-split X into per-slave flat buffers (pre-timing)
     * ============================================================== */
    SlaveTask tasks[64];
    memset(tasks, 0, sizeof(tasks));

    for (int s = 0; s < t; s++) {
        int sub_rows  = row_counts[s];
        int row_start = row_starts[s];

        tasks[s].slave_idx = s;
        tasks[s].sub_rows  = sub_rows;
        tasks[s].cols      = n;
        tasks[s].recv_buf  = NULL;
        tasks[s].result    = -1;
        strncpy(tasks[s].ip, slaves[s].ip, sizeof(tasks[s].ip) - 1);
        tasks[s].port = slaves[s].port;

        size_t bytes = (size_t)sub_rows * n * sizeof(double);
        tasks[s].send_buf = malloc(bytes);
        if (!tasks[s].send_buf) { perror("malloc send_buf"); exit(1); }

        /* flatten: row-major copy */
        for (int r = 0; r < sub_rows; r++)
            memcpy(tasks[s].send_buf + (size_t)r * n,
                   X[row_start + r],
                   n * sizeof(double));

        printf("[MASTER] Pre-split: slave %d rows [%d, %d)  — %d x %d\n",
               s, row_start, row_start + sub_rows, sub_rows, n);
    }

    /* ==============================================================
     * PHASE 2 — Parallel exchange: start clock, spawn, join, stop
     * ============================================================== */
    struct timespec t_before, t_after;
    clock_gettime(CLOCK_MONOTONIC, &t_before);

    pthread_t threads[64];
    for (int s = 0; s < t; s++) {
        if (pthread_create(&threads[s], NULL, slave_exchange_thread, &tasks[s]) != 0) {
            perror("pthread_create"); exit(1);
        }
    }

    for (int s = 0; s < t; s++)
        pthread_join(threads[s], NULL);

    clock_gettime(CLOCK_MONOTONIC, &t_after);

    /* ==============================================================
     * Rebuild T from per-slave flat recv_bufs
     * ============================================================== */
    for (int s = 0; s < t; s++) {
        if (tasks[s].result != 0 || !tasks[s].recv_buf) {
            fprintf(stderr, "[MASTER] slave %d failed — T may be incomplete.\n", s);
            continue;
        }
        int sub_rows  = row_counts[s];
        int row_start = row_starts[s];

        for (int r = 0; r < sub_rows; r++)
            memcpy(T[row_start + r],
                   tasks[s].recv_buf + (size_t)r * n,
                   n * sizeof(double));
    }

    double elapsed = (t_after.tv_sec  - t_before.tv_sec) +
                     (t_after.tv_nsec - t_before.tv_nsec) * 1e-9;

    /* ---- display results ---- */
    print_matrix(T, n, n, "Matrix T (master, fully rebuilt Z-Score)");

    printf("\n[MASTER] Parallel exchange summary:\n");
    for (int s = 0; s < t; s++)
        printf("         slave %d → %s\n",
               s, tasks[s].result == 0 ? "OK" : "FAILED");

    printf("[MASTER] time_elapsed (parallel distribute X → rebuild T) = %.6f seconds\n",
           elapsed);

    /* ---- cleanup ---- */
    for (int s = 0; s < t; s++) {
        free(tasks[s].send_buf);
        free(tasks[s].recv_buf);   /* NULL-safe */
    }
    free_matrix(X, n);
    free_matrix(T, n);
}

/* ================================================================== */
/*  SLAVE — unchanged in behaviour                                     */
/*                                                                     */
/*  Protocol:                                                          */
/*    MASTER --> SLAVE :  hdr[sub_rows, n]  +  X submatrix data       */
/*    (time_before)                                                    */
/*    compute T = zscore(X_sub)                                       */
/*    (time_after)  --> report time_elapsed                           */
/*    SLAVE  --> MASTER:  hdr[sub_rows, n]  +  T submatrix data       */
/*                                                                     */
/*  Timing: ONLY the Z-Score computation kernel.                      */
/* ================================================================== */
static void run_slave(int port)
{
    PeerInfo peers[64];
    int total = read_config(peers, 64);
    char master_ip[64] = "unknown";
    for (int i = 0; i < total; i++)
        if (strcmp(peers[i].role, "master") == 0) {
            strncpy(master_ip, peers[i].ip, sizeof(master_ip) - 1);
            break;
        }
    printf("[SLAVE  port=%d] Master IP from config: %s\n", port, master_ip);

    int sfd = make_server_socket(port);
    printf("[SLAVE  port=%d] Listening...\n", port);

    struct sockaddr_in cli_addr;
    socklen_t cli_len = sizeof(cli_addr);
    int cfd = accept(sfd, (struct sockaddr *)&cli_addr, &cli_len);
    if (cfd < 0) { perror("accept"); exit(1); }

    char peer_ip[INET_ADDRSTRLEN];
    inet_ntop(AF_INET, &cli_addr.sin_addr, peer_ip, sizeof(peer_ip));
    printf("[SLAVE  port=%d] Master connected from %s\n", port, peer_ip);

    /* ---- receive X submatrix ---- */
    int32_t hdr[2];
    recv_all(cfd, hdr, sizeof(hdr));
    int sub_rows = ntohl(hdr[0]);
    int cols     = ntohl(hdr[1]);
    printf("[SLAVE  port=%d] Receiving X submatrix %d x %d\n", port, sub_rows, cols);

    double **X_sub = alloc_matrix(sub_rows, cols);
    for (int r = 0; r < sub_rows; r++)
        recv_all(cfd, X_sub[r], cols * sizeof(double));

    printf("[SLAVE  port=%d] X submatrix received.\n", port);

    char lbl_x[80];
    snprintf(lbl_x, sizeof(lbl_x), "X submatrix received (slave port=%d)", port);
    print_matrix(X_sub, sub_rows, cols, lbl_x);

    /* ---- allocate T submatrix ---- */
    double **T_sub = alloc_matrix(sub_rows, cols);

    /* ================================================================
     * time_before — ONLY the computation
     * ================================================================ */
    struct timespec t_before, t_after;
    clock_gettime(CLOCK_MONOTONIC, &t_before);

    zscore_normalize(X_sub, T_sub, sub_rows, cols);

    /* ================================================================
     * time_after — computation done, before sending T
     * ================================================================ */
    clock_gettime(CLOCK_MONOTONIC, &t_after);

    double elapsed = (t_after.tv_sec  - t_before.tv_sec) +
                     (t_after.tv_nsec - t_before.tv_nsec) * 1e-9;

    char lbl_t[80];
    snprintf(lbl_t, sizeof(lbl_t), "T submatrix (Z-Score, slave port=%d)", port);
    print_matrix(T_sub, sub_rows, cols, lbl_t);

    printf("\n[SLAVE  port=%d] time_elapsed (Z-Score computation ONLY) = %.6f seconds\n",
           port, elapsed);

    /* ---- send T submatrix back to master ---- */
    int32_t hdr_out[2] = { htonl(sub_rows), htonl(cols) };
    send_all(cfd, hdr_out, sizeof(hdr_out));

    for (int r = 0; r < sub_rows; r++)
        send_all(cfd, T_sub[r], cols * sizeof(double));

    printf("[SLAVE  port=%d] T submatrix sent back to master.\n", port);

    free_matrix(X_sub, sub_rows);
    free_matrix(T_sub, sub_rows);
    close(cfd);
    close(sfd);
}

/* ================================================================== */
/*  main                                                               */
/* ================================================================== */
int main(int argc, char *argv[])
{
    if (argc != 4) {
        fprintf(stderr,
            "Usage: %s <n> <p> <s>\n"
            "  n : matrix size (n x n)\n"
            "  p : port number\n"
            "  s : role  0=master  1=slave\n",
            argv[0]);
        return 1;
    }

    int n    = atoi(argv[1]);
    int port = atoi(argv[2]);
    int role = atoi(argv[3]);

    if (n <= 0 || port <= 0 || (role != 0 && role != 1)) {
        fprintf(stderr, "Invalid arguments.\n");
        return 1;
    }

    if (role == 0) {
        run_master(n);
    } else {
        pin_process_to_core(port);
        run_slave(port);
    }

    return 0;
}