// Accept n inputs from user, port number p, and role s (0=master, 1=slave)
// Read config.txt to find slave IPs and ports (for master) or master IP (for slave)
// Master creates random n x n matrix M, computes submatrix ranges for each slave
// Master pre-splits all submatrices, then sends to ALL slaves simultaneously (parallel)
// Slave time_before, receive submatrix, time_after
// Display received submatrix and elapsed time on slave side
// Show/output elapsed time

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

// ---------------------------------------------------------------------------
// Peer info
// ---------------------------------------------------------------------------
typedef struct {
    char role[16];   // "master" or "slave"
    char ip[64];
    int  port;
} PeerInfo;

// ---------------------------------------------------------------------------
// Per-slave send task (passed to each thread)
// ---------------------------------------------------------------------------
typedef struct {
    char   ip[64];
    int    port;
    int    slave_idx;
    int    sub_rows;   // number of rows in this sub-matrix
    int    cols;       // number of columns (= n)
    double *buf;       // flat pre-split buffer: sub_rows * cols doubles
    int    result;     // 0 = success, -1 = failure (written by thread)
} SendTask;

// ---------------------------------------------------------------------------
// Matrix helpers
// ---------------------------------------------------------------------------
static double **alloc_matrix(int rows, int cols) {
    double **M = malloc(rows * sizeof(double *));
    for (int i = 0; i < rows; i++)
        M[i] = malloc(cols * sizeof(double));
    return M;
}

static void free_matrix(double **M, int rows) {
    for (int i = 0; i < rows; i++) free(M[i]);
    free(M);
}

static void gen_rand_matrix(double **M, int rows, int cols) {
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            M[i][j] = (double)(rand() % 100 + 1);   // non-zero positive
}

static void print_matrix(double **M, int rows, int cols, const char *label) {
    printf("\n--- %s (%d x %d) ---\n", label, rows, cols);
    int rlim = (rows < 10) ? rows : 10;
    int clim = (cols < 10) ? cols : 10;
    for (int i = 0; i < rlim; i++) {
        for (int j = 0; j < clim; j++)
            printf("%6.1f ", M[i][j]);
        if (cols > clim) printf("...");
        printf("\n");
    }
    if (rows > rlim) printf("  ... (%d more rows)\n", rows - rlim);
}

// ---------------------------------------------------------------------------
// Config reader
// ---------------------------------------------------------------------------
static int read_config(PeerInfo peers[], int max) {
    FILE *f = fopen("config.txt", "r");
    if (!f) { perror("fopen config"); exit(1); }
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

// ---------------------------------------------------------------------------
// Socket helpers
// ---------------------------------------------------------------------------
static int send_all(int fd, const void *buf, size_t len) {
    const char *p = buf;
    size_t sent = 0;
    while (sent < len) {
        ssize_t n = send(fd, p + sent, len - sent, 0);
        if (n <= 0) return -1;
        sent += n;
    }
    return 0;
}

static int recv_all(int fd, void *buf, size_t len) {
    char *p = buf;
    size_t got = 0;
    while (got < len) {
        ssize_t n = recv(fd, p + got, len - got, 0);
        if (n <= 0) return -1;
        got += n;
    }
    return 0;
}

static int make_server_socket(int port) {
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

static int connect_to(const char *ip, int port) {
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
        sleep(1);
    }
    fprintf(stderr, "Could not connect to %s:%d\n", ip, port);
    exit(1);
}

// ---------------------------------------------------------------------------
// Thread worker: connects to one slave and sends its pre-split submatrix
// ---------------------------------------------------------------------------
static void *send_submatrix_thread(void *arg) {
    SendTask *task = (SendTask *)arg;
    task->result   = -1;   // pessimistic default

    int fd = connect_to(task->ip, task->port);
    printf("[MASTER] [thread slave=%d] Connected to %s:%d\n",
           task->slave_idx, task->ip, task->port);

    // Send header: sub_rows, cols
    int32_t hdr[2] = { htonl(task->sub_rows), htonl(task->cols) };
    if (send_all(fd, hdr, sizeof(hdr)) < 0) {
        fprintf(stderr, "[thread slave=%d] send header failed\n", task->slave_idx);
        close(fd);
        return NULL;
    }

    // Send the entire flat submatrix buffer in one shot
    size_t total_bytes = (size_t)task->sub_rows * task->cols * sizeof(double);
    if (send_all(fd, task->buf, total_bytes) < 0) {
        fprintf(stderr, "[thread slave=%d] send data failed\n", task->slave_idx);
        close(fd);
        return NULL;
    }

    // Wait for ACK
    char ack_buf[4] = {0};
    if (recv_all(fd, ack_buf, 3) < 0) {
        fprintf(stderr, "[thread slave=%d] recv ACK failed\n", task->slave_idx);
        close(fd);
        return NULL;
    }

    printf("[MASTER] [thread slave=%d] ACK received: \"%s\"\n",
           task->slave_idx, ack_buf);

    close(fd);
    task->result = 0;
    return NULL;
}

// ---------------------------------------------------------------------------
// Master logic
//   Phase 1 – pre-split: copy each slave's rows into a flat heap buffer
//   Phase 2 – parallel send: spawn one thread per slave, join them all
// ---------------------------------------------------------------------------
static void run_master(int n) {
    // Read config to find slaves
    PeerInfo peers[64];
    int total = read_config(peers, 64);

    PeerInfo slaves[64];
    int t = 0;
    for (int i = 0; i < total; i++)
        if (strcmp(peers[i].role, "slave") == 0)
            slaves[t++] = peers[i];

    if (t == 0) { fprintf(stderr, "No slaves in config.\n"); exit(1); }
    printf("[MASTER] Found %d slave(s) in config.\n", t);

    // Create random matrix M of size n x n
    srand(time(NULL));
    double **M = alloc_matrix(n, n);
    gen_rand_matrix(M, n, n);
    print_matrix(M, n, n, "Full matrix M (master)");

    int base_rows = n / t;

    // ------------------------------------------------------------------
    // PHASE 1: Pre-split — build one flat buffer per slave BEFORE timing
    // ------------------------------------------------------------------
    SendTask tasks[64];
    int row_start = 0;

    for (int s = 0; s < t; s++) {
        int sub_rows = (s == t - 1) ? (n - row_start) : base_rows;

        tasks[s].slave_idx = s;
        tasks[s].sub_rows  = sub_rows;
        tasks[s].cols      = n;
        tasks[s].result    = -1;
        strncpy(tasks[s].ip,   slaves[s].ip,   sizeof(tasks[s].ip)   - 1);
        tasks[s].port = slaves[s].port;

        // Flatten the submatrix rows into a contiguous buffer
        tasks[s].buf = malloc((size_t)sub_rows * n * sizeof(double));
        if (!tasks[s].buf) { perror("malloc submatrix buf"); exit(1); }

        for (int r = 0; r < sub_rows; r++)
            memcpy(tasks[s].buf + (size_t)r * n,
                   M[row_start + r],
                   n * sizeof(double));

        printf("[MASTER] Pre-split: slave %d gets rows [%d, %d) — %d rows x %d cols\n",
               s, row_start, row_start + sub_rows, sub_rows, n);

        row_start += sub_rows;
    }

    // ------------------------------------------------------------------
    // PHASE 2: Parallel send — start the clock, spawn threads, join all
    // ------------------------------------------------------------------
    struct timespec t_before, t_after;
    clock_gettime(CLOCK_MONOTONIC, &t_before);

    pthread_t threads[64];
    for (int s = 0; s < t; s++) {
        if (pthread_create(&threads[s], NULL, send_submatrix_thread, &tasks[s]) != 0) {
            perror("pthread_create"); exit(1);
        }
    }

    // Wait for all threads to finish
    for (int s = 0; s < t; s++)
        pthread_join(threads[s], NULL);

    clock_gettime(CLOCK_MONOTONIC, &t_after);

    // ------------------------------------------------------------------
    // Report results
    // ------------------------------------------------------------------
    double elapsed = (t_after.tv_sec  - t_before.tv_sec) +
                     (t_after.tv_nsec - t_before.tv_nsec) * 1e-9;

    printf("\n[MASTER] All parallel sends complete.\n");
    for (int s = 0; s < t; s++) {
        printf("[MASTER]   slave %d → %s\n",
               s, tasks[s].result == 0 ? "OK" : "FAILED");
        free(tasks[s].buf);
    }
    printf("[MASTER] time_elapsed (parallel send only) = %.6f seconds\n", elapsed);

    free_matrix(M, n);
}

// ---------------------------------------------------------------------------
// Slave logic — unchanged in behaviour
// ---------------------------------------------------------------------------
static void run_slave(int port) {
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

    // time_before
    struct timespec t_before, t_after;
    clock_gettime(CLOCK_MONOTONIC, &t_before);

    // Receive header: sub_rows, cols
    int32_t hdr[2];
    recv_all(cfd, hdr, sizeof(hdr));
    int sub_rows = ntohl(hdr[0]);
    int cols     = ntohl(hdr[1]);
    printf("[SLAVE  port=%d] Expecting submatrix %d x %d\n", port, sub_rows, cols);

    // Receive submatrix
    double **sub = alloc_matrix(sub_rows, cols);
    for (int r = 0; r < sub_rows; r++)
        recv_all(cfd, sub[r], cols * sizeof(double));

    // Send ACK
    send_all(cfd, "ack", 3);

    // time_after
    clock_gettime(CLOCK_MONOTONIC, &t_after);

    double elapsed = (t_after.tv_sec  - t_before.tv_sec) +
                     (t_after.tv_nsec - t_before.tv_nsec) * 1e-9;

    char label[64];
    snprintf(label, sizeof(label), "Submatrix received (slave port=%d)", port);
    print_matrix(sub, sub_rows, cols, label);
    printf("\n[SLAVE  port=%d] time_elapsed = %.6f seconds\n", port, elapsed);

    free_matrix(sub, sub_rows);
    close(cfd);
    close(sfd);
}

// ---------------------------------------------------------------------------
// Pin process to a CPU core (optional, for better isolation)
// ---------------------------------------------------------------------------
static void pin_process_to_core(int core_id) {
    int num_cores = (int)sysconf(_SC_NPROCESSORS_ONLN);
    core_id = (core_id % num_cores) - 1;
    if (core_id < 0) core_id = 0;   // guard against negative index

    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(core_id, &cpuset);

    if (sched_setaffinity(0, sizeof(cpu_set_t), &cpuset) != 0)
        perror("sched_setaffinity");
}

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------
int main(int argc, char *argv[]) {
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