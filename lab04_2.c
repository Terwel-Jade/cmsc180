/*
 * lab04.c — Master/Slave Matrix Distribution via Sockets
 *
 * Usage:
 *   ./lab04 <n> <p> <s>
 *     n  : size of the square matrix (n x n)
 *     p  : port number this instance listens/sends on
 *     s  : status — 0 = master, 1 = slave
 *
 * Configuration file: config.txt
 *   Format (one entry per line):
 *     <role> <ip> <port>
 *   Roles:
 *     master  <ip>  <port>
 *     slave   <ip>  <port>
 *   Example:
 *     master 127.0.0.1 5000
 *     slave  127.0.0.1 5001
 *     slave  127.0.0.1 5002
 *     slave  127.0.0.1 5003
 */

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

/* ──────────────────────────── constants ──────────────────────────── */
#define CONFIG_FILE   "config.txt"
#define MAX_PEERS     64
#define ACK_MSG       "ack"
#define ACK_LEN       3

/* ──────────────────────────── types ──────────────────────────── */
typedef struct {
    char role[16];   /* "master" or "slave" */
    char ip[64];
    int  port;
    char user[64];   /* SSH username, e.g. "pi" or "student" (optional) */
    char path[256];  /* remote path to lab04 binary dir        (optional) */
} PeerInfo;

/* ──────────────────────────── matrix helpers ──────────────────────────── */
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
            M[i][j] = (double)(rand() % 100 + 1);   /* non-zero positive */
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

/* ──────────────────────────── config reader ──────────────────────────── */
static int read_config(PeerInfo peers[], int max) {
    FILE *f = fopen(CONFIG_FILE, "r");
    if (!f) { perror("fopen config"); exit(1); }
    int count = 0;
    char line[512];
    while (count < max && fgets(line, sizeof(line), f)) {
        /* skip blank lines and comment lines starting with '#' */
        char *p = line;
        while (*p == ' ' || *p == '\t') p++;
        if (*p == '#' || *p == '\n' || *p == '\r' || *p == '\0')
            continue;

        /* Initialize optional fields */
        peers[count].user[0] = '\0';
        peers[count].path[0] = '\0';

        /* Try extended format: role ip port user remotepath */
        int matched = sscanf(p, "%15s %63s %d %63s %255s",
                             peers[count].role,
                             peers[count].ip,
                             &peers[count].port,
                             peers[count].user,
                             peers[count].path);
        if (matched >= 3)
            count++;
    }
    fclose(f);
    return count;
}

/* ──────────────────────────── socket helpers ──────────────────────────── */

/* Send exactly len bytes */
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

/* Receive exactly len bytes */
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

/* Create a listening server socket on given port */
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
    if (listen(sfd, MAX_PEERS) < 0) { perror("listen"); exit(1); }
    return sfd;
}

/* Connect to a remote peer, retry a few times */
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
        sleep(1);   /* slave may not be listening yet */
    }
    fprintf(stderr, "Could not connect to %s:%d\n", ip, port);
    exit(1);
}

/* ══════════════════════════════════════════════════════════════
 *  SSH REMOTE LAUNCH (optional)
 *  Launches slave instances on remote PCs via SSH.
 *  Requires:
 *    - SSH key-based auth set up (no password prompt)
 *    - lab04 binary + config.txt already copied to remote path
 *  config.txt extended format for SSH:
 *    slave <ip> <port> <user> <remote_path>
 *  Example:
 *    slave 192.168.1.11 5001 student /home/student/lab04
 * ══════════════════════════════════════════════════════════════ */
static void ssh_launch_slaves(PeerInfo slaves[], int t, int n) {
    for (int i = 0; i < t; i++) {
        /* Skip if no SSH info provided */
        if (slaves[i].user[0] == '\0' || slaves[i].path[0] == '\0') {
            printf("[MASTER] No SSH info for slave %d — assuming manually started.\n", i);
            continue;
        }

        /* Build SSH command:
           ssh -o StrictHostKeyChecking=no <user>@<ip>
               "cd <path> && nohup ./lab04 <n> <port> 1 > slave_<port>.log 2>&1 &"
        */
        char cmd[512];
        snprintf(cmd, sizeof(cmd),
            "ssh -o StrictHostKeyChecking=no %s@%s "
            "\"cd %s && nohup ./lab04 %d %d 1 > slave_%d.log 2>&1 &\"",
            slaves[i].user, slaves[i].ip,
            slaves[i].path,
            n, slaves[i].port, slaves[i].port);

        printf("[MASTER] Launching slave %d via SSH: %s\n", i, cmd);
        int ret = system(cmd);
        if (ret != 0)
            fprintf(stderr, "[MASTER] Warning: SSH launch for slave %d may have failed.\n", i);
    }

    /* Give slaves a moment to start listening */
    printf("[MASTER] Waiting 2s for slaves to start...\n");
    sleep(2);
}

/* ══════════════════════════════════════════════════════════════
 *  MASTER LOGIC
 * ══════════════════════════════════════════════════════════════ */
static void run_master(int n) {
    /* 1. Read config — collect slaves */
    PeerInfo peers[MAX_PEERS];
    int total = read_config(peers, MAX_PEERS);

    PeerInfo slaves[MAX_PEERS];
    int t = 0;
    for (int i = 0; i < total; i++)
        if (strcmp(peers[i].role, "slave") == 0)
            slaves[t++] = peers[i];

    if (t == 0) { fprintf(stderr, "No slaves in config.\n"); exit(1); }
    printf("[MASTER] Found %d slave(s) in config.\n", t);

    /* Optionally launch slaves via SSH if config includes user+path */
    ssh_launch_slaves(slaves, t, n);

    /* 2. Create random n×n matrix M */
    srand(time(NULL));
    double **M = alloc_matrix(n, n);
    gen_rand_matrix(M, n, n);
    print_matrix(M, n, n, "Full matrix M (master)");

    /* 3. Compute sub-matrix row range for each slave
          rows per slave = n/t  (last slave gets remainder) */
    int base_rows = n / t;

    /* 4. time_before */
    struct timespec t_before, t_after;
    clock_gettime(CLOCK_MONOTONIC, &t_before);

    /* 5. Distribute submatrices */
    int row_start = 0;
    for (int s = 0; s < t; s++) {
        int sub_rows = (s == t - 1) ? (n - row_start) : base_rows;

        /* Connect to slave */
        int fd = connect_to(slaves[s].ip, slaves[s].port);
        printf("[MASTER] Connected to slave %d at %s:%d\n",
               s, slaves[s].ip, slaves[s].port);

        /* Send dimensions first: sub_rows, n */
        int32_t hdr[2] = { htonl(sub_rows), htonl(n) };
        send_all(fd, hdr, sizeof(hdr));

        /* Send sub-matrix row by row */
        for (int r = row_start; r < row_start + sub_rows; r++)
            send_all(fd, M[r], n * sizeof(double));

        /* Wait for ACK */
        char ack_buf[ACK_LEN + 1] = {0};
        recv_all(fd, ack_buf, ACK_LEN);
        printf("[MASTER] Received ACK from slave %d: \"%s\"\n", s, ack_buf);

        close(fd);
        row_start += sub_rows;
    }

    /* 6. time_after */
    clock_gettime(CLOCK_MONOTONIC, &t_after);

    double elapsed = (t_after.tv_sec  - t_before.tv_sec) +
                     (t_after.tv_nsec - t_before.tv_nsec) * 1e-9;
    printf("\n[MASTER] time_elapsed = %.6f seconds\n", elapsed);

    free_matrix(M, n);
}

/* ══════════════════════════════════════════════════════════════
 *  SLAVE LOGIC
 * ══════════════════════════════════════════════════════════════ */
static void run_slave(int port) {
    /* Read config to find master IP (for display/verification) */
    PeerInfo peers[MAX_PEERS];
    int total = read_config(peers, MAX_PEERS);
    char master_ip[64] = "unknown";
    for (int i = 0; i < total; i++)
        if (strcmp(peers[i].role, "master") == 0) {
            strncpy(master_ip, peers[i].ip, sizeof(master_ip) - 1);
            break;
        }
    printf("[SLAVE  port=%d] Master IP from config: %s\n", port, master_ip);

    /* Listen on our port */
    int sfd = make_server_socket(port);
    printf("[SLAVE  port=%d] Listening...\n", port);

    /* Accept connection from master */
    struct sockaddr_in cli_addr;
    socklen_t cli_len = sizeof(cli_addr);
    int cfd = accept(sfd, (struct sockaddr *)&cli_addr, &cli_len);
    if (cfd < 0) { perror("accept"); exit(1); }

    char peer_ip[INET_ADDRSTRLEN];
    inet_ntop(AF_INET, &cli_addr.sin_addr, peer_ip, sizeof(peer_ip));
    printf("[SLAVE  port=%d] Master connected from %s\n", port, peer_ip);

    /* time_before */
    struct timespec t_before, t_after;
    clock_gettime(CLOCK_MONOTONIC, &t_before);

    /* Receive header: sub_rows, cols */
    int32_t hdr[2];
    recv_all(cfd, hdr, sizeof(hdr));
    int sub_rows = ntohl(hdr[0]);
    int cols     = ntohl(hdr[1]);
    printf("[SLAVE  port=%d] Expecting submatrix %d x %d\n", port, sub_rows, cols);

    /* Receive submatrix */
    double **sub = alloc_matrix(sub_rows, cols);
    for (int r = 0; r < sub_rows; r++)
        recv_all(cfd, sub[r], cols * sizeof(double));

    /* Send ACK */
    send_all(cfd, ACK_MSG, ACK_LEN);

    /* time_after */
    clock_gettime(CLOCK_MONOTONIC, &t_after);

    double elapsed = (t_after.tv_sec  - t_before.tv_sec) +
                     (t_after.tv_nsec - t_before.tv_nsec) * 1e-9;

    /* Display received submatrix for verification */
    char label[64];
    snprintf(label, sizeof(label), "Submatrix received (slave port=%d)", port);
    print_matrix(sub, sub_rows, cols, label);
    printf("\n[SLAVE  port=%d] time_elapsed = %.6f seconds\n", port, elapsed);

    free_matrix(sub, sub_rows);
    close(cfd);
    close(sfd);
}

/* ──────────────────────────── core affinity ──────────────────────────── */
#include <sched.h>

static void pin_process_to_core(int core_id) {
    int num_cores = (int)sysconf(_SC_NPROCESSORS_ONLN);
    core_id = core_id % num_cores;
    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(core_id, &cpuset);
    if (sched_setaffinity(0, sizeof(cpu_set_t), &cpuset) != 0)
        perror("sched_setaffinity");
}

/* ──────────────────────────── main ──────────────────────────── */
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

    if (role == 0)
        run_master(n);
    else {
        /* Pin slave process to a core derived from its port number */
        pin_process_to_core(port);
        run_slave(port);
    }

    return 0;
}