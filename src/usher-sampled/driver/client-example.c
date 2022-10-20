#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <stdio.h>
#include <sys/un.h>
//Connect to socket and return a stdlib file handle corresponding to the connection
static FILE*  make_f(char* path){
    int sock_fd=socket(AF_UNIX, SOCK_STREAM, 0);
    if (sock_fd<0) {
        perror("cannot create socket");
        exit(EXIT_FAILURE);
    }
    struct sockaddr_un addr;
    addr.sun_family=AF_UNIX;
    strncpy(addr.sun_path, path, 108);
    int ret=connect(sock_fd, (struct sockaddr *)&addr, sizeof(addr));
    if (ret<0) {
        perror("cannot connect");
        exit(EXIT_FAILURE);
    }

    FILE* ret_f=fdopen(sock_fd, "a+");
    if (!ret_f) {
        perror("cannot create stdio file handle");
        exit(EXIT_FAILURE);
    }
    return ret_f;
}

//send arguments over the socket
static void send_args(FILE* f, int argc, char** argv){
    for (int idx=1; idx<argc; idx++) {
        fprintf(f, "%s\n",argv[idx]);
    }
    fputc('\n',f);
}

#define ArraySize(arr) ( sizeof(arr) / sizeof (*arr) )

int main (int argc, char** argv){
    if (argc != 5) {
        fprintf(stderr, "usage: client-example socket_file samples.vcf tree.pb out_dir\n");
        exit(1);
    }
    FILE* fh=make_f(argv[1]);
    char *cmd[] = { "ignored", "-v", argv[2], "-i", argv[3], "-d", argv[4],
                "-k", "50", "-K", "50", "-u","--no-ignore-prefix","user_"};
    send_args(fh, ArraySize(cmd), cmd);
    char* line=NULL;
    size_t capacity=0;
    while (getline(&line, &capacity, fh)>0) {
        if (line[0]==4) {
            break;
        }
        puts(line);
    }
    free(line);
}
