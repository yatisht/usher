#include <cstdio>
#include <unistd.h>
void wait_debug() {
    volatile int cont=0;
    while (cont==0) {
        fprintf(stderr, "a");
        sleep(1);
    }
}