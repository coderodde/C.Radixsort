#include "com_github_coderodde_util_radixsorts.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifdef _WIN32
#include <windows.h>
#elif defined (__unix__) || (defined (__APPLE__) && defined (__MACH__))
#include <sys/time.h>
#else
#error "Unsupported platform."
#endif

static unsigned* get_random_array(size_t length) {
    size_t i;
    unsigned* array = malloc(length * sizeof(unsigned));
    
    srand(time(NULL));

    for (i = 0; i < length; i++) {
        array[i] = rand() | rand() << 15 | rand() << 30;
    }

    return array;
}

static int is_sorted(unsigned* array, size_t length) {
    size_t i;

    for (i = 0; i < length - 1; i++) {
        if (array[i] > array[i + 1]) {
            return 0;
        }
    }

    return 1;
}

static size_t millis() {
#ifdef _WIN32
    return (size_t) GetTickCount();
#elif defined (__unix__) || (defined (__APPLE__) && defined (__MACH__))
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (size_t)(tv.tv_sec * 1000 + tv.tv_usec / 1000);
#endif
}

int ucmp(const void* pa, const void* pb) {
    unsigned ua = *(unsigned*) pa;
    unsigned ub = *(unsigned*) pb;
    
    if (ua < ub) {
        return -1;
    }

    if (ua > ub) {
        return 1;
    }

    return 0;
}

static size_t LENGTH = 10 * 1000 * 1000;

int main() {
    unsigned* array1;
    unsigned* array2;
    unsigned* array3;
    unsigned* array4;
    size_t t1;
    size_t t2;
    int cmp1;
    int cmp2;
    int cmp3;

    printf("Number of sorting threads: %ld\n", get_number_of_cpus());
    printf("Number of keys to sort: %ld\n", LENGTH);

    t1 = millis();
    
    array1 = get_random_array(LENGTH);
    array2 = malloc(LENGTH * sizeof(unsigned));
    array3 = malloc(LENGTH * sizeof(unsigned));
    array4 = malloc(LENGTH * sizeof(unsigned));

    memcpy(array2, array1, LENGTH * sizeof(unsigned));
    memcpy(array3, array1, LENGTH * sizeof(unsigned));
    memcpy(array4, array1, LENGTH * sizeof(unsigned));

    t2 = millis();
    
    printf("Created the arrays in %ld milliseconds.\n", t2 - t1);

    t1 = millis();
    qsort(array1, LENGTH, sizeof(unsigned), ucmp);
    t2 = millis();

    printf("qsort in %ld milliseconds.\n", t2 - t1);

    t1 = millis();
    bitwise_radix_sort(array2, LENGTH);
    t2 = millis();

    printf("bitwise_radix_sort in %ld milliseconds.\n", t2 - t1);

    t1 = millis();
    radix_sort(array3, LENGTH);
    t2 = millis();

    printf("radix_sort in %ld milliseconds.\n", t2 - t1);

    t1 = millis();
    parallel_radix_sort(array4, LENGTH);
    t2 = millis();

    printf("parallel_radix_sort in %ld milliseconds.\n", t2 - t1);

    cmp1 = memcmp(array1, array2, sizeof(unsigned) * LENGTH);
    cmp2 = memcmp(array1, array3, sizeof(unsigned) * LENGTH);
    cmp3 = memcmp(array1, array4, sizeof(unsigned) * LENGTH);

    puts("");
    
    printf("Algorithms agree: %d\n", (cmp1 == 0 && cmp2 == 0 && cmp3 == 0));

    puts("");

    printf("array1 is sorted: %d\n", is_sorted(array1, LENGTH));
    printf("array2 is sorted: %d\n", is_sorted(array2, LENGTH));
    printf("array3 is sorted: %d\n", is_sorted(array3, LENGTH));
    printf("array4 is sorted: %d\n", is_sorted(array4, LENGTH));

    free(array1);
    free(array2);
    free(array3);
    free(array4);
}

