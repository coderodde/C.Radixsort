#ifndef COM_GITHUB_CODERODDE_UTIL_RADIXSORTS_H
#define COM_GITHUB_CODERODDE_UTIL_RADIXSORTS_H

#include <stdlib.h>

void bitwise_radix_sort(unsigned* data, size_t length);
void radix_sort(unsigned* data, size_t length);
void parallel_radix_sort(unsigned* data, size_t length);
size_t get_number_of_cpus();

#endif /* End of COM_GITHUB_CODERODDE_UTIL_RADIXSORTS_H. */
