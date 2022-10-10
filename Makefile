all: main.c com_github_coderodde_util_radixsorts.h com_github_coderodde_util_radixsorts.c
	gcc -ansi -pedantic -Wall -Werror -Wno-error=unused-command-line-argument -lpthread -o radixsort.demo -O3 -fmax-errors=1 main.c com_github_coderodde_util_radixsorts.c

