#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#ifdef _WIN32
#include <windows.h>
#elif defined (__unix__) || (defined (__APPLE__) && defined (__MACH__))
#include <limits.h>
#include <pthread.h>
#include <sys/time.h>
#include <unistd.h>
#else
#error "Unsupported platform."
#endif

#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define BUCKETS 256

static const size_t BITS_PER_BUCKET = 8;
static const size_t BUCKET_MASK = 0xff;
static const size_t MERGESORT_THRESHOLD = 4096;
static const size_t INSERTION_SORT_THRESHOLD = 16;
static const size_t THREAD_THRESHOLD = 65536;

/******************************************************************************
* Array list data structure.                                                  *
******************************************************************************/
typedef struct {
    void** data;
    size_t size;
} array_t;

static void array_t_init(array_t* array, size_t capacity) {
    array->size = 0;
    array->data = malloc(capacity * sizeof(void*));
}

static void array_t_add(array_t* array, void* datum) {
    array->data[array->size++] = datum;
}

static void* array_t_get(array_t* array, size_t index) {
    return array->data[index];
}

static void array_t_shuffle(array_t* array) {
    size_t i;
    size_t j;
    void* temp;

    srand(time(NULL));

    for (i = 0; i != array->size - 1; ++i) {
        j = i + rand() % (array->size - i);
        temp = array->data[i];
        array->data[i] = array->data[j];
        array->data[j] = temp;
    }
}

static size_t array_t_size(array_t* array) {
    return array->size;
}

static void array_t_destruct(array_t* array) {
    free(array->data);
}

/******************************************************************************
* Thread-specific data structures.                                            *
******************************************************************************/
typedef struct {
    size_t local_bucket_size_map[BUCKETS];
    unsigned* source;
    size_t recursion_deph;
    size_t from_index;
    size_t to_index;
} bucket_size_counter_thread_data;

typedef struct {
    unsigned* source;
    unsigned* target;
    size_t* start_index_map;
    size_t* processed_map;
    size_t recursion_depth;
    size_t from_index;
    size_t to_index;
} bucket_inserter_thread_data;

typedef struct {
    unsigned* source;
    unsigned* target;
    size_t threads;
    size_t recursion_depth;
    size_t from_index;
    size_t to_index;
} task;
/******************************************************************************
* End of data structures.                                                     *
******************************************************************************/

size_t get_number_of_cpus() {
#ifdef _WIN32
    SYSTEM_INFO system_info;
    GetSystemInfo(&system_info);
    return (size_t) system_info.dwNumberOfProcessors;
#elif defined (__unix__) || (defined (__APPLE__) && defined (__MACH__))
    return (size_t) sysconf(_SC_NPROCESSORS_ONLN);
#endif
}

static size_t get_bucket_index(unsigned datum, size_t recursion_depth) {
    size_t bit_shift = CHAR_BIT * sizeof(unsigned) -
        (recursion_depth + 1) * BITS_PER_BUCKET;

    return (((size_t) datum) >> bit_shift) & BUCKET_MASK;
}

static void parallel_radix_sort_impl(unsigned* source,
                                     unsigned* target,
                                     size_t threads,
                                     size_t recursion_depth,
                                     size_t from_index,
                                     size_t to_index);

static void radix_sort_impl_no_threads(unsigned* source,
                                       unsigned* target,
                                       size_t recursion_depth,
                                       size_t from_index,
                                       size_t to_index);

static void process_bucket_size_counter_thread(
    bucket_size_counter_thread_data* data) {

    size_t i;

    memset(data->local_bucket_size_map, 0, BUCKETS * sizeof(size_t));

    for (i = data->from_index; i != data->to_index; ++i) {
        data->local_bucket_size_map[
            get_bucket_index(
                data->source[i],
                data->recursion_deph)]++;
    }
}

static void process_bucket_inserter_thread(bucket_inserter_thread_data* data) {
    size_t bucket_index;
    size_t i;
    unsigned datum;

    for (i = data->from_index; i != data->to_index; ++i) {
        datum = data->source[i];
        bucket_index = get_bucket_index(datum, data->recursion_depth);
        data->target[data->start_index_map[bucket_index] + 
                     data->processed_map[bucket_index]++] = datum;
    }
}

static void process_sorter_thread(array_t* data) {
    size_t i;
    task* t;

    for (i = 0; i != array_t_size(data); ++i) {
        t = array_t_get(data, i);

        if (t->threads > 1) {
            parallel_radix_sort_impl(t->source,
                                     t->target,
                                     t->threads,
                                     t->recursion_depth,
                                     t->from_index,
                                     t->to_index);
        } else {
            radix_sort_impl_no_threads(t->source,
                                       t->target,
                                       t->recursion_depth,
                                       t->from_index,
                                       t->to_index);
        }
    }
}

#ifdef _WIN32

static DWORD WINAPI count_bucket_sizes_thread_func_win(LPVOID parameter) {
    bucket_size_counter_thread_data* thread_data =
        (bucket_size_counter_thread_data*) parameter;

    process_bucket_size_counter_thread(thread_data);
    return 0;
}

#elif defined (__unix__) || (defined (__APPLE__) && defined (__MACH__))

static void* count_bucket_sizes_thread_func_pthreads(void* parameter) {
    bucket_size_counter_thread_data* thread_data =
        (bucket_size_counter_thread_data*) parameter;

    process_bucket_size_counter_thread(thread_data);
    return NULL;
}

#endif

#ifdef _WIN32

static DWORD WINAPI insert_to_buckets_thread_func_win(LPVOID parameter) {
    bucket_inserter_thread_data* thread_data =
        (bucket_inserter_thread_data*) parameter;

    process_bucket_inserter_thread(thread_data);
    return 0;
}

#elif defined (__unix__) || (defined (__APPLE__) && defined (__MACH__))

static void* insert_to_buckets_thread_func_pthreads(void* parameter) {
    bucket_inserter_thread_data* thread_data =
        (bucket_inserter_thread_data*) parameter;

    process_bucket_inserter_thread(thread_data);
    return NULL;
}

#endif

#ifdef _WIN32

static DWORD WINAPI sort_buckets_thread_func_win(LPVOID parameter) {
    array_t* thread_data = (array_t*) parameter;
    process_sorter_thread(thread_data);
    return 0;
}

#elif defined(__unix__) || (defined(__APPLE__) && defined(__MACH__))

static void* sort_buckets_thread_func_pthreads(void* parameter) {
    array_t* thread_data = (array_t*) parameter;
    process_sorter_thread(thread_data);
    return NULL;
}

#endif

static void insertion_sort(unsigned* data, size_t length) {
    size_t i;
    signed long j;
    unsigned datum;

    for (i = 1; i != length; ++i) {
        datum = data[i];
        j = i - 1;

        while (j >= 0 && data[j] > datum) {
            data[j + 1] = data[j];
            --j;
        }

        data[j + 1] = datum;
    }
}

static void merge(unsigned* source,
                  unsigned* target,
                  size_t left_index,
                  size_t left_bound,
                  size_t right_bound) {

    size_t right_index = left_bound;
    size_t target_index = left_index;

    while (left_index < left_bound && right_index < right_bound) {
        target[target_index++] = source[left_index] < source[right_index] ?
            source[left_index++] :
            source[right_index++];
    }

    memcpy(target + target_index,
           source + left_index,
           sizeof(unsigned) * (left_bound - left_index));

    memcpy(target + target_index,
           source + right_index,
           sizeof(unsigned) * (right_bound - right_index));
}

static void radix_sort_mergesort(unsigned* source,
                                 unsigned* target,
                                 size_t recursion_depth,
                                 size_t from_index,
                                 size_t to_index) {

    unsigned* s;
    unsigned* t;
    unsigned* temp;

    int even;

    size_t i;
    size_t left_bound;
    size_t left_index;
    size_t offset = from_index;
    size_t passes = 0;
    size_t range_length;
    size_t right_bound;
    size_t runs;
    size_t run_index;
    size_t run_width;

    range_length = to_index - from_index;
    s = source;
    t = target;

    runs = range_length / INSERTION_SORT_THRESHOLD;

    for (i = 0; i != runs; ++i) {
        insertion_sort(source + offset, INSERTION_SORT_THRESHOLD);
        offset += INSERTION_SORT_THRESHOLD;
    }

    if (range_length % INSERTION_SORT_THRESHOLD != 0) {
        /* Sort the rightmost run that is smaller than */
	/* INSERTION_SORT_THRESHOLD.                   */
        insertion_sort(source + offset, to_index - offset);
        runs++;
    }

    run_width = INSERTION_SORT_THRESHOLD;

    while (runs != 1) {
        passes++;
        run_index = 0;

        for (; run_index < runs - 1; run_index += 2) {
            left_index = from_index + run_index * run_width;
            left_bound = left_index + run_width;
            right_bound = MIN(left_bound + run_width, to_index);

            merge(s,
                  t,
                  left_index,
                  left_bound,
                  right_bound);
        }

        if (run_index < runs) {
            memcpy(t + from_index + run_index * run_width,
                   s + from_index + run_index * run_width,
                   sizeof(unsigned) * (range_length - run_index * run_width));
        }

        runs = (runs / 2) + (runs % 2 == 0 ? 0 : 1);
        temp = s;
        s = t;
        t = temp;
        run_width *= 2;
    }

    even = (passes % 2 == 0) ? 1 : 0;

    if (recursion_depth % 2 == 1) {
        if (even == 1) {
            memcpy(target + from_index, /* Destination */
                   source + from_index,    /* Source      */
                   sizeof(unsigned) * (to_index - from_index));
        }
    }
    else {
        /* Here, recursion_depth % 2 == 0 holds: */
        if (even == 0) {
            memcpy(source + from_index, /* Destination */
                   target + from_index,    /* Source      */
                   sizeof(unsigned) * (to_index - from_index));
        }
    }
}

static void radix_sort_impl_no_threads(unsigned* source,
                                       unsigned* target,
                                       size_t recursion_depth,
                                       size_t from_index,
                                       size_t to_index) {
    size_t bucket_key;
    size_t i;
    size_t bucket_size_map[BUCKETS];
    size_t processed_map[BUCKETS];
    size_t range_length;
    size_t start_index_map[BUCKETS];
    unsigned datum;

    range_length = to_index - from_index;

    if (range_length <= MERGESORT_THRESHOLD) {
        radix_sort_mergesort(source,
                             target,
                             recursion_depth,
                             from_index,
                             to_index);
        return;
    }

    memset(bucket_size_map, 0, BUCKETS * sizeof(size_t));
    memset(start_index_map, 0, BUCKETS * sizeof(size_t));
    memset(processed_map, 0, BUCKETS * sizeof(size_t));

    /* Compute the size of each bucket: */
    for (i = from_index; i != to_index; i++) {
        bucket_size_map[get_bucket_index(source[i], recursion_depth)]++;
    }

    /* Initialize thee start index map: */
    start_index_map[0] = from_index;

    for (i = 1; i != BUCKETS; ++i) {
        start_index_map[i] = start_index_map[i - 1]
                           + bucket_size_map[i - 1];
    }

    /* Insert the data from 'source' into their */
    /* respective position in 'target': */
    for (i = from_index; i != to_index; ++i) {
        datum = source[i];
        bucket_key = get_bucket_index(datum, recursion_depth);
        target[start_index_map[bucket_key] + 
               processed_map[bucket_key]++] = datum;
    }

    if (recursion_depth == sizeof(unsigned) - 1) {
        memcpy(source + from_index, /* Destination */
               target + from_index, /* Source      */
               sizeof(unsigned) * (to_index - from_index));

        /* There is nowhere to recur, return. */
        return;
    }

    for (i = 0; i != BUCKETS; ++i) {
        if (bucket_size_map[i] != 0) {
            radix_sort_impl_no_threads(target,
                                       source,
                                       recursion_depth + 1,
                                       start_index_map[i],
                                       start_index_map[i] + 
                                           bucket_size_map[i]);
        }
    }
}

void radix_sort(unsigned* data, size_t length) {
    unsigned* buffer;

    if (length < 2) {
        return;
    }

    buffer = malloc(sizeof(unsigned) * length);
    radix_sort_impl_no_threads(data, buffer, 0, 0, length);
    free(buffer);
}

static void parallel_radix_sort_impl(unsigned* source,
                                     unsigned* target,
                                     size_t threads,
                                     size_t recursion_depth,
                                     size_t from_index,
                                     size_t to_index) {
    size_t bucket_key;
    size_t f;
    size_t i;
    size_t idx;
    size_t j;
    size_t list_index;
    size_t number_of_nonempty_buckets;
    size_t optimal_subrange_length;
    size_t packed;
    size_t range_length;
    size_t spawn_degree;
    size_t start;
    size_t subrange_length;
    size_t sz;
    size_t sz2;
    size_t tmp;
    size_t* partial_bucket_size_map;
    size_t* thread_count_map;
    size_t bucket_size_map[BUCKETS] = { 0 };
    size_t start_index_map[BUCKETS];
    size_t** processed_map;
    bucket_size_counter_thread_data* bucket_size_counter_threads_data;
    bucket_inserter_thread_data* bucket_inserter_threads_data;
    array_t array_of_task_arrays;
    array_t bucket_index_list_array;
    array_t non_empty_bucket_indices;
    array_t* arr2;
    task* t;

#ifdef _WIN32
    HANDLE windows_thread_handle;
    HANDLE* win_thread_handles;
#elif defined (__unix__) || (defined (__APPLE__) && defined (__MACH__))
    pthread_t pthread_handle;
    pthread_t* unix_thread_ids;
#endif

    range_length = to_index - from_index;

    if (range_length <= MERGESORT_THRESHOLD) {
        radix_sort_mergesort(source,
                             target,
                             recursion_depth,
                             from_index,
                             to_index);
        return;
    }

    if (threads < 2) {
        radix_sort_impl_no_threads(source,
                                   target,
                                   recursion_depth,
                                   from_index,
                                   to_index);
        return;
    }

    bucket_size_counter_threads_data =
        malloc(threads * sizeof(*bucket_size_counter_threads_data));

    start = from_index;
    subrange_length = range_length / threads;

#ifdef _WIN32
    win_thread_handles = malloc(threads * sizeof(HANDLE));
#elif defined(__unix__) || (defined(__APPLE__) && defined(__MACH__))
    unix_thread_ids = malloc(threads * sizeof(pthread_t));
#endif

    for (i = 0; i != threads - 1; ++i) {
        bucket_size_counter_threads_data[i].source = source;
        bucket_size_counter_threads_data[i].recursion_deph = recursion_depth;
        bucket_size_counter_threads_data[i].from_index = start;
        bucket_size_counter_threads_data[i].to_index = start += subrange_length;

        memset(&(bucket_size_counter_threads_data[i]
            .local_bucket_size_map),
            0,
            BUCKETS * sizeof(size_t));

#ifdef _WIN32

        win_thread_handles[i] =
            CreateThread(NULL,
                         0,
                         count_bucket_sizes_thread_func_win,
                         &bucket_size_counter_threads_data[i],
                         0,
                         NULL);

#elif defined (__unix__) || (defined (__APPLE__) && defined (__MACH__))

        pthread_create(&pthread_handle,
                       NULL,
                       count_bucket_sizes_thread_func_pthreads,
                       &bucket_size_counter_threads_data[i]);

        unix_thread_ids[i] = pthread_handle;

#endif
    }

    /* Process the rightmost bucket in THIS thread. No need to spawn */
    /* any more.                                                     */
    bucket_size_counter_threads_data[threads - 1].source = source;
    bucket_size_counter_threads_data[threads - 1].recursion_deph =
        recursion_depth;

    bucket_size_counter_threads_data[threads - 1].from_index = start;
    bucket_size_counter_threads_data[threads - 1].to_index = to_index;

    memset(&(bucket_size_counter_threads_data[threads - 1]
        .local_bucket_size_map),
        0,
        BUCKETS * sizeof(size_t));

    /* Run the rightmost thread routine in THIS thread. */
    /* No need to span another thread:                  */
    process_bucket_size_counter_thread(
        &bucket_size_counter_threads_data[threads - 1]);

    /* Wait for all the bucket counters: */
    for (i = 0; i != threads - 1; ++i) {
#ifdef _WIN32
        WaitForSingleObject(win_thread_handles[i], INFINITE);
#elif defined (__unix__) || (defined (__APPLE__) && defined (__MACH__))
        pthread_join(unix_thread_ids[i], NULL);
#endif
    }

    /* Build the global bucket size map for the entire sorting range: */
    for (i = 0; i != threads; ++i) {
        for (j = 0; j != BUCKETS; ++j) {
            bucket_size_map[j] +=
                bucket_size_counter_threads_data[i].local_bucket_size_map[j];
        }
    }

    number_of_nonempty_buckets = 0;

    for (i = 0; i != BUCKETS; ++i) {
        if (bucket_size_map[i] != 0) {
            number_of_nonempty_buckets++;
        }
    }

    spawn_degree = MIN(number_of_nonempty_buckets, threads);

    /* Prepare the starting indices of each bucket: */
    start_index_map[0] = from_index;

    for (i = 1; i != BUCKETS; ++i) {
        start_index_map[i] = start_index_map[i - 1]
                           + bucket_size_map[i - 1];
    }

    processed_map = malloc(spawn_degree * sizeof(size_t*));

    for (i = 0; i != spawn_degree; ++i) {
        processed_map[i] = calloc(BUCKETS, sizeof(size_t));
    }

    /* Make the preprocessed_map of each thread independent of the other. */
    for (i = 1; i != spawn_degree; ++i) {
        partial_bucket_size_map =
            (bucket_size_counter_threads_data[i - 1].local_bucket_size_map);

        for (j = 0; j != BUCKETS; ++j) {
            processed_map[i][j] = processed_map[i - 1][j]
                                + partial_bucket_size_map[j];
        }
    }

    start = from_index;

    bucket_inserter_threads_data =
        malloc(spawn_degree * sizeof(bucket_inserter_thread_data));

    for (i = 0; i != spawn_degree - 1; ++i) {
        bucket_inserter_threads_data[i].start_index_map = start_index_map;
        bucket_inserter_threads_data[i].processed_map = processed_map[i];
        bucket_inserter_threads_data[i].source = source;
        bucket_inserter_threads_data[i].target = target;
        bucket_inserter_threads_data[i].recursion_depth = recursion_depth;
        bucket_inserter_threads_data[i].from_index = start;
        bucket_inserter_threads_data[i].to_index = start += subrange_length;

#ifdef _WIN32

        win_thread_handles[i] =
            CreateThread(NULL,
                         0,
                         insert_to_buckets_thread_func_win,
                         &bucket_inserter_threads_data[i],
                         0,
                         NULL);

#elif defined (__unix__) || (defined (__APPLE__) && defined (__MACH__))

        pthread_create(&pthread_handle,
                       NULL,
                       insert_to_buckets_thread_func_pthreads,
                       &bucket_inserter_threads_data[i]);

        unix_thread_ids[i] = pthread_handle;

#endif
    }

    /* Process the rightmost bucket in THIS thread. No need to spawn */
    /* any more.                                                     */
    bucket_inserter_threads_data[spawn_degree - 1].start_index_map =
        start_index_map;

    bucket_inserter_threads_data[spawn_degree - 1].processed_map =
        processed_map[spawn_degree - 1];

    bucket_inserter_threads_data[spawn_degree - 1].source = source;
    bucket_inserter_threads_data[spawn_degree - 1].target = target;
    bucket_inserter_threads_data[spawn_degree - 1].recursion_depth =
        recursion_depth;

    bucket_inserter_threads_data[spawn_degree - 1].from_index = start;
    bucket_inserter_threads_data[spawn_degree - 1].to_index = to_index;

    process_bucket_inserter_thread(
        &bucket_inserter_threads_data[spawn_degree - 1]);

    /* Wait for all the bucket inserters: */
    for (i = 0; i != spawn_degree - 1; ++i) {
#ifdef _WIN32
        WaitForSingleObject(win_thread_handles[i], INFINITE);
#elif defined (__unix__) || (defined (__APPLE__) && defined (__MACH__))
        pthread_join(unix_thread_ids[i], NULL);
#endif
    }

    free(bucket_size_counter_threads_data);
    free(bucket_inserter_threads_data);

    for (i = 0; i != spawn_degree; ++i) {
        free(processed_map[i]);
    }

    free(processed_map);

    if (recursion_depth == sizeof(unsigned) - 1) {
        /* Nowhere to recur. */
        return;
    }

    array_t_init(&bucket_index_list_array, spawn_degree);

    for (i = 0; i != spawn_degree; ++i) {
        array_t* bucket_key_array = malloc(sizeof(array_t));
        array_t_init(bucket_key_array, number_of_nonempty_buckets);
        array_t_add(&bucket_index_list_array, bucket_key_array);
    }

    thread_count_map = calloc(spawn_degree, sizeof(size_t));

    for (i = 0; i != spawn_degree; ++i) {
        thread_count_map[i] = threads / spawn_degree;
    }

    for (i = 0; i != threads % spawn_degree; ++i) {
        ++thread_count_map[i];
    }

    array_t_init(&non_empty_bucket_indices, number_of_nonempty_buckets);

    for (bucket_key = 0; bucket_key != BUCKETS; ++bucket_key) {
        if (bucket_size_map[bucket_key] != 0) {
            array_t_add(&non_empty_bucket_indices, (void*) bucket_key);
        }
    }

    array_t_shuffle(&non_empty_bucket_indices);

    f = 0;
    j = 0;
    list_index = 0;
    optimal_subrange_length = range_length / spawn_degree;
    packed = 0;
    sz = array_t_size(&non_empty_bucket_indices);

    while (j != sz) {
        size_t bucket_key =
            (size_t) array_t_get(&non_empty_bucket_indices, j++);

        tmp = bucket_size_map[bucket_key];
        packed += tmp;

        if (packed >= optimal_subrange_length ||
            j == array_t_size(&non_empty_bucket_indices)) {

            packed = 0;

            for (i = f; i != j; ++i) {
                size_t bucket_key =
                    (size_t) array_t_get(&non_empty_bucket_indices, i);

                array_t* arr = array_t_get(&bucket_index_list_array,
                                           list_index);

                array_t_add(arr, (void*) bucket_key);
            }

            ++list_index;
            f = j;
        }
    }

    array_t_init(&array_of_task_arrays, spawn_degree);

    for (i = 0; i != spawn_degree; ++i) {
        array_t* task_array = malloc(sizeof(array_t));
        array_t_init(task_array, BUCKETS);
        arr2 = (array_t*) array_t_get(&bucket_index_list_array, i);
        sz = array_t_size(arr2);

        for (idx = 0; idx != sz; ++idx) {
            bucket_key = (size_t) array_t_get(arr2, idx);

            t = malloc(sizeof(task));

            t->source = target;
            t->target = source;
            t->threads = thread_count_map[i];
            t->recursion_depth = recursion_depth + 1;
            t->from_index = start_index_map[bucket_key];
            t->to_index = start_index_map[bucket_key]
                        + bucket_size_map[bucket_key];

            array_t_add(task_array, t);
        }

        array_t_add(&array_of_task_arrays, task_array);
    }

    for (i = 0; i != spawn_degree - 1; ++i) {
        array_t* task_array = array_t_get(&array_of_task_arrays, i);

#ifdef _WIN32

        win_thread_handles[i] =
            CreateThread(NULL,
                         0,
                         sort_buckets_thread_func_win,
                         task_array,
                         0,
                         NULL);

#elif defined(__unix__) || (defined(__APPLE__) && defined(__MACH__))

        pthread_create(&pthread_handle,
                       NULL,
                       sort_buckets_thread_func_pthreads,
                       task_array);

        unix_thread_ids[i] = pthread_handle;

#endif
    }

    /* Sort the rightmost thread in THIS thread. */
    /* No need to spawn one more thread.         */
    process_sorter_thread(
        array_t_get(
            &array_of_task_arrays,
            spawn_degree - 1));

    for (i = 0; i != spawn_degree - 1; ++i) {

#ifdef _WIN32
        WaitForSingleObject(win_thread_handles[i], INFINITE);
#elif defined (__unix__) || (defined (__APPLE__) && defined (__MACH__))
        pthread_join(unix_thread_ids[i], NULL);
#endif

    }

#ifdef _WIN32
    free(win_thread_handles);
#elif defined (__unix__) || (defined (__APPLE__) && defined (__MACH__))
    free(unix_thread_ids);
#endif

    sz = array_t_size(&array_of_task_arrays);

    for (i = 0; i != sz; ++i) {
        array_t* task_array = array_t_get(&array_of_task_arrays, i);

        sz2 = array_t_size(task_array);

        for (j = 0; j != sz2; ++j) {
            free(array_t_get(task_array, j));
        }

        array_t_destruct(task_array);
        free(task_array);
    }

    sz = array_t_size(&bucket_index_list_array);

    for (i = 0; i != sz; ++i) {
        array_t* array = array_t_get(&bucket_index_list_array, i);
        array_t_destruct(array);
        free(array);
    }

    free(thread_count_map);
    array_t_destruct(&array_of_task_arrays);
    array_t_destruct(&bucket_index_list_array);
    array_t_destruct(&non_empty_bucket_indices);
}

static void bitwise_radix_sort_impl(unsigned* data,
                                    size_t bucket_length,
                                    size_t bit_index) {

    size_t size_of_left_bucket;
    size_t size_of_right_bucket;

    unsigned bit_is_on;
    unsigned datum;
    unsigned mask;
    unsigned temp;

    if (bucket_length < 2) {
        /* Trivially sorted. */
        return;
    }

    size_of_left_bucket = 0;
    size_of_right_bucket = 0;
    mask = 1U << bit_index;

    /* Bucketize the current range: */
    while (size_of_left_bucket + size_of_right_bucket < bucket_length) {
        datum = data[size_of_left_bucket];
        bit_is_on = datum & mask;

        if (bit_is_on) {
            /* Kick the datum to the right 1-bucket: */
            temp = data[bucket_length - size_of_right_bucket - 1];
            data[bucket_length - size_of_right_bucket - 1] = datum;
            data[size_of_left_bucket] = temp;
            size_of_right_bucket++;
        }
        else {
            /* Omit the datum: */
            size_of_left_bucket++;
        }
    }

    /* Any bits to proceed? */
    if (bit_index > 0) {
        /* Sort the 0-bucket of this recursion level: */
        bitwise_radix_sort_impl(data,
                                size_of_left_bucket,
                                bit_index - 1);

        /* Sort the 1-bucket of this recursion level: */
        bitwise_radix_sort_impl(data + size_of_left_bucket,
                                size_of_right_bucket,
                                bit_index - 1);
    }
}

void bitwise_radix_sort(unsigned* data, size_t length) {
    bitwise_radix_sort_impl(data,
                            length,
                            sizeof(unsigned) * CHAR_BIT - 1);
}

void parallel_radix_sort(unsigned* data, size_t length) {
    unsigned* buffer;
    size_t threads;

    if (length < 2) {
        return;
    }

    buffer = malloc(sizeof(unsigned) * length);
    threads = get_number_of_cpus();
    threads = MIN(threads, length / THREAD_THRESHOLD);
    parallel_radix_sort_impl(data, buffer, threads, 0, 0, length);
    free(buffer);
}

