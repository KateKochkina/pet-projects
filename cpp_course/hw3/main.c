#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <errno.h>
#include <string.h>

#define TRUE 1
#define FALSE 0


char * delete_linebreak(const char *str) {
    if (!str) {
        return NULL;
    }

    char *new_str = malloc(strlen(str) + 1);
    if (!new_str) {
        return NULL;
    }

    char *linebreak = strrchr(str, '\n');
    if (!linebreak) {
        linebreak = strrchr(str, '\r');
    }

    if (!linebreak) {
        snprintf(new_str, strlen(str) + 1, "%s", str);
        return new_str;
    }

    snprintf(new_str, linebreak - str + 1, "%s", str);

    linebreak = strrchr(new_str, '\r');
    if (!linebreak) {
        return new_str;
    }

    char *new_new_str = malloc(linebreak - new_str + 1);
    if (!new_new_str) {
        free(new_str);
        return NULL;
    }

    snprintf(new_new_str, linebreak - new_str + 1, "%s", new_str);

    free(new_str);
    return new_new_str;
}


char * delete_spaces(const char *str) {
    if (!str) {
        return NULL;
    }

    int count_spaces = 0;
    for (size_t i = 0; i < strlen(str); ++i) {
        if (str[i] == ' ' || str[i] == '\t') {
            count_spaces += 1;
        } else {
            break;
        }
    }

    char *new_str = malloc(strlen(str) - count_spaces + 1);
    if (!new_str) {
        return NULL;
    }

    snprintf(new_str, strlen(str) - count_spaces + 1, "%s", str + count_spaces);
    return new_str;
}


char * str_proces(const char *str) {
    if (!str) {
        return NULL;
    }

    char *new_str = delete_spaces(str);
    if (!new_str) {
        return NULL;
    }

    char *new_new_str = delete_linebreak(new_str);
    if (!new_new_str) {
        free(new_str);
        return NULL;
    }

    free(new_str);
    return new_new_str;
}


char * get_boundary_if_quotes(const char *boundary) {
    if (!boundary) {
        return NULL;
    }

    char *copy_boundary = malloc(strlen(boundary) + 1);
    if (!copy_boundary) {
        return NULL;
    }

    snprintf(copy_boundary, strlen(boundary) + 1, "-%s", boundary + 1);
    char *quotes = strchr(copy_boundary, '"');

    char *new_boundary = malloc(quotes - copy_boundary + 2);
    if (!new_boundary) {
        free(copy_boundary);
        return NULL;
    }

    snprintf(new_boundary, quotes - copy_boundary + 2, "-%s", copy_boundary);

    free(copy_boundary);
    return new_boundary;
}


char * delete_boundary_spaces(const char *boundary) {
    if (!boundary) {
        return NULL;
    }

    char *boundary_copy = NULL;
    char *new_boundary = NULL;
    if (asprintf(&new_boundary, "%s", boundary) < 0) {
        return NULL;
    }

    while (new_boundary[strlen(new_boundary) - 1] == ' ') {
        if (asprintf(&boundary_copy, "%s", new_boundary) < 0) {
            free(new_boundary);
            return NULL;
        }

        free(new_boundary);
        new_boundary = malloc(strlen(boundary_copy));
        if (!new_boundary) {
            free(boundary_copy);
            return NULL;
        }

        snprintf(new_boundary, strlen(boundary_copy), "%s", boundary_copy);
    }

    free(boundary_copy);
    return new_boundary;
}


char * get_boundary_no_quotes(const char *boundary) {
    if (!boundary) {
        return NULL;
    }

    char *new_boundary = malloc(strlen(boundary) + 3);
    if (!new_boundary) {
        return NULL;
    }

    snprintf(new_boundary, strlen(boundary) + 3, "--%s", boundary);

    char *new_new_boundary = delete_boundary_spaces(new_boundary);
    if (!new_new_boundary) {
        free(new_boundary);
        return NULL;
    }

    free(new_boundary);
    return new_new_boundary;
}


char * delete_equal(const char *boundary) {
    if (!boundary) {
        return 0;
    }

    char *new_boundary = NULL;

    if (boundary[0] == '=') {
        char *non_boundary = malloc(strlen(boundary));
        if (!non_boundary) {
            return NULL;
        }

        snprintf(non_boundary, strlen(boundary), "%s", boundary + 1);

        new_boundary = str_proces(non_boundary);
        if (!new_boundary) {
            free(non_boundary);
            return NULL;
        }

        free(non_boundary);
        return new_boundary;
    }

    if (asprintf(&new_boundary, "%s", boundary) < 0) {
        return NULL;
    }
    return new_boundary;
}


char * boundary_proces(const char *boundary) {
    if (!boundary) {
        return NULL;
    }

    char *non_boundary = delete_equal(boundary);
    if (!non_boundary) {
        return NULL;
    }

    char *new_boundary = NULL;
    if (!strncmp(boundary, "\"", 1)) {
        new_boundary = get_boundary_if_quotes(non_boundary);
        if (!new_boundary) {
            free(non_boundary);
            return NULL;
        }
    } else {
        new_boundary = get_boundary_no_quotes(non_boundary);
        if (!new_boundary) {
            free(non_boundary);
            return NULL;
        }
    }

    free(non_boundary);
    return new_boundary;
}


char * get_boundary_str(const char *cont_type_str) {
    if (!cont_type_str) {
        return NULL;
    }

    const int BOUNDARY_LENGTH = 9;
    char *pointer = strcasestr(cont_type_str, "boundary");

    char *non_boundary = NULL;
    if (asprintf(&non_boundary, "%s", pointer + BOUNDARY_LENGTH) < 0) {
        return NULL;
    }

    char *boundary = boundary_proces(non_boundary);
    if (!boundary) {
        free(non_boundary);
        return NULL;
    }

    free(non_boundary);
    return boundary;
}


char * add_next_str(const char *str, FILE *f, char *next_str) {
    if (!str || !next_str) {
        return NULL;
    }

    const int STR_LENGTH = 1000;

    char *str_with_next = str_proces(str);
    if (!str_with_next) {
        return NULL;
    }

    char *str_copy = NULL;
    char *str_next_copy = NULL;
    char *non_str_next_copy = NULL;

    while (!strncmp(next_str, " ", 1) || !strncmp(next_str, "\t", 1)) {
        if (asprintf(&non_str_next_copy, "%s", next_str) < 0) {
            free(str_with_next);
            return NULL;
        }

        str_next_copy = str_proces(non_str_next_copy);
        if (!str_next_copy) {
            free(str_with_next);
            free(non_str_next_copy);
            return NULL;
        }
        free(non_str_next_copy);

        if (asprintf(&str_copy, "%s", str_with_next) < 0) {
            free(str_with_next);
            free(str_next_copy);
            return NULL;
        }

        free(str_with_next);
        if (asprintf(&str_with_next, "%s %s", str_copy, str_next_copy) < 0) {
            free(str_copy);
            free(str_next_copy);
            return NULL;
        }

        free(str_copy);
        free(str_next_copy);
        fgets(next_str, STR_LENGTH, f);
    }

    return str_with_next;
}


char * get_cont_type_str(FILE *f, const char *str_from_file, char *next_str) {
    if (!str_from_file || !next_str) {
        return NULL;
    }

    const int CONT_TYPE_LENGTH = 13;

    char *str;
    if (asprintf(&str, "%s", str_from_file + CONT_TYPE_LENGTH) < 0) {
        return NULL;
    }

    char *str_with_next = add_next_str(str, f, next_str);
    if (!str_with_next) {
        free(str);
        return NULL;
    }

    free(str);
    return str_with_next;
}


char * get_output_str(FILE *f, const char *str_from_file, char *next_str, const int CAPTION_LENGTH) {
    if (!str_from_file || !next_str) {
        return NULL;
    }

    char *str = NULL;
    if (asprintf(&str, "%s", str_from_file + CAPTION_LENGTH) < 0) {
        return NULL;
    }

    char *str_with_next = add_next_str(str, f, next_str);
    if (!str_with_next) {
        free(str);
        return NULL;
    }

    char *output_str;
    if (asprintf(&output_str, "%s", str_with_next) < 0) {
        free(str);
        free(str_with_next);
        return NULL;
    }

    free(str);
    free(str_with_next);
    return output_str;
}


int head_eml_proces(FILE *f, const char *str_from_file, char *next_str, char **boundary,
                    char **from_str, char **to_str, char **date_str, size_t *check_cont_type) {
    if (!str_from_file || !next_str) {
        return -1;
    }

    const int FROM_LEGTH = 5;
    const int TO_LEGTH = 3;
    const int DATE_LEGTH = 5;
    const int CONT_TYPE_LENGTH = 13;

    char *cont_type_str = NULL;

    if (!strncmp(str_from_file, "From:", FROM_LEGTH)) {
        *from_str = get_output_str(f, str_from_file, next_str, FROM_LEGTH);
        if (!*from_str) {
            return -1;
        }
    }

    if (!strncmp(str_from_file, "To:", TO_LEGTH)) {
        *to_str = get_output_str(f, str_from_file, next_str, TO_LEGTH);
        if (!*to_str) {
            return -1;
        }
    }

    if (!strncmp(str_from_file, "Date:", DATE_LEGTH)) {
        *date_str = get_output_str(f, str_from_file, next_str, DATE_LEGTH);
        if (!*date_str) {
            return -1;
        }
    }

    if (!strncmp(str_from_file, "Content-Type:", CONT_TYPE_LENGTH)) {
        cont_type_str = get_cont_type_str(f, str_from_file, next_str);
        if (!cont_type_str) {
            return -1;
        }

        if (strcasestr(cont_type_str, "multipart")) {
            *check_cont_type = TRUE;

            *boundary = get_boundary_str(cont_type_str);
            if (!*boundary) {
                free(cont_type_str);
                return -1;
            }
        }
    }

    free(cont_type_str);
    return 0;
}


int body_eml_proces(const char *str_from_file, const char *boundary, size_t *n_parts) {
    if (!str_from_file || !n_parts) {
        return -1;
    }

    char *str = NULL;
    if (asprintf(&str, "%s", str_from_file) < 0) {
        return -1;
    }

    char *new_str = str_proces(str);
    if (!new_str) {
        free(str);
        return -1;
    }

    if (!strcmp(new_str, boundary)) {
        *n_parts += 1;
    }

    free(str);
    free(new_str);
    return 0;
}


int file_proces(FILE *f, long file_size, char **from_str, char **to_str, char **date_str, size_t *n_parts) {
    const int STR_LENGTH = 1000;

    unsigned short check_parts = FALSE;
    unsigned short check_body = FALSE;
    size_t check_cont_type = FALSE;
    char *boundary = NULL;

    char *str_from_file = malloc(STR_LENGTH);
    if (!str_from_file) {
        return -1;
    }

    char *next_str = malloc(STR_LENGTH);
    if (!next_str) {
        free(str_from_file);
        return -1;
    }
    fgets(next_str, STR_LENGTH, f);

    while (ftell(f) != file_size) {
        snprintf(str_from_file, strlen(next_str) + 1, "%s", next_str);
        fgets(next_str, STR_LENGTH, f);

        if (!strncmp(str_from_file, "\n", 1) || !strncmp(str_from_file, "\r", 1)) {
            check_parts = TRUE;

            if (check_parts == TRUE && strncmp(next_str, "\n", 1) && strncmp(next_str, "\r", 1)) {
                check_body = TRUE;
            }
        }

        if (check_parts == FALSE) {
            if (head_eml_proces(f, str_from_file, next_str, &boundary,
                            from_str, to_str, date_str, &check_cont_type) != 0) {
                free(str_from_file);
                free(next_str);
                free(boundary);
                return -1;
            }
        } else if (boundary) {
            if (body_eml_proces(str_from_file, boundary, n_parts) != 0) {
                free(str_from_file);
                free(next_str);
                free(boundary);
                return -1;
            }
        }
    }

    if (check_cont_type == FALSE) {
        *n_parts = 1;
    }
    if (check_body == FALSE) {
        *n_parts = 0;
    }

    if (!*from_str) {
        *from_str = strdup("");
    }
    if (!*to_str) {
        *to_str = strdup("");
    }
    if (!*date_str) {
        *date_str = strdup("");
    }

    free(str_from_file);
    free(next_str);
    free(boundary);
    return 0;
}


int main(int argc, const char **argv) {
    if (argc != 2) {
        return -1;
    }

    const char *path_to_eml = argv[1];
    FILE *f = fopen(path_to_eml, "r");
    if (!f) {
        printf("Open file failed: %s\n", strerror(errno));
        return -1;
    }

    struct stat st;
    if (fstat(fileno(f), &st) != 0) {
        printf("fstat failed: %s\n", strerror(errno));
        fclose(f);
        return -1;
	}

    char *ptr = mmap(NULL,
                     st.st_size,
                     PROT_READ,
                     MAP_PRIVATE | MAP_FILE,
                     fileno(f),
                     0);
    if (ptr == MAP_FAILED) {
        printf("mmap failed: %s\n", strerror(errno));
        fclose(f);
        return -1;
    }

    char *from_str = NULL;
    char *to_str = NULL;
    char *date_str = NULL;
    size_t n_parts = 0;

    if (file_proces(f, st.st_size, &from_str, &to_str, &date_str, &n_parts) != 0) {
        printf("Allocation failed/ operation with NULL-pointer: %s\n", strerror(errno));
        free(from_str);
        free(to_str);
        free(date_str);
        munmap(ptr, st.st_size);
        fclose(f);
        return -1;
    }

    fprintf(stdout, "%s|%s|%s|%zu\n", from_str, to_str, date_str, n_parts);

    free(from_str);
    free(to_str);
    free(date_str);
    munmap(ptr, st.st_size);
    fclose(f);
    return 0;
}
