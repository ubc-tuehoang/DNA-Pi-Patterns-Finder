#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <pthread.h>

// Target DNA sequence to find
// Protein-DNA-Ion: PDB 7RCE
// Source: https://alphafoldserver.com/example/examplefold_pdb_7rce
#define PI_DIGITS 10000
#define TARGET_SEQUENCE "GGGGGCATGCAGATC"

// Chudnovsky algorithm for Pi computation
void compute_pi(char *pi_str, int digits) {
    int k, C = 426880 * 10005;
    int M = 1, L = 13591409, X = 1, K = 6, S = L;
    char *pi_buf = malloc(digits + 1);

    for (int i = 0; i < digits; i++) {
        int a = (M * L) / X;
        int b = a / (i + 1);
        pi_buf[i] = (b % 10) + '0';

        S += (a * (i % 2 == 0 ? 1 : -1));
        M = (k * (K + k)) / (i + 1);
        L += 545140134;
        X *= 640320;
        K += 12;
    }
    pi_buf[digits] = '\0';

    // Copy to pi_str
    strncpy(pi_str, pi_buf, digits);
    free(pi_buf);
}

// KMP Pattern Matching Algorithm
int KMP_search(char *text, char *pattern) {
    int n = strlen(text);
    int m = strlen(pattern);
    int *lps = malloc(m * sizeof(int));
    int i = 0, j = 0;

    // Preprocessing the pattern to create the longest prefix-suffix (LPS) array
    lps[0] = 0;
    for (i = 1, j = 0; i < m; i++) {
        while (j > 0 && pattern[i] != pattern[j]) {
            j = lps[j - 1];
        }
        if (pattern[i] == pattern[j]) {
            j++;
        }
        lps[i] = j;
    }

    // Search the pattern in the text
    for (i = 0, j = 0; i < n; i++) {
        while (j > 0 && text[i] != pattern[j]) {
            j = lps[j - 1];
        }
        if (text[i] == pattern[j]) {
            j++;
        }
        if (j == m) {
            free(lps);
            return i - m + 1; // Return index of the start of the match
        }
    }
    free(lps);
    return -1; // Pattern not found
}

// Function to map DNA bases to Pi sequence
void map_pi_to_dna(char *pi_str, char *dna_sequence, int *base_map, int pi_length) {
    int dna_index = 0;
    for (int i = 0; i < pi_length; i++) {
        char digit = pi_str[i];
        if (digit >= '0' && digit <= '9') {
            int num = digit - '0';  
            if (num == base_map[0]) {
                dna_sequence[dna_index++] = 'A';
            } else if (num == base_map[1]) {
                dna_sequence[dna_index++] = 'T';
            } else if (num == base_map[2]) {
                dna_sequence[dna_index++] = 'C';
            } else if (num == base_map[3]) {
                dna_sequence[dna_index++] = 'G';
            }
        }
    }
    dna_sequence[dna_index] = '\0';  
}

// Thread function for base map processing
void* process_base_map(void* arg) {
    int* base_map = (int*)arg;
    char pi_str[PI_DIGITS + 1];
    char dna_sequence[PI_DIGITS + 1];

    // Generate Pi sequence using Chudnovsky algorithm
    compute_pi(pi_str, PI_DIGITS);

    // Map Pi digits to DNA sequence
    map_pi_to_dna(pi_str, dna_sequence, base_map, PI_DIGITS);

    // Search for the target sequence in the DNA sequence
    int found_index = KMP_search(dna_sequence, TARGET_SEQUENCE);
    if (found_index != -1) {
        printf("Target sequence '%s' found at position %d with base map A:%d, T:%d, C:%d, G:%d.\n", 
               TARGET_SEQUENCE, found_index, base_map[0], base_map[1], base_map[2], base_map[3]);
    }

    return NULL;
}

int main() {
    pthread_t threads[10000];  // Attempting up to 10,000 threads for the base combinations
    int base_combinations[10][4]; // Each base can have 10 possible values (0-9)

    // Generate base combinations (0-9 for each base)
    int index = 0;
    for (int a = 0; a < 10; a++) {
        for (int t = 0; t < 10; t++) {
            for (int c = 0; c < 10; c++) {
                for (int g = 0; g < 10; g++) {
                    base_combinations[index][0] = a;
                    base_combinations[index][1] = t;
                    base_combinations[index][2] = c;
                    base_combinations[index][3] = g;
                    index++;
                }
            }
        }
    }

    // Create threads for processing each base combination
    for (int i = 0; i < index; i++) {
        pthread_create(&threads[i], NULL, process_base_map, (void*)base_combinations[i]);
    }

    // Join threads
    for (int i = 0; i < index; i++) {
        pthread_join(threads[i], NULL);
    }

    return 0;
}

