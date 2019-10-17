#include "samtools.pysam.h"

/*  test/split/test_filter_header_rg.c -- split test cases.

    Copyright (C) 2014 Genome Research Ltd.

    Author: Martin O. Pollard <mp15@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#include <config.h>

#include "../../bam_split.c"
#include "../test.h"
#include <unistd.h>

void setup_test_1(bam_hdr_t** hdr_in)
{
    *hdr_in = bam_hdr_init();
    const char *test1 =
    "@HD\tVN:1.4\n"
    "@SQ\tSN:blah\n"
    "@RG\tID:fish\n";
    (*hdr_in)->text = strdup(test1);
    (*hdr_in)->l_text = strlen(test1);
}

bool check_test_1(const bam_hdr_t* hdr) {
    char test1_res[200];
    snprintf(test1_res, 199,
    "@HD\tVN:1.4\n"
    "@SQ\tSN:blah\n"
    "@PG\tID:samtools\tPN:samtools\tVN:%s\tCL:test_filter_header_rg foo bar baz\n", samtools_version());

    if (strcmp(hdr->text, test1_res)) {
        return false;
    }
    return true;
}

void setup_test_2(bam_hdr_t** hdr_in)
{
    *hdr_in = bam_hdr_init();
    const char *test2 =
    "@HD\tVN:1.4\n"
    "@SQ\tSN:blah\n"
    "@RG\tID:fish\n";
    (*hdr_in)->text = strdup(test2);
    (*hdr_in)->l_text = strlen(test2);
}

bool check_test_2(const bam_hdr_t* hdr) {
    char test2_res[200];
    snprintf(test2_res, 199,
    "@HD\tVN:1.4\n"
    "@SQ\tSN:blah\n"
    "@RG\tID:fish\n"
    "@PG\tID:samtools\tPN:samtools\tVN:%s\tCL:test_filter_header_rg foo bar baz\n", samtools_version());

    if (strcmp(hdr->text, test2_res)) {
        return false;
    }
    return true;
}

int samtools_test_filter_header_rg_main(int argc, char *argv[])
{
    // test state
    const int NUM_TESTS = 2;
    int verbose = 0;
    int success = 0;
    int failure = 0;

    int getopt_char;
    char *test_argv[] = { "test_filter_header_rg", "foo\tbar", "baz" };
    char *arg_list = stringify_argv(3, test_argv);
    while ((getopt_char = getopt(argc, argv, "v")) != -1) {
        switch (getopt_char) {
            case 'v':
                ++verbose;
                break;
            default:
                fprintf(samtools_stdout, 
                       "usage: test_filter_header_rg [-v]\n\n"
                       " -v verbose output\n"
                       );
                break;
        }
    }


    // Setup samtools_stderr redirect
    kstring_t res = { 0, 0, NULL };
    FILE* orig_samtools_stderr = fdopen(dup(STDERR_FILENO), "a"); // Save samtools_stderr
    char* tempfname = (optind < argc)? argv[optind] : "test_count_rg.tmp";
    FILE* check = NULL;

    // setup
    if (verbose) fprintf(samtools_stdout, "BEGIN test 1\n");  // test eliminating a tag that isn't there
    bam_hdr_t* hdr1;
    const char* id_to_keep_1 = "1#2.3";
    setup_test_1(&hdr1);
    if (verbose > 0) {
        fprintf(samtools_stdout, "hdr1\n");
        dump_hdr(hdr1);
    }
    if (verbose) fprintf(samtools_stdout, "RUN test 1\n");

    // test
    xfreopen(tempfname, "w", samtools_stderr); // Redirect samtools_stderr to pipe
    bool result_1 = filter_header_rg(hdr1, id_to_keep_1, arg_list);
    fclose(samtools_stderr);

    if (verbose) fprintf(samtools_stdout, "END RUN test 1\n");
    if (verbose > 0) {
        fprintf(samtools_stdout, "hdr1\n");
        dump_hdr(hdr1);
    }

    // check result
    res.l = 0;
    check = fopen(tempfname, "r");
    if ( result_1
        && check_test_1(hdr1)
        && kgetline(&res, (kgets_func *)fgets, check) < 0
        && (feof(check) || res.l == 0)) {
        ++success;
    } else {
        ++failure;
        if (verbose) fprintf(samtools_stdout, "FAIL test 1\n");
    }
    fclose(check);

    // teardown
    bam_hdr_destroy(hdr1);
    if (verbose) fprintf(samtools_stdout, "END test 1\n");

    if (verbose) fprintf(samtools_stdout, "BEGIN test 2\n");  // test eliminating a tag that is there
    bam_hdr_t* hdr2;
    const char* id_to_keep_2 = "fish";
    setup_test_2(&hdr2);
    if (verbose > 0) {
        fprintf(samtools_stdout, "hdr2\n");
        dump_hdr(hdr2);
    }
    if (verbose) fprintf(samtools_stdout, "RUN test 2\n");

    // test
    xfreopen(tempfname, "w", samtools_stderr); // Redirect samtools_stderr to pipe
    bool result_2 = filter_header_rg(hdr2, id_to_keep_2, arg_list);
    fclose(samtools_stderr);

    if (verbose) fprintf(samtools_stdout, "END RUN test 2\n");
    if (verbose > 0) {
        fprintf(samtools_stdout, "hdr2\n");
        dump_hdr(hdr2);
    }

    // check result
    res.l = 0;
    check = fopen(tempfname, "r");
    if ( result_2
        && check_test_2(hdr2)
        && kgetline(&res, (kgets_func *)fgets, check) < 0
        && (feof(check) || res.l == 0)) {
        ++success;
    } else {
        ++failure;
        if (verbose) fprintf(samtools_stdout, "FAIL test 2\n");
    }
    fclose(check);

    // teardown
    bam_hdr_destroy(hdr2);
    if (verbose) fprintf(samtools_stdout, "END test 2\n");


    // Cleanup
    free(res.s);
    free(arg_list);
    remove(tempfname);
    if (failure > 0)
        fprintf(orig_samtools_stderr, "%d failures %d successes\n", failure, success);
    fclose(orig_samtools_stderr);

    return (success == NUM_TESTS)? EXIT_SUCCESS : EXIT_FAILURE;
}