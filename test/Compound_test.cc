/*
 * Copyright(c) 2019 Netflix, Inc.
 * SPDX - License - Identifier: BSD - 2 - Clause - Patent
 */

/******************************************************************************
 * @Compound_test.cc
 *
 * @brief Unit test for compound functions:
 * - av1_wedge_sse_from_residuals_avx2
 * @author Tahani
 *
 ******************************************************************************/
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <limits.h>
#include <new>
// workaround to eliminate the compiling warning on linux
// The macro will conflict with definition in gtest.h
#ifdef __USE_GNU
#undef __USE_GNU  // defined in EbThreads.h
#endif
#ifdef _GNU_SOURCE
#undef _GNU_SOURCE  // defined in EbThreads.h
#endif


#include "EbDefinitions.h"
#include "EbUnitTestUtility.h"
#include "EbUtility.h"
#include "random.h"
#include "util.h"

#include "aom_dsp_rtcd.h"

namespace {

typedef uint64_t (*av1_wedge_sse_from_residuals_func)(const int16_t *r1,
                                                      const int16_t *d,
                                                      const uint8_t *m, int N);
// test av1_wedge_sse_from_residuals_avx2
class WedgeTest 
      : public ::testing::TestWithParam<av1_wedge_sse_from_residuals_func> {
  public:
    WedgeTest() : func_(GetParam()) {
    }

    ~WedgeTest();

    void SetUp() override {
        residual = (int16_t *)(eb_aom_memalign(32, 2 * MAX_SB_SQUARE));
        diff     = (int16_t *)(eb_aom_memalign(32, 2 * MAX_SB_SQUARE));
        mask     = (uint8_t *)(eb_aom_memalign(16, 2 * MAX_SB_SQUARE));

    }

    void TearDown() override {
        if (residual)
            eb_aom_free(residual);
        if (diff)
            eb_aom_free(diff);
    }

  protected:
    void run_test();

    void init_data() {
        eb_buf_random_s16(residual, MAX_SB_SQUARE);
        eb_buf_random_s16(diff, MAX_SB_SQUARE);
        eb_buf_random_u8_with_max(mask, 1, 64);
        N = MAX((int32_t)(64 * (rand() % 256)), MAX_SB_SQUARE);
    }

    void check_output(uint64_t sse_ref, uint64_t sse_test, int N) {
      
        if (sse_ref != sse_test)
        {
            ASSERT_EQ(sse_ref, sse_test)
                << " SSE Error at N : " << N 
                << " reference value : " << sse_ref 
                << " test value : " << sse_test;
               
        } 
        }


    av1_wedge_sse_from_residuals_func func_;
    int16_t *residual;

    int16_t *diff;

    uint8_t *mask;

    int32_t N;
    uint64_t sse_ref, sse_test;
};


WedgeTest::~WedgeTest() {
}

 void WedgeTest::run_test() {
    init_data();
    uint8_t number_runs = 10;
    uint8_t counter;
    for (int counter = 0; counter < number_runs; counter++) {
        sse_ref = av1_wedge_sse_from_residuals_c(residual, diff, mask, N);
        sse_test = func_(residual, diff, mask, N);

        check_output(sse_ref, sse_test, N);
    }
}
TEST_P(WedgeTest, WedgeTest) {
    run_test();
};

INSTANTIATE_TEST_CASE_P(WEDGE, WedgeTest,
                        ::testing::Values(av1_wedge_sse_from_residuals_avx2));



}  // namespace
