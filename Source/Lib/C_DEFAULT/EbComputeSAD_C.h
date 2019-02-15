/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/
#ifndef EbComputeSAD_C_h
#define EbComputeSAD_C_h

#include "EbDefinitions.h"
#ifdef __cplusplus
extern "C" {
#endif

    uint32_t FastLoop_NxMSadKernel(
        uint8_t  *src,                            // input parameter, source samples Ptr
        uint32_t  src_stride,                      // input parameter, source stride
        uint8_t  *ref,                            // input parameter, reference samples Ptr
        uint32_t  refStride,                      // input parameter, reference stride
        uint32_t  height,                         // input parameter, block height (M)
        uint32_t  width);                         // input parameter, block width (N)

    uint32_t CombinedAveragingSAD(
        uint8_t  *src,
        uint32_t  src_stride,
        uint8_t  *ref1,
        uint32_t  ref1Stride,
        uint8_t  *ref2,
        uint32_t  ref2Stride,
        uint32_t  height,
        uint32_t  width);

    void SadLoopKernel(
        uint8_t  *src,                            // input parameter, source samples Ptr
        uint32_t  src_stride,                      // input parameter, source stride
        uint8_t  *ref,                            // input parameter, reference samples Ptr
        uint32_t  refStride,                      // input parameter, reference stride
        uint32_t  height,                         // input parameter, block height (M)
        uint32_t  width,                          // input parameter, block width (N)
        uint64_t *bestSad,
        int16_t *xSearchCenter,
        int16_t *ySearchCenter,
        uint32_t  srcStrideRaw,                   // input parameter, source stride (no line skipping)
        int16_t search_area_width,
        int16_t search_area_height);
#if USE_SAD_ME
    void GetEightHorizontalSearchPointResults_8x8_16x16_PU(
        uint8_t   *src,
        uint32_t   src_stride,
        uint8_t   *ref,
        uint32_t   refStride,
        uint32_t  *p_best_sad8x8,
        uint32_t  *p_best_mv8x8,
        uint32_t  *p_best_sad16x16,
        uint32_t  *p_best_mv16x16,
        uint32_t   mv,
        uint16_t  *p_sad16x16);

    void GetEightHorizontalSearchPointResults_32x32_64x64(
        uint16_t  *p_sad16x16,
        uint32_t  *p_best_sad32x32,
        uint32_t  *p_best_sad64x64,
        uint32_t  *p_best_mv32x32,
        uint32_t  *p_best_mv64x64,
        uint32_t   mv);
#endif

#ifdef __cplusplus
}
#endif
#endif // EbComputeSAD_C_h