/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include "EbComputeSAD_C.h"
#include "EbUtility.h"

/*******************************************
* CombinedAveragingSAD
*
*******************************************/
uint32_t CombinedAveragingSAD(
    uint8_t  *src,
    uint32_t  src_stride,
    uint8_t  *ref1,
    uint32_t  ref1Stride,
    uint8_t  *ref2,
    uint32_t  ref2Stride,
    uint32_t  height,
    uint32_t  width)
{
    uint32_t x, y;
    uint32_t sad = 0;
    uint8_t avgpel;

    for (y = 0; y < height; y++)
    {
        for (x = 0; x < width; x++)
        {
            avgpel = (ref1[x] + ref2[x] + 1) >> 1;
            sad += EB_ABS_DIFF(src[x], avgpel);
        }
        src += src_stride;
        ref1 += ref1Stride;
        ref2 += ref2Stride;
    }

    return sad;
}

/*******************************************
*   returns NxM Sum of Absolute Differences
Note: moved from picture operators.
keep this function here for profiling
issues.
*******************************************/
uint32_t FastLoop_NxMSadKernel(
    uint8_t  *src,                            // input parameter, source samples Ptr
    uint32_t  src_stride,                      // input parameter, source stride
    uint8_t  *ref,                            // input parameter, reference samples Ptr
    uint32_t  refStride,                      // input parameter, reference stride
    uint32_t  height,                         // input parameter, block height (M)
    uint32_t  width)                          // input parameter, block width (N)
{
    uint32_t x, y;
    uint32_t sad = 0;

    for (y = 0; y < height; y++)
    {
        for (x = 0; x < width; x++)
        {
            sad += EB_ABS_DIFF(src[x], ref[x]);
        }
        src += src_stride;
        ref += refStride;
    }

    return sad;
}

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
    int16_t search_area_height)
{
    int16_t xSearchIndex;
    int16_t ySearchIndex;

    *bestSad = 0xffffff;

    for (ySearchIndex = 0; ySearchIndex < search_area_height; ySearchIndex++)
    {
        for (xSearchIndex = 0; xSearchIndex < search_area_width; xSearchIndex++)
        {
            uint32_t x, y;
            uint32_t sad = 0;

            for (y = 0; y < height; y++)
            {
                for (x = 0; x < width; x++)
                {
                    sad += EB_ABS_DIFF(src[y*src_stride + x], ref[xSearchIndex + y * refStride + x]);
                }

            }

            // Update results
            if (sad < *bestSad)
            {
                *bestSad = sad;
                *xSearchCenter = xSearchIndex;
                *ySearchCenter = ySearchIndex;
            }
        }

        ref += srcStrideRaw;
    }

    return;
}
//compute a 8x4 SAD  
static uint32_t Subsad8x8(
    uint8_t  *src,                            // input parameter, source samples Ptr
    uint32_t  srcStride,                      // input parameter, source stride
    uint8_t  *ref,                            // input parameter, reference samples Ptr
    uint32_t  refStride)                     // input parameter, reference stride
{
    uint32_t x, y;
    uint32_t sadBlock8x4 = 0;

#if USE_SAD_ME
    for (y = 0; y < 8; y++)
    {
        for (x = 0; x < 8; x++)
        {
            sadBlock8x4 += EB_ABS_DIFF(src[y*srcStride + x], ref[y*refStride + x]);
        }

    }
#else
    srcStride = srcStride * 2;
    refStride = refStride * 2;

    for (y = 0; y < 4; y++)
    {
        for (x = 0; x < 8; x++)
        {
            sadBlock8x4 += EB_ABS_DIFF(src[y*srcStride + x], ref[y*refStride + x]);
        }

    }
#endif
    return sadBlock8x4;
}
/*******************************************
* GetEightHorizontalSearchPointResults_8x8_16x16_PU
*******************************************/
void GetEightHorizontalSearchPointResults_8x8_16x16_PU(
    uint8_t   *src,
    uint32_t   src_stride,
    uint8_t   *ref,
    uint32_t   ref_stride,
    uint32_t  *p_best_sad8x8,
    uint32_t  *p_best_mv8x8,
    uint32_t  *p_best_sad16x16,
    uint32_t  *p_best_mv16x16,
    uint32_t   mv,
    uint16_t  *p_sad16x16)
{
    uint32_t xSearchIndex;
    int16_t xMv, yMv;
    uint32_t sad8x8_0, sad8x8_1, sad8x8_2, sad8x8_3;
    uint16_t sad16x16;


    /*
    -------------------------------------   -----------------------------------
    | 8x8_00 | 8x8_01 | 8x8_04 | 8x8_05 |   8x8_16 | 8x8_17 | 8x8_20 | 8x8_21 |
    -------------------------------------   -----------------------------------
    | 8x8_02 | 8x8_03 | 8x8_06 | 8x8_07 |   8x8_18 | 8x8_19 | 8x8_22 | 8x8_23 |
    -----------------------   -----------   ----------------------   ----------
    | 8x8_08 | 8x8_09 | 8x8_12 | 8x8_13 |   8x8_24 | 8x8_25 | 8x8_29 | 8x8_29 |
    ----------------------    -----------   ---------------------    ----------
    | 8x8_10 | 8x8_11 | 8x8_14 | 8x8_15 |   8x8_26 | 8x8_27 | 8x8_30 | 8x8_31 |
    -------------------------------------   -----------------------------------

    -------------------------------------   -----------------------------------
    | 8x8_32 | 8x8_33 | 8x8_36 | 8x8_37 |   8x8_48 | 8x8_49 | 8x8_52 | 8x8_53 |
    -------------------------------------   -----------------------------------
    | 8x8_34 | 8x8_35 | 8x8_38 | 8x8_39 |   8x8_50 | 8x8_51 | 8x8_54 | 8x8_55 |
    -----------------------   -----------   ----------------------   ----------
    | 8x8_40 | 8x8_41 | 8x8_44 | 8x8_45 |   8x8_56 | 8x8_57 | 8x8_60 | 8x8_61 |
    ----------------------    -----------   ---------------------    ----------
    | 8x8_42 | 8x8_43 | 8x8_46 | 8x8_48 |   8x8_58 | 8x8_59 | 8x8_62 | 8x8_63 |
    -------------------------------------   -----------------------------------
    */

    /*
    ----------------------    ----------------------
    |  16x16_0  |  16x16_1  |  16x16_4  |  16x16_5  |
    ----------------------    ----------------------
    |  16x16_2  |  16x16_3  |  16x16_6  |  16x16_7  |
    -----------------------   -----------------------
    |  16x16_8  |  16x16_9  |  16x16_12 |  16x16_13 |
    ----------------------    ----------------------
    |  16x16_10 |  16x16_11 |  16x16_14 |  16x16_15 |
    -----------------------   -----------------------
    */

    for (xSearchIndex = 0; xSearchIndex < 8; xSearchIndex++)
    {
#if USE_SAD_ME
        //8x8_0        
        sad8x8_0 = Subsad8x8(src, src_stride, ref + xSearchIndex, ref_stride);
        if (sad8x8_0 < p_best_sad8x8[0]) {
            p_best_sad8x8[0] = sad8x8_0;
            xMv = _MVXT(mv) + (int16_t)xSearchIndex * 4;
            yMv = _MVYT(mv);
            p_best_mv8x8[0] = ((uint16_t)yMv << 16) | ((uint16_t)xMv);
        }

        //8x8_1        
        sad8x8_1 = Subsad8x8(src + 8, src_stride, ref + xSearchIndex + 8, ref_stride);
        if (sad8x8_1 < p_best_sad8x8[1]) {
            p_best_sad8x8[1] = sad8x8_1;
            xMv = _MVXT(mv) + (int16_t)xSearchIndex * 4;
            yMv = _MVYT(mv);
            p_best_mv8x8[1] = ((uint16_t)yMv << 16) | ((uint16_t)xMv);
        }

        //8x8_2        
        sad8x8_2 = Subsad8x8(src + 8 * src_stride, src_stride, ref + xSearchIndex + 8 * ref_stride, ref_stride);
        if (sad8x8_2 < p_best_sad8x8[2]) {
            p_best_sad8x8[2] = sad8x8_2;
            xMv = _MVXT(mv) + (int16_t)xSearchIndex * 4;
            yMv = _MVYT(mv);
            p_best_mv8x8[2] = ((uint16_t)yMv << 16) | ((uint16_t)xMv);
        }

        //8x8_3        
        sad8x8_3 = Subsad8x8(src + 8 + 8 * src_stride, src_stride, ref + 8 + 8 * ref_stride + xSearchIndex, ref_stride);
        if (sad8x8_3 < p_best_sad8x8[3]) {
            p_best_sad8x8[3] = sad8x8_3;
            xMv = _MVXT(mv) + (int16_t)xSearchIndex * 4;
            yMv = _MVYT(mv);
            p_best_mv8x8[3] = ((uint16_t)yMv << 16) | ((uint16_t)xMv);
        }


        //16x16
        sad16x16 = (uint16_t)(sad8x8_0 + sad8x8_1 + sad8x8_2 + sad8x8_3);
        p_sad16x16[xSearchIndex] = sad16x16;  //store the intermediate 16x16 SAD for 32x32 in subsampled form.
        if ((uint32_t)(sad16x16) < p_best_sad16x16[0]) {
            p_best_sad16x16[0] = sad16x16;
            xMv = _MVXT(mv) + (int16_t)xSearchIndex * 4;
            yMv = _MVYT(mv);
            p_best_mv16x16[0] = ((uint16_t)yMv << 16) | ((uint16_t)xMv);
        }
#else
        //8x8_0        
        sad8x8_0 = Subsad8x8(src, src_stride, ref + xSearchIndex, ref_stride);
        if (2 * sad8x8_0 < p_best_sad8x8[0]) {
            p_best_sad8x8[0] = 2 * sad8x8_0;
            xMv = _MVXT(mv) + (int16_t)xSearchIndex * 4;
            yMv = _MVYT(mv);
            p_best_mv8x8[0] = ((uint16_t)yMv << 16) | ((uint16_t)xMv);
        }

        //8x8_1        
        sad8x8_1 = Subsad8x8(src + 8, src_stride, ref + xSearchIndex + 8, ref_stride);
        if (2 * sad8x8_1 < p_best_sad8x8[1]) {
            p_best_sad8x8[1] = 2 * sad8x8_1;
            xMv = _MVXT(mv) + (int16_t)xSearchIndex * 4;
            yMv = _MVYT(mv);
            p_best_mv8x8[1] = ((uint16_t)yMv << 16) | ((uint16_t)xMv);
        }

        //8x8_2        
        sad8x8_2 = Subsad8x8(src + 8 * src_stride, src_stride, ref + xSearchIndex + 8 * ref_stride, ref_stride);
        if (2 * sad8x8_2 < p_best_sad8x8[2]) {
            p_best_sad8x8[2] = 2 * sad8x8_2;
            xMv = _MVXT(mv) + (int16_t)xSearchIndex * 4;
            yMv = _MVYT(mv);
            p_best_mv8x8[2] = ((uint16_t)yMv << 16) | ((uint16_t)xMv);
        }

        //8x8_3        
        sad8x8_3 = Subsad8x8(src + 8 + 8 * src_stride, src_stride, ref + 8 + 8 * ref_stride + xSearchIndex, ref_stride);
        if (2 * sad8x8_3 < p_best_sad8x8[3]) {
            p_best_sad8x8[3] = 2 * sad8x8_3;
            xMv = _MVXT(mv) + (int16_t)xSearchIndex * 4;
            yMv = _MVYT(mv);
            p_best_mv8x8[3] = ((uint16_t)yMv << 16) | ((uint16_t)xMv);
        }


        //16x16
        sad16x16 = (uint16_t)(sad8x8_0 + sad8x8_1 + sad8x8_2 + sad8x8_3);
        p_sad16x16[xSearchIndex] = sad16x16;  //store the intermediate 16x16 SAD for 32x32 in subsampled form.
        if ((uint32_t)(2 * sad16x16) < p_best_sad16x16[0]) {
            p_best_sad16x16[0] = 2 * sad16x16;
            xMv = _MVXT(mv) + (int16_t)xSearchIndex * 4;
            yMv = _MVYT(mv);
            p_best_mv16x16[0] = ((uint16_t)yMv << 16) | ((uint16_t)xMv);
        }
#endif

    }
}


/*******************************************
Calcualte SAD for 32x32,64x64 from 16x16
and check if there is improvement, if yes keep
the best SAD+MV
*******************************************/
void GetEightHorizontalSearchPointResults_32x32_64x64(
    uint16_t  *p_sad16x16,
    uint32_t  *p_best_sad32x32,
    uint32_t  *p_best_sad64x64,
    uint32_t  *p_best_mv32x32,
    uint32_t  *p_best_mv64x64,
    uint32_t   mv)
{
    int16_t xMv, yMv;
    uint32_t sad32x32_0, sad32x32_1, sad32x32_2, sad32x32_3, sad64x64;
    uint32_t xSearchIndex;

    /*--------------------
    |  32x32_0  |  32x32_1
    ----------------------
    |  32x32_2  |  32x32_3
    ----------------------*/


    /*  data ordering in pSad16x16 buffer

    Search    Search            Search
    Point 0   Point 1           Point 7
    ---------------------------------------
    16x16_0    |    x    |    x    | ...... |    x    |
    ---------------------------------------
    16x16_1    |    x    |    x    | ...... |    x    |

    16x16_n    |    x    |    x    | ...... |    x    |

    ---------------------------------------
    16x16_15   |    x    |    x    | ...... |    x    |
    ---------------------------------------
    */



    for (xSearchIndex = 0; xSearchIndex < 8; xSearchIndex++)
    {
#if USE_SAD_ME
        //32x32_0
        sad32x32_0 = p_sad16x16[0 * 8 + xSearchIndex] + p_sad16x16[1 * 8 + xSearchIndex] + p_sad16x16[2 * 8 + xSearchIndex] + p_sad16x16[3 * 8 + xSearchIndex];

        if (sad32x32_0 < p_best_sad32x32[0]) {
            p_best_sad32x32[0] = sad32x32_0;
            xMv = _MVXT(mv) + (int16_t)xSearchIndex * 4;
            yMv = _MVYT(mv);
            p_best_mv32x32[0] = ((uint16_t)yMv << 16) | ((uint16_t)xMv);
        }

        //32x32_1
        sad32x32_1 = p_sad16x16[4 * 8 + xSearchIndex] + p_sad16x16[5 * 8 + xSearchIndex] + p_sad16x16[6 * 8 + xSearchIndex] + p_sad16x16[7 * 8 + xSearchIndex];

        if (sad32x32_1 < p_best_sad32x32[1]) {
            p_best_sad32x32[1] = sad32x32_1;
            xMv = _MVXT(mv) + (int16_t)xSearchIndex * 4;
            yMv = _MVYT(mv);
            p_best_mv32x32[1] = ((uint16_t)yMv << 16) | ((uint16_t)xMv);
        }

        //32x32_2
        sad32x32_2 = p_sad16x16[8 * 8 + xSearchIndex] + p_sad16x16[9 * 8 + xSearchIndex] + p_sad16x16[10 * 8 + xSearchIndex] + p_sad16x16[11 * 8 + xSearchIndex];

        if (sad32x32_2 < p_best_sad32x32[2]) {
            p_best_sad32x32[2] = sad32x32_2;
            xMv = _MVXT(mv) + (int16_t)xSearchIndex * 4;
            yMv = _MVYT(mv);
            p_best_mv32x32[2] = ((uint16_t)yMv << 16) | ((uint16_t)xMv);
        }


        //32x32_3
        sad32x32_3 = p_sad16x16[12 * 8 + xSearchIndex] + p_sad16x16[13 * 8 + xSearchIndex] + p_sad16x16[14 * 8 + xSearchIndex] + p_sad16x16[15 * 8 + xSearchIndex];

        if (sad32x32_3 < p_best_sad32x32[3]) {
            p_best_sad32x32[3] = sad32x32_3;
            xMv = _MVXT(mv) + (int16_t)xSearchIndex * 4;
            yMv = _MVYT(mv);
            p_best_mv32x32[3] = ((uint16_t)yMv << 16) | ((uint16_t)xMv);
        }


        //64x64
        sad64x64 = sad32x32_0 + sad32x32_1 + sad32x32_2 + sad32x32_3;
        if (sad64x64 < p_best_sad64x64[0]) {
            p_best_sad64x64[0] = sad64x64;
            xMv = _MVXT(mv) + (int16_t)xSearchIndex * 4;
            yMv = _MVYT(mv);
            p_best_mv64x64[0] = ((uint16_t)yMv << 16) | ((uint16_t)xMv);

        }
#else
        //32x32_0
        sad32x32_0 = p_sad16x16[0 * 8 + xSearchIndex] + p_sad16x16[1 * 8 + xSearchIndex] + p_sad16x16[2 * 8 + xSearchIndex] + p_sad16x16[3 * 8 + xSearchIndex];

        if (2 * sad32x32_0 < p_best_sad32x32[0]) {
            p_best_sad32x32[0] = 2 * sad32x32_0;
            xMv = _MVXT(mv) + (int16_t)xSearchIndex * 4;
            yMv = _MVYT(mv);
            p_best_mv32x32[0] = ((uint16_t)yMv << 16) | ((uint16_t)xMv);
        }

        //32x32_1
        sad32x32_1 = p_sad16x16[4 * 8 + xSearchIndex] + p_sad16x16[5 * 8 + xSearchIndex] + p_sad16x16[6 * 8 + xSearchIndex] + p_sad16x16[7 * 8 + xSearchIndex];

        if (2 * sad32x32_1 < p_best_sad32x32[1]) {
            p_best_sad32x32[1] = 2 * sad32x32_1;
            xMv = _MVXT(mv) + (int16_t)xSearchIndex * 4;
            yMv = _MVYT(mv);
            p_best_mv32x32[1] = ((uint16_t)yMv << 16) | ((uint16_t)xMv);
        }

        //32x32_2
        sad32x32_2 = p_sad16x16[8 * 8 + xSearchIndex] + p_sad16x16[9 * 8 + xSearchIndex] + p_sad16x16[10 * 8 + xSearchIndex] + p_sad16x16[11 * 8 + xSearchIndex];

        if (2 * sad32x32_2 < p_best_sad32x32[2]) {
            p_best_sad32x32[2] = 2 * sad32x32_2;
            xMv = _MVXT(mv) + (int16_t)xSearchIndex * 4;
            yMv = _MVYT(mv);
            p_best_mv32x32[2] = ((uint16_t)yMv << 16) | ((uint16_t)xMv);
        }


        //32x32_3
        sad32x32_3 = p_sad16x16[12 * 8 + xSearchIndex] + p_sad16x16[13 * 8 + xSearchIndex] + p_sad16x16[14 * 8 + xSearchIndex] + p_sad16x16[15 * 8 + xSearchIndex];

        if (2 * sad32x32_3 < p_best_sad32x32[3]) {
            p_best_sad32x32[3] = 2 * sad32x32_3;
            xMv = _MVXT(mv) + (int16_t)xSearchIndex * 4;
            yMv = _MVYT(mv);
            p_best_mv32x32[3] = ((uint16_t)yMv << 16) | ((uint16_t)xMv);
        }


        //64x64
        sad64x64 = sad32x32_0 + sad32x32_1 + sad32x32_2 + sad32x32_3;
        if (2 * sad64x64 < p_best_sad64x64[0]) {
            p_best_sad64x64[0] = 2 * sad64x64;
            xMv = _MVXT(mv) + (int16_t)xSearchIndex * 4;
            yMv = _MVYT(mv);
            p_best_mv64x64[0] = ((uint16_t)yMv << 16) | ((uint16_t)xMv);

        }

#endif

    }

}