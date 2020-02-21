/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbPictureOperators_AVX2_h
#define EbPictureOperators_AVX2_h

#include <immintrin.h>
#include "EbDefinitions.h"
#include "EbPictureOperators_SSE2.h"

#ifdef __cplusplus
extern "C" {
#endif

extern void eb_enc_msb_pack2d_avx2_intrin_al(uint8_t *in8_bit_buffer, uint32_t in8_stride,
                                             uint8_t *inn_bit_buffer, uint16_t *out16_bit_buffer,
                                             uint32_t inn_stride, uint32_t out_stride,
                                             uint32_t width, uint32_t height);

extern void compressed_packmsb_avx2_intrin(uint8_t *in8_bit_buffer, uint32_t in8_stride,
                                           uint8_t *inn_bit_buffer, uint16_t *out16_bit_buffer,
                                           uint32_t inn_stride, uint32_t out_stride, uint32_t width,
                                           uint32_t height);

void c_pack_avx2_intrin(const uint8_t *inn_bit_buffer, uint32_t inn_stride,
                        uint8_t *in_compn_bit_buffer, uint32_t out_stride, uint8_t *local_cache,
                        uint32_t width, uint32_t height);

void unpack_avg_avx2_intrin(uint16_t *ref16_l0, uint32_t ref_l0_stride, uint16_t *ref16_l1,
                            uint32_t ref_l1_stride, uint8_t *dst_ptr, uint32_t dst_stride,
                            uint32_t width, uint32_t height);

void unpack_avg_safe_sub_avx2_intrin(uint16_t *ref16_l0, uint32_t ref_l0_stride, uint16_t *ref16_l1,
                                     uint32_t ref_l1_stride, uint8_t *dst_ptr, uint32_t dst_stride,
                                     EbBool sub_pred, uint32_t width, uint32_t height);
void full_distortion_kernel_cbf_zero32_bits_avx2(int32_t *coeff, uint32_t coeff_stride,
                                                 int32_t *recon_coeff, uint32_t recon_coeff_stride,
                                                 uint64_t distortion_result[DIST_CALC_TOTAL],
                                                 uint32_t area_width, uint32_t area_height);

void full_distortion_kernel32_bits_avx2(int32_t *coeff, uint32_t coeff_stride, int32_t *recon_coeff,
                                        uint32_t recon_coeff_stride,
                                        uint64_t distortion_result[DIST_CALC_TOTAL],
                                        uint32_t area_width, uint32_t area_height);

uint64_t spatial_full_distortion_kernel4x_n_avx2_intrin(uint8_t *input, uint32_t input_offset,
                                                        uint32_t input_stride, uint8_t *recon,
                                                        uint32_t recon_offset,
                                                        uint32_t recon_stride, uint32_t area_width,
                                                        uint32_t area_height);

uint64_t spatial_full_distortion_kernel8x_n_avx2_intrin(uint8_t *input, uint32_t input_offset,
                                                        uint32_t input_stride, uint8_t *recon,
                                                        uint32_t recon_offset,
                                                        uint32_t recon_stride, uint32_t area_width,
                                                        uint32_t area_height);

uint64_t spatial_full_distortion_kernel16x_n_avx2_intrin(uint8_t *input, uint32_t input_offset,
                                                         uint32_t input_stride, uint8_t *recon,
                                                         uint32_t recon_offset,
                                                         uint32_t recon_stride, uint32_t area_width,
                                                         uint32_t area_height);

uint64_t spatial_full_distortion_kernel32x_n_avx2_intrin(uint8_t *input, uint32_t input_offset,
                                                         uint32_t input_stride, uint8_t *recon,
                                                         uint32_t recon_offset,
                                                         uint32_t recon_stride, uint32_t area_width,
                                                         uint32_t area_height);

uint64_t spatial_full_distortion_kernel64x_n_avx2_intrin(uint8_t *input, uint32_t input_offset,
                                                         uint32_t input_stride, uint8_t *recon,
                                                         uint32_t recon_offset,
                                                         uint32_t recon_stride, uint32_t area_width,
                                                         uint32_t area_height);

uint64_t spatial_full_distortion_kernel128x_n_avx2_intrin(
    uint8_t *input, uint32_t input_offset, uint32_t input_stride, uint8_t *recon,
    uint32_t recon_offset, uint32_t recon_stride, uint32_t area_width, uint32_t area_height);

void convert_8bit_to_16bit_avx2(uint8_t *src, uint32_t src_stride, uint16_t *dst,
                                uint32_t dst_stride, uint32_t width, uint32_t height);

#ifdef __cplusplus
}
#endif

#endif // EbPictureOperators_AVX2_h
