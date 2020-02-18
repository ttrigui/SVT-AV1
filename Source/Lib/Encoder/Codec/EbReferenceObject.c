/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include <stdlib.h>
#include <string.h>

#include "EbThreads.h"
#include "EbReferenceObject.h"
#include "EbPictureBufferDesc.h"

void initialize_samples_neighboring_reference_picture16_bit(EbByte   recon_samples_buffer_ptr,
                                                            uint16_t stride, uint16_t recon_width,
                                                            uint16_t recon_height,
                                                            uint16_t left_padding,
                                                            uint16_t top_padding) {
    uint16_t *recon_samples_ptr;
    uint16_t  sample_count;

    // 1. zero out the top row
    recon_samples_ptr =
        (uint16_t *)recon_samples_buffer_ptr + (top_padding - 1) * stride + left_padding - 1;
    EB_MEMSET((uint8_t *)recon_samples_ptr, 0, sizeof(uint16_t) * (1 + recon_width + 1));

    // 2. zero out the bottom row
    recon_samples_ptr = (uint16_t *)recon_samples_buffer_ptr +
                        (top_padding + recon_height) * stride + left_padding - 1;
    EB_MEMSET((uint8_t *)recon_samples_ptr, 0, sizeof(uint16_t) * (1 + recon_width + 1));

    // 3. zero out the left column
    recon_samples_ptr =
        (uint16_t *)recon_samples_buffer_ptr + top_padding * stride + left_padding - 1;
    for (sample_count = 0; sample_count < recon_height; sample_count++)
        recon_samples_ptr[sample_count * stride] = 0;
    // 4. zero out the right column
    recon_samples_ptr =
        (uint16_t *)recon_samples_buffer_ptr + top_padding * stride + left_padding + recon_width;
    for (sample_count = 0; sample_count < recon_height; sample_count++)
        recon_samples_ptr[sample_count * stride] = 0;
}

void initialize_samples_neighboring_reference_picture_8bit(EbByte   recon_samples_buffer_ptr,
                                                           uint16_t stride, uint16_t recon_width,
                                                           uint16_t recon_height,
                                                           uint16_t left_padding,
                                                           uint16_t top_padding) {
    uint8_t *recon_samples_ptr;
    uint16_t sample_count;

    // 1. zero out the top row
    recon_samples_ptr = recon_samples_buffer_ptr + (top_padding - 1) * stride + left_padding - 1;
    EB_MEMSET(recon_samples_ptr, 0, sizeof(uint8_t) * (1 + recon_width + 1));

    // 2. zero out the bottom row
    recon_samples_ptr =
        recon_samples_buffer_ptr + (top_padding + recon_height) * stride + left_padding - 1;
    EB_MEMSET(recon_samples_ptr, 0, sizeof(uint8_t) * (1 + recon_width + 1));

    // 3. zero out the left column
    recon_samples_ptr = recon_samples_buffer_ptr + top_padding * stride + left_padding - 1;
    for (sample_count = 0; sample_count < recon_height; sample_count++)
        recon_samples_ptr[sample_count * stride] = 0;
    // 4. zero out the right column
    recon_samples_ptr =
        recon_samples_buffer_ptr + top_padding * stride + left_padding + recon_width;
    for (sample_count = 0; sample_count < recon_height; sample_count++)
        recon_samples_ptr[sample_count * stride] = 0;
}

void initialize_samples_neighboring_reference_picture(
    EbReferenceObject *          reference_object,
    EbPictureBufferDescInitData *picture_buffer_desc_init_data_ptr, EbBitDepthEnum bit_depth) {
    if (bit_depth == EB_10BIT) {
        initialize_samples_neighboring_reference_picture16_bit(
            reference_object->reference_picture16bit->buffer_y,
            reference_object->reference_picture16bit->stride_y,
            reference_object->reference_picture16bit->width,
            reference_object->reference_picture16bit->height,
            picture_buffer_desc_init_data_ptr->left_padding,
            picture_buffer_desc_init_data_ptr->top_padding);

        initialize_samples_neighboring_reference_picture16_bit(
            reference_object->reference_picture16bit->buffer_cb,
            reference_object->reference_picture16bit->stride_cb,
            reference_object->reference_picture16bit->width >> 1,
            reference_object->reference_picture16bit->height >> 1,
            picture_buffer_desc_init_data_ptr->left_padding >> 1,
            picture_buffer_desc_init_data_ptr->top_padding >> 1);

        initialize_samples_neighboring_reference_picture16_bit(
            reference_object->reference_picture16bit->buffer_cr,
            reference_object->reference_picture16bit->stride_cr,
            reference_object->reference_picture16bit->width >> 1,
            reference_object->reference_picture16bit->height >> 1,
            picture_buffer_desc_init_data_ptr->left_padding >> 1,
            picture_buffer_desc_init_data_ptr->top_padding >> 1);
    } else {
        initialize_samples_neighboring_reference_picture_8bit(
            reference_object->reference_picture->buffer_y,
            reference_object->reference_picture->stride_y,
            reference_object->reference_picture->width,
            reference_object->reference_picture->height,
            picture_buffer_desc_init_data_ptr->left_padding,
            picture_buffer_desc_init_data_ptr->top_padding);

        initialize_samples_neighboring_reference_picture_8bit(
            reference_object->reference_picture->buffer_cb,
            reference_object->reference_picture->stride_cb,
            reference_object->reference_picture->width >> 1,
            reference_object->reference_picture->height >> 1,
            picture_buffer_desc_init_data_ptr->left_padding >> 1,
            picture_buffer_desc_init_data_ptr->top_padding >> 1);

        initialize_samples_neighboring_reference_picture_8bit(
            reference_object->reference_picture->buffer_cr,
            reference_object->reference_picture->stride_cr,
            reference_object->reference_picture->width >> 1,
            reference_object->reference_picture->height >> 1,
            picture_buffer_desc_init_data_ptr->left_padding >> 1,
            picture_buffer_desc_init_data_ptr->top_padding >> 1);
    }
}

static void eb_reference_object_dctor(EbPtr p) {
    EbReferenceObject *obj = (EbReferenceObject *)p;
    EB_DELETE(obj->reference_picture16bit);
    EB_DELETE(obj->reference_picture);
    EB_FREE_ALIGNED_ARRAY(obj->mvs);
    EB_DESTROY_MUTEX(obj->referenced_area_mutex);
}

/*****************************************
 * eb_picture_buffer_desc_ctor
 *  Initializes the Buffer Descriptor's
 *  values that are fixed for the life of
 *  the descriptor.
 *****************************************/
EbErrorType eb_reference_object_ctor(EbReferenceObject *reference_object,
                                     EbPtr              object_init_data_ptr) {
    EbPictureBufferDescInitData *picture_buffer_desc_init_data_ptr =
        (EbPictureBufferDescInitData *)object_init_data_ptr;
    EbPictureBufferDescInitData picture_buffer_desc_init_data_16bit_ptr =
        *picture_buffer_desc_init_data_ptr;

    reference_object->dctor = eb_reference_object_dctor;
    //TODO:12bit
    if (picture_buffer_desc_init_data_16bit_ptr.bit_depth == EB_10BIT) {
        // Hsan: set split_mode to 0 to construct the packed reference buffer (used @ EP)
        picture_buffer_desc_init_data_16bit_ptr.split_mode = EB_FALSE;
        EB_NEW(reference_object->reference_picture16bit,
               eb_picture_buffer_desc_ctor,
               (EbPtr)&picture_buffer_desc_init_data_16bit_ptr);

        initialize_samples_neighboring_reference_picture(
            reference_object,
            &picture_buffer_desc_init_data_16bit_ptr,
            picture_buffer_desc_init_data_16bit_ptr.bit_depth);

        // Hsan: set split_mode to 1 to construct the unpacked reference buffer (used @ MD)
        picture_buffer_desc_init_data_16bit_ptr.split_mode = EB_TRUE;
        EB_NEW(reference_object->reference_picture,
               eb_picture_buffer_desc_ctor,
               (EbPtr)&picture_buffer_desc_init_data_16bit_ptr);
    } else {
        // Hsan: set split_mode to 0 to as 8BIT input
        picture_buffer_desc_init_data_ptr->split_mode = EB_FALSE;
        EB_NEW(reference_object->reference_picture,
               eb_picture_buffer_desc_ctor,
               (EbPtr)picture_buffer_desc_init_data_ptr);

        initialize_samples_neighboring_reference_picture(
            reference_object,
            picture_buffer_desc_init_data_ptr,
            picture_buffer_desc_init_data_16bit_ptr.bit_depth);
    }
#if ENCDEC_16BIT
    picture_buffer_desc_init_data_16bit_ptr.split_mode = EB_FALSE;
    picture_buffer_desc_init_data_16bit_ptr.bit_depth  = EB_10BIT;
    EB_NEW(reference_object->reference_picture16bit,
           eb_picture_buffer_desc_ctor,
           (EbPtr)&picture_buffer_desc_init_data_16bit_ptr);
    picture_buffer_desc_init_data_16bit_ptr.bit_depth = EB_8BIT;
#endif
    if (picture_buffer_desc_init_data_ptr->mfmv) {
        //MFMV map is 8x8 based.
        uint32_t  mi_rows  = reference_object->reference_picture->height >> MI_SIZE_LOG2;
        uint32_t  mi_cols  = reference_object->reference_picture->width >> MI_SIZE_LOG2;
        const int mem_size = ((mi_rows + 1) >> 1) * ((mi_cols + 1) >> 1);
        EB_CALLOC_ALIGNED_ARRAY(reference_object->mvs, mem_size);
    }
    memset(&reference_object->film_grain_params, 0, sizeof(reference_object->film_grain_params));
    EB_CREATE_MUTEX(reference_object->referenced_area_mutex);
    return EB_ErrorNone;
}

EbErrorType eb_reference_object_creator(EbPtr *object_dbl_ptr, EbPtr object_init_data_ptr) {
    EbReferenceObject *obj;

    *object_dbl_ptr = NULL;
    EB_NEW(obj, eb_reference_object_ctor, object_init_data_ptr);
    *object_dbl_ptr = obj;

    return EB_ErrorNone;
}

static void eb_pa_reference_object_dctor(EbPtr p) {
    EbPaReferenceObject *obj = (EbPaReferenceObject *)p;
    EB_DELETE(obj->input_padded_picture_ptr);
    EB_DELETE(obj->quarter_decimated_picture_ptr);
    EB_DELETE(obj->sixteenth_decimated_picture_ptr);
    EB_DELETE(obj->quarter_filtered_picture_ptr);
    EB_DELETE(obj->sixteenth_filtered_picture_ptr);
}

/*****************************************
 * eb_pa_reference_object_ctor
 *  Initializes the Buffer Descriptor's
 *  values that are fixed for the life of
 *  the descriptor.
 *****************************************/
EbErrorType eb_pa_reference_object_ctor(EbPaReferenceObject *pa_ref_obj_,
                                        EbPtr                object_init_data_ptr) {
    EbPictureBufferDescInitData *picture_buffer_desc_init_data_ptr =
        (EbPictureBufferDescInitData *)object_init_data_ptr;

    pa_ref_obj_->dctor = eb_pa_reference_object_dctor;

    // Reference picture constructor
    EB_NEW(pa_ref_obj_->input_padded_picture_ptr,
           eb_picture_buffer_desc_ctor,
           (EbPtr)picture_buffer_desc_init_data_ptr);
    // Quarter Decim reference picture constructor
    EB_NEW(pa_ref_obj_->quarter_decimated_picture_ptr,
           eb_picture_buffer_desc_ctor,
           (EbPtr)(picture_buffer_desc_init_data_ptr + 1));
    EB_NEW(pa_ref_obj_->sixteenth_decimated_picture_ptr,
           eb_picture_buffer_desc_ctor,
           (EbPtr)(picture_buffer_desc_init_data_ptr + 2));
    // Quarter Filtered reference picture constructor
    if ((picture_buffer_desc_init_data_ptr + 1)->down_sampled_filtered) {
        EB_NEW(pa_ref_obj_->quarter_filtered_picture_ptr,
               eb_picture_buffer_desc_ctor,
               (EbPtr)(picture_buffer_desc_init_data_ptr + 1));
    }
    // Sixteenth Filtered reference picture constructor
    if ((picture_buffer_desc_init_data_ptr + 2)->down_sampled_filtered) {
        EB_NEW(pa_ref_obj_->sixteenth_filtered_picture_ptr,
               eb_picture_buffer_desc_ctor,
               (EbPtr)(picture_buffer_desc_init_data_ptr + 2));
    }

    return EB_ErrorNone;
}

EbErrorType eb_pa_reference_object_creator(EbPtr *object_dbl_ptr, EbPtr object_init_data_ptr) {
    EbPaReferenceObject *obj;

    *object_dbl_ptr = NULL;
    EB_NEW(obj, eb_pa_reference_object_ctor, object_init_data_ptr);
    *object_dbl_ptr = obj;

    return EB_ErrorNone;
}
