/*BSD 2-Clause License
* Copyright(c) 2019, Pekka Astola
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met :
*
* 1. Redistributions of source code must retain the above copyright notice, this
* list of conditions and the following disclaimer.
*
* 2. Redistributions in binary form must reproduce the above copyright notice,
* this list of conditions and the following disclaimer in the documentation
* and/or other materials provided with the distribution.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED.IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
* FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
* DAMAGES(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
*     SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
*     CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
*     OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
*     OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef RESIDUALJP2_HH
#define RESIDUALJP2_HH

#include <cstdint>

using std::int32_t;
using std::uint32_t;

using std::int16_t;
using std::uint16_t;

using std::int8_t;
using std::uint8_t;

#include "view.hh"

enum YUV_FORMAT {
    YUV444,
    YUV420,
    YUV400};

uint16_t *cropImage_for_HM(
    const uint16_t *input_image,
    const uint32_t nr,
    const uint32_t nc,
    const uint32_t ncomp,
    const uint32_t HORP,
    const uint32_t VERP);

std::vector<uint16_t> upscale(
    const std::vector<uint16_t> &input,
    const int32_t nr,
    const int32_t nc,
    const int32_t rz);

std::vector<std::vector<uint16_t>> convertYUV400seqTo444(
    const std::vector<std::vector<uint16_t>> &YUV400,
    const int32_t nr,
    const int32_t nc,
    const int nframes);

std::vector<std::vector<uint16_t>> convertYUVseqTo444(
    const char *inputYUV,
    YUV_FORMAT input_yuv,
    const int32_t nr,
    const int32_t nc,
    const int nframes);

std::vector<std::vector<uint16_t>> convertYUV420seqTo444(
    const std::vector<std::vector<uint16_t>> &YUV420,
    const int32_t nr,
    const int32_t nc,
    const int nframes);

std::vector<std::vector<uint16_t>>  readYUV444_seq_from_disk(
    const char *input_444,
    const int32_t nframes,
    const int32_t nr,
    const int32_t nc);

std::vector<std::vector<uint16_t>>  readYUV420_seq_from_disk(
    const char *input_420,
    const int32_t nframes,
    const int32_t nr,
    const int32_t nc);

std::vector<std::vector<uint16_t>> readYUV400_seq_from_disk(
    const char *input_400,
    const int32_t nframes,
    const int32_t nr,
    const int32_t nc);

std::vector<uint16_t> padArrayUint16_t_for_HM(
    const uint16_t *input_image,
    const uint32_t nr,
    const uint32_t nc,
    const uint32_t ncomp,
    const uint32_t HORP,
    const uint32_t VERP);

void writeYUV444_seq_to_disk(
    const std::vector< std::vector<uint16_t>> &YUV_444_SEQ,
    const char *output_444);

int32_t decodeHM(
    const char *input_hevc,
    const char *outputYUV);

long encodeHM(
    const char *input444,
    const char *output_hevc,
    YUV_FORMAT yuvformat,
    const int32_t QP,
    const int32_t nframes,
    const int32_t nr,
    const int32_t nc,
    const char *outputYUV);

void getJP2Header(
    uint8_t *JP2, 
    uint8_t *&header,
    int32_t JP2Size,
    int32_t &headerSize);

int32_t getJP2DictionaryIndex(
    uint8_t *JP2header,
    int32_t headerSize,
    std::vector<std::vector<unsigned char>> JP2_dict);

void readResidualFromDisk(
    const char *jp2_residual_path_jp2,
    int32_t &n_bytes_residual,
    FILE *input_LF,
    std::vector<std::vector<unsigned char>> &JP2_dict);

void updateJP2Dictionary(
    std::vector<std::vector<unsigned char>> &JP2_dict,
    uint8_t *header,
    int32_t headerSize);

void writeResidualToDisk(
    const char *jp2_residual_path_jp2,
    FILE *output_LF_file, 
    int32_t &n_bytes_residual,
    std::vector<std::vector<unsigned char>> &JP2_dict);

char *kakadu_oparams(
    const double rate,
    const std::string colorspace);

char *kakadu_cparams(
    const double *cweights,
    const int32_t ncomp);

int32_t encodeKakadu(
    const char *ppm_pgm_input_path,
    const char *kdu_compress_path,
    const char *jp2_output_path,
    const char *encoding_parameters,
    const double rate);

int32_t decodeKakadu(
    const char *ppm_pgm_output_path,
    const char *kdu_expand_path,
    const char *jp2_input_path);

double* get_residual(
    const uint16_t *original,
    const uint16_t *prediction,
    const int32_t nr,
    const int32_t nc,
    const int32_t ncomp);

uint16_t* quantize_residual(
    const double *residual,
    const int32_t nr,
    const int32_t nc,
    const int32_t ncomp,
    const int32_t bpc,
    const int32_t Q_i,
    const int32_t offset_i);

double* dequantize_residual(
    const uint16_t *qresidual,
    const int32_t nr,
    const int32_t nc,
    const int32_t ncomp,
    const int32_t bpc,
    const int32_t Q_i,
    const int32_t offset_i);

uint16_t *apply_residual(
    const uint16_t *prediction,
    const double *residual,
    const int32_t nr,
    const int32_t nc,
    const int32_t ncomp,
    const int32_t bpc);

uint16_t* decode_residual_JP2(
    const char *ppm_pgm_output_path,
    const char *kdu_expand_path,
    const char *jp2_input_path);

void encode_residual_JP2(
    const char *ppm_pgm_input_path,
    const char *kdu_compress_path,
    const char *jp2_output_path,
    const char *encoding_parameters,
    const double rate);

#endif
