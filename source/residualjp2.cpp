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

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include "residualjp2.hh"
#include "ycbcr.hh"
#include "ppm.hh"
#include "fileaux.hh"
#include "clip.hh"
#include "medianfilter.hh"

#define USE_JP2_DICTIONARY 1

std::vector<uint16_t> padArrayUint16_t_for_HM(
    const uint16_t *input_image,
    const uint32_t nr,
    const uint32_t nc,
    const uint32_t ncomp,
    const uint32_t HORP,
    const uint32_t VERP) {

    int32_t nr_padded = nr + VERP;
    int32_t nc_padded = nc + HORP;

    std::vector<uint16_t> input_image_padded(nr_padded*nc_padded*ncomp,0);

    for (uint32_t icomp = 0; icomp < ncomp; icomp++){
        for (uint32_t ir = 0; ir < nr_padded - VERP; ir++) {
            for (uint32_t ic = 0; ic < nc_padded - HORP; ic++) {

                int32_t offset = 
                    ir + nr * ic + nr*nc*icomp;

                int32_t offset_padded = 
                    ir + nr_padded * ic + nr_padded*nc_padded*icomp;

                *(input_image_padded.data() + offset_padded) =
                    *(input_image + offset);

            }
        }
    }

    /*bottom borders*/
    for (uint32_t icomp = 0; icomp < ncomp; icomp++) {
        for (uint32_t ir = nr_padded - VERP; ir < nr_padded; ir++) {
            for (uint32_t ic = HORP; ic < nc_padded - HORP; ic++) {

                int32_t offset_from = nr_padded - VERP - 1 + nr_padded * ic + nr_padded*nc_padded*icomp;
                int32_t offset_to = ir + nr_padded * ic + nr_padded*nc_padded*icomp;

                *(input_image_padded.data() + offset_to) =
                    *(input_image_padded.data() + offset_from);

            }
        }
    }

    /*right borders*/
    for (uint32_t icomp = 0; icomp < ncomp; icomp++) {
        for (uint32_t ir = 0; ir < nr_padded; ir++) {
            for (uint32_t ic = nc_padded - HORP; ic < nc_padded; ic++) {

                int32_t offset_from = ir + nr_padded * (nc_padded - HORP - 1) + nr_padded*nc_padded*icomp;
                int32_t offset_to = ir + nr_padded * ic + nr_padded*nc_padded*icomp;

                *(input_image_padded.data() + offset_to) =
                    *(input_image_padded.data() + offset_from);

            }
        }
    }

    return input_image_padded;

}

void writeYUV444_seq_to_disk(
    const std::vector< std::vector<uint16_t>> &YUV_444_SEQ,
    const char *output_444) {

    aux_ensure_directory(output_444);

    FILE *output_444_file = fopen(output_444, "wb");

    for (uint32_t fr = 0; fr < YUV_444_SEQ.size(); fr++) {

        std::vector<uint16_t> Y;
        std::vector<uint16_t> U;
        std::vector<uint16_t> V;

        int32_t nt = YUV_444_SEQ.at(fr).size();
        int32_t np = nt / 3; /*we assume that input format is 3-channel*/

        /*the image will be transposed, 
        since we use column major and YUV format is row major*/
        for (int32_t ii = 0; ii < np; ii++) {

            Y.push_back(YUV_444_SEQ.at(fr).at(ii));
            U.push_back(YUV_444_SEQ.at(fr).at(ii+np));
            V.push_back(YUV_444_SEQ.at(fr).at(ii+2*np));

        }

        fwrite(Y.data(), sizeof(uint16_t), np, output_444_file);
        fwrite(U.data(), sizeof(uint16_t), np, output_444_file);
        fwrite(V.data(), sizeof(uint16_t), np, output_444_file);

    }

    fclose(output_444_file);

}

std::vector<uint16_t> readYUV444_seq_from_disk(
    const char *input_444,
    const int32_t nframes,
    const int32_t nr,
    const int32_t nc,
    const int32_t ncomp) {

    FILE *input_444_file = fopen(input_444, "rb");

    const int32_t np = nframes*nr*nc*ncomp;

    std::vector<uint16_t> YUV444_seq( np,0);

    fread(
        YUV444_seq.data(),
        sizeof(uint16_t),
        np,
        input_444_file);

    fclose(input_444_file);

    return YUV444_seq;

}

std::vector<uint16_t> readYUV420_seq_from_disk(
    const char *input_420,
    const int32_t nframes,
    const int32_t nr,
    const int32_t nc,
    const int32_t ncomp) {

    FILE *input_420_file = fopen(input_420, "rb");

    const int32_t np = nframes*(nr*nc + nr / 2 * nc / 2);

    std::vector<uint16_t> YUV420_seq(np, 0);

    fread(
        YUV420_seq.data(),
        sizeof(uint16_t),
        np,
        input_420_file);

    fclose(input_420_file);

    return YUV420_seq;

}

std::vector<uint16_t> readYUV400_seq_from_disk(
    const char *input_400,
    const int32_t nframes,
    const int32_t nr,
    const int32_t nc,
    const int32_t ncomp) {

    FILE *input_400_file = fopen(input_420, "rb");

    const int32_t np = nframes*(nr*nc);

    std::vector<uint16_t> YUV400_seq(np, 0);

    fread(
        YUV400_seq.data(),
        sizeof(uint16_t),
        np,
        input_400_file);

    fclose(input_400_file);

    return YUV400_seq;

}

std::vector<std::vector<uint16_t>> convertYUVseq(
    YUV_FORMAT input_yuv,
    YUV_FORMAT input_yuv,
    const std::vector<uint16_t> inputSEQ,
    const int32_t nr,
    const int32_t nc,
    const int nframes) {



}


int32_t encodeHM(
    const char *input444,
    const char *output_hevc,
    YUV_FORMAT yuvformat,
    const int32_t QP,
    const int32_t nframes,
    const int32_t nr,
    const int32_t nc,
    const char *outputYUV) {


    aux_ensure_directory(input444);
    aux_ensure_directory(output_hevc);
    aux_ensure_directory(outputYUV);

    const char *input_cfg = 
        "C:/Local/astolap/Data/JPEG_PLENO_2019/SAN_DIEGO/Matlab/hm_inter.cfg";

    const char *hm_encoder =
        "C:/Local/astolap/Data/JPEG_PLENO_2019/BRUSSELS/HEVC-HM/bin/vc2015/x64/Release/TAppEncoder.exe";

    std::string yuvformatstr;

    if (yuvformat == YUV400) {
        yuvformatstr = "400";
    }
    if (yuvformat == YUV420) {
        yuvformatstr = "420";
    }
    if (yuvformat == YUV444) {
        yuvformatstr = "444";
    }

    char hm_call[2048];

    sprintf(hm_call,
        "%s"
        " -c %s"
        " -i %s"
        " -o %s"
        " -fr %d"
        " -wdt %d"
        " -hgt %d"
        " -b %s"
        " --FramesToBeEncoded=%d"
        " --QP=%d"
        " --ChromaFormatIDC=%s",
        hm_encoder,
        input_cfg,
        input444,
        outputYUV,
        1,
        nc,
        nr,
        output_hevc,
        nframes,
        QP,
        yuvformatstr.c_str());

    return system_1(hm_call);

}


int32_t decodeHM(
    const char *input_hevc,
    const char *outputYUV) {

    aux_ensure_directory(outputYUV);

    const char *hm_encoder =
        "C:/Local/astolap/Data/JPEG_PLENO_2019/BRUSSELS/HEVC-HM/bin/vc2015/x64/Release/TAppDecoder.exe";

    char hm_call[2048];

    sprintf(hm_call,
        "%s"
        " -b %s"
        " -o %s",
        hm_encoder,
        input_hevc,
        outputYUV);

    return system_1(hm_call);

}

void getJP2Header(
    uint8_t *JP2, 
    uint8_t *&header, 
    int32_t JP2Size,
    int32_t &headerSize) {

  for (int32_t ii = 0; ii < JP2Size - 1; ii++) {
    if ((uint16_t) (JP2[ii] << 8 | JP2[ii + 1]) == 0xFF90) { /*we have first tile*/
      headerSize = ii + 2;
      header = new uint8_t[headerSize];
      memcpy(header, JP2, headerSize);
      return;
    }
  }
  return;
}

int32_t getJP2DictionaryIndex(
    uint8_t *header, 
    int32_t headerSize,
    std::vector<std::vector<uint8_t>> JP2_dict) {

  for (int32_t ii = 0; ii < JP2_dict.size(); ii++) {

    int32_t L = (int32_t) JP2_dict.at(ii).size();

    if (L == headerSize) {
      if (memcmp(header, &JP2_dict.at(ii)[0], L) == 0) {
        return ii;
      }
    }

  }

  return -1;
}

void updateJP2Dictionary(
    std::vector<std::vector<uint8_t>> &JP2_dict,
    uint8_t *JP2header, 
    int32_t headerSize) {

  std::vector<uint8_t> new_dict_element;

  for (int32_t ii = 0; ii < headerSize; ii++) {
    new_dict_element.push_back(JP2header[ii]);
  }

  JP2_dict.push_back(new_dict_element);

}

void readResidualFromDisk(
    const char *jp2_residual_path_jp2,
    int32_t &n_bytes_residual, 
    FILE *input_LF,
    std::vector<std::vector<uint8_t>> &JP2_dict) {

  int32_t n_bytes_JP2 = 0;
  uint8_t *jp2_residual = 0;

  if (USE_JP2_DICTIONARY) {

    int32_t dict_index = 0, headerSize = 0;
    uint8_t dict_index_char = 0;

    n_bytes_residual += (int32_t) fread(&dict_index_char, sizeof(uint8_t), 1,
                                    input_LF) * sizeof(uint8_t);

    dict_index = (int32_t) dict_index_char;

    if (JP2_dict.size() == 0 || dict_index > (int32_t) JP2_dict.size() - 1) {

      n_bytes_residual += (int32_t) fread(&headerSize, sizeof(int32_t), 1, input_LF)
          * sizeof(int32_t);

      uint8_t *JP2header = new uint8_t[headerSize]();

      n_bytes_residual += (int32_t) fread(JP2header, sizeof(uint8_t),
                                      headerSize, input_LF)
          * sizeof(uint8_t);

      updateJP2Dictionary(JP2_dict, JP2header, headerSize);

      delete[] (JP2header);

    }

    uint8_t *jp2_residual_tmp;

    n_bytes_residual += (int32_t) fread(&n_bytes_JP2, sizeof(int32_t), 1, input_LF)
        * sizeof(int32_t);
    jp2_residual_tmp = new uint8_t[n_bytes_JP2]();
    n_bytes_residual += (int32_t) fread(jp2_residual_tmp, sizeof(uint8_t),
                                    n_bytes_JP2, input_LF);

    headerSize = (int32_t) JP2_dict.at(dict_index).size();
    jp2_residual = new uint8_t[n_bytes_JP2 + headerSize]();

    memcpy(jp2_residual, &JP2_dict.at(dict_index)[0], headerSize);
    memcpy(jp2_residual + headerSize, jp2_residual_tmp, n_bytes_JP2);

    n_bytes_JP2 += headerSize;

    delete[] (jp2_residual_tmp);

  } else {

    n_bytes_residual += (int32_t) fread(&n_bytes_JP2, sizeof(int32_t), 1, input_LF)
        * sizeof(int32_t);
    jp2_residual = new uint8_t[n_bytes_JP2]();
    n_bytes_residual += (int32_t) fread(jp2_residual, sizeof(uint8_t),
                                    n_bytes_JP2, input_LF)
        * sizeof(uint8_t);

  }

  FILE *jp2_res_file;
  jp2_res_file = fopen(jp2_residual_path_jp2, "wb");
  fwrite(jp2_residual, sizeof(uint8_t), n_bytes_JP2, jp2_res_file);
  fclose(jp2_res_file);

  delete[] (jp2_residual);
}

void writeResidualToDisk(
    const char *jp2_residual_path_jp2,
    FILE *output_LF_file, 
    int32_t &n_bytes_residual,
    std::vector<std::vector<uint8_t>> &JP2_dict) {

  int32_t n_bytes_JP2 = aux_GetFileSize(jp2_residual_path_jp2);

  uint8_t *jp2_residual = new uint8_t[n_bytes_JP2]();
  FILE *jp2_color_residual_file = fopen(jp2_residual_path_jp2, "rb");
  fread(jp2_residual, sizeof(uint8_t), n_bytes_JP2,
        jp2_color_residual_file);
  fclose(jp2_color_residual_file);

  if (USE_JP2_DICTIONARY) {
    bool updateDictionary = false;
    /* get header */
    uint8_t *JP2header;
    int32_t headerSize = 0;
    getJP2Header(jp2_residual, JP2header, n_bytes_JP2, headerSize);
    /* get index in dictionary */
    int32_t dict_index = getJP2DictionaryIndex(JP2header, headerSize, JP2_dict);
    /* update dictionary if needed */
    if (dict_index == -1) {
      updateDictionary = true;
      updateJP2Dictionary(JP2_dict, JP2header, headerSize);
      dict_index = (int32_t) JP2_dict.size() - 1;
      //printf("Dictonary update, index=%i\n", dict_index);
    }

    //printf("Using dictionary index:\t%i\n", dict_index);

    delete[] (JP2header);

    /* write index of dictionary to bitstream */
    uint8_t dict_index_char = (uint8_t) dict_index;
    n_bytes_residual += (int32_t) fwrite(&dict_index_char, sizeof(uint8_t), 1,
                                     output_LF_file) * sizeof(uint8_t);

    if (updateDictionary) { /* add update information if necessary */
      n_bytes_residual += (int32_t) fwrite(&headerSize, sizeof(int32_t), 1,
                                       output_LF_file) * sizeof(int32_t);
      n_bytes_residual += (int32_t) fwrite(&JP2_dict.at(dict_index)[0],
                                       sizeof(uint8_t), headerSize,
                                       output_LF_file) * sizeof(uint8_t);
    }

    /* write to codestream without header */

    n_bytes_JP2 = n_bytes_JP2 - headerSize;

    n_bytes_residual += (int32_t) fwrite(&n_bytes_JP2, sizeof(int32_t), 1,
                                     output_LF_file) * sizeof(int32_t);
    n_bytes_residual += (int32_t) fwrite(jp2_residual + headerSize,
                                     sizeof(uint8_t), n_bytes_JP2,
                                     output_LF_file) * sizeof(uint8_t);

  } else {

    n_bytes_residual += (int32_t) fwrite(&n_bytes_JP2, sizeof(int32_t), 1,
                                     output_LF_file) * sizeof(int32_t);
    n_bytes_residual += (int32_t) fwrite(jp2_residual, sizeof(uint8_t),
                                     n_bytes_JP2, output_LF_file)
        * sizeof(uint8_t);

    /*n_bytes_JP2 = n_bytes_JP2 - jp2_headersize;
     n_bytes_residual += (int32_t)fwrite(&n_bytes_JP2, sizeof(int32_t), 1, output_LF_file) * sizeof(int32_t);
     n_bytes_residual += (int32_t)fwrite(jp2_residual + jp2_headersize, sizeof(uint8_t), n_bytes_JP2, output_LF_file) * sizeof(uint8_t);*/

  }

  delete[] (jp2_residual);

}

char *kakadu_oparams(
    const double rate,
    const std::string colorspace) {

    char *oparams = new char[65536]();

    std::string cycc;

    if (colorspace.compare("RGB")==0) {
        cycc = "yes";
    }
    if (colorspace.compare("YCbCr")==0) {
        cycc = "no";
    }

    sprintf(
        oparams, 
        " -num_threads 0"
        " -no_weights"
        " -full"
        " -precise"
        " -no_info"
        " Clevels=6"
        " -rate %2.3f"
        " Cycc=%s ",
        rate,
        cycc.c_str());

    return oparams;

}

char *kakadu_cparams(
    const double *cweights,
    const int32_t ncomp) {

    char *cparams = new char[65536]();

    int32_t s_index = 0;

    for (int32_t icomp = 0; icomp < ncomp; icomp++) {
        s_index += sprintf(
            cparams + s_index, 
            "%s%1d%s%2.2f",
            " Cweight:T0C",
            icomp,
            "=",
            cweights[icomp]);
    }

    return cparams;
}

int32_t encodeKakadu(
    const char *ppm_pgm_input_path,
    const char *kdu_compress_path,
    const char *jp2_output_path,
    const char *encoding_parameters,
    const double rate) {

    aux_ensure_directory(jp2_output_path);

    char kdu_compress_s[1024];

    sprintf(
        kdu_compress_s,
        "\"%s\"%s%s%s%s%s",
        kdu_compress_path,
        " -i ",
        ppm_pgm_input_path,
        " -o ",
        jp2_output_path,
        encoding_parameters);

    return system_1(kdu_compress_s);

}

int32_t decodeKakadu(
    const char *ppm_pgm_output_path,
    const char *kdu_expand_path,
    const char *jp2_input_path) {

    aux_ensure_directory(jp2_input_path);
    aux_ensure_directory(ppm_pgm_output_path);

    char kdu_expand_s[1024];

    sprintf(
        kdu_expand_s,
        "\"%s\"%s%s%s%s", 
        kdu_expand_path, 
        " -precise"
        " -num_threads 0"
        " -i ",
        jp2_input_path, 
        " -o ", 
        ppm_pgm_output_path);

    return system_1(kdu_expand_s);

}

double* get_residual(
    const uint16_t *original,
    const uint16_t *prediction,
    const int32_t nr,
    const int32_t nc,
    const int32_t ncomp) {

    double *residual = new double[nr*nc*ncomp]();

    for (int32_t ii = 0; ii < nr*nc*ncomp; ii++) {
        double o = static_cast<double>(original[ii]);
        double p = static_cast<double>(prediction[ii]);
        residual[ii] = o - p;
    }

    return residual;

}

uint16_t* quantize_residual(
    const double *residual,
    const int32_t nr,
    const int32_t nc,
    const int32_t ncomp,
    const int32_t bpc,
    const int32_t Q_i,
    const int32_t offset_i) {

    double Q = static_cast<double>(Q_i);
    double offset = static_cast<double>(offset_i);

    uint16_t *qresidual = new uint16_t[nr*nc*ncomp]();

    uint16_t maxval = static_cast<uint16_t>( (1 << bpc) -1);

    for (int32_t ii = 0; ii < nr*nc*ncomp; ii++) {
        double qpred = (residual[ii] + offset) / Q;
        qpred = floor(qpred + 0.5);
        qpred = qpred > maxval ? maxval  : qpred;
        qpred = qpred < 0.0 ? 0.0 : qpred;
        qresidual[ii] = static_cast<uint16_t>(qpred);
    }

    return qresidual;
}

double* dequantize_residual(
    const uint16_t *qresidual,
    const int32_t nr,
    const int32_t nc,
    const int32_t ncomp,
    const int32_t bpc,
    const int32_t Q_i,
    const int32_t offset_i) {

    double Q = static_cast<double>(Q_i);
    double offset = static_cast<double>(offset_i);

    double *residual = new double[nr*nc*ncomp]();

    for (int32_t ii = 0; ii < nr*nc*ncomp; ii++) {

        double qres = static_cast<double>(qresidual[ii]);
        residual[ii] = qres * Q - offset;

    }

    return residual;

}

uint16_t *apply_residual(
    const uint16_t *prediction,
    const double *residual,
    const int32_t nr,
    const int32_t nc,
    const int32_t ncomp,
    const int32_t bpc) {

    uint16_t *corrected = new uint16_t[nr*nc*ncomp]();

    uint16_t maxval = static_cast<uint16_t>((1 << bpc) - 1);

    for (int32_t ii = 0; ii < nr*nc*ncomp; ii++) {

        double pred = static_cast<double>(prediction[ii]);
        double corrected_d = pred + residual[ii];

        corrected_d = floor(corrected_d + 0.5);

        corrected_d = corrected_d > maxval ? maxval : corrected_d;
        corrected_d = corrected_d < 0.0 ? 0.0 : corrected_d;

        corrected[ii] = static_cast<uint16_t>(corrected_d);
       
    }

    return corrected;
}

uint16_t* decode_residual_JP2(
    const char *ppm_pgm_output_path,
    const char *kdu_expand_path,
    const char *jp2_input_path) {

    uint16_t *residual_image_decoded = nullptr;

    decodeKakadu(
        ppm_pgm_output_path,
        kdu_expand_path,
        jp2_input_path);

    int32_t nr1, nc1, ncomp1;

    aux_read16PGMPPM(
        ppm_pgm_output_path,
        nc1,
        nr1,
        ncomp1,
        residual_image_decoded);

    return residual_image_decoded;
}

void encode_residual_JP2(
    const char *ppm_pgm_input_path,
    const char *kdu_compress_path,
    const char *jp2_output_path,
    const char *encoding_parameters,
    const double rate) {

    int32_t status = encodeKakadu(
        ppm_pgm_input_path,
        kdu_compress_path,
        jp2_output_path,
        encoding_parameters,
        rate);

}