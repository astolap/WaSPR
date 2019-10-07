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

#include "codestream.hh"
#include "view.hh"
#include "minconf.hh"
#include "fileaux.hh"

#include <iostream>
#include <vector>

viewParametersConstruct::viewParametersConstruct(
    view *LF,
    const int32_t nviews,
    const std::string gzippath,
    const std::string rawfile,
    const std::string mode) : LF(LF), nviews(nviews), gzippath(gzippath), rawfile(rawfile), mode(mode) {

    if (!mode.compare("encode")) {

        convertLFtoBytes();

        FILE *outputfile;
        outputfile = fopen(rawfile.c_str(), "wb");
        fwrite(rawbytes.data(), sizeof(uint8_t), rawbytes.size(), outputfile);
        fclose(outputfile);

        char gzipcall[2048];
        sprintf(gzipcall,
            "%s -f -9 %s", gzippath.c_str(), rawfile.c_str());
        int32_t status = system_1(gzipcall);

    }
    else if (!mode.compare("decode")) {

        char gzipcall[2048];
        sprintf(gzipcall,
            "%s -f -d %s", gzippath.c_str(), (rawfile + ".gz").c_str());
        int32_t status = system_1(gzipcall);

        long nbytes = aux_GetFileSize(rawfile);

        rawbytes = std::vector<uint8_t>(nbytes, 0);

        FILE *inputfile;
        inputfile = fopen(rawfile.c_str(), "rb");
        fread(rawbytes.data(), sizeof(uint8_t), rawbytes.size(), inputfile);
        fclose(inputfile);

        convertBytesToLF();

    }
    else {
        exit(0);
    }


}

viewParametersConstruct::~viewParametersConstruct() {}

void viewParametersConstruct::addMconfs() {
    for (int32_t ii = 0; ii < nviews; ii++) {

        view *SAI = LF + ii;

        minimal_config mconf = makeMinimalConfig(SAI);

        addRawBytes(mconf);

    }
}

void viewParametersConstruct::getMconfs() {
    for (int32_t ii = 0; ii < nviews; ii++) {

        view *SAI = LF + ii;

        minimal_config mconf;

        getRawBytes(mconf);

        setup_form_minimal_config(&mconf, SAI);
    }
}


void viewParametersConstruct::getX()
{
    for (int32_t ii = 0; ii < nviews; ii++) {

        view *SAI = LF + ii;
        if (SAI->has_x_displacement) {
            getRawBytes(SAI->x);
        }
    }
}

void viewParametersConstruct::addX()
{
    for (int32_t ii = 0; ii < nviews; ii++) {

        view *SAI = LF + ii;
        if (SAI->has_x_displacement) {
            addRawBytes(SAI->x);
        }
    }
}

void viewParametersConstruct::getY()
{
    for (int32_t ii = 0; ii < nviews; ii++) {

        view *SAI = LF + ii;

        if (SAI->has_y_displacement) {
            getRawBytes(SAI->y);
        }
    }
}


void viewParametersConstruct::addY()
{
    for (int32_t ii = 0; ii < nviews; ii++) {

        view *SAI = LF + ii;

        if (SAI->has_y_displacement) {
            addRawBytes(SAI->y);
        }
    }
}

void viewParametersConstruct::getNRefs() {
    for (int32_t ii = 0; ii < nviews; ii++) {
        view *SAI = LF + ii;
        if (SAI->has_color_references) {
            getRawBytes(SAI->n_references);
        }
    }
}


void viewParametersConstruct::addNRefs() {
    for (int32_t ii = 0; ii < nviews; ii++) {
        view *SAI = LF + ii;
        if (SAI->has_color_references) {
            addRawBytes(SAI->n_references);
        }
    }
}

void viewParametersConstruct::getRefs() {
    for (int32_t ii = 0; ii < nviews; ii++) {
        view *SAI = LF + ii;
        if (SAI->has_color_references) {
            delete[](SAI->references);
            SAI->references = new int32_t[SAI->n_references]();
            for (int32_t ij = 0; ij < SAI->n_references; ij++) {
                uint16_t nid;
                getRawBytes(nid);
                *(SAI->references + ij) = nid;
            }
        }
    }
}

void viewParametersConstruct::addRefs() {
    for (int32_t ii = 0; ii < nviews; ii++) {
        view *SAI = LF + ii;
        if (SAI->has_color_references) {
            for (int32_t ij = 0; ij < SAI->n_references; ij++) {
                uint16_t nid = (uint16_t) *(SAI->references + ij);
                addRawBytes(nid);
            }
        }
    }
}

void viewParametersConstruct::getNNDRefs() {
    for (int32_t ii = 0; ii < nviews; ii++) {
        view *SAI = LF + ii;
        if (SAI->has_depth_references) {
            getRawBytes(SAI->n_depth_references);
        }
    }
}

void viewParametersConstruct::addNNDRefs() {
    for (int32_t ii = 0; ii < nviews; ii++) {
        view *SAI = LF + ii;
        if (SAI->has_depth_references) {
            addRawBytes(SAI->n_depth_references);
        }
    }
}

void viewParametersConstruct::getNDRefs() {
    for (int32_t ii = 0; ii < nviews; ii++) {

        view *SAI = LF + ii;
        if (SAI->has_depth_references) {

            delete[](SAI->depth_references);
            SAI->depth_references = new int32_t[SAI->n_depth_references]();

            for (int32_t ij = 0; ij < SAI->n_depth_references; ij++) {
                uint16_t nid;
                getRawBytes(nid);
                *(SAI->depth_references + ij) = nid;
            }
        }
    }
}

void viewParametersConstruct::addNDRefs() {
    for (int32_t ii = 0; ii < nviews; ii++) {

        view *SAI = LF + ii;
        if (SAI->has_depth_references) {
            for (int32_t ij = 0; ij < SAI->n_depth_references; ij++) {
                uint16_t nid = (uint16_t) *(SAI->depth_references + ij);
                addRawBytes(nid);
            }
        }
    }
}

void viewParametersConstruct::getLSW() {
    for (int32_t ii = 0; ii < nviews; ii++) {

        view *SAI = LF + ii;

        if (SAI->has_color_references) {
            if (SAI->mmode == 0) {


                int32_t MMM = (1 << SAI->n_references);
                int32_t N_LS = SAI->ncomp*((MMM*SAI->n_references) / 2);
                int32_t N_LS_R = SAI->nc_merge*((MMM*SAI->n_references) / 2);

                delete[](SAI->merge_weights);
                SAI->merge_weights = new int16_t[N_LS]();

                for (int32_t ij = 0; ij < N_LS_R; ij++)
                {
                    uint16_t tmp;
                    getRawBytes(tmp);
                    *(SAI->merge_weights + ij) = tmp;
                }
            }
        }
    }
}

void viewParametersConstruct::addLSW() {
    for (int32_t ii = 0; ii < nviews; ii++) {

        view *SAI = LF + ii;

        if (SAI->has_color_references) {
            if (SAI->mmode == 0) {

                int32_t MMM = (1 << SAI->n_references);
                int32_t N_LS_R = SAI->nc_merge*((MMM*SAI->n_references) / 2);

                for (int32_t ij = 0; ij < N_LS_R; ij++)
                {
                    addRawBytes(*(SAI->merge_weights + ij));
                }
            }
        }
    }
}

void viewParametersConstruct::getSTD() {
    for (int32_t ii = 0; ii < nviews; ii++) {

        view *SAI = LF + ii;

        if (SAI->has_color_references) {
            if (SAI->mmode == 1) {
                getRawBytes(SAI->stdd);
            }
        }
    }
}

void viewParametersConstruct::addSTD() {
    for (int32_t ii = 0; ii < nviews; ii++) {

        view *SAI = LF + ii;

        if (SAI->has_color_references) {
            if (SAI->mmode == 1) {
                addRawBytes(SAI->stdd);
            }
        }
    }
}

void viewParametersConstruct::getSPP() {
    for (int32_t ii = 0; ii < nviews; ii++) {

        view *SAI = LF + ii;

        if (SAI->use_global_sparse) {

            uint8_t tmpNNt;
            uint8_t tmpMs;

            getRawBytes(tmpNNt);
            getRawBytes(tmpMs);
            getRawBytes(SAI->number_of_sp_filters);

            SAI->NNt = tmpNNt;
            SAI->Ms = tmpMs;
        }
    }
}

void viewParametersConstruct::addSPP() {
    for (int32_t ii = 0; ii < nviews; ii++) {

        view *SAI = LF + ii;

        if (SAI->use_global_sparse) {

            uint8_t tmpNNt = (uint8_t)SAI->NNt;
            uint8_t tmpMs = (uint8_t)SAI->Ms;

            addRawBytes(tmpNNt);
            addRawBytes(tmpMs);
            addRawBytes(SAI->number_of_sp_filters);

        }
    }
}

void viewParametersConstruct::getSPW() {
    for (int32_t ii = 0; ii < nviews; ii++) {
        view *SAI = LF + ii;
        if (SAI->use_global_sparse) {

            SAI->sparse_filters.clear();

            for (int32_t ee = 0; ee < SAI->number_of_sp_filters; ee++) {

                spfilter tmpsp;
                SAI->sparse_filters.push_back(tmpsp);

                SAI->sparse_filters.at(ee).Ms = SAI->Ms;
                SAI->sparse_filters.at(ee).NNt = SAI->NNt;

                SAI->sparse_filters.at(ee).MT = SAI->SP_B>0 ?
                    (SAI->n_references + 1)*(SAI->NNt * 2 + 1) * (SAI->NNt * 2 + 1) + 1 :
                    (SAI->NNt * 2 + 1) * (SAI->NNt * 2 + 1) + 1;

                SAI->sparse_filters.at(ee).quantized_filter_coefficients.clear();

                for (int32_t ii = 0; ii < SAI->Ms; ii++) {
                    SAI->sparse_filters.at(ee).quantized_filter_coefficients.push_back(0);
                }
            }

            for (int32_t ee = 0; ee < SAI->number_of_sp_filters; ee++) {
                for (int32_t ij = 0; ij < SAI->Ms; ij++) {
                    getRawBytes(SAI->sparse_filters.at(ee).quantized_filter_coefficients.at(ij));
                }
            }
        }
    }
}
void viewParametersConstruct::addSPW() {
    for (int32_t ii = 0; ii < nviews; ii++) {

        view *SAI = LF + ii;

        if (SAI->use_global_sparse) {
            for (int32_t ee = 0; ee < SAI->number_of_sp_filters; ee++) {
                for (int32_t ij = 0; ij < SAI->Ms; ij++) {
                    addRawBytes(SAI->sparse_filters.at(ee).quantized_filter_coefficients.at(ij));
                }
            }
        }
    }
}

void viewParametersConstruct::getSPM() {
    for (int32_t ii = 0; ii < nviews; ii++) {

        view *SAI = LF + ii;

        int32_t Nsp = SAI->SP_B > 0 ?
            (SAI->n_references + 1)*(SAI->NNt * 2 + 1) * (SAI->NNt * 2 + 1) + 1 :
            (SAI->NNt * 2 + 1) * (SAI->NNt * 2 + 1) + 1;

        int32_t sp_mask_nbytes = (Nsp % 8) ? Nsp / 8 + 1 : Nsp / 8;

        uint8_t *sparsemask = new uint8_t[sp_mask_nbytes*SAI->number_of_sp_filters]();

        for (int32_t ij = 0; ij < sp_mask_nbytes*SAI->number_of_sp_filters; ij++) {
            getRawBytes(sparsemask[ij]);
        }

        for (int32_t ee = 0; ee < SAI->number_of_sp_filters; ee++) {

            SAI->sparse_filters.at(ee).regressor_indexes.clear();

            for (int32_t ii = 0; ii < SAI->Ms; ii++) {
                SAI->sparse_filters.at(ee).regressor_indexes.push_back(0);
            }
        }

        for (int32_t ee = 0; ee < SAI->number_of_sp_filters; ee++) {

            uint32_t ik = 0;

            for (int32_t ij = 0; ij < Nsp; ij++) {

                uint32_t q = ij / 8;

                uint8_t *sparse_mask_byte = &sparsemask[sp_mask_nbytes*ee + q];

                if (*sparse_mask_byte & (1 << (ij - q * 8))) {
                    SAI->sparse_filters.at(ee).regressor_indexes.at(ik) = ij;
                    ik++;
                }

            }
        }

        delete[](sparsemask);

    }
}


void viewParametersConstruct::addSPM() {
    for (int32_t ii = 0; ii < nviews; ii++) {

        view *SAI = LF + ii;

        int32_t Nsp = SAI->SP_B > 0 ?
            (SAI->n_references + 1)*(SAI->NNt * 2 + 1) * (SAI->NNt * 2 + 1) + 1 :
            (SAI->NNt * 2 + 1) * (SAI->NNt * 2 + 1) + 1;

        int32_t sp_mask_nbytes = (Nsp % 8) ? Nsp / 8 + 1 : Nsp / 8;

        uint8_t *sparsemask = new uint8_t[sp_mask_nbytes*SAI->number_of_sp_filters]();

        for (int32_t ee = 0; ee < SAI->number_of_sp_filters; ee++) {

            for (int32_t ij = 0; ij < SAI->Ms; ij++) {

                uint32_t regr_indx =
                    SAI->sparse_filters.at(ee).regressor_indexes.at(ij);

                uint32_t q = regr_indx / 8;

                uint8_t *sparse_mask_byte = &sparsemask[q + sp_mask_nbytes*ee];

                *sparse_mask_byte = *sparse_mask_byte
                    | (1 << (regr_indx - q * 8));

            }
        }

        for (int32_t ij = 0; ij < sp_mask_nbytes*SAI->number_of_sp_filters; ij++) {
            addRawBytes(sparsemask[ij]);
        }

        delete[](sparsemask);
    }
}

void viewParametersConstruct::convertLFtoBytes() {

    addMconfs();

    addX();
    addY();

    addNRefs();
    addNNDRefs();

    addRefs();
    addNDRefs();

    addLSW();
    addSTD();

    addSPP();
    addSPW();
    addSPM();

}

void viewParametersConstruct::convertBytesToLF() {

    dec_iter = rawbytes.data();

    getMconfs();

    getX();
    getY();

    getNRefs();
    getNNDRefs();

    getRefs();
    getNDRefs();

    getLSW();
    getSTD();

    getSPP();
    getSPW();
    getSPM();

}

void viewHeaderToCodestream(
    int32_t &n_bytes_prediction,
    view *SAI,
    FILE *output_LF_file) {

    minimal_config mconf = makeMinimalConfig(SAI);

    //printf("size of minimal_config %i bytes\n", (int32_t)sizeof(minimal_config));

    n_bytes_prediction += (int32_t)fwrite(
        &mconf,
        sizeof(minimal_config),
        1,
        output_LF_file) * sizeof(minimal_config);

    /* lets see what else needs to be written to bitstream */

    if (SAI->has_x_displacement) {
        n_bytes_prediction += (int32_t)fwrite(&SAI->x, sizeof(float), 1,
            output_LF_file) * sizeof(float);
    }

    if (SAI->has_y_displacement) {
        n_bytes_prediction += (int32_t)fwrite(&SAI->y, sizeof(float), 1,
            output_LF_file) * sizeof(float);
    }

    if (SAI->has_color_references) {
        unsigned char tmpNREF = (unsigned char)SAI->n_references;
        n_bytes_prediction += (int32_t)fwrite(&tmpNREF, sizeof(unsigned char), 1,
            output_LF_file) * sizeof(unsigned char);
        for (int32_t ij = 0; ij < SAI->n_references; ij++) {
            uint16_t nid = (uint16_t) *(SAI->references + ij);
            n_bytes_prediction += (int32_t)fwrite(&nid, sizeof(uint16_t), 1,
                output_LF_file)
                * sizeof(uint16_t);
        }
    }

    if (SAI->has_depth_references) {
        unsigned char tmpNDREF = (unsigned char)SAI->n_depth_references;
        n_bytes_prediction += (int32_t)fwrite(&tmpNDREF, sizeof(unsigned char), 1,
            output_LF_file) * sizeof(unsigned char);
        for (int32_t ij = 0; ij < SAI->n_depth_references; ij++) {
            uint16_t nid = (uint16_t) *(SAI->depth_references + ij);
            n_bytes_prediction += (int32_t)fwrite(&nid, sizeof(uint16_t), 1,
                output_LF_file)
                * sizeof(uint16_t);
        }
    }

    int32_t MMM = (1 << SAI->n_references);
    int32_t N_LS = SAI->ncomp*((MMM*SAI->n_references) / 2);

    if (SAI->has_color_references) {
        if (SAI->mmode == 0) {
            /* use LS merging weights */

            int32_t N_LS_R = SAI->nc_merge*((MMM*SAI->n_references) / 2);

            n_bytes_prediction += (int32_t)fwrite(
                SAI->merge_weights,
                sizeof(int16_t),
                N_LS_R,
                output_LF_file)
                * sizeof(int16_t);
        }
        if (SAI->mmode == 1) {
            /* use standard deviation */
            n_bytes_prediction += (int32_t)fwrite(
                &SAI->stdd,
                sizeof(float),
                1,
                output_LF_file) * sizeof(float);
        }
    }

    if (SAI->use_global_sparse) {

        unsigned char tmpNNt = (unsigned char)SAI->NNt;
        unsigned char tmpMs = (unsigned char)SAI->Ms;

        n_bytes_prediction += (int32_t)fwrite(
            &SAI->number_of_sp_filters,
            sizeof(uint16_t),
            1,
            output_LF_file) * sizeof(uint16_t);

        n_bytes_prediction += (int32_t)fwrite(
            &tmpNNt,
            sizeof(unsigned char),
            1,
            output_LF_file) * sizeof(unsigned char);

        n_bytes_prediction += (int32_t)fwrite(
            &tmpMs,
            sizeof(unsigned char),
            1,
            output_LF_file) * sizeof(unsigned char);

        for (int32_t ee = 0; ee < SAI->number_of_sp_filters; ee++) {
            n_bytes_prediction += (int32_t)fwrite(
                &SAI->sparse_filters.at(ee).quantized_filter_coefficients[0],
                sizeof(int16_t),
                SAI->Ms,
                output_LF_file)
                * sizeof(int16_t);
        }

        int32_t Nsp = SAI->SP_B>0 ?
            (SAI->n_references + 1)*(SAI->NNt * 2 + 1) * (SAI->NNt * 2 + 1) + 1 :
            (SAI->NNt * 2 + 1) * (SAI->NNt * 2 + 1) + 1;

        int32_t sp_mask_nbytes = (Nsp % 8) ? Nsp / 8 + 1 : Nsp / 8;

        uint8_t *sparsemask = new uint8_t[sp_mask_nbytes*SAI->number_of_sp_filters]();

        for (int32_t ee = 0; ee < SAI->number_of_sp_filters; ee++) {

            for (int32_t ij = 0; ij < SAI->Ms; ij++) {

                uint32_t regr_indx =
                    SAI->sparse_filters.at(ee).regressor_indexes.at(ij);

                uint32_t q = regr_indx / 8;

                uint8_t *sparse_mask_byte = &sparsemask[q + sp_mask_nbytes*ee];

                *sparse_mask_byte = *sparse_mask_byte
                    | (1 << (regr_indx - q * 8));

            }
        }

        n_bytes_prediction += (int32_t)fwrite(
            sparsemask,
            sizeof(uint8_t),
            sp_mask_nbytes*SAI->number_of_sp_filters,
            output_LF_file)
            * sizeof(uint8_t);

        delete[](sparsemask);

    }

    return;

}

void codestreamToViewHeader(
    int32_t &n_bytes_prediction,
    view *SAI,
    FILE *input_LF) {

    minimal_config mconf;

    n_bytes_prediction += (int32_t)fread(
        &mconf,
        sizeof(minimal_config),
        1,
        input_LF)
        * sizeof(minimal_config);

    setup_form_minimal_config(&mconf, SAI);

    if (SAI->has_x_displacement) {
        n_bytes_prediction += (int32_t)fread(&SAI->x, sizeof(float), 1, input_LF)
            * sizeof(float);
    }

    if (SAI->has_y_displacement) {
        n_bytes_prediction += (int32_t)fread(&SAI->y, sizeof(float), 1, input_LF)
            * sizeof(float);
    }

    if (SAI->has_color_references) {

        unsigned char tmpNREF = 0;

        n_bytes_prediction += (int32_t)fread(&tmpNREF, sizeof(unsigned char), 1,
            input_LF) * sizeof(unsigned char);

        SAI->n_references = tmpNREF;

        SAI->references = new int32_t[SAI->n_references]();
        for (int32_t ij = 0; ij < SAI->n_references; ij++) {
            uint16_t nid;
            n_bytes_prediction += (int32_t)fread(&nid, sizeof(uint16_t), 1,
                input_LF) * sizeof(uint16_t);
            *(SAI->references + ij) = (int32_t)nid;
        }
    }

    if (SAI->has_depth_references) {

        unsigned char tmpNDREF = 0;

        n_bytes_prediction += (int32_t)fread(&tmpNDREF, sizeof(unsigned char), 1,
            input_LF) * sizeof(unsigned char);

        SAI->n_depth_references = tmpNDREF;

        SAI->depth_references = new int32_t[SAI->n_depth_references]();
        for (int32_t ij = 0; ij < SAI->n_depth_references; ij++) {
            uint16_t nid;
            n_bytes_prediction += (int32_t)fread(&nid, sizeof(uint16_t), 1,
                input_LF) * sizeof(uint16_t);
            *(SAI->depth_references + ij) = (int32_t)nid;
        }
    }

    SAI->NB = (1 << SAI->n_references) * SAI->n_references;

    int32_t MMM = (1 << SAI->n_references);

    if (SAI->has_color_references) {
        if (SAI->mmode == 0) {

            int32_t N_LS = SAI->ncomp*((MMM*SAI->n_references) / 2);

            SAI->merge_weights = new int16_t[N_LS]();
            /* use LS merging weights */

            int32_t N_LS_R = SAI->nc_merge*((MMM*SAI->n_references) / 2);

            n_bytes_prediction += (int32_t)fread(
                SAI->merge_weights,
                sizeof(int16_t),
                N_LS_R,
                input_LF)
                * sizeof(int16_t);
        }
        if (SAI->mmode == 1) {
            /* use standard deviation */
            n_bytes_prediction += (int32_t)fread(
                &SAI->stdd,
                sizeof(float),
                1,
                input_LF) * sizeof(float);
        }
    }

    if (SAI->use_global_sparse) {

        unsigned char tmpNNt = 0;
        unsigned char tmpMs = 0;

        n_bytes_prediction += (int32_t)fread(
            &SAI->number_of_sp_filters,
            sizeof(uint16_t),
            1,
            input_LF) * sizeof(uint16_t);

        n_bytes_prediction += (int32_t)fread(
            &tmpNNt,
            sizeof(unsigned char),
            1,
            input_LF) * sizeof(unsigned char);

        n_bytes_prediction += (int32_t)fread(
            &tmpMs,
            sizeof(unsigned char),
            1,
            input_LF) * sizeof(unsigned char);

        SAI->NNt = (int32_t)tmpNNt;
        SAI->Ms = (int32_t)tmpMs;

        SAI->sparse_filters.clear();

        for (int32_t ee = 0; ee < SAI->number_of_sp_filters; ee++) {

            spfilter tmpsp;
            SAI->sparse_filters.push_back(tmpsp);

            SAI->sparse_filters.at(ee).Ms = SAI->Ms;
            SAI->sparse_filters.at(ee).NNt = SAI->NNt;

            SAI->sparse_filters.at(ee).MT = SAI->SP_B>0 ? 
                (SAI->n_references + 1)*(SAI->NNt * 2 + 1) * (SAI->NNt * 2 + 1) + 1 :
                (SAI->NNt * 2 + 1) * (SAI->NNt * 2 + 1) + 1;

            SAI->sparse_filters.at(ee).quantized_filter_coefficients.clear();

            for (int32_t ii = 0; ii < SAI->Ms; ii++) {
                SAI->sparse_filters.at(ee).quantized_filter_coefficients.push_back(0);
            }
        }

        for (int32_t ee = 0; ee < SAI->number_of_sp_filters; ee++)
        {
            n_bytes_prediction += (int32_t)fread(
                &SAI->sparse_filters.at(ee).quantized_filter_coefficients[0],
                sizeof(int16_t),
                SAI->Ms,
                input_LF)
            * sizeof(int16_t);

        }

        int32_t Nsp = SAI->sparse_filters.at(0).MT;// (SAI->NNt * 2 + 1)* (SAI->NNt * 2 + 1) + 1;
        int32_t sp_mask_nbytes = (Nsp % 8) ? Nsp / 8 + 1 : Nsp / 8;

        uint8_t *sparsemask = new uint8_t[sp_mask_nbytes*SAI->number_of_sp_filters]();

        n_bytes_prediction += (int32_t)fread(
            sparsemask,
            sizeof(uint8_t),
            sp_mask_nbytes*SAI->number_of_sp_filters,
            input_LF)
            * sizeof(uint8_t);

        for (int32_t ee = 0; ee < SAI->number_of_sp_filters; ee++) {

            SAI->sparse_filters.at(ee).regressor_indexes.clear();

            for (int32_t ii = 0; ii < Nsp; ii++) {
                SAI->sparse_filters.at(ee).regressor_indexes.push_back(0);
            }
        }

        for (int32_t ee = 0; ee < SAI->number_of_sp_filters; ee++) {

            uint32_t ik = 0;

            for (int32_t ij = 0; ij < Nsp; ij++) {

                uint32_t q = ij / 8;

                uint8_t *sparse_mask_byte = &sparsemask[sp_mask_nbytes*ee + q];

                if (*sparse_mask_byte & (1 << (ij - q * 8))) {
                    SAI->sparse_filters.at(ee).regressor_indexes.at(ik) = ij;
                    ik++;
                }

            }
        }

        delete[](sparsemask);

    }

    return;

}