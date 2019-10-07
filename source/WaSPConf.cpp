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

#include <string.h>
#include "WaSPConf.hh"

WaSPConfig::WaSPConfig(int argc, char *argv[], const char *type) {

    print_intro();

    if (!strcmp(type,"encoder") && !parseCommandLine_encoder(argc, argv)) {
        print_encoder_help();
        exit(0);
    }

    if (!strcmp(type, "decoder") && !parseCommandLine_decoder(argc, argv)) {
        print_decoder_help();
        exit(0);
    }

}

WaSPConfig::~WaSPConfig() {

}

bool WaSPConfig::parseCommandLine_decoder(int argc, char *argv[]) {

    for (int32_t ii = 1; ii < argc - 1; ii += 2) {

        if (!strcmp(argv[ii], "-i")) {
            WaSP_setup.input_directory = std::string(argv[ii + 1]);
        }

        else if (!strcmp(argv[ii], "--input")) {
            WaSP_setup.input_directory = std::string(argv[ii + 1]);
        }

        else if (!strcmp(argv[ii], "-o")) {
            WaSP_setup.output_directory = std::string(argv[ii + 1]);
        }

        else if (!strcmp(argv[ii], "--output")) {
            WaSP_setup.output_directory = std::string(argv[ii + 1]);
        }

        else if (!strcmp(argv[ii], "-k")) {
            WaSP_setup.wasp_kakadu_directory = std::string(argv[ii + 1]);
        }

        else if (!strcmp(argv[ii], "--kakadu")) {
            WaSP_setup.wasp_kakadu_directory = std::string(argv[ii + 1]);
        }

        else if (!strcmp(argv[ii], "-t")) {
            WaSP_setup.hm_decoder = std::string(argv[ii + 1]);
        }

        else if (!strcmp(argv[ii], "--TAppDecoder")) {
            WaSP_setup.hm_decoder = std::string(argv[ii + 1]);
        }

        else if (!strcmp(argv[ii], "--kvazaar-path")) {
            WaSP_setup.kvazaarpath = std::string(argv[ii + 1]);

        }

        else if (!strcmp(argv[ii], "--gzip-path")) {
            WaSP_setup.gzipath = std::string(argv[ii + 1]);

        }

        else {
            return false;
        }

    }

    if (WaSP_setup.input_directory.length() == 0) {
        printf("\n Input directory not set\n");
        return false;
    }

    if (WaSP_setup.output_directory.length() == 0) {
        printf("\n Output directory not set\n");
        return false;
    }

    if (WaSP_setup.wasp_kakadu_directory.length() == 0) {
        printf("\n Kakadu directory not set\n");
        return false;
    }

    if (WaSP_setup.hm_decoder.length() == 0) {
        printf("\n Path to TAppDecoder needs to be defined\n");
        return false;
    }

    WaSP_setup.stats_file = WaSP_setup.output_directory + "/stats.json";

    return true;

}

void WaSPConfig::print_encoder_help() {
    printf("\n\tUsage: wasp-encoder"
        "\n\t--input [INPUT DIRECTORY .PPM/.PGM]"
        "\n\t--output [OUTPUT DIRECTORY .LF/.PPM/.PGM]"
        "\n\t--config [JSON CONFIG]"
        "\n\t--kakadu [KAKADU BINARY DIRECTORY]"
        "\n\t--TAppEncoder [Path to TAppEncoder executable]"
        "\n\t--TAppDecoder [Path to TAppDecoder executable]"
        "\n\t--HEVCcfg [Path to TAppEncoder config file]"
        "\n\t--kvazaar-path [path to Kvazaar binary]"
        "\n\t--gzip-path [path to gzip binary]"
        "\n\t--sparse_subsampling [Subsampling factor when solving sparse filter,"
        " needs to be integer >0. Values 2 or 4 will increase encoder speed with some loss in PSNR.]\n\n");
    return;
}

void WaSPConfig::print_decoder_help() {
    printf("\n\tUsage: wasp-decoder"
        "\n\t--input [INPUT .LF]"
        "\n\t--output [OUTPUT DIRECTORY .PPM/.PGM]"
        "\n\t--kakadu [KAKADU BINARY DIRECTORY]"
        "\n\t--TAppDecoder [Path to TAppDecoder executable]"
        "\n\t--kvazaar-path [path to Kvazaar binary]"
        "\n\t--gzip-path [path to gzip binary]\n\n");
    return;
}

void WaSPConfig::print_intro() {
    printf("\n\t--------------------------------------------------------\n");
    printf(
        "\n\tWaSP - Warping and Sparse Prediction"
        "\n\tAuthor: Pekka Astola (pekka.astola@tuni.fi), 2018-2019"
        "\n\tWeb page: https://github.com/astolap/WaSP"
        "\n\tHome page of author: http://www.cs.tut.fi/~astolap/"
        "\n\tPublication: P. Astola and I. Tabus,"
        "\n\t\t\tWaSP: Hierarchical Warping, Merging, and Sparse Prediction for Light Field Image Compression,"
        "\n\t\t\t2018 7th European Workshop on Visual Information Processing (EUVIP), Tampere, 2018, pp. 1-6."
        "\n\n\tBuild date: %s %s\n\n",__DATE__,__TIME__);
    printf("\n\t--------------------------------------------------------\n");
    return;
}

bool WaSPConfig::parseCommandLine_encoder(int argc, char *argv[]) {

    WaSP_setup.sparse_subsampling = 1;

    for (int32_t ii = 1; ii < argc-1; ii+=2) {

        if (!strcmp(argv[ii], "-c")) {
            WaSP_setup.config_file = std::string(argv[ii + 1]);
        }

        else if (!strcmp(argv[ii], "--config")) {
            WaSP_setup.config_file = std::string(argv[ii + 1]);
        }

        else if (!strcmp(argv[ii], "-i")) {
            WaSP_setup.input_directory = std::string(argv[ii+1]);
        }

        else if (!strcmp(argv[ii], "--input")) {
            WaSP_setup.input_directory = std::string(argv[ii + 1]);
        }

        else if (!strcmp(argv[ii], "-o")) {
            WaSP_setup.output_directory = std::string(argv[ii+1]);
        }

        else if (!strcmp(argv[ii], "--output")) {
            WaSP_setup.output_directory = std::string(argv[ii + 1]);
        }

        else if (!strcmp(argv[ii], "-k")) {
            WaSP_setup.wasp_kakadu_directory = std::string(argv[ii+1]);
        }

        else if (!strcmp(argv[ii], "--kakadu")) {
            WaSP_setup.wasp_kakadu_directory = std::string(argv[ii + 1]);

        }

        else if (!strcmp(argv[ii], "-s")) {
            WaSP_setup.sparse_subsampling = atoi(argv[ii + 1]);
        }

        else if (!strcmp(argv[ii], "--sparse_subsampling")) {
            WaSP_setup.sparse_subsampling = atoi(argv[ii + 1]);

        }

        else if (!strcmp(argv[ii], "-t")) {
            WaSP_setup.hm_encoder = std::string(argv[ii + 1]);

        }

        else if (!strcmp(argv[ii], "--TAppEncoder")) {
            WaSP_setup.hm_encoder = std::string(argv[ii + 1]);

        }

        else if (!strcmp(argv[ii], "-d")) {
            WaSP_setup.hm_decoder = std::string(argv[ii + 1]);

        }

        else if (!strcmp(argv[ii], "--TAppDecoder")) {
            WaSP_setup.hm_decoder = std::string(argv[ii + 1]);

        }

        else if (!strcmp(argv[ii], "-x")) {
            WaSP_setup.hm_cfg = std::string(argv[ii + 1]);

        }

        else if (!strcmp(argv[ii], "--HEVCcfg")) {
            WaSP_setup.hm_cfg = std::string(argv[ii + 1]);

        }

        else if (!strcmp(argv[ii], "--kvazaar-path")) {
            WaSP_setup.kvazaarpath = std::string(argv[ii + 1]);

        }

        else if (!strcmp(argv[ii], "--gzip-path")) {
            WaSP_setup.gzipath = std::string(argv[ii + 1]);

        }

        else {
            return false;
        }

    }


    if (WaSP_setup.input_directory.length() == 0) {
        printf("\n Input directory not set\n");
        return false;
    }

    if (WaSP_setup.output_directory.length() == 0) {
        printf("\n Output directory not set\n");
        return false;
    }

    if (WaSP_setup.config_file.length() == 0) {
        printf("\n Config file (.json) not set\n");
        return false;
    }

    if (WaSP_setup.wasp_kakadu_directory.length() == 0) {
        printf("\n Kakadu directory not set\n");
        return false;
    }

    if (WaSP_setup.hm_encoder.length() < 1) {
        printf("\n Path to TAppEncoder needs to be defined\n");
        return false;
    }

    if (WaSP_setup.hm_decoder.length() < 1) {
        printf("\n Path to TAppDecoder needs to be defined\n");
        return false;
    }

    if (WaSP_setup.hm_cfg.length() < 1) {
        printf("\n Path to HM .cfg needs to be defined\n");
        return false;
    }

    if (WaSP_setup.sparse_subsampling < 1) {
        printf("\n Sub sampling factor needs to be >= 1\n");
        return false;
    }

    WaSP_setup.stats_file = WaSP_setup.output_directory + "/stats.json";

    return true;

}