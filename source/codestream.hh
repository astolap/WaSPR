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

#ifndef CODESTREAM_HH
#define CODESTREAM_HH

#include "view.hh"
#include "minconf.hh"
#include <iostream>

#include <cstdint>

using std::int32_t;
using std::uint32_t;

using std::int16_t;
using std::uint16_t;

using std::int8_t;
using std::uint8_t;

class viewParametersConstruct {

private:

    std::vector<uint8_t> rawbytes; /*stores LF parameters as bytes*/

    const std::string gzippath;
    const std::string rawfile;
    const std::string mode;

    view *LF;
    int32_t nviews;

    uint8_t *dec_iter;

    void convertLFtoBytes(); /*from LF to rawbytes*/
    void convertBytesToLF(); /*from rawbytes to LF*/

    template<class T1>
    void addRawBytes(T1 datai)
    {
        std::vector<uint8_t> bytesi(sizeof(T1), 0);
        memcpy(bytesi.data(), &datai, sizeof(T1));
        rawbytes.insert(rawbytes.end(), bytesi.begin(), bytesi.end());
    }

    template<class T1>
    void getRawBytes(T1 &dataout)
    {
        memcpy(&dataout, dec_iter, sizeof(T1));
        dec_iter += sizeof(T1);
    }

    void addSPM();
    void addSPW();
    void addSPP();
    void addSTD();
    void addLSW();

    void addNRefs();
    void addNNDRefs();

    void addNDRefs();
    void addRefs();
    void addY();
    void addX();
    void addMconfs();

    void getSPM();
    void getSPW();
    void getSPP();
    void getSTD();
    void getLSW();

    void getNRefs();
    void getNNDRefs();

    void getNDRefs();
    void getRefs();
    void getY();
    void getX();
    void getMconfs();


public:

    viewParametersConstruct(
        view *LF,
        const int32_t nviews,
        const std::string gzippath,
        const std::string rawfile,
        const std::string mode);

    ~viewParametersConstruct();

};

void viewHeaderToCodestream(
    int32_t &n_bytes_prediction, 
    view *SAI,
    FILE *output_LF_file);

void codestreamToViewHeader(
    int32_t &n_bytes_prediction, 
    view *SAI, 
    FILE *input_LF);

#endif
