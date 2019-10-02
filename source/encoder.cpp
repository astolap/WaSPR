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

#include <fstream>
#include <iomanip>
#include <string>

#include "encoder.hh"
#include "ppm.hh"
#include "fileaux.hh"
#include "codestream.hh"
#include "residual.hh"
#include "clip.hh"
#include "json.hh"
#include "view.hh"
#include "merging.hh"
#include "medianfilter.hh"
#include "inpainting.hh"
#include "predictdepth.hh"
#include "ycbcr.hh"
#include "psnr.hh"
#include "warping.hh"
#include "bitdepth.hh"
#include "segmentation.hh"

encoder::encoder(const WaSPsetup encoder_setup)
{

    setup = encoder_setup;

    load_config_json(setup.config_file);

}

encoder::~encoder() {
}

void encoder::encode() {

    generate_normalized_disparity();
    generate_texture();
    write_bitstream();

}

void encoder::write_statsfile() {

    nlohmann::json conf_out;

    conf_out["type"] = "WaSP";

    conf_out["n_seg_iterations"] = n_seg_iterations;

    vector<nlohmann::json::object_t> views;

    for (int32_t ii = 0; ii < n_views_total; ii++) {

        view *SAI = LF + ii;

        nlohmann::json view_configuration;

        view_configuration["column_index"] = SAI->c;
        view_configuration["row_index"] = SAI->r;
        view_configuration["index"] = SAI->i_order;
        view_configuration["level"] = SAI->level;
        view_configuration["finalQP"] = SAI->finalQP;

        view_configuration["NNt"] = SAI->NNt;
        view_configuration["Ms"] = SAI->Ms;

        view_configuration["num_of_sparse_filters"] = SAI->sparse_filters.size();

        view_configuration["QP_range"] = SAI->QP_range;
        view_configuration["bpp_range"] = SAI->bpp_range;

        view_configuration["real_bpp_texture"] = SAI->real_rate_texture;
        view_configuration["real_bpp_normdisp"] = SAI->real_rate_normpdisp;

        for (int32_t ij = 0; ij < SAI->sparse_filters.size(); ij++) {
            view_configuration[std::string("sp_qcoeffs_"+std::to_string(ij)).c_str()] = 
                SAI->sparse_filters.at(ij).quantized_filter_coefficients;
            view_configuration[std::string("sp_regr_indices_" + std::to_string(ij)).c_str()] =
                SAI->sparse_filters.at(ij).regressor_indexes;
        }

        views.push_back(view_configuration);

    }

    conf_out["views"] = views;

    std::ofstream file(setup.stats_file);

    file << std::setw(2) << conf_out << std::endl;

}

void encoder::write_config(string config_json_file_out) {

    nlohmann::json conf_out;

    conf_out["type"] = "WaSP";

    conf_out["number_of_views_to_be_encoded"] = n_views_total;
    conf_out["enable_fixed_weight_parameter_search"] = std_search_s;

    conf_out["image_height"] = nr;
    conf_out["image_width"] = nc;

    conf_out["colorspace"] = colorspace_LF;

    vector<nlohmann::json::object_t> views;

    for (int32_t ii = 0; ii < n_views_total; ii++) {

        view *SAI = LF + ii;

        nlohmann::json view_configuration;

        view_configuration["id"] = SAI->i_order;

        view_configuration["level"] = SAI->level;

        view_configuration["row_index"] = SAI->r;
        view_configuration["column_index"] = SAI->c;

        view_configuration["cweight_search"] = SAI->cweight_search;

        view_configuration["rate_texture"] = SAI->residual_rate_color,
            view_configuration["rate_inverse_depth"] = SAI->residual_rate_depth;

        view_configuration["view_merging_mode"] = SAI->mmode;

        view_configuration["horizontal_camera_center_position"] = SAI->x;
        view_configuration["vertical_camera_center_position"] = SAI->y;

        view_configuration["sparse_filter_order"] = SAI->Ms;
        view_configuration["sparse_filter_neighborhood_size"] = SAI->NNt;

        view_configuration["fixed_merging_weight_parameter"] = SAI->stdd;

        view_configuration["minimum_inverse_depth"] = SAI->min_inv_d;



        view_configuration["number_of_texture_references"] = SAI->n_references;
        view_configuration["number_of_inverse_depth_references"] =
            SAI->n_depth_references;

        std::vector<int32_t> sai_texture_references;
        for (int32_t ij = 0; ij < SAI->n_references; ij++) {
            sai_texture_references.push_back(SAI->references[ij]);
        }
        view_configuration["texture_reference_indices"] =
            sai_texture_references;

        std::vector<int32_t> sai_inverse_depth_references;
        for (int32_t ij = 0; ij < SAI->n_depth_references; ij++) {
            sai_inverse_depth_references.push_back(SAI->depth_references[ij]);
        }
        view_configuration["inverse_depth_reference_indices"] =
            sai_inverse_depth_references;


        views.push_back(view_configuration);

    }

    conf_out["views"] = views;

    std::ofstream file(config_json_file_out);
    file << std::setw(2) << conf_out << std::endl;
}

void encoder::load_config_json(string config_json_file) {

    ifstream ifs(config_json_file);
    nlohmann::json conf = nlohmann::json::parse(ifs);

    n_views_total = conf["number_of_views_to_be_encoded"].get<int32_t>();
    std_search_s = conf["enable_fixed_weight_parameter_search"].get<int32_t>();

    STD_SEARCH = std_search_s > 0 ? true : false;

    nr = conf["image_height"].get<int32_t>();
    nc = conf["image_width"].get<int32_t>();

    nc_sparse = conf["nc_sparse"].get<int32_t>();
    nc_merge = conf["nc_merge"].get<int32_t>();
    nc_color_ref = conf["nc_color_ref"].get<int32_t>();
    SP_B = conf["SP_B"].get<int32_t>();

    n_seg_iterations = conf["n_seg_iterations"].get<int32_t>();

    vector<nlohmann::json::object_t> conf_views =
        conf["views"].get<vector<nlohmann::json::object_t>>();

    colorspace_LF = conf["colorspace"].get<std::string>();

    LF = new view[n_views_total]();
    for (nlohmann::json::object_t view_configuration : conf_views) {

        int32_t ii = view_configuration["id"].get<int32_t>();
        view *SAI = LF + ii;

        SAI->nr = nr;
        SAI->nc = nc;

        SAI->colorspace = colorspace_LF;

        initView(SAI);
        SAI->i_order = ii;

        SAI->ncomp = 3;

        SAI->nc_sparse = nc_sparse;
        SAI->nc_merge = nc_merge;

        SAI->SP_B = SP_B;

        SAI->r = view_configuration["row_index"].get<int32_t>();
        SAI->c = view_configuration["column_index"].get<int32_t>();

        SAI->cweight_search =
            view_configuration["cweight_search"].get<int32_t>() > 0 ? true : false;

        SAI->residual_rate_color =
            view_configuration["rate_texture"].get<float>();
        SAI->residual_rate_depth =
            view_configuration["rate_inverse_depth"].get<float>();

        SAI->mmode = view_configuration["view_merging_mode"].get<unsigned char>();

        SAI->x = view_configuration["horizontal_camera_center_position"].get<float>();
        SAI->y = view_configuration["vertical_camera_center_position"].get<float>();

        SAI->Ms = view_configuration["sparse_filter_order"].get<int32_t>();

        SAI->NNt = view_configuration["sparse_filter_neighborhood_size"].get<int32_t>();

        SAI->stdd = view_configuration["fixed_merging_weight_parameter"].get<float>();

        SAI->min_inv_d = view_configuration["minimum_inverse_depth"].get<int32_t>();
        SAI->n_references =
            view_configuration["number_of_texture_references"].get<int32_t>();
        SAI->n_depth_references =
            view_configuration["number_of_inverse_depth_references"].get<int32_t>();

        SAI->level = view_configuration["level"].get<int32_t>();

        SAI->preset_QP = view_configuration["QP"].get<int32_t>();

        if (abs(SAI->x) > 0.0001) {
            SAI->has_x_displacement = true;
        }

        if (abs(SAI->y) > 0.0001) {
            SAI->has_y_displacement = true;
        }

        if (SAI->n_references > 0) {

            SAI->has_color_references = true;

            SAI->references = new int32_t[SAI->n_references]();

            vector<int32_t> texture_references =
                view_configuration["texture_reference_indices"].get<vector<int32_t>>();

            memcpy(
                SAI->references,
                &texture_references[0],
                sizeof(int32_t)*SAI->n_references);

        }

        if (SAI->n_depth_references > 0) {

            SAI->has_depth_references = true;

            SAI->depth_references = new int32_t[SAI->n_depth_references]();

            vector<int32_t> depth_references =
                view_configuration["inverse_depth_reference_indices"].get<vector<int32_t>>();

            memcpy(
                SAI->depth_references,
                &depth_references[0],
                sizeof(int32_t)*SAI->n_depth_references);

        }

        setPaths(
            SAI,
            setup.input_directory.c_str(),
            setup.output_directory.c_str());

    }
}

void encoder::forward_warp_texture_references(
    view *LF,
    view *SAI,
    uint16_t **warped_texture_views,
    uint16_t **warped_depth_views,
    float **DispTargs) {

    for (int32_t ij = 0; ij < SAI->n_references; ij++) {

        view *ref_view = LF + SAI->references[ij];

        int32_t tmp_w, tmp_r, tmp_ncomp;

        aux_read16PGMPPM(
            ref_view->path_out_pgm,
            tmp_w,
            tmp_r,
            tmp_ncomp,
            ref_view->depth);

        aux_read16PGMPPM(
            ref_view->path_internal_colorspace_out_ppm,
            tmp_w,
            tmp_r,
            tmp_ncomp,
            ref_view->color);

        /* FORWARD warp color */
        warpView0_to_View1(
            ref_view,
            SAI,
            warped_texture_views[ij],
            warped_depth_views[ij],
            DispTargs[ij]);

        delete[](ref_view->depth);
        delete[](ref_view->color);

        ref_view->depth = nullptr;
        ref_view->color = nullptr;

        if (SAVE_PARTIAL_WARPED_VIEWS) {

            char tmp_str[1024];

            sprintf(tmp_str, "%s/%03d_%03d_warped_to_%03d_%03d.ppm",
                setup.output_directory.c_str(), (ref_view)->c, (ref_view)->r,
                SAI->c, SAI->r);

            aux_write16PGMPPM(tmp_str, SAI->nc, SAI->nr, 3, warped_texture_views[ij]);

            sprintf(tmp_str, "%s/%03d_%03d_warped_to_%03d_%03d.pgm",
                setup.output_directory.c_str(), (ref_view)->c, (ref_view)->r,
                SAI->c, SAI->r);

            aux_write16PGMPPM(tmp_str, SAI->nc, SAI->nr, 1, warped_depth_views[ij]);

        }

    }


}

void encoder::merge_texture_views(
    view *SAI,
    view *LF,
    uint16_t **warped_texture_views,
    float **DispTargs) {

    initViewW(SAI, DispTargs);

    if (SAI->mmode == 0) {

        int32_t nr1, nc1, ncomp1;

        uint16_t *original_color_view = read_input_ppm(
            SAI->path_input_ppm,
            nr1,
            nc1,
            ncomp1,
            10,
            SAI->colorspace);

        for (int32_t icomp = 0; icomp < SAI->nc_merge; icomp++) {
            getViewMergingLSWeights_icomp(
                SAI,
                warped_texture_views,
                original_color_view,
                icomp);
        }

        delete[](original_color_view);

        for (int32_t icomp = 0; icomp < SAI->nc_merge; icomp++) {
            mergeWarped_N_icomp(
                warped_texture_views,
                SAI,
                icomp);
        }

        /* hole filling for texture*/
        for (int32_t icomp = 0; icomp < SAI->nc_merge; icomp++) {
            uint32_t nholes = holefilling(
                SAI->color + icomp*SAI->nr*SAI->nc,
                SAI->nr,
                SAI->nc,
                (uint16_t)0,
                SAI->seg_vp);
        }


    }

    if (SAI->mmode == 1) {

        /* we don't use LS weights but something derived on geometric distance in view array*/
        for (int32_t icomp = 0; icomp < SAI->nc_merge; icomp++) {
            getGeomWeight_icomp(
                SAI,
                LF,
                icomp);
        }

        /* merge color with prediction */
        for (int32_t icomp = 0; icomp < SAI->nc_merge; icomp++) {
            mergeWarped_N_icomp(
                warped_texture_views,
                SAI,
                icomp);
        }

        /* hole filling for texture*/
        for (int32_t icomp = 0; icomp < nc_merge; icomp++) {
            uint32_t nholes = holefilling(
                SAI->color + icomp*SAI->nr*SAI->nc,
                SAI->nr,
                SAI->nc,
                (uint16_t)0,
                SAI->seg_vp);
        }

    }

    if (SAI->mmode == 2) {

        /*merge with median operator*/
        mergeMedian_N(warped_texture_views, DispTargs, SAI, nc_merge);

        /* hole filling for texture*/
        for (int32_t icomp = 0; icomp < nc_merge; icomp++) {
            uint32_t nholes = holefilling(
                SAI->color + icomp*SAI->nr*SAI->nc,
                SAI->nr,
                SAI->nc,
                (uint16_t)0,
                SAI->seg_vp);
        }
    }

    delete[](SAI->seg_vp);
    SAI->seg_vp = nullptr;
    delete[](SAI->bmask);
    SAI->bmask = nullptr;
}

void encoder::generate_normalized_disparity() {

    /*the term inverse depth is used interchangeably with normalized disparity */

    int32_t maxh = get_highest_level(LF, n_views_total);

    /*we DO have to do the last level*/
    for (int32_t hlevel = 1; hlevel <= maxh; hlevel++) {

        std::vector< int32_t > view_indices;

        for (int32_t ii = 0; ii < n_views_total; ii++) {
            if ((LF + ii)->level == hlevel) {
                view_indices.push_back(ii);
            }
        }

        /*ascending order of view index at a particular level*/
        sort(view_indices.begin(), view_indices.end());

        for (int32_t ii = 0; ii < view_indices.size(); ii++) {

            view *SAI = LF + view_indices.at(ii);

            SAI->depth = new uint16_t[SAI->nr * SAI->nc]();

            if (hlevel == 1 && SAI->residual_rate_depth > 0) { /*intra coding of inverse depth*/

                bool depth_file_exist = false;

                int32_t nc1, nr1, ncomp1;

                delete[](SAI->depth);
                SAI->depth = nullptr;

                SAI->depth_file_exist = aux_read16PGMPPM(
                    SAI->path_input_pgm,
                    nc1,
                    nr1,
                    ncomp1,
                    SAI->depth);

                if (SAI->depth_file_exist) {

                    /* ------------------------------
                    INVERSE DEPTH ENCODING STARTS
                    -------------------------------*/

                    printf("Encoding normalized disparity for view %03d_%03d\n", SAI->c, SAI->r);

                    /*write inverse depth to .pgm (path SAI->path_out_pgm) */
                    aux_write16PGMPPM(
                        SAI->path_out_pgm,
                        SAI->nc,
                        SAI->nr,
                        1,
                        SAI->depth);

                    delete[](SAI->depth);
                    SAI->depth = nullptr;

                    char *oparams = kakadu_oparams(
                        SAI->residual_rate_depth,
                        "YCbCr"); /*cycc shouldn't matter for single-channel*/

                    char *encoding_parameters = new char[65535]();
                    sprintf(
                        encoding_parameters,
                        "%s",
                        oparams);

                    encodeKakadu(
                        SAI->path_out_pgm,
                        (setup.wasp_kakadu_directory + "/kdu_compress").c_str(),
                        SAI->jp2_residual_depth_path_jp2,
                        encoding_parameters,
                        SAI->residual_rate_depth);

                    delete[](encoding_parameters);
                    /* ------------------------------
                    INVERSE DEPTH ENCODING ENDS
                    ------------------------------*/

                    /* ------------------------------
                    INVERSE DEPTH DECODING STARTS
                    ------------------------------*/

                    double bytesndisp = aux_GetFileSize(SAI->jp2_residual_depth_path_jp2);
                    double bppndisp = 
                        bytesndisp * 8.0 
                        / static_cast<double>(SAI->nr) 
                        / static_cast<double>(SAI->nc);

                    SAI->real_rate_normpdisp = bppndisp;

                    printf("Decoding normalized disparity for view %03d_%03d\n", SAI->c, SAI->r);

                    decodeKakadu(
                        SAI->path_out_pgm,
                        (setup.wasp_kakadu_directory + "/kdu_expand").c_str(),
                        SAI->jp2_residual_depth_path_jp2);

                    /*------------------------------
                    INVERSE DEPTH DECODING ENDS
                    ------------------------------*/

                    SAI->has_depth_residual = true;

                    aux_read16PGMPPM(
                        SAI->path_out_pgm,
                        nc1,
                        nr1,
                        ncomp1,
                        SAI->depth);

                    if (MEDFILT_DEPTH) {

                        uint16_t *filtered_depth = medfilt2D(
                            SAI->depth,
                            3,
                            SAI->nr,
                            SAI->nc);

                        memcpy(
                            SAI->depth,
                            filtered_depth,
                            sizeof(uint16_t) * SAI->nr * SAI->nc);

                        delete[](filtered_depth);

                    }

                    //aux_write16PGMPPM(
                    //    SAI->path_out_pgm,
                    //    nc1,
                    //    nr1,
                    //    ncomp1,
                    //    SAI->depth);

                }
            }
            else { /*prediction only*/

                printf("Predicting normalized disparity for view %03d_%03d\n", SAI->c, SAI->r);

                /*color array needed in warping function ...*/
                SAI->color = new uint16_t[SAI->nr * SAI->nc * 3]();
                /*SAI->depth = new uint16_t[SAI->nr * SAI->nc]();*/

                WaSP_predict_depth(SAI, LF);

                if (MEDFILT_DEPTH) {

                    uint16_t *filtered_depth = medfilt2D(
                        SAI->depth,
                        3,
                        SAI->nr,
                        SAI->nc);

                    memcpy(
                        SAI->depth,
                        filtered_depth,
                        sizeof(uint16_t) * SAI->nr * SAI->nc);

                    delete[](filtered_depth);

                }

                delete[](SAI->color);
                SAI->color = nullptr;

            }

            aux_write16PGMPPM(
                SAI->path_out_pgm,
                SAI->nc,
                SAI->nr,
                1,
                SAI->depth);

            delete[](SAI->depth);
            SAI->depth = nullptr;

        }
    }
}

void encoder::generate_texture() {

    //FILE *tmp;
    //tmp = fopen("C:/Temp/coeffs.data", "wb");
    //int32_t qvalc = (1<<BIT_DEPTH_SPARSE);
    //fwrite(&qvalc, sizeof(int32_t), 1, tmp);
    //fclose(tmp);

    maxh = get_highest_level(LF, n_views_total);

    for (int32_t hlevel = 1; hlevel <= maxh; hlevel++) {

        printf("\n\tProcessing of hierarchical level: %d\n\n", hlevel);

        const int32_t bpc = 10;

        int32_t Q = 1;
        int32_t offset = 0;

        if (hlevel > 1) {
            Q = 2;

            offset = (1 << bpc) - 1; /* 10bit images currently */
        }

        std::vector< int32_t > view_indices;

        for (int32_t ii = 0; ii < n_views_total; ii++) {
            if ((LF + ii)->level == hlevel) {
                view_indices.push_back(ii);
            }
        }

        /*ascending order of view index at level=hlevel*/
        sort(view_indices.begin(), view_indices.end());

        /* predict (i.e., warp and merge) all views at level=hlevel*/
        for (int32_t ii = 0; ii < view_indices.size(); ii++) {

            view *SAI = LF + view_indices.at(ii);

            printf("Encoding view %03d_%03d\n", SAI->c, SAI->r);

            SAI->color = new uint16_t[SAI->nr * SAI->nc * 3]();

            if (SAI->n_references > 0) {

                printf("View prediction for view %03d_%03d\n", SAI->c, SAI->r);

                uint16_t **warped_texture_views = nullptr;
                uint16_t **warped_depth_views = nullptr;
                float **DispTargs = nullptr;

                init_warping_arrays(
                    SAI->n_references,
                    warped_texture_views,
                    warped_depth_views,
                    DispTargs,
                    SAI->nr,
                    SAI->nc,
                    SAI->ncomp);

                forward_warp_texture_references(
                    LF,
                    SAI,
                    warped_texture_views,
                    warped_depth_views,
                    DispTargs);

                merge_texture_views(
                    SAI,
                    LF,
                    warped_texture_views,
                    DispTargs);

                clean_warping_arrays(
                    SAI->n_references,
                    warped_texture_views,
                    warped_depth_views,
                    DispTargs);

                if (SAI->Ms > 0 && SAI->NNt > 0) {
        
                    /* OBTAIN SEGMENTATION*/
                    segmentation seg = makeSegmentation(SAI, n_seg_iterations);

                    uint16_t *original_color_view = read_input_ppm(
                        SAI->path_input_ppm,
                        SAI->nr,
                        SAI->nc,
                        SAI->ncomp,
                        bpc,
                        SAI->colorspace);

                    SAI->sparse_filters.clear();

                    /*READ DECODED REFERENCE VIEWS*/
                    for (int ikr = 0; ikr < SAI->n_references; ikr++) {

                        view *ref_view = LF + SAI->references[ikr];

                        int32_t tmp_w, tmp_r, tmp_ncomp;

                        aux_read16PGMPPM(
                            ref_view->path_internal_colorspace_out_ppm,
                            tmp_w,
                            tmp_r,
                            tmp_ncomp,
                            ref_view->color);

                    }

                    for (int32_t icomp = 0; icomp < SAI->nc_sparse; icomp++) {

                        std::vector<std::vector<uint16_t>> padded_regressors;

                        padded_regressors.push_back(
                            padArrayUint16_t_vec(
                                SAI->color + SAI->nr*SAI->nc*icomp,
                                SAI->nr,
                                SAI->nc,
                                SAI->NNt));

                        for (int ikr = 0; ikr < SAI->n_references; ikr++) {

                            view *ref_view = LF + SAI->references[ikr];

                            padded_regressors.push_back(
                                padArrayUint16_t_vec(
                                    ref_view->color + SAI->nr*SAI->nc*icomp,
                                    SAI->nr,
                                    SAI->nc,
                                    SAI->NNt));

                        }

                        uint16_t *padded_icomp_orig =
                            padArrayUint16_t(
                                original_color_view + SAI->nr*SAI->nc*icomp,
                                SAI->nr,
                                SAI->nc,
                                SAI->NNt);

                        for (int ir = 1;
                            ir <= seg.number_of_regions;
                            ir++)
                        {

                            SAI->sparse_filters.push_back(getGlobalSparseFilter_vec_reg(
                                padded_icomp_orig,
                                padded_regressors,
                                seg.seg,
                                ir,
                                SAI->nr + 2 * SAI->NNt,
                                SAI->nc + 2 * SAI->NNt,
                                SAI->NNt,
                                SAI->Ms,
                                SPARSE_BIAS_TERM));
                        }

                        delete[](padded_icomp_orig);

                        /* exit(0);*/
                    }

                    /* APPLY FILTER */

                    uint16_t *sp_filtered_image_padded =
                        new uint16_t[(SAI->nr + 2 * SAI->NNt)*(SAI->nc + 2 * SAI->NNt)*SAI->ncomp]();

                    std::vector<uint16_t> sp_filtered_image(
                        SAI->color, 
                        SAI->color+ SAI->nr*SAI->nc*SAI->ncomp );

                    int ee = 0;

                    for (int32_t icomp = 0; icomp < nc_sparse; icomp++) {

                        std::vector<std::vector<uint16_t>> padded_regressors;

                        padded_regressors.push_back(
                            padArrayUint16_t_vec(
                                SAI->color + SAI->nr*SAI->nc*icomp,
                                SAI->nr,
                                SAI->nc,
                                SAI->NNt));

                        if (SP_B) {
                            for (int ikr = 0; ikr < SAI->n_references; ikr++) {

                                view *ref_view = LF + SAI->references[ikr];

                                padded_regressors.push_back(
                                    padArrayUint16_t_vec(
                                        ref_view->color + SAI->nr*SAI->nc*icomp,
                                        SAI->nr,
                                        SAI->nc,
                                        SAI->NNt));

                            }
                        }

                        std::vector<double> filtered_icomp( (SAI->nr + 2 * SAI->NNt)*(SAI->nc + 2 * SAI->NNt), 0);

                        for (int ir = 1;
                            ir <= seg.number_of_regions;
                            ir++)
                        {

                            quantize_and_reorder_spfilter(
                                SAI->sparse_filters.at(ee));

                            dequantize_and_reorder_spfilter(
                                SAI->sparse_filters.at(ee));

                            applyGlobalSparseFilter_vec_reg(
                                padded_regressors,
                                seg.seg,
                                ir,
                                SAI->nr + 2 * SAI->NNt,
                                SAI->nc + 2 * SAI->NNt,
                                SAI->Ms,
                                SAI->NNt,
                                SPARSE_BIAS_TERM,
                                SAI->sparse_filters.at(ee++).filter_coefficients,
                                filtered_icomp);

                        }

                        for (int32_t iij = 0; iij < (SAI->nr + 2 * SAI->NNt)*(SAI->nc + 2 * SAI->NNt); iij++) {

                            double mmax = static_cast<double>((1 << BIT_DEPTH) - 1);

                            filtered_icomp[iij] =
                                clip(filtered_icomp[iij], 0.0, mmax);

                            sp_filtered_image_padded[iij + (SAI->nr + 2 * SAI->NNt)*(SAI->nc + 2 * SAI->NNt)*icomp] =
                                static_cast<uint16_t>(floor(filtered_icomp[iij] + 0.5));

                        }

                        uint16_t *cropped_icomp =
                            cropImage(sp_filtered_image_padded + (SAI->nr + 2 * SAI->NNt)*(SAI->nc + 2 * SAI->NNt)*icomp,
                            (SAI->nr + 2 * SAI->NNt),
                                (SAI->nc + 2 * SAI->NNt),
                                SAI->NNt);

                        memcpy(
                            sp_filtered_image.data() + SAI->nr*SAI->nc*icomp,
                            cropped_icomp,
                            sizeof(uint16_t)*SAI->nr*SAI->nc);

                        delete[](cropped_icomp);

                    }

                    SAI->number_of_sp_filters = SAI->sparse_filters.size();

                    /* CLEAN */

                    for (int ikr = 0; ikr < SAI->n_references; ikr++) {

                        view *ref_view = LF + SAI->references[ikr];

                        delete[](ref_view->color);
                        ref_view->color = nullptr;

                    }

                    delete[](sp_filtered_image_padded);

                    double psnr_without_sparse = PSNR(
                        original_color_view,
                        SAI->color,
                        SAI->nr,
                        SAI->nc,
                        nc_sparse,
                        (1 << 10) - 1);

                    double psnr_with_sparse = PSNR(
                        original_color_view,
                        sp_filtered_image.data(),
                        SAI->nr,
                        SAI->nc,
                        nc_sparse,
                        (1 << bpc) - 1);

                    if (psnr_with_sparse > psnr_without_sparse) {

                        memcpy(
                            SAI->color,
                            sp_filtered_image.data(),
                            sizeof(uint16_t)*SAI->nr*SAI->nc*SAI->ncomp);

                        SAI->use_global_sparse = true;

                    }

                    delete[](original_color_view);

                }

            }

            /* write raw prediction to .ppm */

            aux_write16PGMPPM(
                SAI->path_raw_prediction_at_encoder_ppm,
                SAI->nc,
                SAI->nr,
                3,
                SAI->color);

            delete[](SAI->color);
            SAI->color = nullptr;

        }

        /* get residue for all views at level=hlevel,
        AFTER THIS LOOP,
        YOU CAN FIND ALL RESIDUAL IMAGES IN directories "outputdir/residual/RAW/<level>"
        */
        for (int32_t ii = 0; ii < view_indices.size(); ii++) {

            view *SAI = LF + view_indices.at(ii);

            if (SAI->residual_rate_color > 0) {

                printf("Obtaining texture residual for view %03d_%03d\n", SAI->c, SAI->r);

                aux_read16PGMPPM(
                    SAI->path_raw_prediction_at_encoder_ppm,
                    SAI->nc,
                    SAI->nr,
                    SAI->ncomp,
                    SAI->color);

                uint16_t *original_color_view = read_input_ppm(
                    SAI->path_input_ppm,
                    SAI->nr,
                    SAI->nc,
                    SAI->ncomp,
                    bpc,
                    SAI->colorspace);

                double *residual_image_double = get_residual(
                    original_color_view,
                    SAI->color,
                    SAI->nr,
                    SAI->nc,
                    3);

                delete[](original_color_view);

                SAI->residual_image = quantize_residual(
                    residual_image_double,
                    SAI->nr,
                    SAI->nc,
                    SAI->ncomp,
                    bpc,
                    Q,
                    offset);

                delete[](residual_image_double);

                /* write raw quantized residual to .ppm */

                aux_write16PGMPPM(
                    SAI->path_raw_texture_residual_at_encoder_ppm,
                    SAI->nc,
                    SAI->nr,
                    SAI->ncomp,
                    SAI->residual_image);

                delete[](SAI->color);
                SAI->color = nullptr;
                delete[](SAI->residual_image);
                SAI->residual_image = nullptr;

            }

        }

        /* encode residual images using HEVC for all views at level=hlevel
        which have texture residual rate >0,
        here we can substitute HEVC with MuLE etc*/

        /*make scan order "serpent" in vector "hevc_i_order" */
        if ((LF + view_indices.at(0))->residual_rate_color > 0) {

            std::vector<int32_t> hevc_i_order =
                getScanOrder(LF, view_indices);

            std::vector< std::vector<uint16_t>> YUV_444_SEQ;

            /*padding to mincusize*/

            const int32_t mincusize = 8;

            const int32_t VERP = (mincusize - LF->nr%mincusize);
            const int32_t HORP = (mincusize - LF->nc%mincusize);

            int32_t nr1 = LF->nr + VERP;
            int32_t nc1 = LF->nc + HORP;

            for (int32_t ii = 0; ii < view_indices.size(); ii++) {

                view *SAI = LF + hevc_i_order.at(ii);

                /* ------------------------------
                TEXTURE RESIDUAL ENCODING STARTS
                ------------------------------*/

                printf("Encoding texture residual for view %03d_%03d\n", SAI->c, SAI->r);

                aux_read16PGMPPM(
                    SAI->path_raw_texture_residual_at_encoder_ppm,
                    SAI->nc,
                    SAI->nr,
                    SAI->ncomp,
                    SAI->residual_image);

                std::vector<uint16_t> paddedi = padArrayUint16_t_for_HM(
                    SAI->residual_image,
                    SAI->nr,
                    SAI->nc,
                    SAI->ncomp,
                    HORP,
                    VERP);

                YUV_444_SEQ.push_back(paddedi);

                delete[](SAI->residual_image);
                SAI->residual_image = nullptr;

            }

            view *SAI0 = LF + hevc_i_order.at(0);

            writeYUV444_seq_to_disk(
                YUV_444_SEQ,
                SAI0->encoder_raw_output_444);

            /* ------------------------------
            TEXTURE RESIDUAL ENCODING STARTS
            ------------------------------*/

            /* encode HM, YUV444 -> .hevc (any YUV format) */

            int32_t QPfinal = SAI0->preset_QP;

            long(*hevc_encoder)(
                const char *,
                const char *,
                YUV_FORMAT,
                const int32_t,
                const int32_t,
                const int32_t,
                const int32_t,
                const char *,
                const char *,
                const char *);

            if (USE_KVAZAAR) {

                hevc_encoder = &encodeKVAZAAR;

                /*super cumbersome YUV transform to suite kvazaar ...*/

                if ( YUVTYPE==YUV420 && nc_color_ref>1)
                {

                    std::vector<std::vector<uint16_t>> yuv420_seq = convertYUVseqTo420(
                        SAI0->encoder_raw_output_444,
                        YUV444,
                        nr1,
                        nc1,
                        YUV_444_SEQ.size());

                    writeYUV420_seq_to_disk(
                        yuv420_seq,
                        SAI0->encoder_raw_output_444);

                }
                else if( nc_color_ref<2 || (SAI0->level>1 && YUVTYPE==YUV400) )
                {

                    std::vector<std::vector<uint16_t>> yuv400_seq = convertYUVseqTo400(
                        SAI0->encoder_raw_output_444,
                        YUV444,
                        nr1,
                        nc1,
                        YUV_444_SEQ.size());

                    writeYUV400_seq_to_disk(
                        yuv400_seq,
                        SAI0->encoder_raw_output_444);

                }

            }
            else
            {
                hevc_encoder = &encodeHM;
            }

            if (SAI0->preset_QP < 0) {

                std::vector<double> bpps;
                std::vector<int32_t> QPs;

                int32_t QPstep = 3;

                for (int32_t QP = 0; QP <= 51; QP += QPstep) {

                    //for (int32_t QP = 25; QP = 25; QP=25 ) {
                    long bytes_hevc = hevc_encoder(
                        SAI0->encoder_raw_output_444,
                        SAI0->hevc_texture,
                        hlevel > 1 ? YUVTYPE: (nc_color_ref>1 ? YUV420 : YUV400),
                        QP,
                        YUV_444_SEQ.size(),
                        nc1,
                        nr1,
                        SAI0->decoder_raw_output_YUV,
                        setup.hm_encoder.c_str(),
                        setup.hm_cfg.c_str()); /*transpose for nr,nc*/

                    double bpphevc =
                        double(bytes_hevc * 8) / double((LF->nr*LF->nc*view_indices.size()));

                    printf("\nQP=%d\tbpp=\t%f\n", QP, bpphevc);

                    bpps.push_back(bpphevc);
                    QPs.push_back(QP);

                    if (*(bpps.end() - 1) < SAI0->residual_rate_color) {
                        break;
                    }

                }

                SAI0->bpp_range = std::vector<double>(bpps);
                SAI0->QP_range = std::vector<int32_t>(QPs);

                /*use first one*/
                if (bpps.size() == 1) {
                    QPfinal = QPs.at(0);
                }
                /*use last one*/
                else if (*(bpps.end() - 1) > SAI0->residual_rate_color) {
                    QPfinal = *(QPs.end() - 1);
                }
                /*interpolate*/
                else {

                    double dx = *(bpps.end() - 1) - *(bpps.end() - 2);
                    double dy = QPstep;// *(QPs.end() - 1) - *(QPs.end() - 2);

                    double diffx = SAI0->residual_rate_color - *(bpps.end() - 2);

                    QPfinal = round(double(*(QPs.end() - 2) + diffx*(dy / dx)));

                }
            }

            long bytes_hevc = hevc_encoder(
                SAI0->encoder_raw_output_444,
                SAI0->hevc_texture,
                hlevel>1 ? YUVTYPE : (nc_color_ref>1 ? YUV420 : YUV400),
                QPfinal,
                YUV_444_SEQ.size(),
                nc1,
                nr1,
                SAI0->decoder_raw_output_YUV,
                setup.hm_encoder.c_str(),
                setup.hm_cfg.c_str()); /*transpose for nr,nc*/

            double bpphevc =
                double(bytes_hevc * 8) / double((LF->nr*LF->nc*view_indices.size()));

            printf("\nFinal QP=%d\tbpp=\t%f\n", QPfinal, bpphevc);

            for (int32_t ii = 0; ii < view_indices.size(); ii++) {

                view *SAI = LF + hevc_i_order.at(ii);

                SAI->finalQP = QPfinal;

                SAI->real_rate_texture = bpphevc;

                SAI->QP_range = SAI0->QP_range;
                SAI->bpp_range = SAI0->bpp_range;

            }

            /* ------------------------------
            TEXTURE RESIDUAL ENCODING ENDS
            ------------------------------*/


            /* ------------------------------
            TEXTURE RESIDUAL DECODING STARTS
            ------------------------------*/

            /* decode HM, .hevc (any YUV format) -> (any YUV format)  */

            int32_t status = decodeHM(
                SAI0->hevc_texture,
                SAI0->decoder_raw_output_YUV,
                setup.hm_decoder.c_str());

            /*convert (any YUV format) -> YUV444 */

            std::vector<std::vector<uint16_t>> YUV444_dec = convertYUVseqTo444(
                SAI0->decoder_raw_output_YUV,
                hlevel>1 ? YUVTYPE : (nc_color_ref>1 ? YUV420 : YUV400),
                nr1,
                nc1,
                hevc_i_order.size());

            /* write PPM back to correct places, YUV444 -> .ppm */

            for (int32_t ii = 0; ii < view_indices.size(); ii++) {

                view *SAI = LF + hevc_i_order.at(ii);

                uint16_t *cropped = cropImage_for_HM(
                    YUV444_dec.at(ii).data(),
                    nr1,
                    nc1,
                    SAI->ncomp,
                    HORP,
                    VERP);

                aux_write16PGMPPM(
                    SAI->path_raw_texture_residual_at_decoder_ppm,
                    SAI->nc,
                    SAI->nr,
                    SAI->ncomp,
                    cropped);

                SAI->has_color_residual = true;

                delete[](cropped);

            }

            /* ------------------------------
            TEXTURE RESIDUAL DECODING ENDS
            ------------------------------*/
        }



        for (int32_t ii = 0; ii < view_indices.size(); ii++) {

            view *SAI = LF + view_indices.at(ii);

            aux_read16PGMPPM(
                SAI->path_raw_prediction_at_encoder_ppm,
                SAI->nc,
                SAI->nr,
                SAI->ncomp,
                SAI->color);

            if (SAI->has_color_residual) {

                printf("Decoding texture residual for view %03d_%03d\n", SAI->c, SAI->r);

                uint16_t *decoded_residual_image;

                aux_read16PGMPPM(
                    SAI->path_raw_texture_residual_at_decoder_ppm,
                    SAI->nc,
                    SAI->nr,
                    SAI->ncomp,
                    decoded_residual_image);

                double *residual = dequantize_residual(
                    decoded_residual_image,
                    SAI->nr,
                    SAI->nc,
                    SAI->ncomp,
                    bpc,
                    Q,
                    offset);

                uint16_t *corrected = apply_residual(
                    SAI->color,
                    residual,
                    SAI->nr,
                    SAI->nc,
                    SAI->ncomp,
                    bpc);

                /* update SAI->color to contain
                corrected (i.e., prediction + residual) version*/
                memcpy(
                    SAI->color,
                    corrected,
                    sizeof(uint16_t)*SAI->nr*SAI->nc*SAI->ncomp);

                delete[](residual);
                delete[](corrected);
                delete[](decoded_residual_image);

            }

            /*WRITE result to disk*/

            /*writing .ppm in internal colorspace*/
            aux_write16PGMPPM(
                SAI->path_internal_colorspace_out_ppm,
                SAI->nc,
                SAI->nr,
                SAI->ncomp,
                SAI->color);

            /*writing .ppm in output colorspace,
            If we only encode luminance, we have luminance as the first
            component of the .ppm file !
            */
            write_output_ppm(
                SAI->color,
                SAI->path_out_ppm,
                SAI->nr,
                SAI->nc,
                nc_color_ref,
                bpc,
                SAI->colorspace);

            delete[](SAI->color);
            SAI->color = nullptr;

        }

    }

}

void encoder::write_bitstream() {

    printf("Writing encoder statistics file\n");

    write_statsfile();

    printf("Writing header information to codestream\n");

    uint8_t colorspace_enumerator;

    if (colorspace_LF.compare("RGB") == 0) {
        colorspace_enumerator = 0;
    }

    if (colorspace_LF.compare("YCbCr") == 0) {
        colorspace_enumerator = 1;
    }

    char path_out_LF_data[1024];
    sprintf(
        path_out_LF_data,
        "%s/%s",
        setup.output_directory.c_str(),
        "output.LF");

    FILE* output_LF_file = fopen(path_out_LF_data, "wb");

    int32_t n_bytes_prediction = 0;
    int32_t n_bytes_residual = 0;

    n_bytes_prediction += (int32_t)fwrite(
        &n_views_total,
        sizeof(int32_t),
        1,
        output_LF_file);
    n_bytes_prediction += (int32_t)fwrite(
        &LF->nr,
        sizeof(int32_t),
        1,
        output_LF_file)
        * sizeof(int32_t);  // needed only once per LF
    n_bytes_prediction += (int32_t)fwrite(
        &LF->nc, sizeof(int32_t),
        1,
        output_LF_file)
        * sizeof(int32_t);  //
    n_bytes_prediction += (int32_t)fwrite(
        &LF->min_inv_d,
        sizeof(uint16_t),
        1,
        output_LF_file) * sizeof(uint16_t);
    n_bytes_prediction += (int32_t)fwrite(
        &colorspace_enumerator,
        sizeof(uint8_t),
        1,
        output_LF_file) * sizeof(uint8_t);
    n_bytes_prediction += (int32_t)fwrite(
        &maxh,
        sizeof(int32_t),
        1,
        output_LF_file) * sizeof(int32_t);

    n_bytes_prediction += (int32_t)fwrite(
        &nc_sparse,
        sizeof(uint8_t),
        1,
        output_LF_file) * sizeof(int32_t);

    n_bytes_prediction += (int32_t)fwrite(
        &nc_merge,
        sizeof(uint8_t),
        1,
        output_LF_file) * sizeof(int32_t);

    n_bytes_prediction += (int32_t)fwrite(
        &SP_B,
        sizeof(uint8_t),
        1,
        output_LF_file) * sizeof(int32_t);

    n_bytes_prediction += (int32_t)fwrite(
        &nc_color_ref,
        sizeof(uint8_t),
        1,
        output_LF_file) * sizeof(int32_t);

    n_bytes_prediction += (int32_t)fwrite(
        &n_seg_iterations,
        sizeof(uint8_t),
        1,
        output_LF_file) * sizeof(int32_t);

    std::vector< bool > levels_already_written_to_codestream(maxh, 0);

    for (int32_t ii = 0; ii < n_views_total; ii++) {

        view *SAI = LF + ii;

        printf("Writing codestream for view %03d_%03d\n", SAI->c, SAI->r);

        viewHeaderToCodestream(
            n_bytes_prediction,
            SAI,
            output_LF_file);

        if (SAI->has_color_residual) {
            if (!levels_already_written_to_codestream.at(SAI->level - 1)) {
                writeResidualToDisk(
                    SAI->hevc_texture,
                    output_LF_file,
                    n_bytes_residual,
                    JP2_dict);
            }
            levels_already_written_to_codestream.at(SAI->level - 1) = true;
        }

        if (SAI->has_depth_residual) {
            writeResidualToDisk(
                SAI->jp2_residual_depth_path_jp2,
                output_LF_file,
                n_bytes_residual,
                JP2_dict);
        }

    }

    fclose(output_LF_file);
}
