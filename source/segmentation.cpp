#include "segmentation.hh"
#include "medianfilter.hh"
#include "ppm.hh"

#include <algorithm>
#include <iterator>

#define SAVE_SEGMENTATION true

segmentation makeSegmentation(
    view* SAI, 
    const int32_t n_seg_iterations)
{

    segmentation seg;

    if (n_seg_iterations > 0) {

        int32_t tmp_w, tmp_r, tmp_ncomp;
        aux_read16PGMPPM(
            SAI->path_out_pgm,
            tmp_w,
            tmp_r,
            tmp_ncomp,
            SAI->depth);

        uint16_t *img16_padded =
            padArrayUint16_t(
                SAI->depth,
                SAI->nr,
                SAI->nc,
                SAI->NNt);

        delete[](SAI->depth);
        SAI->depth = nullptr;

        std::vector<uint16_t> img16_padded_v(
            img16_padded,
            img16_padded + (SAI->nr + 2 * SAI->NNt)*(SAI->nc + 2 * SAI->NNt));

        seg = normdispsegmentation(
            img16_padded_v,
            n_seg_iterations,
            SAI->nr + 2 * SAI->NNt,
            SAI->nc + 2 * SAI->NNt);

        std::vector<uint16_t> seg16(seg.seg.begin(), seg.seg.end());

        if (SAVE_SEGMENTATION)
        {
            aux_write16PGMPPM(SAI->path_output_seg,
                SAI->nc + 2 * SAI->NNt,
                SAI->nr + 2 * SAI->NNt,
                1,
                &seg16[0]);
        }
    }
    else {
        seg.seg = std::vector<int32_t>((SAI->nr + 2 * SAI->NNt)*(SAI->nc + 2 * SAI->NNt), 1);
        seg.number_of_regions = 1;
    }
    return seg;
}

segmentation normdispsegmentation(
    const std::vector<uint16_t> &img,
    const uint32_t iterations,
    const int32_t nr,
    const int32_t nc) {

    /*init segmentation with zeros*/
    std::vector<uint16_t> prev_iter(img.size(), 0);

    for (uint32_t iter = 0; iter < iterations; iter++) {

        std::vector<uint16_t> uniquevals(prev_iter.begin(), prev_iter.end());
        sort(uniquevals.begin(), uniquevals.end());
        std::vector<uint16_t>::iterator it;
        it = std::unique(uniquevals.begin(), uniquevals.end());
        /*unique values in previuos iteration*/
        uniquevals.resize(std::distance(uniquevals.begin(), it));

        std::vector<uint16_t> cur_iter(img.size(), 0);

        for (uint32_t ij = 0; ij < uniquevals.size(); ij++) {

            /*max value in segmentation*/
            uint16_t maxval = 
                *std::max_element(cur_iter.begin(), cur_iter.end());

            std::vector<uint16_t> vals;
            std::vector<uint32_t> indices;

            /*get median for current region*/
            for (uint32_t ii = 0; ii < prev_iter.size(); ii++) {

                if (prev_iter.at(ii) == uniquevals.at(ij)) {
                    vals.push_back(img.at(ii));
                    indices.push_back(ii);
                }

            }
            uint16_t medianval = getMedian(vals);
            
            /*split current region based on median*/
            for (uint32_t ii = 0; ii < indices.size(); ii++) {

                if (img.at(indices.at(ii)) >= medianval) {
                    cur_iter.at(indices.at(ii)) = maxval + 2;
                }
                else {
                    cur_iter.at(indices.at(ii)) = maxval + 1;
                }

            }

        }

        /*save current segmentation for next iteration*/
        prev_iter = cur_iter;

    }

    uint16_t maxval = *std::max_element(prev_iter.begin(), prev_iter.end());

    //prev_iter = std::vector<uint16_t>(prev_iter.size(), 1);

    segmentation seg;
    seg.seg = std::vector<int32_t>(prev_iter.begin(), prev_iter.end());
    seg.number_of_regions = maxval; /*zero is not a region*/

    return seg;

}