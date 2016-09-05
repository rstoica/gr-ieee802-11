/*
 * Copyright (C) 2016 Andrei Stoica <r.stoica@jacobs-university.de>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef INCLUDED_IEEE802_11_EQUALIZER_SOF_H
#define INCLUDED_IEEE802_11_EQUALIZER_SOF_H

#include "base.h"
#include "utils.h"
#include "viterbi_decoder.h"
#include "ieee802-11/constellations.h"
#include <vector>

namespace gr {
namespace ieee802_11 {
namespace equalizer {

class sof: public base {
public:
    virtual void equalize(gr_complex *in, int n, gr_complex *symbols, uint8_t *bits, boost::shared_ptr<gr::digital::constellation> mod);
    double get_snr();
    void set_ofdm_frame_params(const Encoding ofdm_enc, const int frame_len);

private:
    gr_complex d_H[64];
    int d_frame_length;
    Encoding d_encoding;
//    viterbi_decoder d_decoder;
    const bool d_debug = true;

    const double d_snr_dbths = 30.0;
    const double d_snr_dbthg = 30.0;
    double d_snr_dbthm;

    double d_snr_db;

    double conf_lvl;
    double freq_corr_val[6];
    double time_corr_val;
    bool correction_flag;
    double sigma;
    std::vector<float> ftaps;
    //double corr_dbth;

    static const float freq_corr_lo_cutoff;
    static const float freq_corr_hi_cutoff;

    constellation_bpsk::sptr d_bpsk;
    constellation_qpsk::sptr d_qpsk;
    constellation_16qam::sptr d_16qam;
    constellation_64qam::sptr d_64qam;

    void decode_curr_ofdm_symbol(uint8_t *symbols_in, uint8_t *symbols_out);
    void compute_corr_dbth();
    void compute_conf_lvl();
    void compute_freq_corr(gr_complex *interm_ch_resp);
    void compute_time_corr(gr_complex *interm_ch_resp);
    void set_correction_flag();
    float map_cutoff_len_to_sigma(uint8_t cutoff_val);
    void filter_channel_estimate(gr_complex *ch_resp);
    void reestimate_channel_response(gr_complex *interm_ch_resp);
};

} /* namespace channel_estimation */
} /* namespace ieee802_11 */
} /* namespace gr */

#endif /* INCLUDED_IEEE802_11_EQUALIZER_STA_H */
