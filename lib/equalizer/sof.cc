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

#include "sof.h"

using namespace gr::ieee802_11::equalizer;

const float sof::freq_corr_hi_cutoff = 0.9;
const float sof::freq_corr_lo_cutoff = 0.7;

void sof::equalize (gr_complex *in, int n, gr_complex *symbols, uint8_t *bits, boost::shared_ptr<gr::digital::constellation> mod) {

    dout << "SOF: symb #" << n << std::endl;
    if (n == 0) {
        d_encoding = BPSK_1_2;
        std::memcpy(d_H, in, 64 * sizeof(gr_complex));
    } else if (n == 1) {
        double signal = 0;
        double noise = 0;
        for (int i = 0; i < 64 ; ++i) {
            noise += std::pow(std::abs(d_H[i] - in[i]), 2);
            signal += std::pow(std::abs(d_H[i] + in[i]), 2);
        }
        d_snr_db = 10 * std::log10(signal / noise - 1);

        for (int i = 0; i < 64; ++i) {
            d_H[i] += in[i];
            d_H[i] /= LONG[i] * gr_complex(2, 0);
        }
        filter_channel_estimate(d_H);
    } else {
        gr_complex H_new[64];
        gr_complex p = POLARITY[(n-2) % 127];

        uint8_t bits_demapped[48];   // current symbol bits to decode
        uint8_t err_flags[48] = {0};


        if (n == 2) {
            compute_conf_lvl();
            dout << "conf_lvl: " << conf_lvl << std::endl;
        }


        if (n == 3) {
            compute_corr_dbth();
            set_correction_flag();
            dout << "flag ? corr_dbth: " << correction_flag << " " << d_snr_dbthm << std::endl;
        }

        H_new[11] = in[11] *  p;
        H_new[25] = in[25] *  p;
        H_new[39] = in[39] *  p;
        H_new[53] = in[53] * -p;

        int c = 0;
        for (int i = 6; i <= 58 ; ++i) {
            if( (i == 11) || (i == 25) || (i == 32) || (i == 39) || (i == 53) ) {
                continue;
            } else {
                symbols[c] = in[i] / d_H[i];
                bits_demapped[c] = mod->decision_maker(&symbols[c]);
                ++c;
            }
        }

//        dout << "Before decoding..." << std::endl;
        decode_curr_ofdm_symbol(bits_demapped, bits);
//        dout << "After decoding..." << std::endl;

        // TODO-AS: Remove the line below when decoding of the current equalized OFDM symbol is fixed
        std::memcpy(bits, bits_demapped, 48 * sizeof(uint8_t));

        for (int i = 0; i < 48; ++i) {
            if (bits_demapped[i] != bits[i]) {
                err_flags[i] = 1;
            }
        }

        c = 0;
        for (int i = 6; i <= 58 ; ++i) {
            if( (i == 11) || (i == 25) || (i == 32) || (i == 39) || (i == 53) ) {
                continue;
            } else {
                gr_complex point;
                mod->map_to_points(bits[c], &point);
                H_new[i] = in[i] / point;
//                dout << i << ": " << in[i] << d_H[i] << H_new[i] << point << "<->" << unsigned(bits[c]) << "<<->>" << symbols[c] << std::endl;
                c++;
            }
        }

        filter_channel_estimate(H_new);
        reestimate_channel_response(H_new);
        if (correction_flag) {
            for (int i = 6; i <= 58; ++i) {
                if ((i == 11) || (i == 25) || (i == 32) || (i == 39) || (i == 53)) {
                    continue;
                } else {
                    if (err_flags[i]) {
                        d_H[i] = H_new[i];
                    }
                }
            }
        }
    }
}

void sof::decode_curr_ofdm_symbol(uint8_t *symbols_in, uint8_t *symbols_out) {
    // 1 OFDM symbol has a max of 6 bpsc -> 48 * 6 = 288 encoded bits
    uint8_t rx_bits_in[288] = {0};
    uint8_t rx_bits_deinter[288] = {0};
    uint8_t rx_bits_dec[288] = {0};
    uint8_t *rx_bits_dec_ptr;
    uint8_t rx_bits_enc[288] = {0};
    uint8_t rx_bits_pct[288] = {0};
    uint8_t rx_bits_int[288] = {0};

    dout << "Curr OFDM symb enc: " << d_encoding << std::endl;

    ofdm_param ofdm(d_encoding);

    // split symbols into individual bits first
    for(int i = 0; i < 48; i++) {
        for(int k = 0; k < ofdm.n_bpsc; k++) {
            rx_bits_in[i * ofdm.n_bpsc + k] = !!(symbols_in[i] & (1 << k));
            dout << "bit #" << i * ofdm.n_bpsc + k << ": " << unsigned(rx_bits_in[i * ofdm.n_bpsc + k]) << std::endl;
        }
    }

    // setup the OFDM symbol properties & propagate info to frame params
    frame_param curr_symb_param(ofdm, std::floor(ofdm.n_dbps/8));
    // 1 OFDM symbol only
    curr_symb_param.n_sym = 1;
    // total data bits is just the n_dbps from the OFDM parameters
    curr_symb_param.n_data_bits = ofdm.n_dbps;
    // padding number is 0 -> do not care about the content
//    curr_symb_param.n_pad = 0;
    curr_symb_param.n_pad = curr_symb_param.n_data_bits - 8 * curr_symb_param.psdu_size;
    // number of encoded bites is just the n_cbps from the OFDM params
    curr_symb_param.n_encoded_bits = ofdm.n_cbps;

    viterbi_decoder vdecoder;

    // TODO-AS: Fix the decoding and recoding of the bits and packing into symbols.
    // This seems not work somewhere in between decoding, coding and puncturing.
    // I suspect there may either be something wrong with the frame_param I set up, or
    // with decoding and getting the results - do I copy the right bits into the end
    // decoding buffer ?


    // first deinterleave the coded received bits
    interleave((char *)rx_bits_in, (char *)rx_bits_deinter, curr_symb_param, ofdm, true);
    // next perform error correction decoding
    rx_bits_dec_ptr = vdecoder.decode(&ofdm, &curr_symb_param, rx_bits_deinter);
    std::memcpy(rx_bits_dec, rx_bits_dec_ptr, ofdm.n_dbps * sizeof(uint8_t));
    // given corrected errors of the channel reencode at given rate
    convolutional_encoding((char *)rx_bits_dec, (char *)rx_bits_enc, curr_symb_param);
    puncturing ((char *)rx_bits_enc, (char *)rx_bits_pct, curr_symb_param, ofdm);
    // reinterleave the bitstream obtained
    interleave((char *)rx_bits_pct, (char *)rx_bits_int, curr_symb_param, ofdm);
    // split the bits into subcarrier symbols
    split_symbols((char *)rx_bits_int, (char *)symbols_out, curr_symb_param, ofdm);
}

void sof::set_ofdm_frame_params(const Encoding ofdm_enc, const int frame_len) {
    d_encoding = ofdm_enc;
    d_frame_length = frame_len;
}

double sof::get_snr() {
    return d_snr_db;
}

void sof::compute_corr_dbth() {
    switch(d_encoding) {
        case BPSK_1_2:
        case BPSK_3_4:
            // 10*log10(5^2*(2/2)^2)
            d_snr_dbthm = 13.9794;
            break;
        case QPSK_1_2:
        case QPSK_3_4:
            // 10*log10(5^2*(2*sqrt(2)/2)^2)
            d_snr_dbthm = 16.9897;
            break;
        case QAM16_1_2:
        case QAM16_3_4:
            // 10*log10(5^2*(2*sqrt(10)/2)^2)
            d_snr_dbthm = 23.9794;
            break;
        case QAM64_2_3:
        case QAM64_3_4:
            // 10 * log10(5^2*(2*sqrt(42)/2)^2)
            d_snr_dbthm = 30.2119;
            break;
        default:
            throw std::runtime_error("Unknown OFDM modulation");
    }
}

void sof::compute_conf_lvl() {
    if (d_snr_db <= 0.0) {
        conf_lvl = 1.0;
    } else if (d_snr_db <= d_snr_dbthg) {
        conf_lvl = 1.0 + d_snr_db / d_snr_dbthg;
    } else {
        conf_lvl = 2.0;
    }
}

void sof::compute_freq_corr(gr_complex *interm_ch_resp) {
    double res_idx0 = 0.0;

    for (int i = 6; i <= 58; ++i) {
        res_idx0 += std::pow(std::abs(interm_ch_resp[i]),2);
    }

    freq_corr_val[0] = 1.0;

    for (int s = 1; s <= 5; ++s) {
        gr_complexd res_curr(0, 0);
        for (int i = 6; i <= 58-s; ++i) {
            res_curr += (interm_ch_resp[i] * std::conj(interm_ch_resp[i+s]));
        }
        freq_corr_val[s] = (double)std::abs(res_curr) / res_idx0;
//        dout << "freq_corr " << s << ": " << freq_corr_val[s] << std::endl;
    }
}

void sof::compute_time_corr(gr_complex *interm_ch_resp) {
    gr_complexd res_num(0, 0);
    double res_den = 0.0;

    for (int i = 6; i <= 58; ++i) {
//        dout << "tcc " << i << ": " << d_H[i] << interm_ch_resp[i] << std::endl;
        res_num += (d_H[i] * std::conj(interm_ch_resp[i]));
        res_den += (std::pow(std::abs(d_H[i]), 2) +
                    std::pow(std::abs(interm_ch_resp[i]), 2));
    }

    time_corr_val = std::abs(res_num) / res_den;

    dout << "time corr: " << time_corr_val << std::endl;
}

void sof::set_correction_flag()
{
    correction_flag = false;
    if (d_snr_db > d_snr_dbthm) {
        correction_flag = true;
    }
}

float sof::map_cutoff_len_to_sigma(uint8_t cutoff_val) {
    cutoff_val = cutoff_val * 2 + 1;
    switch(cutoff_val) {
        case 1:
            return 0.46;
        case 3:
            return 0.93;
        case 5:
            return 1.39;
        case 7:
            return 1.86;
        case 9:
            return 2.32;
//        case 11:
//            return 2.79;
        default:
//            return 2.79;
            return 2.32;
    }
}


void sof::filter_channel_estimate(gr_complex *ch_resp) {
    gr_complex ch_resp_replica[64];

    // make sure DC is interpolation of close by subcarriers
    ch_resp[32] = gr_complex(0.5, 0) * (ch_resp[31] + ch_resp[33]);

    compute_freq_corr(ch_resp);
    std::memcpy(ch_resp_replica, ch_resp, 64 * sizeof(gr_complex));

    uint8_t cutoff_lo = 0, cutoff_hi = 0;

    for (uint8_t i = 1; i < 6; ++i) {
        if ((freq_corr_val[i] > freq_corr_lo_cutoff) && (i > cutoff_lo)) {
            cutoff_lo = i;
            dout << "cutoff_lo" << unsigned(cutoff_lo) << std::endl;
        }
        if ((freq_corr_val[i] > freq_corr_hi_cutoff) && (i > cutoff_hi)) {
            cutoff_hi = i;
            dout << "cutoff_hi" << unsigned(cutoff_hi) << std::endl;
        }
    }

    float a, b;
    b = map_cutoff_len_to_sigma(cutoff_lo);
    a = map_cutoff_len_to_sigma(cutoff_hi) - b;

    if (d_snr_db <= 0) {
        sigma = b;
    } else if (d_snr_db <= d_snr_dbths) {
        sigma = a * (d_snr_db / d_snr_dbths) + b;
    } else {
        // no filtering, so just return - nothing to do here
        return;
    }

    dout << "sigma: " << sigma << std::endl;

    sigma *= sigma;
    float truncate_th = 0.01, taps_sum = 0.0;
//    uint8_t flen = 11;
    uint8_t flen = 9;
    int8_t start_idx = -flen / 2;
    ftaps.clear();

    for (int i = 0; i < flen; ++i) {
        float exparg = (float)start_idx + (float)i;
        exparg *= exparg / sigma;
        exparg = std::exp(-exparg);
        if (exparg >= truncate_th) {
            ftaps.push_back(exparg);
            taps_sum += exparg;
        }
    }

    flen = ftaps.size();
    for (int i = 0 ; i < flen ; ++i) {
        ftaps[i] /= taps_sum;
//        dout << "ftap " << i << ": " << ftaps[i] << std::endl;
    }

    int offset = (flen - 1) / 2;
    gr_complex scaling;
    gr_complex sum;
    for (int i = 6 ; i <= 58 ; ++i) {
        sum = gr_complex(0, 0);
        scaling = gr_complex(0, 0);
        for (int j = -offset ; j <= offset ; ++j) {
            int curr_idx = i+j;
            if ((curr_idx >= 6) && (curr_idx<=58)) {
                scaling += gr_complex(ftaps[offset + j], 0);
                sum += ch_resp_replica[i + j] * ftaps[offset + j];
            }
        }
//        dout << "i: " << i << " scaling: " << scaling << std::endl;
        ch_resp[i] = sum / scaling;
//        dout << "i: " << ch_resp_replica[i] << " vs " << ch_resp[i] << std::endl;
    }
}

void sof::reestimate_channel_response (gr_complex *interm_ch_resp) {
    compute_time_corr(interm_ch_resp);
    gr_complex fact(conf_lvl * time_corr_val, 0);
    dout << "update rate: " << fact << std::endl;
    for (int i = 6 ; i <= 58; ++i) {
        d_H[i] = d_H[i] + fact * (interm_ch_resp[i] - d_H[i]);
    }
}
