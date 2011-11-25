#include "fdmdv.h"
#include <stdlib.h>
#include "complex.h"

void fdmdv::fdmdv_write_qpsk(int sym, int carrier, double* buffer, int len)
{
        double delta;
        double  ival, qval, shapeA;
        complex symbol;

				// Sanity check
				sym &= 3;

        switch (sym) {
        case 0:
                symbol = complex (1.0, 1.0); 
                break;
        case 1:
                symbol = complex (-1.0, 1.0);
                break;
        case 2:
                symbol = complex (1.0, -1.0);
                break;
        case 3:
                symbol = complex (-1.0, -1.0);
                break;
        }

        delta = 2.0 * M_PI * (get_txfreq_woffset() + CARRIER_BW * carrier) / samplerate;
        for (int i = 0; i < symbollen; i++) {

                shapeA = tx_shape[i];

                ival = shapeA * symbol.real();
                qval = shapeA * symbol.imag();

                buffer[i] += ival * cos(phaseacc) + qval * sin(phaseacc);

                phaseacc += delta;
                if (phaseacc > M_PI)
                        phaseacc -= 2.0 * M_PI;
        }
}

void fdmdv::fdmdv_write_bpsk(int sym, int carrier, double* buffer, int len)
{
        double delta;
        double  ival, qval, shapeA;
        complex symbol;

        if(sym == FDMDV_FRAME_1) bpsk_bit = ~bpsk_bit & 1;

        switch (bpsk_bit) {
        case 0:
                symbol = complex ( 0.0, -1.0 );
                break;
        case 1:
                symbol = complex ( 0.0, 1.0);
                break;
        }

        delta = 2.0 * M_PI * (get_txfreq_woffset() + CARRIER_BW * carrier) / samplerate;
        for (int i = 0; i < symbollen; i++) {

                shapeA = tx_shape[i];

                ival = shapeA * symbol.real();
                qval = shapeA * symbol.imag();

                buffer[i] += ival * cos(phaseacc) + qval * sin(phaseacc);

                phaseacc += delta;
                if (phaseacc > M_PI)
                        phaseacc -= 2.0 * M_PI;
        }
}
