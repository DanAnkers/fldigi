// ----------------------------------------------------------------------------
// cw.cxx  --  morse code modem
//
// Copyright (C) 2006
//		Dave Freese, W1HKJ
//
// This file is part of fldigi.  Adapted from code contained in gmfsk source code 
// distribution.
//  gmfsk Copyright (C) 2001, 2002, 2003
//  Tomi Manninen (oh2bns@sral.fi)
//  Copyright (C) 2004
//  Lawrence Glaister (ve7it@shaw.ca)
//
// fldigi is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// fldigi is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with fldigi; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
// ----------------------------------------------------------------------------


#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <sys/time.h>

#include "cw.h"
#include "CWdialog.h"
#include "misc.h"
//#include "modeIO.h"
#include "configuration.h"

void cw::tx_init(cSound *sc)
{
	scard = sc;
	phaseacc = 0;
	lastsym = 0;
}

void cw::rx_init()
{
	cw_receive_state = RS_IDLE;	
	smpl_ctr = 0;				
	cw_rr_current = 0;			
	agc_peak = 0;	
	digiscope->mode(Digiscope::SCOPE);
	put_MODEstatus(mode);
}

void cw::init()
{
	modem::init();
	set_cwXmtWPM(sldrCWxmtWPM->value());
	rx_init();
}

cw::~cw() {
//	if (KeyLine) {
//		delete KeyLine;
//		KeyLine = (modeIO *)0;
//	}
	if (cwfilter) delete cwfilter;
	if (bitfilter) delete bitfilter;
	if (trackingfilter) delete trackingfilter;
}


cw::cw() : morse(), modem()
{
	double lp;

	mode = MODE_CW;
	freqlock = false;
	cw_speed = progdefaults.CWspeed;
	frequency = 800;
	tx_frequency = 800;
	bandwidth = progdefaults.CWbandwidth;
	
	samplerate = CWSampleRate;
	fragmentsize = CWMaxSymLen;

	cw_send_speed = cw_speed;
	cw_receive_speed = cw_speed;
//	cw_noise_spike_threshold = INITIAL_NOISE_THRESHOLD;
	cw_adaptive_receive_threshold = 2 * DOT_MAGIC / cw_speed;
	cw_noise_spike_threshold = cw_adaptive_receive_threshold / 4;
	
	memset(rx_rep_buf, 0, sizeof(rx_rep_buf));

// block of variables that get updated each time speed changes
	pipesize = 512;
	cwTrack = true;
	phaseacc = 0.0;
	pipeptr = 0;
	agc_peak = 1.0;

	lp = 0.5 * bandwidth / samplerate;
	cwfilter = new C_FIR_filter();
	cwfilter->init_lowpass (CW_FIRLEN, DEC_RATIO, lp);
	
	bitfilter = new Cmovavg(8);
	trackingfilter = new Cmovavg(TRACKING_FILTER_SIZE);

	sync_parameters();
	wf->Bandwidth ((int)bandwidth);

}

// sync_parameters()
// Synchronize the dot, dash, end of element, end of character, and end
// of word timings and ranges to new values of Morse speed, or receive tolerance.
 
void cw::sync_parameters()
{
	int lowerwpm, upperwpm;
// check if user changed the tracking or the cw default speed
	if (cw_send_speed != progdefaults.CWspeed ||
		cwTrack != progdefaults.CWtrack) {

		cw_send_speed = progdefaults.CWspeed;
		cw_send_dot_length = DOT_MAGIC / cw_send_speed;
		cw_send_dash_length = 3 * cw_send_dot_length;
		symbollen = (int)(1.0 * samplerate * cw_send_dot_length / USECS_PER_SEC);

		trackingfilter->reset();
		cw_adaptive_receive_threshold = (long int)trackingfilter->run(2 * cw_send_dot_length);
		
		put_cwRcvWPM(cw_send_speed);
	} else {
		cw_send_speed = progdefaults.CWspeed;
		cw_send_dot_length = DOT_MAGIC / cw_send_speed;
		cw_send_dash_length = 3 * cw_send_dot_length;
		symbollen = (int)(1.0 * samplerate * cw_send_dot_length / USECS_PER_SEC);
	}
	cwTrack = progdefaults.CWtrack;

// Receive parameters:
	lowerwpm = cw_send_speed - progdefaults.CWrange;
	upperwpm = cw_send_speed + progdefaults.CWrange;
	if (lowerwpm < CW_MIN_SPEED) lowerwpm = CW_MIN_SPEED;
	cw_lower_limit = 2 * DOT_MAGIC / upperwpm;
	cw_upper_limit = 2 * DOT_MAGIC / lowerwpm;

	if (cwTrack)
		cw_receive_speed = DOT_MAGIC / (cw_adaptive_receive_threshold / 2);
	else {
		cw_receive_speed = cw_send_speed;
		cw_adaptive_receive_threshold = 2 * cw_send_dot_length;
	}
		
// receive routines track speeds, but we put hard limits
// on the speeds here if necessary.
// may not need with new algorithm which limits tracking range
//	if (cw_receive_speed < CW_MIN_SPEED)
//		cw_receive_speed = CW_MIN_SPEED;
//	if (cw_receive_speed > CW_MAX_SPEED)
//		cw_receive_speed = CW_MAX_SPEED;
		
	cw_receive_dot_length = DOT_MAGIC / cw_receive_speed;

//	cw_adaptive_receive_threshold = 2 * cw_receive_dot_length;
	cw_receive_dash_length = 3 * cw_receive_dot_length;

	cw_noise_spike_threshold = cw_receive_dot_length / 4;

// Set the parameters in sync flag.
//	cw_in_sync = true;
}


//=======================================================================
// cw_update_tracking()
// This gets called everytime we have a dot dash sequence or a dash dot
// sequence. Since we have semi validated tone durations, we can try and
// track the cw speed by adjusting the cw_adaptive_receive_threshold variable.
// This is done with moving average filters for both dot & dash.
//=======================================================================

void cw::update_tracking(int idot, int idash)
{
	int dot, dash;
	if (idot > cw_lower_limit && idot < cw_upper_limit)
		dot = idot;
	else
		dot = cw_send_dot_length;
	if (idash > cw_lower_limit && idash < cw_upper_limit)
		dash = idash;
	else
		dash = cw_send_dash_length;
	
	cw_adaptive_receive_threshold = (long int)trackingfilter->run((dash + dot) / 2);
	sync_parameters();
}

//=======================================================================
//update_syncscope()
//Routine called to update the display on the sync scope display.
//For CW this is an o scope pattern that shows the cw data stream.
//=======================================================================
//
void cw::update_syncscope()
{
	int j;

	for (int i = 0; i < pipesize; i++) {
		j = (i + pipeptr) % pipesize;
		scopedata[i] = 0.1 + 0.8 * pipe[j] / agc_peak;
	}
	set_scope(scopedata, pipesize, false);

//	cwRcvWPM = cw_receive_speed;
//	put_cwRcvWPM(cwRcvWPM);
	put_cwRcvWPM(cw_receive_speed);
}



//=====================================================================
// cw_rxprocess()
// Called with a block (512 samples) of audio.
//=======================================================================

int cw::rx_process(double *buf, int len)
{
	complex z;
	double delta;
	double value;
	char *c;

// check if user changed filter bandwidth
	if (bandwidth != progdefaults.CWbandwidth) {
		bandwidth = progdefaults.CWbandwidth;
		cwfilter->init_lowpass (CW_FIRLEN, DEC_RATIO, 0.5 * bandwidth / samplerate);
		wf->Bandwidth ((int)bandwidth);
	}

// compute phase increment expected at our specific rx tone freq 
	delta = 2.0 * M_PI * frequency / samplerate;

	while (len-- > 0) {
		// Mix with the internal NCO 
		z = complex ( *buf * cos(phaseacc), *buf * sin(phaseacc) );
		buf++;
		phaseacc += delta;
		if (phaseacc > M_PI)
			phaseacc -= 2.0 * M_PI;
		if (cwfilter->run ( z, z )) {
		
// update the basic sample counter used for morse timing 
			smpl_ctr += DEC_RATIO;
// demodulate 
			value = z.mag();
			
			value = bitfilter->run(value);
// Compute a variable threshold value for tone 
// detection. Fast attack and slow decay.
			if (value > agc_peak)
				agc_peak = decayavg(agc_peak, value, 10.0);
			else
				agc_peak = decayavg(agc_peak, value, 800.0);

			metric = clamp(agc_peak * 1000.0 , 0.0, 100.0);
			display_metric(metric);
			
// save correlation amplitude value for the sync scope
			pipe[pipeptr] = value;
			pipeptr = (pipeptr + 1) % pipesize;
			if (pipeptr == pipesize - 1)
				update_syncscope();

			if (!squelchon || metric > squelch ) {
// upward trend means tone starting 
				if ((value > 0.66 * agc_peak) && (cw_receive_state != RS_IN_TONE))
					handle_event(CW_KEYDOWN_EVENT, NULL);
// downward trend means tone stopping 
				if ((value < 0.33 * agc_peak) && (cw_receive_state == RS_IN_TONE))
					handle_event(CW_KEYUP_EVENT, NULL);
			}
			if (handle_event(CW_QUERY_EVENT, &c) == CW_SUCCESS) {
				while (*c)
					put_rx_char(*c++);
			}
		}
	}

	return 0;
}

// ---------------------------------------------------------------------- 

// Compare two timestamps, and return the difference between them in usecs.
 
int cw::usec_diff(unsigned int earlier, unsigned int later)
{
// Compare the timestamps.
// At 4 WPM, the dash length is 3*(1200000/4)=900,000 usecs, and
// the word gap is 2,100,000 usecs.
	if (earlier >= later) {
		return 0;
	} else
		return (int) (((double) (later - earlier) * USECS_PER_SEC) / samplerate);
}


//=======================================================================
// handle_event()
//    high level cw decoder... gets called with keyup, keydown, reset and
//    query commands. 
//   Keyup/down influences decoding logic.
//    Reset starts everything out fresh.
//    The query command returns CW_SUCCESS and the character that has 
//    been decoded (may be '*',' ' or [a-z,0-9] or a few others)
//    If there is no data ready, CW_ERROR is returned.
//=======================================================================

int cw::handle_event(int cw_event, char **c)
{
	static int space_sent = true;	// for word space logic
	static int last_element = 0;	// length of last dot/dash
	int element_usec;				// Time difference in usecs

	switch (cw_event) {
	case CW_RESET_EVENT:
		sync_parameters();
		cw_receive_state = RS_IDLE;
		cw_rr_current = 0;			// reset decoding pointer
		smpl_ctr = 0;					// reset audio sample counter
		memset(rx_rep_buf, 0, sizeof(rx_rep_buf));
		break;
	case CW_KEYDOWN_EVENT:
// A receive tone start can only happen while we
// are idle, or in the middle of a character.
		if (cw_receive_state == RS_IN_TONE)
			return CW_ERROR;
// first tone in idle state reset audio sample counter
		if (cw_receive_state == RS_IDLE) {
			smpl_ctr = 0;
			memset(rx_rep_buf, 0, sizeof(rx_rep_buf));
			cw_rr_current = 0;
		}
// save the timestamp
		cw_rr_start_timestamp = smpl_ctr;
// Set state to indicate we are inside a tone.
		cw_receive_state = RS_IN_TONE;
		return CW_ERROR;
		break;
	case CW_KEYUP_EVENT:
// The receive state is expected to be inside a tone.
		if (cw_receive_state != RS_IN_TONE)
			return CW_ERROR;
// Save the current timestamp 
		cw_rr_end_timestamp = smpl_ctr;
		element_usec = usec_diff(cw_rr_start_timestamp, cw_rr_end_timestamp);
							 
// make sure our timing values are up to date
		sync_parameters();
// If the tone length is shorter than any noise cancelling 
// threshold that has been set, then ignore this tone.
		if (cw_noise_spike_threshold > 0
		    && element_usec < cw_noise_spike_threshold)
			return CW_ERROR;

// Set up to track speed on dot-dash or dash-dot pairs for this test to work, we need a dot dash pair or a 
// dash dot pair to validate timing from and force the speed tracking in the right direction. This method 
// is fundamentally different than the method in the unix cw project. Great ideas come from staring at the
// screen long enough!. Its kind of simple really ... when you have no idea how fast or slow the cw is... 
// the only way to get a threshold is by having both code elements and setting the threshold between them
// knowing that one is supposed to be 3 times longer than the other. with straight key code... this gets 
// quite variable, but with most faster cw sent with electronic keyers, this is one relationship that is 
// quite reliable. Lawrence Glaister (ve7it@shaw.ca)
		if (last_element > 0) {
// check for dot dash sequence (current should be 3 x last)
			if ((element_usec > 2 * last_element) &&
			    (element_usec < 4 * last_element)) {
				update_tracking(last_element, element_usec);
			}
// check for dash dot sequence (last should be 3 x current)
			if ((last_element > 2 * element_usec) &&
			    (last_element < 4 * element_usec)) {
				update_tracking(element_usec, last_element);
			}
		}
		last_element = element_usec;
// ok... do we have a dit or a dah?
// a dot is anything shorter than 2 dot times
		if (element_usec <= cw_adaptive_receive_threshold) {
			rx_rep_buf[cw_rr_current++] = CW_DOT_REPRESENTATION;
		} else {
// a dash is anything longer than 2 dot times
			rx_rep_buf[cw_rr_current++] = CW_DASH_REPRESENTATION;
		}
// We just added a representation to the receive buffer.  
// If it's full, then reset everything as it probably noise
		if (cw_rr_current == RECEIVE_CAPACITY - 1) {
			cw_receive_state = RS_IDLE;
			cw_rr_current = 0;	// reset decoding pointer
			smpl_ctr = 0;		// reset audio sample counter
			return CW_ERROR;
		} else
// zero terminate representation
			rx_rep_buf[cw_rr_current] = 0;
// All is well.  Move to the more normal after-tone state.
		cw_receive_state = RS_AFTER_TONE;
		return CW_ERROR;
		break;
	case CW_QUERY_EVENT:
// this should be called quite often (faster than inter-character gap) It looks after timing
// key up intervals and determining when a character, a word space, or an error char '*' should be returned.
// CW_SUCCESS is returned when there is a printable character. Nothing to do if we are in a tone
		if (cw_receive_state == RS_IN_TONE)
			return CW_ERROR;
// in this call we expect a pointer to a char to be valid
		if (c == NULL) {
// else we had no place to put character...
			cw_receive_state = RS_IDLE;
			cw_rr_current = 0;	
// reset decoding pointer
			return CW_ERROR;
		}
// compute length of silence so far
		sync_parameters();
		element_usec = usec_diff(cw_rr_end_timestamp, smpl_ctr);
		
// SHORT time since keyup... nothing to do yet
		if (element_usec < (2 * cw_receive_dot_length))
			return CW_ERROR;
// MEDIUM time since keyup... check for character space
// one shot through this code via receive state logic
		if (element_usec >= (2 * cw_receive_dot_length) &&
		    element_usec <= (4 * cw_receive_dot_length) &&
		    cw_receive_state == RS_AFTER_TONE) {
// Look up the representation 
			*c = morse::rx_lookup(rx_rep_buf);
			if (*c == NULL)
// invalid decode... let user see error
				*c = "*";
			cw_receive_state = RS_IDLE;
			cw_rr_current = 0;	// reset decoding pointer
			space_sent = false;
			return CW_SUCCESS;
		}
// LONG time since keyup... check for a word space
		if ((element_usec > (4 * cw_receive_dot_length)) && !space_sent) {
			*c = " "; 
			space_sent = true;
			return CW_SUCCESS;
		}
// should never get here... catch all
		return CW_ERROR;
		break;
	}
// should never get here... catch all
	return CW_ERROR;
}

//===========================================================================
// cw transmit routines
// Define the amplitude envelop for key down events (32 samples long)      
// this is 1/2 cycle of a raised cosine                                    
// the tables with 32 entries give about 4ms rise and fall times           
// when using 8000 samples/sec. This shaping of the cw pulses is           
// very necssary to avoid having a very wide and clicky cw signal          
// when using the sound card to gen cw. When using the rig key input       
// the shaping is done in the rig hardware, but we want to be able to      
// pick one cw signal out of a cluster and be able to respond on his freq. 
//===========================================================================

#define KNUM 32
// keydown wave shape
double kdshape[KNUM] = {
	0.00240750255310301, 0.00960708477768751,
	0.02152941088003600, 0.03805966253618680,
	0.05903864465505320, 0.08426431851158830, 
	0.11349374748686800, 0.14644543667658500,
	0.18280204383628200, 0.22221343555548300, 
	0.26430005922814900, 0.30865659834558700,
	0.35485587590940700, 0.40245296837259500, 
	0.45098949048925500, 0.49999800980765500,
	0.54900654829266300, 0.59754312772456200, 
	0.64514031509964400, 0.69133972425796200,
	0.73569643038517400, 0.77778325487450100, 
	0.81719487928327800, 0.85355174876454100,
	0.88650372738152000, 0.91573347010241700, 
	0.94095947900139100, 0.96193881423287900,
	0.97846943367117300, 0.99039213868324900, 
	0.99759210729604500, 0.99999999999295900
};

// keyup wave shape
double kushape[KNUM] = {
	0.99999999999295900, 0.99759210729604500, 
	0.99039213868324900, 0.97846943367117300,
	0.96193881423287900, 0.94095947900139100, 
	0.91573347010241700, 0.88650372738152000,
	0.85355174876454100, 0.81719487928327800, 
	0.77778325487450100, 0.73569643038517400,
	0.69133972425796200, 0.64514031509964400, 
	0.59754312772456200, 0.54900654829266300,
	0.49999800980765500, 0.45098949048925500, 
	0.40245296837259500, 0.35485587590940700,
	0.30865659834558700, 0.26430005922814900, 
	0.22221343555548300, 0.18280204383628200,
	0.14644543667658500, 0.11349374748686800, 
	0.08426431851158830, 0.05903864465505320,
	0.03805966253618680, 0.02152941088003600, 
	0.00960708477768751, 0.00240750255310301
};


inline double cw::nco(double freq)
{
	phaseacc += 2.0 * M_PI * freq / samplerate;

	if (phaseacc > M_PI)
		phaseacc -= 2.0 * M_PI;

	return cos(phaseacc);
}

//=====================================================================
// send_symbol()
// Sends a part of a morse character (one dot duration) of either
// sound at the correct freq or silence. Rise and fall time is controlled
// with a raised cosine shape.
//=======================================================================


void cw::send_symbol(int symbol)
{
	double freq;
	int sample = 0, i;
	int currsym = symbol & 1;
	int nextsym = (symbol >> 1) & 1;
//	int symlen100 = (int) (samplerate * 0.012);
	int delta = 0;
	int keydown;
	int keyup;
	int duration = 0;

	freq = tx_frequency;

	if ((currsym == 1 && lastsym == 0) || (currsym == 0 && lastsym == 1))
		delta = (int) (symbollen * (progdefaults.CWweight - 50) / 100.0);
	keydown = symbollen - 2 * KNUM + delta;
	keyup = symbollen - delta;

	if (currsym == 1) {
		if (cw_send_speed <= 100) {
			for (i = 0; i < KNUM; i++, sample++) {
				if (lastsym == 0)
					outbuf[sample] = nco(freq) * kdshape[i];
				else
					outbuf[sample] = nco(freq);
			}
			for (i = 0; i < keydown; i++, sample++) {
				outbuf[sample] = nco(freq);
			}
			for (i = 0; i < KNUM; i++, sample++ ) {
				if (nextsym == 0)
					outbuf[sample] = nco(freq) * kushape[i];
				else
					outbuf[sample] = nco(freq);
			}
		} else {
			for (i = 0; i < KNUM; i += 2, sample++) {
				if (lastsym == 0)
					outbuf[sample] = nco(freq) * kdshape[i];
				else
					outbuf[sample] = nco(freq);
			}
			duration += KNUM / 2;
			for (i = 0; i < keydown + KNUM; i++, sample++) {
				outbuf[sample] = nco(freq);
			}
			for (i = 0; i < KNUM; i += 2, sample++ ) {
				if (nextsym == 0)
					outbuf[sample] = nco(freq) * kushape[i];
				else
					outbuf[sample] = nco(freq);
			}
		}
		duration = keydown + 2 * KNUM;
	} else {
		if (lastsym == 1) {
			for (i = 0; i < keyup; i++)
				outbuf[i] = 0.0;
			duration = keyup;
		} else {
			for (i = 0; i < symbollen; i++)
				outbuf[i] = 0.0;
			duration = symbollen;
		}
	}

//	if (progdefaults.useCWkeylineRTS || progdefaults.useCWkeylineDTR) {
//		if (currsym != lastsym)
//			cw_keyline(currsym);
//		std::cout << currsym; fflush(stdout);
//		}

	ModulateXmtr(outbuf, duration);
	
	lastsym = currsym;
}

//=====================================================================
// send_ch()
// sends a morse character and the space afterwards
//=======================================================================

void cw::send_ch(int ch)
{
	int code;

	sync_parameters();
// handle word space separately (7 dots spacing) 
// last char already had 2 dots of inter-character spacing sent with it 
	if ((ch == ' ') || (ch == '\n')) {
		send_symbol(0);
		send_symbol(0);
		send_symbol(0);
		send_symbol(0);
		send_symbol(0);
		put_echo_char(ch);
		return;
	}

// convert character code to a morse representation 
	if ((ch < 256) && (ch >= 0))
		code = tx_lookup(ch); //cw_tx_lookup(ch);
	else
		code = 0x4; // two dot spaces
// loop sending out binary bits of cw character 
	while (code > 1) {
		send_symbol(code);// & 1);
		code = code >> 1;
	}
	if (ch != 0)
		put_echo_char(ch);
}

//=====================================================================
// cw_txprocess()
// Read charcters from screen and send them out the sound card.
// This is called repeatedly from a thread during tx.
//=======================================================================

int cw::tx_process()
{
	int c;
	c = get_tx_char();
	if (c == 0x03 || stopflag) {			
		send_symbol(0);
		stopflag = false;
			return -1;
	}
	if (c != 0)
		send_ch(c);
	else
		send_symbol(0);

	return 0;
}

/*
void cw::cw_keyup()
{
	if (progdefaults.useCWkeylineRTS)
		KeyLine->clearRTS();
	else if (progdefaults.useCWkeylineDTR)
		KeyLine->clearDTR();
}

void cw::cw_keydown()
{
	if (progdefaults.useCWkeylineRTS)
		KeyLine->setRTS();
	else if (progdefaults.useCWkeylineDTR)
		KeyLine->setDTR();
}

void cw::cw_keyline(int symbol)
{
	if (symbol)
		cw_keydown();
	else
		cw_keyup();
}
*/
