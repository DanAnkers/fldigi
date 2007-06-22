/* 
 * Waterfall Spectrum Analyzer Widget
 * Copyright (C) 2006 Dave Freese, W1HKJ
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Library General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Library General Public License for more details.
 *
 * You should have received a copy of the GNU Library General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 * Please report all bugs and problems to "w1hkj@w1hkj.com".
 */

#ifndef _WF_H
#define _WF_H

#include <math.h>
#include <stdio.h>
#include <string.h>

#include "config.h"

#include "complex.h"
#include "fft.h"
#include "sound.h"

#include <FL/Fl.H>
#include <FL/fl_draw.H>
#include <FL/Fl_Widget.H>
#include <FL/Fl_Repeat_Button.H>
#include <FL/Fl_Light_Button.H>
#include <FL/Fl_Window.H>
#include <FL/Fl_Group.H>
#include <FL/Fl_Box.H>
#include <FL/Fl_Counter.H>
#include <FL/Fl_Choice.H>
#include <FL/Enumerations.H>

/*
#ifdef HAVE_DFFTW_H
#  include <dfftw.h>
#endif
#ifdef HAVE_FFTW_H
#  include <fftw.h>
#endif
*/

// recommended minimum size for the control is width = 504, height = 104;
// the actual waterfall will be width -4 (bezel size) and
//                              height - 4 - 24 (bezel, text, scale & marker)

//#define FFT_LEN     2048
#define FFT_LEN		4096
#define SC_SMPLRATE	8000

#define BEZEL		 2
#define WFTEXT		10
#define WFSCALE     10
#define WFMARKER     6
#define BTN_HEIGHT	20

#define bwColor		35
#define bwFFT		35
#define bwX1		25
#define bgX14		75
#define bwMov		20
#define bgMov		60
#define cwCnt		120
#define cwRef		60
#define cwMode		85
#define bwQsy		45
#define bwRate		45
#define bwXmtLock	45
#define bwRev		45
#define bwXmtRcv	60
#define wSpace		4
 
#define fftabs(a,b) sqrt((a)*(a) + (b)*(b))

struct RGB {
	uchar R;
	uchar G;
	uchar B;
};

struct RGBint {
	int R;
	int G;
	int B;
};

struct RGBI {
	uchar R;
	uchar G;
	uchar B;
	uchar I;
};

extern 	RGBI	mag2RGBI[256];
extern	RGB		palette[9];

class WFdisp : public Fl_Widget {
public:
enum WFmode {
	WATERFALL,
	SPECTRUM,
	SCOPE
};

#define MAG_1 1
#define MAG_2 2
#define MAG_4 3

enum WFspeed {FAST = 1, NORMAL = 2, SLOW = 4};

	WFdisp (int x, int y, int w, int h, char *lbl = 0);
	~WFdisp ();
	int wfmag();
	void Mode(WFmode M) {
		mode = M;
	}
	WFmode Mode() {
		return mode;
	}
	int cursorFreq(int xpos) {
		return (offset + step * xpos);
	}
	void DispColor(bool Y) {
		dispcolor = Y;
	}
	bool DispColor() {
		return dispcolor;
	}
	void Ampspan(double AmpSpn) {
		ampspan = (int)AmpSpn;
		update_fft_db();
	}
	double Ampspan() {
		return ampspan;
	}
	void Reflevel(double RefLev) {
		reflevel = (int)RefLev;
		update_fft_db();
	}
	double Reflevel() {
		return reflevel;
	}
	void Bandwidth (int bw) {
		bandwidth = bw;
		makeMarker();
	}
	int  Bandwidth () {
		return bandwidth;
	}
	void Overload(int ovr) {
		if (overload == ovr) return;
		overload = ovr;
	}
	WFspeed Speed() { return wfspeed;}
	void Speed(WFspeed rate) { wfspeed = rate;}
	
	void initmaps();
	void draw();
//	void resize (int, int, int, int);
	void update_sigmap();
	void update_waterfall();
	void checkoffset();
	void slew(int);
	void movetocenter();
	void carrier(int cf);
	int  carrier();
	void makeMarker();
	void process_analog(double *sig, int len);
	void processFFT();
	void sig_data( double *sig, int len );
	void rfcarrier(long long f) { 
		rfc = f;
	}
	void USB(bool b) { 
		usb = b;
	}
	bool USB() {return usb;};
	long long rfcarrier() { return rfc;};
	
	void useBands(bool b) { usebands = b;};
	
	void updateMarker() { 
		drawMarker();};
	int peakFreq(int f0, int delta);
	double powerDensity(double f0, double bw);
	void setPrefilter(int v) {
		switch (v) {
			case 0: wfft->setWindow(FFT_NONE); break;
			case 1: wfft->setWindow(FFT_BLACKMAN); break;
			case 2: wfft->setWindow(FFT_HAMMING); break;
			case 3: wfft->setWindow(FFT_HANNING); break;
			case 4: wfft->setWindow(FFT_TRIANGULAR); break;
		}
	}
	void setcolors();
	double dFreq() {return dfreq;}
	void redrawCursor();
	void defaultColors();
	
private:
	int disp_width;
	int image_width;
	int scale_width;
	int RGBwidth;
	int RGBsize;
	int image_height;
	int image_area;
	int sig_image_area;
	int	mag;
	int magset;
	WFmode	mode;
	bool	overload;
	bool	usb;
	long long	rfc;
	bool	usebands;
	int		offset;
	int		sigoffset;
	int		step;
	int		carrierfreq;
	int		bandwidth;
	int		wfspdcnt;
	int		dispcnt;
	int 	ampspan;
	int 	reflevel;
	double	dfreq;
	bool	centercarrier;
	bool	dispcolor;
	bool	cursormoved;
	WFspeed	wfspeed;
	int		srate;
	RGBI	*fft_img;
//	RGBI	mag2RGBI[256];
	RGB		*markerimage;
	RGB		RGBmarker;
	RGB		RGBcursor;
	double	*fftout;
	uchar	*scaleimage;
	uchar	*fft_sig_img;
	uchar	*sig_img;
	uchar	*scline;
	
	short int	*fft_hist;
	short int	*fft_db;
	double	 	*circbuff;
	double	*pwr;
	Cfft	*wfft;

	int checkMag();
	void checkWidth();
	void initMarkers();
	void makeScale();
	void drawScale();
	void drawMarker();

	int	 log2disp(int v);
	void update_fft_db();
	void drawcolorWF();
	void drawgrayWF();
	void drawspectrum();
	void drawsignal();
protected:
public:
	bool	wantcursor;
	int		cursorpos;
};

class waterfall: public Fl_Group {
	friend void x1_cb(Fl_Widget *w, void* v);
	friend void bwclr_cb(Fl_Widget *w, void * v);
//	friend void slew_cb(Fl_Widget *w, void * v);
	friend void slew_left(Fl_Widget *w, void * v);
	friend void slew_right(Fl_Widget *w, void * v);
	friend void center_cb(Fl_Widget *w, void *v);
	friend void carrier_cb(Fl_Widget *w, void *v);
	friend void mode_cb(Fl_Widget *w, void *v);
	friend void reflevel_cb(Fl_Widget *w, void *v);
	friend void ampspan_cb(Fl_Widget *w, void *v);
	friend void qsy_cb(Fl_Widget *w, void *v);
	friend void rate_cb(Fl_Widget *w, void *v);
public:
	waterfall(int x, int y, int w, int h, char *lbl= 0);
	~waterfall(){};
	void opmode();
	void sig_data(double *sig, int len){
		wfdisp->sig_data(sig, len);
	}
	void Overload(bool ovr) { 
		wfdisp->Overload(ovr);
	}
	int carrier() {
		return wfdisp->carrier();
	}
	void carrier(int f);
	void rfcarrier(long long cf);
	long long rfcarrier();
	void set_XmtRcvBtn(bool val);
	void USB(bool b);
	bool USB();
	void Reverse( bool v) { reverse = v;}
	bool Reverse() { return reverse;}
	void Bandwidth(int bw)
	{
		wfdisp->Bandwidth(bw);
	}
	int peakFreq(int f0, int delta)
	{
		return (wfdisp->peakFreq(f0, delta));
	}
	double powerDensity(double f0, double bw)
	{
		return (wfdisp->powerDensity(f0,bw));
	}

	void movetocenter() { wfdisp->movetocenter();}
	
	void setPrefilter(int v) {wfdisp->setPrefilter(v);}
	
	void setcolors() { wfdisp->setcolors(); }
	void setRefLevel();
	void setAmpSpan();
	double dFreq() { return wfdisp->dFreq();}
	
	void setQSY(bool on) {
		if (on)
			qsy->activate();
		else
			qsy->deactivate();
		wfdisp->useBands(!on);
	}
	
	
	int handle(int event);
/*
*/
	Fl_Button	*btnRev;
private:
	Fl_Box		*bezel;
	WFdisp		*wfdisp;
	Fl_Button	*bwclr;
	Fl_Button	*mode;
	Fl_Button	*x1;
	Fl_Button	*left; 
	Fl_Button	*center;
	Fl_Button	*right;
	Fl_Counter	*wfcarrier;
	Fl_Counter	*wfRefLevel;
	Fl_Counter	*wfAmpSpan;
	Fl_Button	*qsy;
	Fl_Button	*wfrate;
	Fl_Light_Button	*xmtrcv;
	Fl_Light_Button *xmtlock;
	int			buttonrow;
	bool	reverse;
protected:
};

#endif