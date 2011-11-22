// ----------------------------------------------------------------------------
// Digital Modem Program for the Fast Light Toolkit
//
// Copyright 2006-2010, Dave Freese, W1HKJ
// Copyright 2007-2010, Stelios Bounanos, M0GLD
//
// This file is part of fldigi.
//
// Fldigi is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Fldigi is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with fldigi.  If not, see <http://www.gnu.org/licenses/>.
//
// Please report all bugs and problems to fldigi-devel@lists.berlios.de.
// ----------------------------------------------------------------------------

#include <config.h>

#include <iostream>
#include <iomanip>
#include <sstream>
#include <cstdlib>
#include <getopt.h>
#include <sys/types.h>

#if !defined(__WOE32__) && !defined(__APPLE__)
#  include <sys/ipc.h>
#  include <sys/msg.h>
#endif

#ifdef __MINGW32__
#  include "compat.h"
#endif

#include <sys/stat.h>

#if HAVE_SYS_UTSNAME_H
#  include <sys/utsname.h>
#endif

#include <unistd.h>
#ifndef __MINGW32__
#  include <dirent.h>
#endif
#include <exception>
#include <signal.h>
#include <locale.h>

#include <FL/Fl.H>
#include <FL/Enumerations.H>
#include <FL/Fl_Window.H>
#include <FL/Fl_Shared_Image.H>
#include <FL/x.H>
#ifdef __MINGW32__
#  define dirent fl_dirent_no_thanks
#endif
#include <FL/filename.H>
#ifdef __MINGW32__
#  undef dirent
#endif

#include "gettext.h"
#include "main.h"
#include "waterfall.h"
#include "trx.h"
#include "soundconf.h"
#include "fl_digi.h"
#include "rigio.h"
#include "globals.h"
#include "confdialog.h"
#include "colorsfonts.h"
#include "configuration.h"
#include "macros.h"
#include "status.h"
#include "fileselect.h"
#include "timeops.h"
#include "debug.h"
#include "pskrep.h"
#include "notify.h"
#include "logbook.h"
#include "dxcc.h"
#include "newinstall.h"
#include "Viewer.h"

#if USE_HAMLIB
	#include "rigclass.h"
#endif
#include "rigsupport.h"

#include "log.h"

#include "qrunner.h"
#include "stacktrace.h"

#if USE_XMLRPC
	#include "xmlrpc.h"
#endif

#if BENCHMARK_MODE
	#include "benchmark.h"
#endif

#include "icons.h"


using namespace std;

string appname;

string scDevice[2];
string vscDevice[2];

string HomeDir = "";
string RigsDir = "";
string ScriptsDir = "";
string PalettesDir = "";
string LogsDir = "";
string PicsDir = "";
string HelpDir = "";
string MacrosDir = "";
string WrapDir = "";
string TalkDir = "";
string TempDir = "";
string PskMailDir = "";

string NBEMS_dir = "";
string ARQ_dir = "";
string ARQ_files_dir = "";
string ARQ_recv_dir = "";
string ARQ_send = "";
string WRAP_dir = "";
string WRAP_recv_dir = "";
string WRAP_send_dir = "";
string WRAP_auto_dir = "";
string ICS_dir = "";
string ICS_msg_dir = "";
string ICS_tmp_dir = "";

string FLMSG_dir = "";
string FLMSG_dir_default = "";
string FLMSG_WRAP_dir = "";
string FLMSG_WRAP_recv_dir = "";
string FLMSG_WRAP_send_dir = "";
string FLMSG_WRAP_auto_dir = "";
string FLMSG_ICS_dir = "";
string FLMSG_ICS_msg_dir = "";
string FLMSG_ICS_tmp_dir = "";

string PskMailFile;
string ArqFilename;
string xmlfname;

PTT		*push2talk = (PTT *)0;
#if USE_HAMLIB
Rig		*xcvr = (Rig *)0;
#endif

bool tlfio = false;
cLogfile	*logfile = 0;
cLogfile	*Maillogfile = (cLogfile *)0;
FILE	*server;
FILE	*client;
bool	mailserver = false, mailclient = false, arqmode = false;
static bool show_cpucheck = false;
static bool iconified = false;

RXMSGSTRUC rxmsgst;
int		rxmsgid = -1;

TXMSGSTRUC txmsgst;
int txmsgid = -1;

string option_help, version_text, build_text;

qrunner *cbq[NUM_QRUNNER_THREADS];

void arqchecks(void);
void generate_option_help(void);
int parse_args(int argc, char **argv, int& idx);
void generate_version_text(void);
void debug_exec(char** argv);
void set_platform_ui(void);
double speed_test(int converter, unsigned repeat);
static void setup_signal_handlers(void);
static void checkdirectories(void);

static void arg_error(const char* name, const char* arg, bool missing) noreturn__;

// TODO: find out why fldigi crashes on OS X if the wizard window is
// shown before fldigi_main.
#ifndef __APPLE__
#  define SHOW_WIZARD_BEFORE_MAIN_WINDOW 1
#else
#  define SHOW_WIZARD_BEFORE_MAIN_WINDOW 0
#endif

int main(int argc, char ** argv)
{
	appname = argv[0];
	debug_exec(argv);
	CREATE_THREAD_ID(); // only call this once
	SET_THREAD_ID(FLMAIN_TID);

	for (int i = 0; i < NUM_QRUNNER_THREADS; i++) {
		cbq[i] = new qrunner;
		cbq[i]->attach();
	}

	set_unexpected(handle_unexpected);
	set_terminate(diediedie);
	setup_signal_handlers();

#ifndef ENABLE_NLS
	setlocale(LC_TIME, "");
#endif

	{
		char dirbuf[FL_PATH_MAX + 1];
#ifdef __WOE32__
		fl_filename_expand(dirbuf, sizeof(dirbuf) - 1, "$USERPROFILE/fldigi.files/");
		HomeDir = dirbuf;
		fl_filename_expand(dirbuf, sizeof(dirbuf) - 1, "$USERPROFILE/NBEMS.files/");
		NBEMS_dir = dirbuf;
		fl_filename_expand(dirbuf, sizeof(dirbuf) - 1, "$USERPROFILE/");
		PskMailDir = dirbuf;
		FLMSG_dir_default = "$USERPROFILE/NBEMS.files/";
#else
		fl_filename_expand(dirbuf, sizeof(dirbuf) - 1, "$HOME/.fldigi/");
		HomeDir = dirbuf;
		fl_filename_expand(dirbuf, sizeof(dirbuf) - 1, "$HOME/.nbems/");
		NBEMS_dir = dirbuf;
		fl_filename_expand(dirbuf, sizeof(dirbuf) - 1, "$HOME/");
		PskMailDir = dirbuf;
		FLMSG_dir_default = "$HOME/.nbems/";
#endif
	}

	set_platform_ui();

	generate_option_help();
	generate_version_text();
	int arg_idx;
	if (Fl::args(argc, argv, arg_idx, parse_args) != argc)
		arg_error(argv[0], NULL, false);
	if (main_window_title.empty())
		main_window_title = PACKAGE_TARNAME;

	{
		char dirbuf[FL_PATH_MAX + 1];
		if (FLMSG_dir_default[FLMSG_dir_default.length()-1] != '/')
			FLMSG_dir_default += '/';
		fl_filename_expand(dirbuf, sizeof(dirbuf) - 1, FLMSG_dir_default.c_str());
		FLMSG_dir = dirbuf;
	}
	checkdirectories();
	check_nbems_dirs();

	try {
		debug::start(string(HomeDir).append("status_log.txt").c_str());
		time_t t = time(NULL);
		LOG(debug::QUIET_LEVEL, debug::LOG_OTHER, _("%s log started on %s"), PACKAGE_STRING, ctime(&t));
		LOG_THREAD_ID();
	}
	catch (const char* error) {
		cerr << error << '\n';
		debug::stop();
	}

	bool have_config = progdefaults.readDefaultsXML();

	xmlfname = HomeDir;
	xmlfname.append("rig.xml");

#if !defined(__WOE32__) && !defined(__APPLE__)
   	txmsgid = msgget( (key_t) progdefaults.tx_msgid, 0666 );
#else
	txmsgid = -1;
#endif

	checkTLF();


	Fl::lock();  // start the gui thread!!
	Fl::visual(FL_RGB); // insure 24 bit color operation

	fl_register_images();
	Fl::set_fonts(0);

	Fl::scheme(progdefaults.ui_scheme.c_str());
	progdefaults.initFonts();

	dxcc_open(string(HomeDir).append("cty.dat").c_str());
	qsl_open(string(HomeDir).append("lotw1.txt").c_str(), QSL_LOTW);
	if (!qsl_open(string(HomeDir).append("eqsl.txt").c_str(), QSL_EQSL))
		qsl_open(string(HomeDir).append("AGMemberList.txt").c_str(), QSL_EQSL);

	progStatus.loadLastState();
	create_fl_digi_main(argc, argv);

	if (!have_config || show_cpucheck) {
		double speed = speed_test(SRC_SINC_FASTEST, 8);

		if (speed > 150.0) {      // fast
			progdefaults.slowcpu = false;
			progdefaults.sample_converter = SRC_SINC_BEST_QUALITY;
		}
		else if (speed > 60.0) {  // ok
			progdefaults.slowcpu = false;
			progdefaults.sample_converter = SRC_SINC_MEDIUM_QUALITY;
		}
		else if (speed > 15.0) { // slow
			progdefaults.slowcpu = true;
			progdefaults.sample_converter = SRC_SINC_FASTEST;
		}
		else {                   // recycle me
			progdefaults.slowcpu = true;
			progdefaults.sample_converter = SRC_LINEAR;
		}

		LOG_INFO("CPU speed factor=%f: setting slowcpu=%s, sample_converter=\"%s\"", speed,
			 progdefaults.slowcpu ? "true" : "false",
			 src_get_name(progdefaults.sample_converter));
	}
	if (progdefaults.XmlRigFilename.empty())
		progdefaults.XmlRigFilename = xmlfname;

	if (progStatus.LOGenabled == true) {
    	Date tdy;
	    string lfname = HomeDir;
	    lfname.append("fldigi");
	    lfname.append(tdy.szDate(2));
	    lfname.append(".log");
	    logfile = new cLogfile(lfname);
	    logfile->log_to_file_start();
	}

#if BENCHMARK_MODE
	return setup_benchmark();
#endif

	FSEL::create();

	make_colorsfonts();
	setTabColors();

//	start_logbook();

	progdefaults.testCommPorts();

	macros.loadDefault();

#if USE_HAMLIB
	xcvr = new Rig();
#endif

	push2talk = new PTT();

	progdefaults.setDefaults();

	atexit(sound_close);
	sound_init();

	progdefaults.initInterface();
	trx_start();

#if SHOW_WIZARD_BEFORE_MAIN_WINDOW
	if (!have_config) {
		show_wizard(argc, argv);
		Fl_Window* w;
		while ((w = Fl::first_window()) && w->visible())
			Fl::wait();
	}
#endif

	dlgViewer = createViewer();
	create_logbook_dialogs();
	connect_to_log_server();


// OS X will prevent the main window from being resized if we change its
// size *after* it has been shown. With some X11 window managers, OTOH,
// the main window will not be restored at its exact saved position if
// we move it *after* it has been shown.
#ifndef __APPLE__
	fl_digi_main->show(argc, argv);
	progStatus.initLastState();
#else
	progStatus.initLastState();
	fl_digi_main->show(argc, argv);
#endif

	if (iconified)
		for (Fl_Window* w = Fl::first_window(); w; w = Fl::next_window(w))
			w->iconize();
	update_main_title();

	arq_init();

#ifdef __WIN32__
	if (progdefaults.auto_talk) open_talker();
#else
	grpTalker->hide();
#endif

#if USE_XMLRPC
	XML_RPC_Server::start(progdefaults.xmlrpc_address.c_str(), progdefaults.xmlrpc_port.c_str());
#endif

	if (progdefaults.usepskrep)
		if (!pskrep_start())
			LOG_ERROR("Could not start PSK reporter: %s", pskrep_error());

	notify_start();
	mode_browser = new Mode_Browser;

#if !SHOW_WIZARD_BEFORE_MAIN_WINDOW
	if (!have_config)
		show_wizard();
#endif
//	connect_to_log_server();

	int ret = Fl::run();

	arq_close();

#if USE_XMLRPC
	XML_RPC_Server::stop();
#endif

	if (progdefaults.usepskrep)
		pskrep_stop();

	for (int i = 0; i < NUM_QRUNNER_THREADS; i++) {
		cbq[i]->detach();
		delete cbq[i];
	}
	FSEL::destroy();
	debug::stop();

	return ret;
}


void generate_option_help(void) {
	ostringstream help;
	help << "Usage:\n"
	     << "    " << PACKAGE_NAME << " [option...]\n\n";

	help << PACKAGE_NAME << " options:\n\n"
	     << "  --config-dir DIRECTORY\n"
	     << "    Look for configuration files in DIRECTORY\n"
	     << "    The default is: " << HomeDir << "\n\n"

#if !defined(__WOE32__) && !defined(__APPLE__)
	     << "  --rx-ipc-key KEY\n"
	     << "    Set the receive message queue key\n"
	     << "    May be given in hex if prefixed with \"0x\"\n"
	     << "    The default is: " << progdefaults.rx_msgid
	     << " or 0x" << hex << progdefaults.rx_msgid << dec << "\n\n"

	     << "  --tx-ipc-key KEY\n"
	     << "    Set the transmit message queue key\n"
	     << "    May be given in hex if prefixed with \"0x\"\n"
	     << "    The default is: " << progdefaults.tx_msgid
	     << " or 0x" << hex << progdefaults.tx_msgid << dec << "\n\n"
#endif

	     << "  --arq-server-address HOSTNAME\n"
	     << "    Set the ARQ TCP server address\n"
	     << "    The default is: " << progdefaults.arq_address << "\n\n"
	     << "  --arq-server-port PORT\n"
	     << "    Set the ARQ TCP server port\n"
	     << "    The default is: " << progdefaults.arq_port << "\n\n"
	     << "  --flmsg-dir DIRECTORY\n"
	     << "    Look for flmsg files in DIRECTORY\n"
	     << "    The default is " << FLMSG_dir_default << "\n\n"
	     << "  --auto-dir DIRECTORY\n"
	     << "    Look for auto-send files in DIRECTORY\n"
	     << "    The default is " << HomeDir << "/autosend" << "\n\n"

#if USE_XMLRPC
	     << "  --xmlrpc-server-address HOSTNAME\n"
	     << "    Set the XML-RPC server address\n"
	     << "    The default is: " << progdefaults.xmlrpc_address << "\n\n"
	     << "  --xmlrpc-server-port PORT\n"
	     << "    Set the XML-RPC server port\n"
	     << "    The default is: " << progdefaults.xmlrpc_port << "\n\n"
	     << "  --xmlrpc-allow REGEX\n"
	     << "    Allow only the methods whose names match REGEX\n\n"
	     << "  --xmlrpc-deny REGEX\n"
	     << "    Allow only the methods whose names don't match REGEX\n\n"
	     << "  --xmlrpc-list\n"
	     << "    List all available methods\n\n"
#endif

#if BENCHMARK_MODE
	     << "  --benchmark-modem ID\n"
	     << "    Specify the modem\n"
	     << "    Default: " << mode_info[benchmark.modem].sname << "\n\n"
	     << "  --benchmark-frequency FREQ\n"
	     << "    Specify the modem frequency\n"
	     << "    Default: " << benchmark.freq << "\n\n"
	     << "  --benchmark-afc BOOLEAN\n"
	     << "    Set modem AFC\n"
	     << "    Default: " << benchmark.afc
	     << " (" << boolalpha << benchmark.afc << noboolalpha << ")\n\n"
	     << "  --benchmark-squelch BOOLEAN\n"
	     << "    Set modem squelch\n"
	     << "    Default: " << benchmark.sql
	     << " (" << boolalpha << benchmark.sql << noboolalpha << ")\n\n"
	     << "  --benchmark-squelch-level LEVEL\n"
	     << "    Set modem squelch level\n"
	     << "    Default: " << benchmark.sqlevel << " (%)\n\n"
	     << "  --benchmark-input INPUT\n"
	     << "    Specify the input\n"
	     << "    Must be a positive integer indicating the number of samples\n"
		"    of silence to generate as the input"
#  if USE_SNDFILE
		", or a filename containing\n"
		"    non-digit characters"
#endif
		"\n\n"

	     << "  --benchmark-output FILE\n"
	     << "    Specify the output data file\n"
	     << "    Default: decoder output is discarded\n\n"
	     << "  --benchmark-src-ratio RATIO\n"
	     << "    Specify the sample rate conversion ratio\n"
	     << "    Default: 1.0 (input is not resampled)\n\n"
	     << "  --benchmark-src-type TYPE\n"
	     << "    Specify the sample rate conversion type\n"
	     << "    Default: " << benchmark.src_type << " (" << src_get_name(benchmark.src_type) << ")\n\n"
#endif

	     << "  --cpu-speed-test\n"
	     << "    Perform the CPU speed test, show results in the event log\n"
	     << "    and possibly change options.\n\n"

	     << "  --noise\n"
	     << "    Unhide controls for noise tests\n\n"

	     << "  --wfall-only\n"
	     << "    Hide all controls but the waterfall\n\n"

	     << "  --debug-level LEVEL\n"
	     << "    Set the event log verbosity\n\n"

	     << "  --version\n"
	     << "    Print version information\n\n"

	     << "  --build-info\n"
	     << "    Print build information\n\n"

	     << "  --help\n"
	     << "    Print this option help\n\n";

// Fl::help looks ugly so we'll write our own

	help << "Standard FLTK options:\n\n"

	     << "   -bg COLOR, -background COLOR\n"
	     << "    Set the background color\n"

	     << "   -bg2 COLOR, -background2 COLOR\n"
	     << "    Set the secondary (text) background color\n\n"

	     << "   -di DISPLAY, -display DISPLAY\n"
	     << "    Set the X display to use DISPLAY,\n"
	     << "    format is ``host:n.n''\n\n"

	     << "   -dn, -dnd or -nodn, -nodnd\n"
	     << "    Enable or disable drag and drop copy and paste in text fields\n\n"

	     << "   -fg COLOR, -foreground COLOR\n"
	     << "    Set the foreground color\n\n"

	     << "   -g GEOMETRY, -geometry GEOMETRY\n"
	     << "    Set the initial window size and position\n"
	     << "    GEOMETRY format is ``WxH+X+Y''\n"
	     << "    ** " << PACKAGE_NAME << " may override this setting **\n\n"

	     << "   -i, -iconic\n"
	     << "    Start " << PACKAGE_NAME << " in iconified state\n\n"

	     << "   -k, -kbd or -nok, -nokbd\n"
	     << "    Enable or disable visible keyboard focus in non-text widgets\n\n"

	     << "   -na CLASSNAME, -name CLASSNAME\n"
	     << "    Set the window class to CLASSNAME\n\n"

	     << "   -ti WINDOWTITLE, -title WINDOWTITLE\n"
	     << "    Set the window title\n\n";

	help << "Additional UI options:\n\n"

	     << "  --font FONT[:SIZE]\n"
	     << "    Set the widget font and (optionally) size\n"
	     << "    The default is: " << Fl::get_font(FL_HELVETICA)
	     << ':' << FL_NORMAL_SIZE << "\n\n"

		;

	option_help = help.str();
}

void exit_cb(void*) { fl_digi_main->do_callback(); }

int parse_args(int argc, char **argv, int& idx)
{
	// Only handle long options
	if (!(strlen(argv[idx]) >= 2 && strncmp(argv[idx], "--", 2) == 0)) {
		// Store the window title. We may need this early in the initialisation
		// process, before FLTK uses it to set the main window title.
		if (main_window_title.empty() && argc > idx &&
		    (!strcmp(argv[idx], "-ti") || !strcmp(argv[idx], "-title")))
			main_window_title = argv[idx + 1];
		else if (!strcmp(argv[idx], "-i") || !strcmp(argv[idx], "-iconic"))
			iconified = true;
		return 0;
	}

        enum { OPT_ZERO,
#ifndef __WOE32__
	       OPT_RX_IPC_KEY, OPT_TX_IPC_KEY,
#endif
	       OPT_CONFIG_DIR,
	       OPT_ARQ_ADDRESS, OPT_ARQ_PORT,
	       OPT_SHOW_CPU_CHECK,
	       OPT_FLMSG_DIR,
	       OPT_AUTOSEND_DIR,

#if USE_XMLRPC
	       OPT_CONFIG_XMLRPC_ADDRESS, OPT_CONFIG_XMLRPC_PORT,
	       OPT_CONFIG_XMLRPC_ALLOW, OPT_CONFIG_XMLRPC_DENY, OPT_CONFIG_XMLRPC_LIST,
#endif

#if BENCHMARK_MODE
	       OPT_BENCHMARK_MODEM, OPT_BENCHMARK_AFC, OPT_BENCHMARK_SQL, OPT_BENCHMARK_SQLEVEL,
	       OPT_BENCHMARK_FREQ, OPT_BENCHMARK_INPUT, OPT_BENCHMARK_OUTPUT,
	       OPT_BENCHMARK_SRC_RATIO, OPT_BENCHMARK_SRC_TYPE,
#endif

               OPT_FONT, OPT_WFALL_HEIGHT,
               OPT_WINDOW_WIDTH, OPT_WINDOW_HEIGHT, OPT_WFALL_ONLY,
#if USE_PORTAUDIO
               OPT_FRAMES_PER_BUFFER,
#endif
	       OPT_NOISE, OPT_DEBUG_LEVEL,
               OPT_EXIT_AFTER,
               OPT_DEPRECATED, OPT_HELP, OPT_VERSION, OPT_BUILD_INFO };

	const char shortopts[] = ":";
	static struct option longopts[] = {
#ifndef __WOE32__
		{ "rx-ipc-key",	   1, 0, OPT_RX_IPC_KEY },
		{ "tx-ipc-key",	   1, 0, OPT_TX_IPC_KEY },
#endif
		{ "config-dir",	   1, 0, OPT_CONFIG_DIR },

		{ "arq-server-address", 1, 0, OPT_ARQ_ADDRESS },
		{ "arq-server-port",    1, 0, OPT_ARQ_PORT },
		{ "flmsg-dir", 1, 0, OPT_FLMSG_DIR },
		{ "auto-dir", 1, 0, OPT_AUTOSEND_DIR },

		{ "cpu-speed-test", 0, 0, OPT_SHOW_CPU_CHECK },

#if USE_XMLRPC
		{ "xmlrpc-server-address", 1, 0, OPT_CONFIG_XMLRPC_ADDRESS },
		{ "xmlrpc-server-port",    1, 0, OPT_CONFIG_XMLRPC_PORT },
		{ "xmlrpc-allow",          1, 0, OPT_CONFIG_XMLRPC_ALLOW },
		{ "xmlrpc-deny",           1, 0, OPT_CONFIG_XMLRPC_DENY },
		{ "xmlrpc-list",           0, 0, OPT_CONFIG_XMLRPC_LIST },
#endif

#if BENCHMARK_MODE
		{ "benchmark-modem", 1, 0, OPT_BENCHMARK_MODEM },
		{ "benchmark-frequency", 1, 0, OPT_BENCHMARK_FREQ },
		{ "benchmark-afc", 1, 0, OPT_BENCHMARK_AFC },
		{ "benchmark-squelch", 1, 0, OPT_BENCHMARK_SQL },
		{ "benchmark-squelch-level", 1, 0, OPT_BENCHMARK_SQLEVEL },
		{ "benchmark-input", 1, 0, OPT_BENCHMARK_INPUT },
		{ "benchmark-output", 1, 0, OPT_BENCHMARK_OUTPUT },
		{ "benchmark-src-ratio", 1, 0, OPT_BENCHMARK_SRC_RATIO },
		{ "benchmark-src-type", 1, 0, OPT_BENCHMARK_SRC_TYPE },
#endif

		{ "font",	   1, 0, OPT_FONT },

		{ "wfall-height",  1, 0, OPT_WFALL_HEIGHT },
		{ "window-width",  1, 0, OPT_WINDOW_WIDTH },
		{ "window-height", 1, 0, OPT_WINDOW_HEIGHT },
		{ "wfall-only",    0, 0, OPT_WFALL_ONLY },
		{ "wo",            0, 0, OPT_WFALL_ONLY },

#if USE_PORTAUDIO
		{ "frames-per-buffer",1, 0, OPT_FRAMES_PER_BUFFER },
#endif
		{ "exit-after",    1, 0, OPT_EXIT_AFTER },

		{ "noise", 0, 0, OPT_NOISE },
		{ "debug-level",   1, 0, OPT_DEBUG_LEVEL },

		{ "help",	   0, 0, OPT_HELP },
		{ "version",	   0, 0, OPT_VERSION },
		{ "build-info",	   0, 0, OPT_BUILD_INFO },
		{ 0 }
	};

	int longindex;
	optind = idx;
	opterr = 0;
	int c = getopt_long(argc, argv, shortopts, longopts, &longindex);

	switch (c) {
		case -1:
			return 0;
		case 0:
			// handle options with non-0 flag here
			return 0;

#if !defined(__WOE32__) && !defined(__APPLE__)
		case OPT_RX_IPC_KEY: case OPT_TX_IPC_KEY:
		{
			errno = 0;
			int key = strtol(optarg, NULL, (strncasecmp(optarg, "0x", 2) ? 10 : 16));
			if (errno || key <= 0)
				cerr << "Hmm, " << key << " doesn't look like a valid IPC key\n";
			if (c == OPT_RX_IPC_KEY)
				progdefaults.rx_msgid = key;
			else
				progdefaults.tx_msgid = key;
		}
			break;
#endif

		case OPT_CONFIG_DIR: {
			char buf[FL_PATH_MAX + 1];
			fl_filename_absolute(buf, sizeof(buf) - 1, optarg);
			HomeDir = buf;
		}
			if (*HomeDir.rbegin() != '/')
			       HomeDir += '/';
			break;

		case OPT_ARQ_ADDRESS:
			progdefaults.arq_address = optarg;
			break;
		case OPT_ARQ_PORT:
			progdefaults.arq_port = optarg;
			break;

		case OPT_FLMSG_DIR:
			FLMSG_dir_default = optarg;
			break;

		case OPT_AUTOSEND_DIR:
			FLMSG_WRAP_auto_dir = optarg;
			break;

#if USE_XMLRPC
		case OPT_CONFIG_XMLRPC_ADDRESS:
			progdefaults.xmlrpc_address = optarg;
			break;
		case OPT_CONFIG_XMLRPC_PORT:
			progdefaults.xmlrpc_port = optarg;
			break;
		case OPT_CONFIG_XMLRPC_ALLOW:
			progdefaults.xmlrpc_allow = optarg;
			break;
		case OPT_CONFIG_XMLRPC_DENY:
			if (!progdefaults.xmlrpc_allow.empty())
				cerr << "W: --" << longopts[longindex].name
				     << " cannot be used together with --"
				     << longopts[OPT_CONFIG_XMLRPC_ALLOW-1].name
				     << " and will be ignored\n";
			else
				progdefaults.xmlrpc_deny = optarg;
			break;
		case OPT_CONFIG_XMLRPC_LIST:
			XML_RPC_Server::list_methods(cout);
			exit(EXIT_SUCCESS);
#endif

#if BENCHMARK_MODE
		case OPT_BENCHMARK_MODEM:
			benchmark.modem = strtol(optarg, NULL, 10);
			if (!(benchmark.modem >= 0 && benchmark.modem < NUM_MODES)) {
				cerr << "Bad modem id\n";
				exit(EXIT_FAILURE);
			}
			break;

		case OPT_BENCHMARK_FREQ:
			benchmark.freq = strtol(optarg, NULL, 10);
			if (benchmark.freq < 0) {
				cerr << "Bad frequency\n";
				exit(EXIT_FAILURE);
			}
			break;

		case OPT_BENCHMARK_AFC:
			benchmark.afc = strtol(optarg, NULL, 10);
			break;

		case OPT_BENCHMARK_SQL:
			benchmark.sql = strtol(optarg, NULL, 10);
			break;

		case OPT_BENCHMARK_SQLEVEL:
			benchmark.sqlevel = strtod(optarg, NULL);
			break;

		case OPT_BENCHMARK_INPUT:
			benchmark.input = optarg;
			break;

		case OPT_BENCHMARK_OUTPUT:
			benchmark.output = optarg;
			break;

		case OPT_BENCHMARK_SRC_RATIO:
			benchmark.src_ratio = strtod(optarg, NULL);
			break;

		case OPT_BENCHMARK_SRC_TYPE:
			benchmark.src_type = strtol(optarg, NULL, 10);
			break;
#endif

		case OPT_FONT:
		{
			char *p;
			if ((p = strchr(optarg, ':'))) {
				*p = '\0';
				FL_NORMAL_SIZE = strtol(p + 1, 0, 10);
			}
		}
			Fl::set_font(FL_HELVETICA, optarg);
			break;

		case OPT_WFALL_HEIGHT:
			progdefaults.wfheight = strtol(optarg, NULL, 10);
			break;

		case OPT_WINDOW_WIDTH:
			WNOM = strtol(optarg, NULL, 10);
			break;

		case OPT_WINDOW_HEIGHT:
			HNOM = strtol(optarg, NULL, 10);
			break;

#if USE_PORTAUDIO
		case OPT_FRAMES_PER_BUFFER:
			progdefaults.PortFramesPerBuffer = strtol(optarg, 0, 10);
			break;
#endif // USE_PORTAUDIO

		case OPT_EXIT_AFTER:
			Fl::add_timeout(strtod(optarg, 0), exit_cb);
			break;

		case OPT_WFALL_ONLY:
			bWF_only = true;
			break;

		case OPT_NOISE:
			withnoise = true;
			break;

		case OPT_SHOW_CPU_CHECK:
			show_cpucheck = true;
			break;

		case OPT_DEBUG_LEVEL:
		{
			int v = strtol(optarg, 0, 10);
			debug::level = (debug::level_e)CLAMP(v, 0, debug::LOG_NLEVELS-1);
		}
			break;

		case OPT_DEPRECATED:
			cerr << "W: the --" << longopts[longindex].name
			     << " option has been deprecated and will be removed in a future version\n";
			break;

		case OPT_HELP:
			cout << option_help;
			exit(EXIT_SUCCESS);

		case OPT_VERSION:
			cout << version_text;
			exit(EXIT_SUCCESS);

		case OPT_BUILD_INFO:
			cout << build_text;
			exit(EXIT_SUCCESS);

		case '?': case ':': default:
			arg_error(argv[0], argv[idx], (c == ':'));
	}

	// Increment idx by the number of args we used and return that number.
	// We must check whether the option argument is in the same argv element
	// as the option name itself, i.e., --opt=arg.
        c = longopts[longindex].has_arg ? 2 : 1;
        if (c == 2) {
                string arg = argv[idx];
                string::size_type p;
                if ((p = arg.rfind(optarg)) != string::npos && arg[p-1] == '=')
                        c = 1;
        }
	idx += c;
	return c;
}

void generate_version_text(void)
{
	version_text.assign(PACKAGE_STRING "\nCopyright (C) 2007-2010 " PACKAGE_AUTHORS ".\n");
	version_text.append(_("License GPLv3+: GNU GPL version 3 or later "
			      "<http://www.gnu.org/licenses/gpl-3.0.html>\n"
			      "This is free software: you are free to change and redistribute it.\n"
			      "There is NO WARRANTY, to the extent permitted by law.\n"));

	ostringstream s;
	s << "Build information:\n";
	s << "  built          : " << BUILD_DATE << " by " << BUILD_USER
	  << '@' << BUILD_HOST << " on " << BUILD_BUILD_PLATFORM
	  << " for " << BUILD_TARGET_PLATFORM << "\n\n"
	  << "  configure flags: " << BUILD_CONFIGURE_ARGS << "\n\n"
	  << "  compiler       : " << BUILD_COMPILER << "\n\n"
	  << "  compiler flags : " << FLDIGI_BUILD_CXXFLAGS << "\n\n"
	  << "  linker flags   : " << FLDIGI_BUILD_LDFLAGS << "\n\n"

	  << "  libraries      : " "FLTK " FLTK_BUILD_VERSION "\n"
	  << "                   " "libsamplerate " << SAMPLERATE_BUILD_VERSION "\n";
#if USE_SNDFILE
	s << "                   " "libsndfile " << SNDFILE_BUILD_VERSION "\n";
#endif
#if USE_PORTAUDIO
	s << "                   " "PortAudio " << PORTAUDIO_BUILD_VERSION "\n";
#endif
#if USE_PULSEAUDIO
	s << "                   " "PulseAudio " << PULSEAUDIO_BUILD_VERSION "\n";
#endif
#if USE_HAMLIB
	s << "                   " "Hamlib " << HAMLIB_BUILD_VERSION "\n";
#endif
#if USE_XMLRPC
	s << "                   " "XMLRPC-C " << XMLRPC_BUILD_VERSION "\n\n";
#endif

	s << "\nRuntime information:\n";
        struct utsname u;
        if (uname(&u) != -1) {
		s << "  system         : " << u.sysname << ' ' << u.nodename
		  << ' ' << u.release << ' ' << u.version << ' ' << u.machine << "\n\n";
	}

	s << "  libraries      : " << src_get_version() << '\n';
#if USE_SNDFILE
	char sndfile_version[32];
	sf_command(NULL, SFC_GET_LIB_VERSION, sndfile_version, sizeof(sndfile_version));
	s << "                   " << sndfile_version << '\n';
#endif
#if USE_PORTAUDIO
	s << "                   " << Pa_GetVersionText() << ' ' << Pa_GetVersion() << '\n';
#endif
#if USE_PULSEAUDIO
	s << "                   " << "Pulseaudio " << pa_get_library_version() << '\n';
#endif
#if USE_HAMLIB
	s << "                   " << hamlib_version << '\n';
#endif

	build_text = s.str();
}

// When debugging is enabled, reexec with malloc debugging hooks enabled, unless
// the env var FLDIGI_NO_EXEC is set, or our parent process is gdb.
void debug_exec(char** argv)
{
#if !defined(NDEBUG) && defined(__GLIBC__)
        if (getenv("FLDIGI_NO_EXEC"))
                return;

	char ppath[32], lname[32];
	ssize_t n;
	snprintf(ppath, sizeof(ppath), "/proc/%u/exe", getppid());
	if ((n = readlink(ppath, lname, sizeof(lname))) > 0) {
		lname[n] = '\0';
		if (strstr(lname, "gdb")) {
                        cerr << "Not using malloc debugging hooks\n";
                        return;
                }
	}

        setenv("FLDIGI_NO_EXEC", "1", 0);
        setenv("MALLOC_CHECK_", "3", 0);
        setenv("MALLOC_PERTURB_", "42", 0);
        if (execvp(*argv, argv) == -1)
                perror("execvp");
#endif
}

void set_platform_ui(void)
{
#if defined(__APPLE__)
       FL_NORMAL_SIZE = 12;
       progdefaults.WaterfallFontsize = 12;
       progdefaults.RxFontsize = 12;
       progdefaults.TxFontsize = 12;
#elif defined(__WOE32__)
       Fl::set_font(FL_HELVETICA, "Tahoma");
       FL_NORMAL_SIZE = 11;
       progdefaults.WaterfallFontnbr = FL_HELVETICA;
       progdefaults.WaterfallFontsize = 12;
       progdefaults.RxFontsize = 12;
       progdefaults.TxFontsize = 12;
#else
       FL_NORMAL_SIZE = 12;
#endif
}

// Convert 1 second of 1-channel silence from IN_RATE Hz to OUT_RATE Hz,
// Repeat test "repeat" times. Return (repeat / elapsed_time),
// the faster-than-realtime factor averaged over "repeat" runs.
// Some figures for SRC_SINC_FASTEST:
// Pentium 4 2.8GHz:     70
// Pentium 3 550MHz:     13
// UltraSparc II 270MHz: 3.5
// Atom N280 1.66GHz:    17.7
#define IN_RATE 48000
#define OUT_RATE 8000
double speed_test(int converter, unsigned repeat)
{
	SRC_DATA src;
	src.src_ratio = (double)OUT_RATE / IN_RATE;
	src.input_frames = IN_RATE;
	src.output_frames = OUT_RATE;
	src.data_in = new float[src.input_frames];
	src.data_out = new float[src.output_frames];

	memset(src.data_in, 0, src.input_frames * sizeof(float));

	// warm up
	src_simple(&src, converter, 1);

	struct timespec t0, t1;
#ifdef _POSIX_MONOTONIC_CLOCK
	clock_gettime(CLOCK_MONOTONIC, &t0);
#else
	clock_gettime(CLOCK_REALTIME, &t0);
#endif
	for (unsigned i = 0; i < repeat; i++)
		src_simple(&src, converter, 1);
#ifdef _POSIX_MONOTONIC_CLOCK
	clock_gettime(CLOCK_MONOTONIC, &t1);
#else
	clock_gettime(CLOCK_REALTIME, &t1);
#endif

	delete [] src.data_in;
	delete [] src.data_out;

	t0 = t1 - t0;
	return repeat / (t0.tv_sec + t0.tv_nsec/1e9);
}

static void setup_signal_handlers(void)
{
#ifndef __WOE32__
	struct sigaction action;
	memset(&action, 0, sizeof(struct sigaction));

	// no child stopped notifications, no zombies
	action.sa_handler = SIG_DFL;
	action.sa_flags = SA_NOCLDSTOP;
#ifdef SA_NOCLDWAIT
	action.sa_flags |= SA_NOCLDWAIT;
#endif
	sigaction(SIGCHLD, &action, NULL);
	action.sa_flags = 0;

	action.sa_handler = handle_signal;
	sigaction(SIGSEGV, &action, NULL);
	sigaction(SIGILL, &action, NULL);
	sigaction(SIGABRT, &action, NULL);
	sigaction(SIGUSR2, &action, NULL);

	action.sa_handler = SIG_IGN;
	sigaction(SIGPIPE, &action, NULL);

	sigemptyset(&action.sa_mask);
	sigaddset(&action.sa_mask, SIGUSR2);
	pthread_sigmask(SIG_BLOCK, &action.sa_mask, NULL);
#else
	signal(SIGSEGV, handle_signal);
	signal(SIGILL, handle_signal);
	signal(SIGABRT, handle_signal);
#endif
}

static void checkdirectories(void)
{
	struct DIRS {
		string& dir;
		const char* suffix;
		void (*new_dir_func)(void);
	};
	DIRS fldigi_dirs[] = {
		{ HomeDir, 0, 0 },
		{ RigsDir, "rigs", 0 },
		{ ScriptsDir, "scripts", 0 },
		{ PalettesDir, "palettes", create_new_palettes },
		{ LogsDir, "logs", 0 },
		{ PicsDir, "images", 0 },
		{ HelpDir, "help", 0 },
		{ MacrosDir, "macros", create_new_macros },
		{ WrapDir, "wrap", 0 },
		{ TalkDir, "talk", 0 },
		{ TempDir, "temp", 0 },
	};

	int r;
	for (size_t i = 0; i < sizeof(fldigi_dirs)/sizeof(*fldigi_dirs); i++) {
		if (fldigi_dirs[i].suffix)
			fldigi_dirs[i].dir.assign(HomeDir).append(fldigi_dirs[i].suffix).append(PATH_SEP);

		if ((r = mkdir(fldigi_dirs[i].dir.c_str(), 0777)) == -1 && errno != EEXIST) {
			cerr << _("Could not make directory") << ' ' << fldigi_dirs[i].dir
			     << ": " << strerror(errno) << '\n';
			exit(EXIT_FAILURE);
		}
		else if (r == 0 && fldigi_dirs[i].new_dir_func)
			fldigi_dirs[i].new_dir_func();
	}

}

bool nbems_dirs_checked = false;

void check_nbems_dirs(void)
{
	if (nbems_dirs_checked) return;

	struct DIRS {
		string& dir;
		const char* suffix;
		void (*new_dir_func)(void);
	};
	DIRS NBEMS_dirs[] = {
		{ NBEMS_dir,     0, 0 },
		{ ARQ_dir,       "ARQ", 0 },
		{ ARQ_files_dir, "ARQ/files", 0 },
		{ ARQ_recv_dir,  "ARQ/recv", 0 },
		{ ARQ_send,      "ARQ/send", 0 },
		{ WRAP_dir,      "WRAP", 0 },
		{ WRAP_recv_dir, "WRAP/recv", 0 },
		{ WRAP_send_dir, "WRAP/send", 0 },
		{ WRAP_auto_dir, "WRAP/auto", 0 },
		{ ICS_dir,       "ICS", 0 },
		{ ICS_msg_dir,   "ICS/messages", 0 },
		{ ICS_tmp_dir,   "ICS/templates", 0 },
	};

	int r;
	for (size_t i = 0; i < sizeof(NBEMS_dirs)/sizeof(*NBEMS_dirs); i++) {
		if (NBEMS_dirs[i].suffix)
			NBEMS_dirs[i].dir.assign(NBEMS_dir).append(NBEMS_dirs[i].suffix).append(PATH_SEP);

		if ((r = mkdir(NBEMS_dirs[i].dir.c_str(), 0777)) == -1 && errno != EEXIST) {
			cerr << _("Could not make directory") << ' ' << NBEMS_dirs[i].dir
			     << ": " << strerror(errno) << '\n';
			exit(EXIT_FAILURE);
		}
		else if (r == 0 && NBEMS_dirs[i].new_dir_func)
			NBEMS_dirs[i].new_dir_func();
	}

	DIRS FLMSG_dirs[] = {
		{ FLMSG_dir,               0, 0 },
		{ FLMSG_WRAP_dir,          "WRAP", 0 },
		{ FLMSG_WRAP_recv_dir,     "WRAP/recv", 0 },
		{ FLMSG_WRAP_send_dir,     "WRAP/send", 0 },
		{ FLMSG_WRAP_auto_dir,     "WRAP/auto", 0 },
		{ FLMSG_ICS_dir,           "ICS", 0 },
		{ FLMSG_ICS_msg_dir,       "ICS/messages", 0 },
		{ FLMSG_ICS_tmp_dir,       "ICS/templates", 0 },
	};

	for (size_t i = 0; i < sizeof(FLMSG_dirs)/sizeof(*FLMSG_dirs); i++) {
		if (FLMSG_dirs[i].dir.empty() && FLMSG_dirs[i].suffix)
			FLMSG_dirs[i].dir.assign(FLMSG_dir).append(FLMSG_dirs[i].suffix).append("/");

		if ((r = mkdir(FLMSG_dirs[i].dir.c_str(), 0777)) == -1 && errno != EEXIST) {
			cerr << _("Could not make directory") << ' ' << FLMSG_dirs[i].dir
			     << ": " << strerror(errno) << '\n';
			exit(EXIT_FAILURE);
		}
		else if (r == 0 && FLMSG_dirs[i].new_dir_func)
			FLMSG_dirs[i].new_dir_func();
	}

	nbems_dirs_checked = true;
}

// Print an error message and exit. If stderr is not a terminal
// show an error dialog instead. On win32 we always use Fl::fatal,
// which displays its own error window.
static void arg_error(const char* name, const char* arg, bool missing)
{
	ostringstream msg;
	msg << name << ": ";
	if (arg && *arg) {
		if (missing)
			msg << "option '" << arg << "' requires an argument\n";
		else
			msg << "unrecognized option '" << arg << "'\n";
	}
	else
		msg << "error while parsing command line\n";

#ifdef __WOE32__
	bool woe32 = true;
#else
	bool woe32 = false;
#endif
	bool tty = isatty(STDERR_FILENO);

	if (tty && !woe32)
		msg << "Try `" << name << " --help' for more information.";
	else
		msg << "See command line help for more information.";
	if (tty || woe32)
		Fl::fatal("%s", msg.str().c_str());
	else {
		fl_message_font(FL_HELVETICA, FL_NORMAL_SIZE);
		fl_alert2("%s", msg.str().c_str());
	}

	exit(EXIT_FAILURE);
}
