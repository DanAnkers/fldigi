// pskmail.cxx
// support for pskmail server/client system

#include <iostream>
#include <string>

#include "main.h"
#include "configuration.h"
#include "fl_digi.h"

using namespace std;

static string mailtext;
string::iterator pText;
static char mailline[1000];

bool pskmail_text_available = false;

void ParseMode(string src)
{
	if (src.find("QPSK31") != string::npos)
		initQPSK31();
	else if (src.find("QPSK63") != string::npos)
		initQPSK63();
	else if (src.find("QPSK125") != string::npos)
		initQPSK125();
	else if (src.find("PSK31") != string::npos)
		initPSK31();
	else if (src.find("PSK63") != string::npos)
		initPSK63();
	else if (src.find("PSK125") != string::npos)
		initPSK125();
	else if (src.find("DOMINOEX4") != string::npos)
		initDOMINOEX4();
	else if (src.find("DOMINOEX5") != string::npos)
		initDOMINOEX5();
	else if (src.find("DOMINOEX8") != string::npos)
		initDOMINOEX8();
	else if (src.find("DOMINOEX11") != string::npos)
		initDOMINOEX11();
	else if (src.find("DOMINOEX16") != string::npos)
		initDOMINOEX16();
	else if (src.find("DOMINOEX22") != string::npos)
		initDOMINOEX22();
	else if (src.find("MFSK8") != string::npos)
		initMFSK8();
	else if (src.find("MFSK16") != string::npos)
		initMFSK16();
	else if (src.find("PTTTUNE") != string::npos)
	{
		int msecs = 100;
		if (src.length() > 7)
			sscanf( src.substr(7, src.length() - 7).c_str(), "%d", &msecs);
		push2talk->set(true);
		MilliSleep(msecs);
		push2talk->set(false);
	}
}

void parse_mailtext()
{
	string strCmdText;
	string strSubCmd;
	unsigned long int idxCmd, idxCmdEnd, idxSubCmd, idxSubCmdEnd;

	idxCmd = mailtext.find("<cmd>");
	idxCmdEnd = mailtext.find("</cmd>");
	
	if ( idxCmd != string::npos && idxCmdEnd != string::npos && idxCmdEnd > idxCmd ) {

		strCmdText = mailtext.substr(idxCmd + 5, idxCmdEnd - idxCmd - 5);
		while ((idxSubCmd = strCmdText.find("<mode>")) != string::npos) {
			idxSubCmdEnd = strCmdText.find("</mode>");
			if (	idxSubCmdEnd != string::npos && 
					idxSubCmdEnd > idxSubCmd ) {
				strSubCmd = strCmdText.substr(idxSubCmd + 6, idxSubCmdEnd - idxSubCmd - 6);
				ParseMode(strSubCmd);
				strCmdText.erase(idxSubCmd, idxSubCmdEnd - idxSubCmd + 7);
			}
		}
		mailtext.erase(idxCmd, idxCmdEnd - idxCmd + 8);
		if (mailtext.length() == 1 && mailtext[0] == '\n')
			mailtext = "";
	}
}

void check_formail() {
	ifstream autofile("gmfsk_autofile");
	if(autofile) {
		mailtext = "";
		while (!autofile.eof()) {
			memset(mailline,0,1000);
			autofile.getline(mailline, 998); // leave space for "\n" and null byte
			mailtext.append(mailline);
			mailtext.append("\n");
		}
		autofile.close();
		std::remove ("gmfsk_autofile");
		
		parse_mailtext();

		if (mailtext.length() > 0) {
			if (mailserver)
				active_modem->set_freq(progdefaults.PSKsweetspot);

			pText = mailtext.begin();
			pskmail_text_available = true;

			active_modem->set_stopflag(false);

			fl_lock(&trx_mutex);
			trx_state = STATE_TX;
			fl_unlock(&trx_mutex);
			wf->set_XmtRcvBtn(true);
		}
	} 
}

void pskmail_loop(void *)
{
	check_formail();
	Fl::repeat_timeout(1.0, pskmail_loop);
}

char pskmail_get_char()
{
	if (pText != mailtext.end())
		return *pText++;

	pskmail_text_available = false;
	return 0x03; // tells psk modem to return to rx
}