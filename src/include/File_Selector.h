//
// "$Id: File_Selector.H 4473 2005-08-08 00:50:02Z mike $"
//
// File_Selector dialog for the Fast Light Tool Kit (FLTK).
//
// Copyright 1998-2005 by Bill Spitzak and others.
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Library General Public
// License as published by the Free Software Foundation; either
// version 2 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA.
//
// Please report all bugs and problems on the following page:
//
//     http://www.fltk.org/str.php
//

// generated by Fast Light User Interface Designer (fluid) version 1.0107

#ifndef File_Selector_H
#define File_Selector_H
#include <FL/Fl.H>
#include <FL/Fl_Double_Window.H>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <FL/Fl_Group.H>
#include <FL/Fl_Choice.H>
#include <FL/Fl_Menu_Button.H>
#include <FL/Fl_Button.H>
#include <FL/Fl_Tile.H>
#include <FL/Fl_File_Browser.H>
#include <FL/Fl_Box.H>
#include <FL/Fl_Check_Button.H>
#include <FL/Fl_File_Input.H>
#include <FL/Fl_Return_Button.H>
#include <FL/fl_ask.H>

class File_Selector {
public:
  enum { SINGLE = 0, MULTI = 1, CREATE = 2, DIRECTORY = 4 };
private:
  void (*callback_)(File_Selector*, void *);
  void *data_;
  char directory_[1024];
  char pattern_[1024];
  char preview_text_[2048];
  int type_;
  void fileListCB();
  void fileNameCB();
  void newdir();
  static void previewCB(File_Selector *fc);
  void showChoiceCB();
  void update_favorites();
  void update_preview();
public:
  File_Selector(const char *d, const char *p, int t, const char *title);
private:
  Fl_Double_Window *window;
  void cb_window_i(Fl_Double_Window*, void*);
  static void cb_window(Fl_Double_Window*, void*);
  Fl_Choice *showChoice;
  void cb_showChoice_i(Fl_Choice*, void*);
  static void cb_showChoice(Fl_Choice*, void*);
  Fl_Menu_Button *favoritesButton;
  void cb_favoritesButton_i(Fl_Menu_Button*, void*);
  static void cb_favoritesButton(Fl_Menu_Button*, void*);
public:
  Fl_Button *newButton;
private:
  void cb_newButton_i(Fl_Button*, void*);
  static void cb_newButton(Fl_Button*, void*);
  void cb__i(Fl_Tile*, void*);
  static void cb_(Fl_Tile*, void*);
  Fl_File_Browser *fileList;
  void cb_fileList_i(Fl_File_Browser*, void*);
  static void cb_fileList(Fl_File_Browser*, void*);
  Fl_Box *previewBox;
public:
  Fl_Check_Button *previewButton;
private:
  void cb_previewButton_i(Fl_Check_Button*, void*);
  static void cb_previewButton(Fl_Check_Button*, void*);
  Fl_File_Input *fileName;
  void cb_fileName_i(Fl_File_Input*, void*);
  static void cb_fileName(Fl_File_Input*, void*);
  Fl_Return_Button *okButton;
  void cb_okButton_i(Fl_Return_Button*, void*);
  static void cb_okButton(Fl_Return_Button*, void*);
  Fl_Button *cancelButton;
  void cb_cancelButton_i(Fl_Button*, void*);
  static void cb_cancelButton(Fl_Button*, void*);
  Fl_Double_Window *favWindow;
  Fl_File_Browser *favList;
  void cb_favList_i(Fl_File_Browser*, void*);
  static void cb_favList(Fl_File_Browser*, void*);
  Fl_Button *favUpButton;
  void cb_favUpButton_i(Fl_Button*, void*);
  static void cb_favUpButton(Fl_Button*, void*);
  Fl_Button *favDeleteButton;
  void cb_favDeleteButton_i(Fl_Button*, void*);
  static void cb_favDeleteButton(Fl_Button*, void*);
  Fl_Button *favDownButton;
  void cb_favDownButton_i(Fl_Button*, void*);
  static void cb_favDownButton(Fl_Button*, void*);
  Fl_Button *favCancelButton;
  void cb_favCancelButton_i(Fl_Button*, void*);
  static void cb_favCancelButton(Fl_Button*, void*);
  Fl_Return_Button *favOkButton;
  void cb_favOkButton_i(Fl_Return_Button*, void*);
  static void cb_favOkButton(Fl_Return_Button*, void*);
public:
  ~File_Selector();
  void callback(void (*cb)(File_Selector *, void *), void *d = 0);
  void color(Fl_Color c);
  Fl_Color color();
  int count();
  void directory(const char *d);
  char * directory();
  void filter(const char *p);
  const char * filter();
  int filter_value();
  void filter_value(int f);
  void hide();
  void iconsize(uchar s);
  uchar iconsize();
  void label(const char *l);
  const char * label();
  void ok_label(const char *l);
  const char * ok_label();
  void preview(int e);
  int preview() const { return previewButton->value(); };
  void rescan();
  void show();
  int shown();
  void textcolor(Fl_Color c);
  Fl_Color textcolor();
  void textfont(uchar f);
  uchar textfont();
  void textsize(uchar s);
  uchar textsize();
  void type(int t);
  int type();
  void * user_data() const;
  void user_data(void *d);
  const char *value(int f = 1);
  void value(const char *filename);
  int visible();
  static const char *add_favorites_label;
  static const char *all_files_label;
  static const char *custom_filter_label;
  static const char *existing_file_label;
  static const char *favorites_label;
  static const char *filename_label;
  static const char *filesystems_label;
  static const char *manage_favorites_label;
  static const char *new_directory_label;
  static const char *new_directory_tooltip;
  static const char *preview_label;
  static const char *save_label;
  static const char *show_label;
  static Fl_File_Sort_F *sort;
};

extern char *						// O - Directory or NULL
Dir_Chooser(const char *message,	// I - Message for titlebar
            const char *fname,		// I - Initial directory name
			int        relative);	// I - 0 for absolute

extern char *						// O - Filename or NULL
File_Select(const char *message,	// I - Message in titlebar
            const char *pat,		// I - Filename pattern
			const char *fname,		// I - Initial filename selection
			int       relative); 	// I - 0 for absolute path

#endif

//
// End of "$Id: File_Selector.H 4473 2005-08-08 00:50:02Z mike $".
//