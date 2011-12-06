// ----------------------------------------------------------------------------
// fdmdvvaricode.cxx  --  FDMDV Varicode
//
// Copyright (C) 2006
//		Dave Freese, W1HKJ
//
// This file is part of fldigi.  Adapted from code contained in gmfsk source code 
// distribution.
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
// ----------------------------------------------------------------------------


#include <config.h>
#include <map>
#include <initializer_list>

#include "fdmdvvaricode.h"
using namespace std;

map<char,char> create_encodetable()
{
	map<char, char> m;
	m[0x1]='/';
	m[0x2]= '#';
	m[0x3]= ' ';
	m[0xF]= '1';
	m[0x7]= '2';
	m[0xB]= '3';
	m[0xD]= '4';
	m[0x5]= '5';
	m[0x9]= '6';
	m[0xE]= '7';
	m[0x6]= '8';
	m[0xA]= '9';
	m[0x3F]='a';
	m[0x1F]='b';
	m[0x2F]='c';
	m[0x37]='d';
	m[0x17]='e';
	m[0x27]='f';
	m[0x3B]='g';
	m[0x1B]='h';
	m[0x2B]='i';
	m[0x3D]='j';
	m[0x1D]='k';
	m[0x2D]='l';
	m[0x35]='m';
	m[0x15]='n';
	m[0x25]='o';
	m[0x39]='p';
	m[0x19]='q';
	m[0x29]='r';
	m[0x3E]='s';
	m[0x1E]='t';
	m[0x2E]='u';
	m[0x36]='v';
	m[0x16]='w';
	m[0x26]='x';
	m[0x3A]='y';
	m[0x1A]='z';
	m[0x2A]='0';
	return m;
}

map<char,char> create_decodetable()
{
	map<char, char> m;
	m[0x1]='/';
	m[0x2]= '#';
	m[0x3]= ' ';
	m[0xF]= '1';
	m[0x7]= '2';
	m[0xB]= '3';
	m[0xD]= '4';
	m[0x5]= '5';
	m[0x9]= '6';
	m[0xE]= '7';
	m[0x6]= '8';
	m[0xA]= '9';
	m[0x3F]='a';
	m[0x1F]='b';
	m[0x2F]='c';
	m[0x37]='d';
	m[0x17]='e';
	m[0x27]='f';
	m[0x3B]='g';
	m[0x1B]='h';
	m[0x2B]='i';
	m[0x3D]='j';
	m[0x1D]='k';
	m[0x2D]='l';
	m[0x35]='m';
	m[0x15]='n';
	m[0x25]='o';
	m[0x39]='p';
	m[0x19]='q';
	m[0x29]='r';
	m[0x3E]='s';
	m[0x1E]='t';
	m[0x2E]='u';
	m[0x36]='v';
	m[0x16]='w';
	m[0x26]='x';
	m[0x3A]='y';
	m[0x1A]='z';
	m[0x2A]='0';
	return m;
}

// FDMDV Varicode for encoding

char fdmdv_varicode_encode(unsigned char c)
{
	map<char, char> varicodetable = create_encodetable();
	return varicodetable[c];
}

int fdmdv_varicode_decode(unsigned int symbol)
{
	map<char, char> varicodetable = create_decodetable();
	return varicodetable[symbol];
}
