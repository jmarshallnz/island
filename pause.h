/********************************************/
/*	pause.h 23rd February 2005				*/
/*	(c) Danny Wilson.						*/
/*	www.danielwilson.me.uk					*/
/********************************************/

#ifndef _MYUTILS_PAUSE_H_
#define _MYUTILS_PAUSE_H_

#ifdef _WIN32

	#include <conio.h>
	#include <stdio.h>

	namespace myutils
	{
		inline void pause()
		{
			printf("\nPress any key\n");
			int ch=-99;
			while (ch==-99)
				ch=_getch();
		}
		inline void silent_pause()
		{
			int ch=-99;
			while (ch==-99)
				ch=_getch();
		}
	};

#else

	namespace myutils
	{
		inline void pause() {}
	};

#endif

#endif