bin_PROGRAMS = fusepipe

fusepipe_CPPFLAGS = \
     -I$(top_srcdir)/extra/

fusepipe_LDADD = \
        $(LDFLAGS)

fusepipe_SOURCES = main.cpp sjm.cpp util.cpp

clean:
	rm -rf .deps Makefile.in Makefile *.o ${bin_PROGRAMS}
