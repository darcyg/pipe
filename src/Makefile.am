bin_PROGRAMS = fusepipe

fusepipe_CPPFLAGS = \
     -I$(top_srcdir)/extra/

fusepipe_LDADD = \
        $(LDFLAGS)

fusepipe_SOURCES = main.cpp genjob.cpp genpipe.cpp job.cpp options.cpp pipeline.cpp task.cpp

clean:
	rm -rf .deps Makefile.in Makefile *.o ${bin_PROGRAMS}
